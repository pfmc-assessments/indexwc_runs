library(dplyr)
library(indexwc)

savedir <- here::here("indices")
# The configuration file is a rda file in the package
configuration_all <- indexwc::configuration

#===============================================================================
# NWFSC Combo
#===============================================================================
benchmark_list <- c(
  "lingcod",
  "redbanded rockfish"
)
update_list <- c(
  "petrale sole",
  "shortspine thornyhead"
)
species_list <- c(benchmark_list, update_list)

configuration_sub <- dplyr::bind_rows(
  configuration_all |>
    dplyr::filter(
      source == "NWFSC.Combo",
      species %in% benchmark_list,
      !formula %in%
        c(
          "catch_weight ~ 0 + fyear*split_mendocino + pass_scaled",
          "catch_weight ~ 0 + fyear*split_mendocino"
        )
    ),
  configuration_all |>
    dplyr::filter(source == "Triennial", species %in% benchmark_list),
  configuration_all |>
    dplyr::filter(
      source == "NWFSC.Combo",
      species %in% update_list,
      used == TRUE
    )
)

for (sp in species_list) {
  configuration_to_run <- configuration_sub |>
    dplyr::filter(species == sp)

  for (run in 5:nrow(configuration_to_run)) {
    my_data <- pull_and_format_data(
      configuration_to_run = configuration_to_run[run, ]
    )

    fit <- run_sdmtmb(
      data = my_data$data_filtered[[1]],
      family = my_data$family,
      formula = my_data$formula,
      n_knots = my_data$knots,
      share_range = my_data$share_range,
      anisotropy = my_data$anisotropy,
      spatiotemporal = list(my_data$spatiotemporal1, my_data$spatiotemporal2)
    )

    diag <- diagnose(
      fit = fit
    )

    index <- calc_index_areas(
      data = fit$data,
      fit = fit,
      boundaries = "Coastwide"
    )

    save_index_outputs(
      fit = fit,
      diagnostics = diag,
      indices = index,
      dir = here::here("indices"),
      overwrite = FALSE
    )
  }
}

#===============================================================================
# Special runs
#===============================================================================
benchmark_list <- c(
  "lingcod"
)
species_list <- c(benchmark_list)

configuration_sub <-
  configuration_all |>
  dplyr::filter(
    source %in% c("NWFSC.Combo", "Triennial"),
    species %in% benchmark_list,
    formula %in%
      c(
        "catch_weight ~ 0 + fyear*split_mendocino + pass_scaled",
        "catch_weight ~ 0 + fyear*split_mendocino"
      )
  ) |>
  dplyr::arrange(source)


for (run in 1:nrow(configuration_sub)) {
  my_data <- pull_and_format_data(
    configuration_to_run = configuration_sub[run, ]
  )
  my_data$data_filtered[[1]] <- my_data$data_filtered[[1]] |>
    dplyr::mutate(split_mendocino = ifelse(latitude > 40.1666667, "N", "S"))

  check <- dplyr::filter(my_data$data_filtered[[1]], catch_weight > 0) |>
    dplyr::group_by(split_mendocino, year) |>
    dplyr::summarise(n = dplyr::n())

  lm <- lm(
    formula = as.formula(configuration_sub[run, "formula"]),
    data = my_data$data_filtered[[1]]
  )
  lm_pos <- lm(
    formula = as.formula(configuration_sub[run, "formula"]),
    data = dplyr::filter(my_data$data_filtered[[1]], catch_weight > 0)
  )
  map_pres <- names(coef(lm))
  map_pres <- factor(map_pres)
  start_pres <- rep(0, length(map_pres))
  map_pos <- names(coef(lm_pos))
  map_pos <- factor(map_pos)
  start_pos <- rep(0, length(map_pos))

  fit <- run_sdmtmb(
    data = my_data$data_filtered[[1]],
    family = my_data$family,
    formula = my_data$formula,
    n_knots = my_data$knots,
    share_range = my_data$share_range,
    anisotropy = my_data$anisotropy,
    spatiotemporal = list(my_data$spatiotemporal1, my_data$spatiotemporal2),
    sdmtmb_control = sdmTMB::sdmTMBcontrol(
      map = list(
        b_j = map_pres,
        b_j2 = map_pos
      ),
      start = list(b_j = start_pres, b_j2 = start_pos),
      newton_loops = 1
    )
  )

  prediction_grid_north <- lookup_grid(
    x = fit$data[["survey_name"]][1],
    max_latitude = fit$ranges$latitude_max,
    min_latitude = 40.167,
    max_longitude = fit$ranges$longitude_max,
    min_longitude = fit$ranges$longitude_min,
    max_depth = abs(fit$ranges$depth_max),
    years = sort(unique(fit$data$year)),
    data = california_current_grid
  )

  prediction_grid_south <- lookup_grid(
    x = fit$data[["survey_name"]][1],
    max_latitude = 40.167,
    min_latitude = fit$ranges$latitude_min,
    max_longitude = fit$ranges$longitude_max,
    min_longitude = fit$ranges$longitude_min,
    max_depth = abs(fit$ranges$depth_max),
    years = sort(unique(fit$data$year)),
    data = california_current_grid
  )

  diag_north <- diagnose(
    fit = fit,
    prediction = prediction_grid_north
  )

  diag_south <- diagnose(
    fit = fit,
    prediction = prediction_grid_south
  )

  index_north <- calc_index_areas(
    data = fit$data,
    fit = fit,
    prediction_grid = prediction_grid_north,
    boundaries = "Coastwide"
  )

  index_south <- calc_index_areas(
    data = fit$data,
    fit = fit,
    prediction_grid = prediction_grid_south,
    boundaries = "Coastwide"
  )

  save_index_outputs(
    fit = fit,
    diagnostics = diag_north,
    indices = index_north,
    dir = here::here("indices", "lingcod-north-split_mendocino"),
    overwrite = TRUE
  )

  save_index_outputs(
    fit = fit,
    diagnostics = diag_south,
    indices = index_south,
    dir = here::here("indices", "lingcod-south-split_mendocino"),
    overwrite = TRUE
  )
}

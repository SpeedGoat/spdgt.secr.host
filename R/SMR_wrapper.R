#' Spatial Mark-Recapture Model Wrapper
#'
#' @description
#' A wrapper function that coordinates the parallel execution of multiple spatial
#' mark-recapture (SMR) models with different combinations of telemetry and
#' mobility settings. The function handles model selection, parallel processing,
#' and result aggregation.
#'
#' @param data List. Contains capture-recapture and spatial data:
#' \itemize{
#'   \item X1: Optional matrix of trap locations
#'   \item locs: Telemetry locations
#'   \item area: Numeric study area size (for density calculation)
#'   \item Other components required by mcmc_SMR
#' }
#' @param input List. Contains model settings:
#' \itemize{
#'   \item model_choices: Vector of model IDs to run (1-4)
#'   \item mobile_center: Will be set based on model choice
#'   \item Other parameters required by mcmc_SMR
#' }
#'
#' @return A tibble with columns:
#' \itemize{
#'   \item model: Character. Model description
#'   \item Parameter: Character. Parameter name
#'   \item Mean: Numeric. Parameter estimate
#' }
#'
#' @export
SMR_wrapper <- function(data, input) {

  model_choices <- unlist(input$model_choices)

  # Mobile cannot be used when traps are not available
  if (is.null(data$X1)) {
    if (any(model_choices %in% c(1, 3))) {
      warning("Mobile cannot be used when traps are not available.
              Using non-mobile model.")
    }

    model_choices <- model_choices[!(model_choices %in% c(1, 3))]
  }

  # Validate inputs
  if (length(model_choices) == 0) {
    stop("At least one model must be selected")
  }

  # Define model options
  all_models <- tibble::tibble(
    ID = 1:4,
    Name = c(
      "Tele and Mobile",
      "Tele and No Mobile",
      "No Tele and Mobile",
      "No Tele and No Mobile"
    )
  )

  # Filter IDs based on selected combinations
  selected_models <- all_models %>%
    dplyr::filter(.data$ID %in% model_choices) %>%
    dplyr::pull(.data$ID)

  my_cluster <- parallel::makeCluster(6, type = "PSOCK")
  #register cluster to be used by %dopar%
  doParallel::registerDoParallel(cl = my_cluster)

  system.time({
    all_out <- foreach::foreach(
      mm = selected_models,
      .packages = c("dplyr", "spdgt.secr.host")
    ) %dopar% {

      if (mm == 1) {
        # Tele, Mobile
        input$mobile_center <- T

        out <- mcmc_SMR(data, input)

        tibble::as_tibble(out$out) %>%
          dplyr::mutate(model = all_models$Name[mm])

      } else if (mm == 2) {
        # Tele, No Mobile
        input$mobile_center <- F

        out <- mcmc_SMR(data, input)

        tibble::as_tibble(out$out) %>%
          dplyr::mutate(model = all_models$Name[mm])
      } else if (mm == 3) {
        # No Tele, Mobile
        data_no_tele <- data
        data_no_tele$locs <- NA

        input$mobile_center <- T

        out <- mcmc_SMR(data_no_tele, input)

        tibble::as_tibble(out$out) %>%
          dplyr::mutate(model = all_models$Name[mm])
      } else if (mm == 4) {
        # No Tele, No Mobile
        data_no_tele <- data
        data_no_tele$locs <- NA

        input$mobile_center <- F

        out <- mcmc_SMR(data_no_tele, input)

        tibble::as_tibble(out$out) %>%
          dplyr::mutate(model = all_models$Name[mm])
      }
    }
  }) # Time stop
  # Stop cluster
  parallel::stopCluster(my_cluster)

  final_out <- dplyr::bind_rows(all_out) %>%
    dplyr::mutate(
      D = .data$N * 100 / .data$area
    ) %>%
    tidyr::pivot_longer(
      cols = -.data$model,
      names_to = "Parameter",
      values_to = "Mean"
    )

  return(final_out)
}


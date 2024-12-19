SMR_wrapper <- function(data, input) {

  model_choices <- unlist(input$model_choices)

  # Mobile cannot be used when traps are not available
  if (is.null(data$X1)) {
    if (any(model_choices %in% c(1, 3))) {
      warning("Mobile cannot be used when traps are not available.
              Using non-mobile model.")
    }

    model_choices <- model_choices %>%
      dplyr::filter(!(ID %in% c(1, 3)))
  }

  # Validate inputs
  if (length(model_choices) == 0) {
    stop("At least one option must be chosen")
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
    dplyr::filter(ID %in% model_choices) %>%
    dplyr::pull(ID)

  my_cluster <- makeCluster(6, type = "PSOCK")
  #register cluster to be used by %dopar%
  doParallel::registerDoParallel(cl = my_cluster)

  system.time({
    all_out <- foreach::foreach(
      model = selected_models,
      .packages = c("dplyr")
    ) %dopar% {

      # Remove this and use spdgt.secr.host package once it's built
      devtools::load_all()

      if (model == 1) {
        # Tele, Mobile
        input$mobile_center <- T

        out <- mcmc_SMR(data, input)

        tibble::as_tibble(out$out) %>%
          dplyr::mutate(model = all_models$Name[model])

      } else if (model == 2) {
        # Tele, No Mobile
        input$mobile_center <- "no move"

        out <- mcmc_SMR(data, input)

        tibble::as_tibble(out$out) %>%
          dplyr::mutate(model = all_models$Name[model])
      } else if (model == 3) {
        # No Tele, Mobile
        data_no_tele <- data
        data_no_tele$locs <- NA

        input$mobile_center <- "move"

        out <- mcmc_SMR(data_no_tele, input)

        tibble::as_tibble(out$out) %>%
          dplyr::mutate(model = all_models$Name[model])
      } else if (model == 4) {
        # No Tele, No Mobile
        data_no_tele <- data
        data_no_tele$locs <- NA

        input$mobile_center <- "no move"

        out <- mcmc_SMR(data_no_tele, input)

        tibble::as_tibble(out$out) %>%
          dplyr::mutate(model = all_models$Name[model])
      }
    }
  }) # Time stop
  # Stop cluster
  stopCluster(my_cluster)

  final_out <- dplyr::bind_rows(all_out) %>%
    dplyr::mutate(
      D = N * 100 / area
    ) %>%
    tidyr::pivot_longer(
      cols = -model,
      names_to = "Parameter",
      values_to = "Mean"
    )

  return(final_out)
}


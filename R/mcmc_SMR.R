# NEED TO INCORPORATE CONVENTIONAL MODEL (OMIT ALL MARKING PROCESSES)


#' Fit the generalized categorical spatial mark resight model. SMR_move fits the activity centers to cameras and snares separately whereas SMR fits them concurrently.
#'
#' Allowing for activity center relocation (mobile_center = "move") means the model
#' incorporates a process for updating activity centers between the marking and
#' sighting phases. If the data are sparse, there is no telemetry data, or if
#' the study period is short, it may better to not allow for activity center
#' relocation (mobile_center = "no move").
#'
#' @param data a data list
#' @param input Model definitions and inputs
#' @description This function fits the generalized categorical spatial mark resight model.
#' Modelling the marking process relaxes the assumption that the distribution of marked
#' individuals across the landscape is spatially uniform.
#'
#' the data list should all elements of that object are necessary.
#' y_mark, y_sight_marked, y_sight_unmarked, G_marked, and G_unmarked are necessary
#' list elements. y_sight_x and G_x for x=unk and marke.noID are necessary if there are samples
#' of unknown marked status or samples from marked samples without individual identities.
#'
#' "niter" the number of MCMC iterations to perform
#' "nburn" the number of MCMC iterations to discard as burnin
#' "nthin" the MCMC thinning interval. Keep every nthin iterations.
#' "M" the level of data augmentation
#' "inits" a list of initial values for lam0_mark,lam0_sight, sigma, gamma, and psi.
#' The list element for gamma is itself a list with ncat elements. See the example below.
#' "obstype" a vector of length two indicating the observation model, "bernoulli" or "poisson", for the
#' marking and sighting process
#' "nswap" an integer indicating how many samples for which the latent identities
#' are updated on each iteration.
#' "proppars" a list of proposal distribution tuning parameters for lam0_mark,
#' lam0_sight, sigma, s, and st, for the activity centers of untelemetered and
#' telemetered individuals, respectively. The tuning parameter should be
#' smaller for individuals with telemetry and increasingly so as the number of
#' locations per individual increases
#' "storeLatent" a logical indicator for whether or not the posteriors of the
#' latent individual identities, z, and s are stored and returned
#' "storeGamma" a logical indicator for whether or not the posteriors for gamma are stored and returned
#' "IDup" a character string indicating whether the latent identity update is
#' done by Gibbs or Metropolis-Hastings, "Gibbs", or "MH". For
#' obstype="bernoulli", only "MH" is available because the full conditional is not known.
#' "tf1" a trap operation vector or matrix for the marking process. If exposure to capture does
#' not vary by individual, tf1 should be a vector of length J1 indicating how many of the K1 occasions
#' each marking location was operational. If exposure to capture varies by individual or by trap and
#' individual, tf1 should be a matrix of dimension M x J1 indicating how many
#' of the K1 occasions individual i was exposed to at trap j. This allows known
#' additions or removals during the marking process to be accounted for.
#' Exposure for the n_marked+1 ... M uncaptured individuals should be the
#' same as the number of occasions each trap was operational. We can't account for unknown
#' additions and removals.
#' "tf2" a trap operation vector or matrix for the sighting process. If exposure to capture does
#' not vary by indiviudal, tf2 should be a vector of length J2 indicating how many of the K2 occasions
#' each sighting location was operational. If exposure to capture varies by individual or by trap and
#' individual, tf2 should be a matrix of dimension M x J2 indicating how many
#' of the K2 occasions individual i was exposed to at trap j. This allows
#' known additions or removals between the marking and sighting processes and
#' during the sighting process to be accounted for. Exposure for
#' the n_marked+1 ... M uncaptured individuals should be the
#' same as the number of occasions each trap was operational. We can't account for unknown
#' additions and removals.
#' "X1" is a matrix of marking coordinates
#' "X2" is a matrix of sighting coordinates, ,
#' "K1", the integer number of marking occasions, and an element "K2",
#' the integer number of sighting occasions are necessary.
#' "IDlist" is a list containing elements ncat and IDcovs. ncat is an integer for the number
#' of categorical identity covariates and IDcovs is a list of length ncat with elements containing the
#' values each categorical identity covariate may take.
#' "locs" is an n_marked x nloc x  2 array of telemetry locations is optional. This array can
#' have missing values if not all individuals have the same number of locations
#' and the entry for individuals with no telemetry should all be missing values (coded NA).
#' "markedS" is required if marking and sighting sessions are interspersed. This is a
#' n_marked x K2 matrix with 0 indicating an individual was not marked on occasion k and 1 if it
#' was.
#' "lam0_mark" is the coefficient determining the probability of marking an individual
#' "lam0_sight" is the coefficient determining the probability of sighting an individual
#' "sigma_d" is the variance for sighting probability distributions
#' "sigma_p" is the variance for marking probability distributions
#' "gamma" is a list of the category level probabilities of the same dimension
#' as IDcovs. The category level probabilities for each covariate must sum to 1
#' "psi" is probability of keeping augmented activity centers
#' "s1" and "s2" are the activity centers for marked and sighted augmented data
#' "proppars$s" is variance for s
#' "proppars$st" is variance for s w/ telemetry data
#'
#' @export
#'
mcmc_SMR <- function(
  data,
  input = list(
    niter = 2400, nburn = 1200, nthin = 5, M = 200, mobile_center = F,
    inits = list(lam0_mark = 0.05, lam0_sight = 0.009, sigma_d = 0.1,
                 sigma_p = 0.1, s1 = 5, s2 = 5, psi = 0.5,
                 gamma = vector("list", 1)),
    obstype = c("bernoulli", "poisson"), nswap = NA,
    proppars = list(lam0_mark = 0.05, lam0_sight = 0.1, sigma_d = 0.02,
                    sigma_p = 0.2, s1 = 0.5, s2 = 0.25, s2t = 0.1),
    max_proppars <- list(lam0_mark = 100, lam0_sight = 100, sigma_d = 100,
                         sigma_p = 100, s1 = 100, s2 = 100, s2t = 100),
    min_proppars <- list(lam0_mark = .001, lam0_sight = .001, sigma_d = .001,
                         sigma_p = .001, s1 = .001, s2 = .001, s2t = .001),
    storeLatent = TRUE, storeGamma = TRUE, IDup = "Gibbs", tf1 = NA, tf2 = NA
  )
){

  # Toggle activity center movement
  mobile_center <- input$mobile_center

  # Retrieve input parms
  niter <- input$niter
  nburn <- input$nburn
  nthin <- input$nthin

  # Data augmentation
  M <- input$M

  # List of initial values for lam0, sigma, gamma, and psi
  inits <- if (is.null(input$inits)) NA else input$inits

  # Character string indicating the observation model, "bernoulli" or "poisson".
  obstype <- input$obstype

  # An integer indicating how many samples for which the latent identities
  # are updated on each iteration (default NA)
  nswap <- if (is.null(input$nswap)) NA else input$nswap

  # List of proposal distribution tuning parameters for lam0, sigma, s, and st
  proppars <- input$proppars
  max_proppars <- input$max_proppars
  min_proppars <- input$min_proppars

  # A logical indicator to store posteriors for the latent individual identities
  # and gamma
  storeLatent <- input$storeLatent
  storeGamma <- input$storeGamma

  # Character string indicating whether the latent identity update is done by
  # Gibbs or Metropolis-Hastings, "Gibbs", or "MH"
  # For obstype = "bernoulli", "MH" must be used
  IDup <- input$IDup
  # priors <- input$priors

  # pull out initial values
  psi <- inits$psi
  lam0_mark <- inits$lam0_mark
  lam0_sight <- inits$lam0_sight
  sigma_d <- inits$sigma_d
  sigma_p <- inits$sigma_p
  gamma <- inits$gamma

  # Check data and inputs
  data <- data_check(data)
  input_check(input)

  # Initializations
  y_mark <- data$y_mark
  y_sight_marked <- data$y_sight_marked
  y_sight_unmarked <- data$y_sight_unmarked
  X1 <- data$X1
  X2 <- data$X2
  J1 <- nrow(X1)
  J2 <- nrow(X2)
  K1 <- data$K1
  K2 <- data$K2
  ncat <- data$IDlist$ncat
  nallele <- data$IDlist$nallele
  IDcovs <- data$IDlist$IDcovs
  buff <- data$buff
  vertices <- data$vertices
  Xall <- rbind(X1, X2)
  n_marked <- data$n_marked
  G_marked <- data$G_marked
  G_marked_noID <- data$G_marked_noID
  tf1 <- data$tf1
  tf2 <- data$tf2
  y_sight_marked_noID <- data$y_sight_marked_noID
  y_sight_unk <- data$y_sight_unk

  useUm <- init_useUM(data, ncat)
  nlevels <- unlist(lapply(IDcovs, length))

  # Define state space with buffer or vertices
  state_space <- define_state_space(data = data, Xall = Xall)

  #### Acceptance rate tuning
  tune_parms <- c("lam0_mark", "lam0_sight", "sigma_p", "sigma_d", "s1", "s2")

  accept_rates <- tibble::tibble(
    acc_iter = rep(1:niter, each = length(proppars)),
    label = rep(names(proppars), niter),
    accept = dplyr::case_when(.data$label %in% tune_parms ~ 0,
                              .default = 0.45)
  )

  tune_check <- 50
  batch_n <- 10

  # Check for unknown marked individuals
  useUnk <- FALSE
  if ("G_unk" %in% names(data)) {
    if (!is.na(data$G_unk[1])) {
      G_unk <- data$G_unk
      if (!is.matrix(G_unk)) {
        G_marked <- matrix(G_unk)
      }
      if (ncol(G_unk) != ncat) {
        stop("G_unk needs ncat number of columns")
      }
      y_sight_unk <- data$y_sight_unk
      useUnk <- TRUE
    }
  }

  # Check for marked no ID individuals
  useMarkednoID <- FALSE
  if ("G_marked_noID" %in% names(data)) {
    if (!is.na(data$G_marked_noID[1])) {
      G_marked_noID <- data$G_marked_noID
      if (!is.matrix(G_marked_noID)) {
        G_marked_noID <- matrix(G_marked_noID)
      }
      if (ncol(G_marked_noID) != ncat) {
        stop("G_marked_noID needs ncat number of columns")
      }
      y_sight_marked_noID <- data$y_sight_marked_noID
      useMarkednoID <- TRUE
    }
  }

  G_data <- combine_genetic_data(
    data,
    ncat,
    nlevels,
    useUnk,
    useMarkednoID
  )

  G_use <- G_data$G_use
  G_marked <- G_data$G_marked
  y_sight_latent <- G_data$y_sight_latent
  ncat <- G_data$ncat
  nlevels <- G_data$nlevels
  status <- G_data$status

  # number of unmarked individuals + uncertainty
  n_samp_latent <- nrow(y_sight_latent)
  if (is.na(nswap)) {
    nswap <- round(n_samp_latent / 2)
    # warning("nswap not specified, using round(n_samp_latent/2)")
  }

  # make constraints for data initialization
  all_constraints <- generate_constraints(y_sight_latent, obstype[2], G_use)
  constraints <- all_constraints$constraints
  binconstraints <- all_constraints$binconstraints

  # marking occasion order constraints
  Kconstraints <- generate_Kconstraints(
    data,
    y_sight_latent,
    n_marked,
    M,
    n_samp_latent,
    K2
  )

  # Build y_sight_true
  y_sight_true <- array(0, dim = c(M, J2, K2))
  y_sight_all <- initialize_capture_histories(
    y_sight_latent, y_sight_true, y_sight_marked, y_mark, status, ID, M,
    n_marked, n_samp_latent, X1, X2, G_marked, G_marked_noID, constraints,
    Kconstraints, useMarkednoID, data$markedS
  )
  y_sight_true <- y_sight_all$y_sight_true
  ID <- y_sight_all$ID

  if (binconstraints) {
    if (any(y_sight_true > 1)) stop("bernoulli data not initialized correctly")
  }

  # Reduce y_mark dimension
  if (!is.null(y_mark)) {
    y_mark <- combine_matrices(y_mark, array(0, dim = c(M - n_marked, J1, K1)))
    y_mark2D <- apply(y_mark, c(1, 2), sum)
  }

  # Check marked guy k constraints
  y_sight_true <- apply(y_sight_true, c(1, 2), sum)
  known_vector <- c(rep(1, n_marked), rep(0, M - n_marked))
  known_vector[(n_marked + 1):M] <- 1 * (rowSums(y_sight_true[(n_marked + 1):M, ]) > 0)

  # Initialize z
  z <- 1 * (known_vector > 0)
  add <- M * (0.5 - sum(z) / M)
  if (add > 0) {
    z[sample(which(z == 0), add)] <- 1 # switch some uncaptured z's to 1.
  }
  unmarked <- c(rep(FALSE, n_marked), rep(TRUE, M - n_marked))

  # Optimize starting locations given where they are trapped.
  s1 <- process_spatial_points(M, state_space, y_sight_true, Xall,
                               y_mark)
  s2 <- s1

  # collapse unmarked data to 2D
  y_sight_latent <- apply(y_sight_latent, c(1, 2), sum)

  # Initialize G_true
  G_true_dat <- processGeneticData(M, ncat, n_marked, G_marked, ID, G_use,
                                   IDcovs, gamma, useUnk, useMarkednoID)
  G_true <- G_true_dat$G_true
  G_use <- G_true_dat$G_use
  Mark_obs <- G_true_dat$Mark_obs
  G_latent <- G_true_dat$G_latent

  #### MCMC initializations
  # some objects to hold the MCMC output
  nstore <- (niter - nburn) / nthin
  if (nburn %% nthin != 0) {
    nstore <- nstore + 1
  }

  out <- matrix(NA, nrow = nstore, ncol = 7)
  dimnames(out) <- list(NULL, c("lam0_mark", "lam0_sight", "sigma_d", "sigma_p", "N", "n.um", "psi"))
  if (storeLatent) {
    s1xout <- s1yout <- s2xout <- s2yout <- zout <- matrix(NA, nrow = nstore, ncol = M)
    IDout <- matrix(NA, nrow = nstore, ncol = length(ID))
  }
  idx <- 1 # for storing output not recorded every iteration
  if (storeGamma) {
    gammaOut <- vector("list", ncat)
    for (i in 1:ncat) {
      gammaOut[[i]] <- matrix(NA, nrow = nstore, ncol = nlevels[i])
      colnames(gammaOut[[i]]) <- paste("Lo", i, "G", 1:nlevels[i], sep = "")
    }
  }
  if (!is.na(data$locs[1])) {
    uselocs <- TRUE
    locs <- data$locs
    telguys <- which(rowSums(!is.na(locs[, , 1])) > 0)
    ll_tel <- matrix(0, nrow = dim(locs)[1], ncol = dim(locs)[2])
    # update starting locations using telemetry data
    for (i in telguys) {
      s2[i, ] <- c(mean(locs[i, , 1], na.rm = TRUE), mean(locs[i, , 2], na.rm = TRUE))
    }
    for (i in telguys) {
      ll_tel[i, ] <- dnorm(locs[i, , 1], s2[i, 1], sigma_d, log = TRUE) +
        dnorm(locs[i, , 2], s2[i, 2], sigma_d, log = TRUE)
    }
    ll_tel_cand <- ll_tel
  } else {
    uselocs <- FALSE
    telguys <- c()
  }

  if (!is.null(y_mark)) {
    if (!any(is.na(tf1))) {
      if (any(tf1 > K1)) {
        stop("Some entries in tf1 are greater than K1.")
      }
      if (is.null(dim(tf1))) {
        if (length(tf1) != J1) {
          stop("2D tf1 vector must be of length J1.")
        }
        K2D1 <- matrix(rep(tf1, M), nrow = M, ncol = J1, byrow = TRUE)
        # warning("Since 1D tf1 entered, assuming all individuals exposed to equal capture")
      } else {
        if (!all(dim(tf1) == c(M, J1))) {
          stop("tf1 must be dim M by J1 if tf1 varies by individual")
        }
        K2D1 <- tf1
        # warning("Since 2D tf1 entered, assuming individual exposure to traps differ")
      }
    } else {
      tf1 <- rep(K1, J1)
      K2D1 <- matrix(rep(tf1, M), nrow = M, ncol = J1, byrow = TRUE)
    }

    # Make sure all K2D1 is >= all y_mark2D
    K2D1[K2D1 < max(y_mark2D)] <- max(y_mark2D)
  }

  if (!any(is.na(tf2))) {
    if (any(tf2 > K2)) {
      stop("Some entries in tf2 are greater than K2.")
    }
    if (is.null(dim(tf2))) {
      if (length(tf2) != J2) {
        stop("tf2 vector must be of length J2.")
      }
      # if (!all(dim(K2D2) == c(M, J2))) {
      #   stop("K2D2 must be dim M by J2 if K2D2 varies by individual")
      # }
      K2D2 <- matrix(rep(tf2, M), nrow = M, ncol = J2, byrow = TRUE)
      # warning("Since 1D tf2 entered, assuming all individuals exposed to equal sighting effort")
    } else {
      K2D2 <- tf2
      # warning("Since 2D tf2 entered, assuming individual exposure to sighting effort differs")
    }
  } else {
    tf2 <- rep(K2, J2)
    K2D2 <- matrix(rep(tf2, M), nrow = M, ncol = J2, byrow = TRUE)
  }

  # Trap
  if (!is.null(y_mark)) {
    D1 <- pairwise_distances(s1, X1)
    lamd_trap <- lam0_mark * exp(-D1 * D1 / (2 * sigma_d * sigma_d))
    ll_y_mark <- array(0, dim = c(M, J1))
    if (obstype[1] == "bernoulli") {
      pd_trap <- 1 - exp(-lamd_trap)
      pd_trap_cand <- pd_trap
      ll_y_mark <- dbinom(y_mark2D, K2D1, pd_trap * z, log = TRUE)
    } else if (obstype[1] == "poisson") {
      ll_y_mark <- dpois(y_mark2D, K2D1 * lamd_trap * z, log = TRUE)
    }
    lamd_trap_cand <- lamd_trap
    ll_y_mark_cand <- ll_y_mark

    if (!is.finite(sum(ll_y_mark))) {
      stop("Trap obs likelihood not finite. Try raising lam0_mark and/or sigma_d inits")
    }
  } else {
    pd_trap <- array(0, dim = c(M, 1))
    pd_trap_cand <- array(0, dim = c(M, 1))
    lamd_trap_cand <- array(0, dim = c(M, 1))
    ll_y_mark_cand <- array(0, dim = c(M, 1))
    D1 <- array(0, dim = c(M, 1))
  }

  # Sight
  D2 <- pairwise_distances(s2, X2)
  lamd_sight <- lam0_sight * exp(-D2 * D2 / (2 * sigma_d * sigma_d))
  ll_y_sight <- array(0, dim = c(M, J2))
  if (obstype[2] == "bernoulli") {
    pd_sight <- 1 - exp(-lamd_sight)
    pd_sight_cand <- pd_sight
    ll_y_sight <- dbinom(y_sight_true, K2D2, pd_sight * z, log = TRUE)
  } else if (obstype[2] == "poisson") {
    pd_sight <- matrix(0, dim(lamd_sight)[1], dim(lamd_sight)[2])
    pd_sight_cand <- pd_sight
    ll_y_sight <- dpois(y_sight_true, K2D2 * lamd_sight * z, log = TRUE)
  }
  lamd_sight_cand <- lamd_sight
  ll_y_sight_cand <- ll_y_sight

  if (!is.finite(sum(ll_y_sight))) {
    stop("Sighting obs likelihood not finite. Try raising lam0_sight and/or sigma_d inits")
  }

  # movement likelihood.
  if (mobile_center & !is.null(y_mark)){
    ll_s2 <- log(dnorm(s2[, 1], s1[, 1], sigma_p) /
                   (pnorm(state_space$xlim[2], s1[, 1], sigma_p) -
                      pnorm(state_space$xlim[1], s1[, 1], sigma_p)))
    ll_s2 <- ll_s2 + log(dnorm(s2[, 2], s1[, 2], sigma_p) /
                           (pnorm(state_space$ylim[2], s1[, 2], sigma_p) -
                              pnorm(state_space$ylim[1], s1[, 2], sigma_p)))
    ll_s2_cand <- ll_s2
  }

  ###################################
  # Begin MCMC algorithm ----
  ###################################
  for (iter in 1:niter) {
    ###################################
    # Update lam0_mark ----
    ###################################
    if (!is.null(y_mark)) {
      lam0_mark_cand <- rnorm(1, lam0_mark, proppars$lam0_mark)
      if (lam0_mark_cand > 0) {
        lamd_trap_cand <- lam0_mark_cand * exp(-D1 * D1 / (2 * sigma_d * sigma_d))
        pd_trap_cand <-
          if (obstype[1] == "bernoulli") 1 - exp(-lamd_trap_cand) else 0

        # update log likelihood for marked individuals
        ll_y_mark_cand <- calculate_ll_bern_pois(
          obstype[1],
          y_mark2D,
          K2D1,
          lamd_trap_cand,
          z
        )

        if (runif(1) < exp(sum(ll_y_mark_cand) - sum(ll_y_mark))) {
          lam0_mark <- lam0_mark_cand
          lamd_trap <- lamd_trap_cand
          pd_trap <- pd_trap_cand
          ll_y_mark <- ll_y_mark_cand

          # Update acceptance rates
          accept_rates$accept[accept_rates$acc_iter == iter &
                                accept_rates$label == "lam0_mark"] <- 1

        }
      }
    } else {
      lam0_mark <- 0
      lamd_trap <- array(0, dim = c(M, 1))
      ll_y_mark <- array(0, dim = c(M, 1))
    }

    ###################################
    # Update lam0_sight ----
    ###################################
    lam0_sight_cand <- rnorm(1, lam0_sight, proppars$lam0_sight)
    if (lam0_sight_cand > 0) {
      lamd_sight_cand <- lam0_sight_cand * exp(-D2 * D2 / (2 * sigma_d * sigma_d))

      pd_sight_cand <- if (obstype[2] == "bernoulli") {
        1 - exp(-lamd_sight_cand)
      } else {
        matrix(0, dim(lamd_sight_cand)[1], dim(lamd_sight_cand)[2])
      }

      # update log likelihood for marked individuals
      ll_y_sight_cand <- calculate_ll_bern_pois(
        obstype[2],
        y_sight_true,
        K2D2,
        lamd_sight_cand,
        z
      )

      if (runif(1) < exp(sum(ll_y_sight_cand) - sum(ll_y_sight))) {
        lam0_sight <- lam0_sight_cand
        lamd_sight <- lamd_sight_cand
        pd_sight <- pd_sight_cand
        ll_y_sight <- ll_y_sight_cand
        llysightsum <- sum(ll_y_sight_cand)

        # Update acceptance rates
        accept_rates$accept[accept_rates$acc_iter == iter &
                              accept_rates$label == "lam0_sight"] <- 1

      }
    }

    ###################################
    # Update sigma_d
    ###################################
    sigma_d_cand <- rnorm(1, sigma_d, proppars$sigma_d)
    if (sigma_d_cand > 0) {
      if (!is.null(y_mark)) {
        # log likelihood for marked
        lamd_trap_cand <- lam0_mark * exp(-D1 * D1 / (2 * sigma_d_cand * sigma_d_cand))
        ll_y_mark_cand <-
          calculate_ll_bern_pois(obstype[1], y_mark2D, K2D1, lamd_trap_cand, z)
        llytrapcandsum <- sum(ll_y_mark_cand)
      } else {
        lamd_trap_cand <- array(0, dim = c(M, 1))
        ll_y_mark_cand <- array(0, dim = c(M, 1))
        llytrapcandsum <- 0
      }

      # log likelihood for sighted
      lamd_sight_cand <- lam0_sight * exp(-D2 * D2 / (2 * sigma_d_cand * sigma_d_cand))
      ll_y_sight_cand <-
        calculate_ll_bern_pois(obstype[2], y_sight_true, K2D2, lamd_sight_cand, z)
      llysightcandsum <- sum(ll_y_sight_cand)

      if (uselocs) {
        for (i in telguys) {
          ll_tel_cand[i, ] <- dnorm(locs[i, , 1], s2[i, 1], sigma_d_cand, log = TRUE) +
            dnorm(locs[i, , 2], s2[i, 2], sigma_d_cand, log = TRUE)
        }
      } else {
        ll_tel_cand <- ll_tel <- 0
      }

      # if (usePriors) {
      #   prior.curr <- dgamma(sigma, priors$sigma[1], priors$sigma[2], log = TRUE)
      #   prior_cand <- dgamma(sigma_cand, priors$sigma[1], priors$sigma[2], log = TRUE)
      # } else {
      #   prior.curr <- prior_cand <- 0
      # }

      if (runif(1) < exp((llytrapcandsum + llysightcandsum + sum(ll_tel_cand, na.rm = TRUE)) -
                         (sum(ll_y_mark) + sum(ll_y_sight) + sum(ll_tel, na.rm = TRUE)))) { # + prior.curr
        sigma_d <- sigma_d_cand
        lamd_trap <- lamd_trap_cand
        lamd_sight <- lamd_sight_cand
        ll_y_mark <- ll_y_mark_cand
        ll_y_sight <- ll_y_sight_cand
        ll_tel <- ll_tel_cand

        # Update acceptance rates
        accept_rates$accept[accept_rates$acc_iter == iter &
                              accept_rates$label == "sigma_d"] <- 1

        if (obstype[1] == "bernoulli") {
          pd_trap <- pd_trap_cand
        }
        if (obstype[2] == "bernoulli") {
          pd_sight <- pd_sight_cand
        }
      }
    }

    ###################################
    # ID update
    ###################################
    if (IDup == "Gibbs") {
      # Update y_sight_true from full conditional canceling out inconsistent combos with constraints.
      up <- sample(1:n_samp_latent, nswap, replace = FALSE)
      for (l in up) {
        nj <- which(y_sight_latent[l, ] > 0)
        # Can only swap if IDcovs match
        idx2 <- which(G_use[l, ] != 0)
        if (length(idx2) > 1) { # multiple loci observed
          possible <- which(z == 1 & apply(G_true[, idx2], 1, function(x) {
            all(x == G_use[l, idx2])
          }))
        } else if (length(idx2) == 1) { # single loci observed
          possible <- which(z == 1 & G_true[, idx2] == G_use[l, idx2])
        } else { # fully latent G_obs
          possible <- which(z == 1) # Can match anyone
        }
        if (!(useUnk | useMarkednoID)) { # mark status exclusions handled through G_true
          if (any(data$markedS == 0 | data$markedS == 2)) {
            possible <- possible[which(Kconstraints[possible, l] == 0)] # k marked status constraints
          } else {
            possible <- possible[possible > n_marked] # Can't swap to a marked guy
          }
        } else {
          if (Mark_obs[l] == 2) { # This is an unmarked sample
            if (any(data$markedS == 0 | data$markedS == 2)) {
              possible <- possible[which(Kconstraints[possible, l] == 0)] # k marked status constraints
            } else {
              possible <- possible[possible > n_marked] # Can't swap to a marked guy
            }
          }
          if (Mark_obs[l] == 1) { # This is a marked sample
            possible <- possible[possible <= n_marked] # Can't swap to an unmarked guy
            if (any(data$markedS == 0 | data$markedS == 2)) {
              possible <- possible[which(Kconstraints[possible, l] == 1)] # k marked status constraints
            }
          }
        }
        if (length(possible) == 0) next
        njprobs <- lamd_sight[, nj]
        njprobs[setdiff(1:M, possible)] <- 0
        njprobs <- njprobs / sum(njprobs)
        newID <- sample(1:M, 1, prob = njprobs)
        if (ID[l] != newID) {
          swapped <- c(ID[l], newID)
          # update y.true
          y_sight_true[ID[l], ] <- y_sight_true[ID[l], ] - y_sight_latent[l, ]
          y_sight_true[newID, ] <- y_sight_true[newID, ] + y_sight_latent[l, ]
          ID[l] <- newID
          if (obstype[2] == "bernoulli") {
            ll_y_sight[swapped, ] <- dbinom(y_sight_true[swapped, ], K2D2[swapped, ], pd_sight[swapped, ], log = TRUE)
          } else {
            ll_y_sight[swapped, ] <- dpois(y_sight_true[swapped, ], K2D2[swapped, ] * lamd_sight[swapped, ], log = TRUE)
          }
        }
      }
    } else {
      up <- sample(1:n_samp_latent, nswap, replace = FALSE)
      y_sight_cand <- y_sight_true
      for (l in up) {
        nj <- which(y_sight_latent[l, ] > 0)
        # Can only swap if IDcovs match
        idx2 <- which(G_use[l, ] != 0)
        if (length(idx2) > 1) { # multiple loci observed
          possible <- which(z == 1 & apply(G_true[, idx2], 1, function(x) {
            all(x == G_use[l, idx2])
          }))
        } else if (length(idx2) == 1) { # single loci observed
          possible <- which(z == 1 & G_true[, idx2] == G_use[l, idx2])
        } else { # fully latent G_obs
          possible <- which(z == 1) # Can match anyone
        }
        if (!(useUnk | useMarkednoID)) { # mark status exclusions handled through G_true
          if (any(data$markedS == 0 | data$markedS == 2)) {
            possible <- possible[which(Kconstraints[possible, l] == 0)] # k marked status constraints
          } else {
            possible <- possible[possible > n_marked] # Can't swap to a marked guy
          }
        } else {
          if (Mark_obs[l] == 2) { # This is an unmarked sample
            if (any(data$markedS == 0 | data$markedS == 2)) {
              possible <- possible[which(Kconstraints[possible, l] == 0)] # k marked status constraints
            } else {
              possible <- possible[possible > n_marked] # Can't swap to a marked guy
            }
          }
          if (Mark_obs[l] == 1) { # This is a marked sample
            possible <- possible[possible <= n_marked] # Can't swap to an unmarked guy
            if (any(data$markedS == 0 | data$markedS == 2)) {
              possible <- possible[which(Kconstraints[possible, l] == 1)] # k marked status constraints
            }
          }
        }
        if (binconstraints) { # can't have a y[i,j,k]>1
          legal <- rep(TRUE, length(possible))
          for (i in 1:length(possible)) {
            check <- which(ID == possible[i]) # Who else is currently assigned this possible new ID?
            if (length(check) > 0) { # if false, no samples assigned to this guy and legal stays true
              if (any(constraints[l, check] == 0)) { # if any members of the possible cluster are inconsistent with sample, illegal move
                legal[i] <- FALSE
              }
            }
          }
          possible <- possible[legal]
        }
        if (length(possible) == 0) next
        njprobs <- lamd_sight[, nj]
        njprobs[setdiff(1:M, possible)] <- 0
        njprobs <- njprobs / sum(njprobs)
        newID <- ID
        newID[l] <- sample(1:M, 1, prob = njprobs)
        if (ID[l] == newID[l]) next

        swapped <- c(ID[l], newID[l]) # order swap.out then swap.in
        propprob <- njprobs[swapped[2]]
        backprob <- njprobs[swapped[1]]
        # focalprob=1/n_samp_latent
        # focalbackprob=1/length(possible)
        # update y.true
        y_sight_cand[ID[l], ] <- y_sight_true[ID[l], ] - y_sight_latent[l, ]
        y_sight_cand[newID[l], ] <- y_sight_true[newID[l], ] + y_sight_latent[l, ]
        focalprob <- (sum(ID == ID[l]) / n_samp_latent) * (y_sight_true[ID[l], nj] / sum(y_sight_true[ID[l], ]))
        focalbackprob <- (sum(newID == newID[l]) / n_samp_latent) * (y_sight_cand[newID[l], nj] / sum(y_sight_cand[newID[l], ]))
        ## update ll.y
        if (obstype[2] == "poisson") {
          if (any(data$markedS == 0 | data$markedS == 2)) {
            ll_y_sight_cand[swapped, ] <- dpois(y_sight_cand[swapped, ], K2[swapped, ] * lamd_sight[swapped, ], log = TRUE)
          } else {
            ll_y_sight_cand[swapped, ] <- dpois(y_sight_cand[swapped, ], K2D2[swapped, ] * lamd_sight[swapped, ], log = TRUE)
          }
        } else {
          if (any(data$markedS == 0 | data$markedS == 2)) {
            ll_y_sight_cand[swapped, ] <- dbinom(y_sight_cand[swapped, ], K2[swapped, ], pd_sight[swapped, ], log = TRUE)
          } else {
            ll_y_sight_cand[swapped, ] <- dbinom(y_sight_cand[swapped, ], K2D2[swapped, ], pd_sight[swapped, ], log = TRUE)
          }
        }
        if (runif(1) < exp(sum(ll_y_sight_cand[swapped, ]) - sum(ll_y_sight[swapped, ])) *
            (backprob / propprob) * (focalbackprob / focalprob)) {
          y_sight_true[swapped, ] <- y_sight_cand[swapped, ]
          ll_y_sight[swapped, ] <- ll_y_sight_cand[swapped, ]
          ID[l] <- newID[l]
        }
      }
    }

    ###################################
    # update known_vector and G_latent
    ###################################
    known_vector[(n_marked + 1):M] <- 1 * (rowSums(y_sight_true[(n_marked + 1):M, ]) > 0)
    G_true.tmp <- matrix(0, nrow = M, ncol = ncat)
    G_true.tmp[1:n_marked, ] <- 1
    for (i in unique(ID[ID > n_marked])) {
      idx2 <- which(ID == i)
      if (length(idx2) == 1) {
        G_true.tmp[i, ] <- G_use[idx2, ]
      } else {
        if (ncol(G_use) > 1) {
          G_true.tmp[i, ] <- apply(G_use[idx2, ], 2, max) # consensus
        } else {
          G_true.tmp[i, ] <- max(G_use[idx2, ]) # consensus
        }
      }
    }

    # update G_true
    G_latent <- G_true.tmp == 0
    for (j in 1:ncat) {
      swap <- G_latent[, j]
      G_true[swap, j] <- sample(IDcovs[[j]], sum(swap), replace = TRUE, prob = gamma[[j]])
    }

    # update genotype frequencies
    for (j in 1:ncat) {
      x <- rep(NA, nlevels[[j]])
      for (k in 1:nlevels[[j]]) {
        x[k] <- sum(G_true[z == 1, j] == k) # genotype freqs in pop
      }
      gam <- rgamma(rep(1, nlevels[[j]]), 1 + x)
      gamma[[j]] <- gam / sum(gam)
    }

    # probability of not being captured in a trap AT ALL by either method
    if (obstype[1] == "poisson") {
      pd_trap <- 1 - exp(-lamd_trap)
    }
    if (obstype[2] == "poisson") {
      pd_sight <- 1 - exp(-lamd_sight)
    }

    if (!is.null(y_mark)) {
      pbar.trap <- (1 - pd_trap)^K2D1
      prob0.trap <- exp(rowSums(log(pbar.trap)))
    } else {
      prob0.trap <- 1
    }
    pbar.sight <- (1 - pd_sight)^K2D2
    prob0.sight <- exp(rowSums(log(pbar.sight)))
    prob0 <- prob0.trap * prob0.sight


    fc <- prob0 * psi / (prob0 * psi + 1 - psi)
    z[known_vector == 0] <- rbinom(sum(known_vector == 0), 1, fc[known_vector == 0])

    if (!is.null(y_mark)) {
      if (obstype[1] == "bernoulli") {
        ll_y_mark <- dbinom(y_mark2D, K2D1, pd_trap * z, log = TRUE)
      } else {
        ll_y_mark <- dpois(y_mark2D, K2D1 * lamd_trap * z, log = TRUE)
      }
    } else {
      ll_y_mark <- array(0, dim = c(M, 1))
    }

    if (obstype[2] == "bernoulli") {
      ll_y_sight <- dbinom(y_sight_true, K2D2, pd_sight * z, log = TRUE)
    } else {
      ll_y_sight <- dpois(y_sight_true, K2D2 * lamd_sight * z, log = TRUE)
    }
    psi <- rbeta(1, 1 + sum(z), 1 + M - sum(z))

    ###################################
    # Update activity centers
    ###################################
    if (mobile_center) {
      s_2_accept <- 0
      s_1_accept <- 0
      # s1 - activity centers for marked individuals
      for (i in 1:M) {
        Scand <- c(rnorm(1, s1[i, 1], proppars$s1), rnorm(1, s1[i, 2], proppars$s1))

        inbox <- point_in_area(Scand, state_space)
        if (inbox) {
          d1tmp <- sqrt((Scand[1] - X1[, 1])^2 + (Scand[2] - X1[, 2])^2)
          lamd_trap_cand[i, ] <- lam0_mark *
            exp(-d1tmp * d1tmp / (2 * sigma_d * sigma_d))

          # update log likelihood for s2
          ll_s2_cand[i] <-
            calculate_log_likelihood(s2[i, 1], Scand[1], state_space$xlim, sigma_p) +
            calculate_log_likelihood(s2[i, 2], Scand[2], state_space$ylim, sigma_p)

          pd_trap_cand[i, ] <-
            if (obstype[1] == "bernoulli") 1 - exp(-lamd_trap_cand[i, ]) else 0

          # update log likelihood for marked individuals
          ll_y_mark_cand[i, ] <- calculate_ll_bern_pois(
            obstype[1],
            y_mark2D[i, ],
            K2D1[i, ],
            lamd_trap_cand[i, ],
            z[i]
            )

          if (runif(1) < exp((sum(ll_y_mark_cand[i, ]) + ll_s2_cand[i]) -
                             (sum(ll_y_mark[i, ]) + ll_s2[i]))) {
            s1[i, ] <- Scand
            D1[i, ] <- d1tmp
            lamd_trap[i, ] <- lamd_trap_cand[i, ]
            pd_trap[i, ] <- pd_trap_cand[i, ]
            ll_y_mark[i, ] <- ll_y_mark_cand[i, ]
            ll_s2[i] <- ll_s2_cand[i]
            s_1_accept <- s_1_accept + 1
          }
        }
      }
    # s2 - activity centers for sighted individuals
      for (i in 1:M) {
        if (i %in% telguys) {
          Scand <- c(rnorm(1, s2[i, 1], proppars$s2t),
                     rnorm(1, s2[i, 2], proppars$s2t))
        } else {
          Scand <- c(rnorm(1, s2[i, 1], proppars$s2),
                     rnorm(1, s2[i, 2], proppars$s2))
        }

        inbox <- point_in_area(Scand, state_space)
        if (inbox) {
          d2tmp <- sqrt((Scand[1] - X2[, 1])^2 + (Scand[2] - X2[, 2])^2)
          lamd_sight_cand[i, ] <- lam0_sight *
            exp(-d2tmp * d2tmp / (2 * sigma_d * sigma_d))

          # movement likelihood
          ll_s2_cand[i] <-
            calculate_log_likelihood(Scand[1], s1[i, 1], state_space$xlim, sigma_p) +
            calculate_log_likelihood(Scand[2], s1[i, 2], state_space$ylim, sigma_p)

          ll_y_sight_cand[i, ] <- calculate_ll_bern_pois(
            obstype[2],
            y_sight_true[i, ],
            K2D2[i, ],
            lamd_sight_cand[i, ],
            z[i]
          )

          if (uselocs & (i %in% telguys)) {
            ll_tel_cand[i, ] <-
              dnorm(locs[i, , 1], Scand[1], sigma_d, log = TRUE) +
              dnorm(locs[i, , 2], Scand[2], sigma_d, log = TRUE)

            sum_ll_tel_cand <- sum(ll_tel_cand[i, ], na.rm = TRUE)
            sum_ll_tel <- sum(ll_tel[i, ], na.rm = TRUE)
          } else {
            sum_ll_tel_cand <- 0
            sum_ll_tel <- 0
          }

          if (runif(1) < exp((sum(ll_y_sight_cand[i, ]) +
                              sum_ll_tel_cand +
                              ll_s2_cand[i]) -
                             (sum(ll_y_sight[i, ]) +
                              sum_ll_tel_cand +
                              ll_s2[i]))
          ) {
            s2[i, ] <- Scand
            D2[i, ] <- d2tmp
            lamd_sight[i, ] <- lamd_sight_cand[i, ]
            pd_sight[i, ] <- pd_sight_cand[i, ]
            ll_y_sight[i, ] <- ll_y_sight_cand[i, ]
            if (i %in% telguys) ll_tel[i, ] <- ll_tel_cand[i, ]
            ll_s2[i] <- ll_s2_cand[i]
            s_2_accept <- s_2_accept + 1
          }
        }
      }
    } else {
      # Activity centers do not move between marking and sighting processes
      s_1_accept <- 0
      s_2_accept <- 0
      for (i in 1:M) {
        if (i %in% telguys) {
          Scand <- c(rnorm(1, s2[i, 1], proppars$s2t),
                     rnorm(1, s2[i, 2], proppars$s2t))
        } else {
          Scand <- c(rnorm(1, s2[i, 1], proppars$s2),
                     rnorm(1, s2[i, 2], proppars$s2))
        }

        inbox <- point_in_area(Scand, state_space)
        if (inbox) {

          if (!is.null(y_mark)) {
            # Marked individuals
            d1tmp <- sqrt((Scand[1] - X1[, 1])^2 + (Scand[2] - X1[, 2])^2)
            lamd_trap_cand[i, ] <- lam0_mark *
              exp(-d1tmp * d1tmp / (2 * sigma_d * sigma_d))

            ll_y_mark_cand[i, ] <- calculate_ll_bern_pois(
              obstype[1],
              y_mark2D[i, ],
              K2D1[i, ],
              lamd_trap_cand[i, ],
              z[i]
            )
          } else {
            d1tmp <- 0
          }

          # Sighted individuals
          d2tmp <- sqrt((Scand[1] - X2[, 1])^2 + (Scand[2] - X2[, 2])^2)
          lamd_sight_cand[i, ] <- lam0_sight *
            exp(-d2tmp * d2tmp / (2 * sigma_d * sigma_d))

          ll_y_sight_cand[i, ] <- calculate_ll_bern_pois(
            obstype[2],
            y_sight_true[i, ],
            K2D2[i, ],
            lamd_sight_cand[i, ],
            z[i]
          )

          if (uselocs & (i %in% telguys)) {
            ll_tel_cand[i, ] <-
              dnorm(locs[i, , 1], Scand[1], sigma_d, log = TRUE) +
              dnorm(locs[i, , 2], Scand[2], sigma_d, log = TRUE)

            sum_ll_tel_cand <- sum(ll_tel_cand[i, ], na.rm = TRUE)
            sum_ll_tel <- sum(ll_tel[i, ], na.rm = TRUE)
          } else {
            sum_ll_tel_cand <- 0
            sum_ll_tel <- 0
          }

          if (runif(1) < exp((sum(ll_y_mark_cand[i, ]) +
                              sum(ll_y_sight_cand[i, ]) +
                              sum_ll_tel_cand) -
                             (sum(ll_y_mark[i, ]) +
                              sum(ll_y_sight[i, ]) +
                              sum_ll_tel))
          ) {
            s1[i, ] <- Scand # not actually fit, but needed for accept rates
            s2[i, ] <- Scand
            D1[i, ] <- d1tmp
            D2[i, ] <- d2tmp
            lamd_trap[i, ] <- lamd_trap_cand[i, ]
            lamd_sight[i, ] <- lamd_sight_cand[i, ]
            ll_y_mark[i, ] <- ll_y_mark_cand[i, ]
            ll_y_sight[i, ] <- ll_y_sight_cand[i, ]
            if (i %in% telguys) ll_tel[i, ] <- ll_tel_cand[i, ]
            s_1_accept <- s_1_accept + 1
            s_2_accept <- s_2_accept + 1
            if (obstype[1] == "bernoulli") {
              pd_trap[i, ] <- pd_trap_cand[i, ]
            }
            if (obstype[2] == "bernoulli") {
              pd_sight[i, ] <- pd_sight_cand[i, ]
            }
          }
        }
      }
    }

    # Update acceptance rates
    accept_rates$accept[accept_rates$acc_iter == iter &
                          accept_rates$label == "s1"] <- sum(s_1_accept) / M
    accept_rates$accept[accept_rates$acc_iter == iter &
                          accept_rates$label == "s2"] <- sum(s_2_accept) / M


    # update sigma_p
    if (mobile_center & !is.null(y_mark)){
      sigma_p_cand <- rnorm(1, sigma_p, proppars$sigma_p)

      if (sigma_p_cand > 0) {
        ll_s2_cand <- log(dnorm(s2[, 1], s1[, 1], sigma_p_cand) /
                            (pnorm(state_space$xlim[2], s1[, 1], sigma_p_cand) -
                               pnorm(state_space$xlim[1], s1[, 1], sigma_p_cand)))
        ll_s2_cand <- ll_s2_cand +
          log(dnorm(s2[, 2], s1[, 2], sigma_p_cand) /
                (pnorm(state_space$ylim[2], s1[, 2], sigma_p_cand) -
                   pnorm(state_space$ylim[1], s1[, 2], sigma_p_cand)))


        if (runif(1) < exp(sum(ll_s2_cand) - sum(ll_s2))) {
          sigma_p <- sigma_p_cand
          ll_s2 <- ll_s2_cand

          # Update acceptance rates
          accept_rates$accept[accept_rates$acc_iter == iter &
                                accept_rates$label == "sigma_p"] <- 1
        }
      }
    }

    # Do we record output on this iteration?
    if (iter > nburn & iter %% nthin == 0) {
      if (storeLatent) {
        s1xout[idx, ] <- s1[, 1]
        s1yout[idx, ] <- s1[, 2]
        s2xout[idx, ] <- s2[, 1]
        s2yout[idx, ] <- s2[, 2]
        zout[idx, ] <- z
        IDout[idx, ] <- ID
      }
      if (storeGamma) {
        for (k in 1:ncat) {
          gammaOut[[k]][idx, ] <- gamma[[k]]
        }
      }
      if (useUnk | useMarkednoID) {
        n <- length(unique(ID[ID > n_marked]))
      } else {
        n <- length(unique(ID))
      }

      out[idx, ] <- c(lam0_mark, lam0_sight, sigma_d, sigma_p, sum(z), n, psi)
      idx <- idx + 1
    }

    # Update tuning parms
    if(iter%%tune_check == 0){
      # Update batch counter and delta
      batch_n <- batch_n + 1
      delta_n <- batch_n^-1

      # Calculate mean acceptance rates
      mean_accept <- accept_rates %>%
        dplyr::group_by(.data$label) %>%
        dplyr::filter(.data$acc_iter > (iter - tune_check) & .data$acc_iter <= iter) %>%
        dplyr::summarise(
          val = mean(.data$accept)
        ) %>%
        dplyr::mutate(label = factor(.data$label, names(proppars))) %>%
        dplyr::arrange(.data$label)

      # Adjust tuning parms
      proppars[mean_accept$val > 0.55] <-
        unlist(proppars[mean_accept$val > 0.55]) /
        (1-abs(mean_accept$val[mean_accept$val > 0.55] - 0.55))
      proppars[mean_accept$val < 0.35] <-
        unlist(proppars[mean_accept$val < 0.35]) *
        (1-abs(mean_accept$val[mean_accept$val < 0.35] - 0.35))

      proppars[which(unlist(proppars) < 0)] <- 1e-5

      # Bound tuning parms
      proppars[proppars > unlist(max_proppars)] <- max_proppars[proppars > unlist(max_proppars)]
      proppars[proppars < unlist(min_proppars)] <- min_proppars[proppars < unlist(min_proppars)]
    }
  } # end of MCMC algorithm
  # CheckID
  if (any(data$markedS == 0 | data$markedS == 2)) { # capture order constraints
    if (storeLatent) {
      if (!(useUnk | useMarkednoID)) {
        Mark_obs <- rep(2, n_samp_latent)
      }
      for (l in 1:length(ID)) {
        cons <- Kconstraints[IDout[, l], l]
        if (any(cons == 1) & Mark_obs[l] == 2) {
          stop("Unmarked samples incorrectly assigned in MCMC!")
        }
        if (any(cons == 0) & Mark_obs[l] == 1) {
          stop("Marked samples incorrectly assigned in MCMC!")
        }
      }
    }
  }

  # Check acceptance rates
  if (any(mean_accept$val > 0.6) || any(mean_accept$val < 0.3)) {
    warning("Acceptance rates out of bounds. Check tuning parameters or increase MCMC iterations")
  }

  if (storeLatent & storeGamma) {
    list(out = out, s1xout = s1xout, s1yout = s1yout, s2xout = s2xout, s2yout = s2yout, zout = zout, IDout = IDout, gammaOut = gammaOut, proppars = proppars, accept = mean_accept)
  } else if (storeLatent & !storeGamma) {
    list(out = out, s1xout = s1xout, s1yout = s1yout, s2xout = s2xout, s2yout = s2yout, zout = zout, IDout = IDout)
  } else if (!storeLatent & storeGamma) {
    list(out = out, gammaOut = gammaOut)
  } else {
    list(out = out)
  }
}

################################################################################
# Helper Functions

#' Process and Initialize Spatial Points for Capture-Recapture Analysis
#'
#' @description
#' Initializes and processes spatial locations for individuals in a capture-recapture
#' study, handling both marked and unmarked individuals. The function generates
#' random starting locations, processes observed capture locations, and ensures all
#' points fall within the defined state space.
#'
#' @param M Integer. Total number of individuals to process.
#' @param state_space List. Defines the study area with components:
#' \itemize{
#'   \item xlim: Vector of length 2 with min, max x-coordinates
#'   \item ylim: Vector of length 2 with min, max y-coordinates
#'   \item useverts: Logical. If TRUE, uses vertex-defined state space
#'   \item vertices: Matrix of vertex coordinates (required if useverts=TRUE)
#' }
#' @param y_sight_true Matrix. Sighting history matrix where rows represent
#'                     individuals and columns represent occasions/locations.
#' @param Xall Matrix. Coordinates of all sampling locations.
#' @param y_mark Matrix. Optional marking history matrix. If provided, will be
#'               combined with sighting data.
#'
#' @return Matrix. An M × 2 matrix of processed spatial coordinates where each row
#'         contains x, y coordinates for an individual.
#'
#' @details
#' The function processes spatial points in several steps:
#'
#' 1. Initial Point Generation:
#'    - Randomly generates M points within xlim/ylim bounds
#'    - Uses uniform distribution for initial placement
#'
#' 2. Capture Location Processing:
#'    - Combines marking and sighting histories if both available
#'    - For individuals with captures/sightings:
#'      * Single capture: Uses exact capture location
#'      * Multiple captures: Uses mean location of captures
#'
#' 3. State Space Validation:
#'    - If using vertices (useverts=TRUE):
#'      * Checks if points fall within polygon
#'      * Regenerates points that fall outside until valid
#'      * Uses point_in_area function for validation
process_spatial_points <- function(M, state_space, y_sight_true, Xall,
                                   y_mark = NULL) {
  # Initialize random points within state space
  s1 <- cbind(runif(M, state_space$xlim[1], state_space$xlim[2]),
              runif(M, state_space$ylim[1], state_space$ylim[2]))


  # Combine mark and sight data if available
  y_all2D <- if (!is.null(y_mark)) {
    cbind(apply(y_mark, c(1, 2), sum), y_sight_true)
  } else {
    y_sight_true
  }

  # Process points with positive captures/sightings
  idx <- which(rowSums(y_all2D) > 0)
  for (i in idx) {
    trps <- matrix(Xall[y_all2D[i, ] > 0, 1:2], ncol = 2, byrow = FALSE)
    if (nrow(trps) > 1) {
      s1[i, ] <- c(mean(trps[, 1]), mean(trps[, 2]))
    } else {
      s1[i, ] <- trps
    }
  }

  # Check and adjust points against state space vertices if needed
  if (state_space$useverts == TRUE) {
    inside <- rep(NA, nrow(s1))

    # Check all points
    for (i in 1:nrow(s1)) {
      inside[i] <- point_in_area(s1[i, ], state_space)
    }

    # Adjust points that fall outside the area
    idx <- which(inside == FALSE)
    if (length(idx) > 0) {
      for (i in 1:length(idx)) {
        while (inside[idx[i]] == FALSE) {
          s1[idx[i], ] <- c(runif(1, state_space$xlim[1], state_space$xlim[2]),
                            runif(1, state_space$ylim[1], state_space$ylim[2]))
          inside[idx[i]] <- point_in_area(s1[idx[i], ], state_space)
        }
      }
    }
  }

  return(s1)
}

#' Process Genetic Data with Marking Status
#'
#' @description
#' Processes genetic data while accounting for marking status and ID information. This function
#' handles genetic markers for both marked and unmarked individuals, consolidates multiple
#' observations per individual, and can impute missing genetic data based on provided probability
#' distributions.
#'
#' @param M Integer. Total number of individuals in the population to be processed.
#' @param ncat Integer. Number of genetic categories or markers being analyzed.
#' @param n_marked Integer. Number of marked individuals in the dataset.
#' @param G_marked Matrix. Genetic data for marked individuals, with dimensions n_marked × ncat.
#' @param ID Vector. Individual identifiers corresponding to the genetic observations.
#' @param G_use Matrix. Observed genetic data to be processed, with rows corresponding to observations.
#' @param IDcovs List. List of possible values for each genetic category.
#' @param gamma List. List of probability vectors for sampling missing genetic values.
#' @param useUnk Logical. If TRUE, handles unknown individuals differently. Default is FALSE.
#' @param useMarkednoID Logical. If TRUE, processes marked individuals without IDs. Default is FALSE.
#'
#' @return A list containing:
#' \itemize{
#'   \item G_true: Matrix of processed genetic data for all individuals
#'   \item G_use: Matrix of observed genetic data (potentially modified)
#'   \item Mark_obs: Vector of marking status (if useUnk or useMarkednoID is TRUE)
#'   \item G_latent: Logical matrix indicating which genetic values were initially missing
#' }
#'
#' @details
#' The function performs several key operations:
#' 1. Initializes a matrix for true genetic values using marked individual data
#' 2. Consolidates multiple observations per individual by taking the maximum value
#' 3. Handles special cases for unknown or unmarked individuals
#' 4. Imputes missing genetic data using provided probability distributions
#' 5. Ensures output consistency regardless of input dimensions
#'
#' When useUnk or useMarkednoID is TRUE, the function treats the last column of genetic
#' data as marking status and processes it separately.
#'
#' @examples
#' \dontrun{
#' # Basic usage with marked individuals only
#' M <- 10
#' ncat <- 3
#' n_marked <- 5
#' G_marked <- matrix(c(1,2,1, 2,1,2, 1,1,1, 2,2,2, 1,2,2), nrow=5, ncol=3, byrow=TRUE)
#' ID <- c(1,1,2,3,4,5)
#' G_use <- matrix(c(1,2,1, 1,2,1, 2,1,2, 1,1,1, 2,2,2, 1,2,2), nrow=6, ncol=3, byrow=TRUE)
#' IDcovs <- list(c(1,2), c(1,2), c(1,2))
#' gamma <- list(c(0.5,0.5), c(0.5,0.5), c(0.5,0.5))
#'
#' result <- processGeneticData(M, ncat, n_marked, G_marked, ID, G_use, IDcovs, gamma)
#'}
processGeneticData <- function(M, ncat, n_marked, G_marked, ID, G_use, IDcovs,
                               gamma, useUnk = FALSE, useMarkednoID = FALSE) {
  # Initialize G_true matrix
  G_true <- matrix(0, nrow = M, ncol = ncat)
  G_true[1:n_marked, ] <- G_marked

  # Process each unique ID
  for (i in unique(ID)) {
    idx <- which(ID == i)
    if (length(idx) == 1) {
      G_true[i, ] <- G_use[idx, ]
    } else {
      if (ncol(G_use) > 1) {
        G_true[i, ] <- apply(G_use[idx, ], 2, max)
      } else {
        G_true[i, ] <- max(G_use[idx, ])
      }
    }
  }

  # Handle unknown or unmarked ID cases
  if (useUnk | useMarkednoID) {
    if (max(ID) < M) {
      G_true[(max(ID) + 1):M, ncol(G_true)] <- 2
    }
    unkguys <- which(G_use[, ncol(G_use)] == 0)
  }

  # Determine which genotypes can be updated
  G_latent <- G_true == 0

  # Update genotypes based on conditions
  if (!(useUnk | useMarkednoID)) {
    for (j in 1:ncat) {
      fix <- G_true[, j] == 0
      G_true[fix, j] <- sample(IDcovs[[j]], sum(fix), replace = TRUE, prob = gamma[[j]])
    }
  } else {
    for (j in 1:(ncat - 1)) {
      fix <- G_true[, j] == 0
      G_true[fix, j] <- sample(IDcovs[[j]], sum(fix), replace = TRUE, prob = gamma[[j]])
    }
    # Split marked status back off
    Mark_obs <- G_use[, ncat]
    ncat <- ncat - 1
    G_use <- G_use[, 1:ncat]
    G_true <- G_true[, 1:ncat]
  }

  # Ensure outputs are matrices
  if (!is.matrix(G_use)) {
    G_use <- matrix(G_use, ncol = 1)
  }
  if (!is.matrix(G_true)) {
    G_true <- matrix(G_true, ncol = 1)
  }

  # Return results as a list
  return(list(
    G_true = G_true,
    G_use = G_use,
    Mark_obs = if(exists("Mark_obs")) Mark_obs else NULL,
    G_latent = G_latent
  ))
}

#' Initialize Unmarked Data Usage Flag
#'
#' @description
#' Initializes and validates the usage of unmarked individual data in capture-recapture studies.
#' The function checks for the presence of unmarked individual data and ensures proper
#' dimensionality of sighting data.
#'
#' @param data List. A data structure containing capture-recapture data with the following possible components:
#' \itemize{
#'   \item G_unmarked: Matrix of genetic data for unmarked individuals
#'   \item y_sight_unmarked: Array of sighting data for unmarked individuals
#' }
#' @param ncat Integer. Number of genetic categories/markers
#'
#' @return Logical. TRUE if unmarked individual data should be used (G_unmarked exists in data),
#' FALSE otherwise.
#'
#' @details
#' The function performs two main operations:
#' 1. Checks if G_unmarked exists in the input data
#' 2. If G_unmarked exists, validates that y_sight_unmarked has 3 dimensions
#'
#' If G_unmarked is not present, the function creates an empty matrix with the
#' appropriate number of columns (ncat) and returns FALSE.
#'
#' @examples
#' \dontrun{
#' # Example with unmarked data
#' data <- list(
#'   G_unmarked = matrix(c(1,2,1, 2,1,2), nrow=2, ncol=3, byrow=TRUE),
#'   y_sight_unmarked = array(0, dim=c(2,3,4))  # 2 individuals, 3 occasions, 4 sites
#' )
#' useUM <- init_useUM(data, ncat = 1)  # Returns TRUE
#'
#' # Example without unmarked data
#' data <- list(
#'   y_sight_marked = matrix(c(0,1,0, 1,0,1), nrow=2, ncol=3)
#' )
#' useUM <- init_useUM(data, ncat = 1)  # Returns FALSE
#'}
init_useUM <- function(data, ncat) {
  if ("G_unmarked" %in% names(data)) {
    if (length(dim(data$y_sight_unmarked)) != 3) {
      stop("dim(y_sight_unmarked) must be 3. Reduced to 2 during initialization")
    }
    useUM <- TRUE
  } else {
    data$G_unmarked <- matrix(0, nrow = 0, ncol = ncat)
    useUM <- FALSE
  }

  return(useUM)
}

#' Initialize Capture Histories for Spatial Capture-Recapture
#'
#' @description
#' Initializes and processes spatial capture-recapture histories by assigning
#' observations to individuals, handling both marked and unmarked individuals,
#' and managing spatial constraints. The function processes observations in two
#' passes: first for unmarked individuals, then for marked individuals with
#' unknown IDs.
#'
#' @param y_sight_latent Array. Latent capture histories from sighting data.
#' @param y_sight_true Array. True capture histories to be initialized.
#' @param y_sight_marked Array. Capture histories for marked individuals.
#' @param y_mark Array. Optional. Marking capture histories.
#' @param status Vector. Status indicator for each capture (1 for marked, 0 for unmarked).
#' @param ID Vector. Individual identifiers to be initialized.
#' @param M Integer. Maximum number of individuals in the population.
#' @param n_marked Integer. Number of marked individuals.
#' @param n_samp_latent Integer. Number of latent samples.
#' @param X1 Matrix. Coordinates of marking locations.
#' @param X2 Matrix. Coordinates of sighting locations.
#' @param G_marked Matrix. Genetic data for marked individuals.
#' @param G_marked_noID Matrix. Genetic data for marked individuals without IDs.
#' @param constraints Matrix. Binary matrix indicating compatible sample pairs.
#' @param Kconstraints Matrix. Additional constraints for sample assignment.
#' @param useMarkednoID Logical. Whether to process marked individuals without IDs. Default is TRUE.
#' @param markedS Vector. Status indicators for marked samples.
#'
#' @return A list containing:
#' \itemize{
#'   \item y_sight_true: Array of initialized true capture histories
#'   \item ID: Vector of assigned individual identifiers
#' }
#'
#' @details
#' The function performs initialization in multiple steps:
#' 1. Validates input arrays
#' 2. Processes unmarked and unknown samples:
#'    - Identifies detection locations
#'    - Finds candidates caught at same locations
#'    - Applies spatial and genetic constraints
#' 3. Processes marked individuals with unknown IDs:
#'    - Calculates mean locations for marked individuals
#'    - Checks genetic compatibility
#'    - Assigns based on spatial proximity
#' 4. Validates final assignments against constraints
#'
initialize_capture_histories <- function(
    y_sight_latent, y_sight_true, y_sight_marked, y_mark = NULL, status, ID, M,
    n_marked, n_samp_latent, X1, X2, G_marked, G_marked_noID, constraints,
    Kconstraints, useMarkednoID = TRUE, markedS) {

  # Input validation
  if (!is.array(y_sight_latent) || !is.array(y_sight_true)) {
    stop("y_sight_latent and y_sight_true must be arrays")
  }

  # Initialize index for new individuals
  y_sight_true[1:n_marked, , ] <- y_sight_marked
  ID <- rep(NA, n_samp_latent)
  idx <- n_marked + 1

  # First pass: assign unmarked and unknown samples
  for (i in 1:n_samp_latent) {
    if (useMarkednoID && status[i] == 1) next

    if (idx > M) {
      stop("Need to raise M to initialize y.true")
    }

    # Find traps where individual was detected
    traps <- which(rowSums(y_sight_latent[i, , ]) > 0)
    y_sight_true2D <- apply(y_sight_true, c(1, 2), sum)

    # Find candidates caught at same traps
    if (length(traps) == 1) {
      cand <- which(y_sight_true2D[, traps] > 0)
    } else {
      cand <- which(rowSums(y_sight_true2D[, traps]) > 0)
    }
    cand <- cand[cand > n_marked]

    if (length(cand) > 0) {
      cand <- cand[1] # if multiple candidates, choose first

      # Check constraints
      cands <- which(ID %in% cand)
      if (all(constraints[i, cands] == 1)) {
        y_sight_true[cand, , ] <- y_sight_true[cand, , ] + y_sight_latent[i, , ]
        ID[i] <- cand
      } else {
        y_sight_true[idx, , ] <- y_sight_latent[i, , ]
        ID[i] <- idx
        idx <- idx + 1
      }
    } else {
      y_sight_true[idx, , ] <- y_sight_latent[i, , ]
      ID[i] <- idx
      idx <- idx + 1
    }
  }

  # Second pass: assign marked unknown ID individuals
  if (useMarkednoID) {
    fix <- which(status == 1)
    meanloc <- matrix(NA, nrow = n_marked, ncol = 2)

    # Calculate mean locations for marked individuals
    for (i in 1:n_marked) {
      locs2 <- matrix(0, nrow = 0, ncol = 2)
      if (!is.null(y_mark)) {
        trap1 <- which(rowSums(y_mark[i, , ]) > 0)
        if (length(trap1) > 0) {
          locs2 <- rbind(locs2, X1[trap1, ])
        }
      }
      trap2 <- which(rowSums(y_sight_marked[i, , ]) > 0)
      if (length(trap2) > 0) {
        locs2 <- rbind(locs2, X2[trap2, ])
      }
      if (nrow(locs2) > 1) {
        meanloc[i, ] <- colMeans(locs2)
      } else if (nrow(locs2) > 0) {
        meanloc[i, ] <- locs2
      }
    }

    # Assign marked unknown ID samples
    for (i in 1:nrow(G_marked_noID)) {
      trap <- which(rowSums(y_sight_latent[i, , ]) > 0)
      compatible <- rep(FALSE, n_marked)

      # Check genetic compatibility
      for (j in 1:n_marked) {
        nonzero1 <- G_marked[j, 1:(ncol(G_marked) - 1)] != 0
        nonzero2 <- G_marked_noID[i, ] != 0
        nonzero <- which(nonzero1 & nonzero2)
        if (all(G_marked[j, nonzero] == G_marked_noID[i, nonzero])) {
          compatible[j] <- TRUE
        }
      }

      if (all(!compatible)) {
        stop(sprintf("No G_marked compatible with G_marked_noID %d", i))
      }

      # Calculate distances and find closest compatible individual
      dists <- sqrt((X2[trap, 1] - meanloc[, 1])^2 + (X2[trap, 2] - meanloc[, 2])^2)
      dists[!compatible] <- Inf
      dists[which(Kconstraints[1:n_marked, fix[i]] == 0)] <- Inf

      if (all(!is.finite(dists))) {
        stop(sprintf("No G_marked compatible with G_marked_noID %d", i))
      }

      ID[fix[i]] <- which.min(dists)
      y_sight_true[ID[fix[i]], , ] <- y_sight_true[ID[fix[i]], , ] +
        y_sight_latent[fix[i], , ]
    }
  }

  # Check assignment consistency with constraints
  checkID <- unique(ID)
  checkID <- checkID[checkID > n_marked]
  for (i in 1:length(checkID)) {
    idx <- which(ID == checkID[i])
    if (!all(constraints[idx, idx] == 1)) {
      stop("ID initialized improperly")
    }
  }

  if (!is.null(y_mark)) {
    if (any(markedS == 0 | markedS == 2)) {
      check <- which(ID <= n_marked)
      for (i in check) {
        if (Kconstraints[ID[i], i] == 0) {
          stop("ID initialized improperly for marked no ID samples")
        }
      }
    }
  }

  return(list(
    y_sight_true = y_sight_true,
    ID = ID
  ))
}

#' Generate Constraint Matrices for Sample Data
#'
#' @description
#' Generates constraint matrices for sample data in capture-recapture studies by evaluating
#' genetic compatibility between samples and handling special cases for Bernoulli-type
#' observations. The function creates constraints that prevent incompatible samples from
#' being assigned to the same individual.
#'
#' @param y_sight_latent Array. Latent capture histories containing detection data.
#' @param obstype Character. Observation type, with special handling for "bernoulli".
#' @param G_use Matrix. Genetic data matrix where rows represent samples and columns
#'              represent genetic markers.
#'
#' @return A list containing:
#' \itemize{
#'   \item constraints: Matrix. Binary constraint matrix where 0 indicates incompatible
#'         sample pairs and 1 indicates compatible pairs.
#'   \item binconstraints: Logical. TRUE if Bernoulli-specific constraints were added,
#'         FALSE otherwise.
#' }
#'
#' @details
#' The function generates constraints in two phases:
#'
#' 1. Genetic Compatibility:
#'    - Creates an n × n matrix where n is the number of samples
#'    - Compares genetic markers between all sample pairs
#'    - Sets constraint to 0 if samples have incompatible genetic markers
#'
#' 2. Bernoulli-specific Constraints (if applicable):
#'    - Identifies samples from same trap and occasion
#'    - Prevents assignment to same individual when detected simultaneously
#'    - Ensures constraint matrix symmetry
#'
#' @examples
#' \dontrun{
#' # Create example data
#' y_sight_latent <- array(0, dim=c(3,2,2))  # 3 samples, 2 occasions, 2 traps
#' y_sight_latent[1,1,1] <- 1
#' y_sight_latent[2,1,2] <- 1
#' y_sight_latent[3,1,1] <- 1
#'
#' # Genetic data with 2 markers
#' G_use <- matrix(c(
#'   1,2,  # Sample 1
#'   1,0,  # Sample 2
#'   2,2   # Sample 3
#' ), ncol=2, byrow=TRUE)
#'
#' # Generate constraints
#' result <- generate_constraints(
#'   y_sight_latent = y_sight_latent,
#'   obstype = "bernoulli",
#'   G_use = G_use
#' )
#'}
generate_constraints <- function(y_sight_latent, obstype, G_use) {

  n_samp_latent <- nrow(y_sight_latent)

  # make constraints for data initialization
  constraints <- matrix(1, nrow = n_samp_latent, ncol = n_samp_latent)
  for (i in 1:n_samp_latent) {
    for (j in 1:n_samp_latent) {
      n1 <- which(G_use[i, ] != 0)
      n2 <- which(G_use[j, ] != 0)
      comp <- n1[which(n1 %in% n2)]
      if (any(G_use[i, comp] != G_use[j, comp])) {
        constraints[i, j] <- 0
      }
    }
  }

  # Initialize return values
  binconstraints <- FALSE
  n_samp_latent <- nrow(y_sight_latent)

  # Only proceed if observation type is bernoulli
  if (obstype == "bernoulli") {
    # Get indices of non-zero elements
    idx <- t(apply(y_sight_latent, 1, function(x) {
      which(x > 0, arr.ind = TRUE)
    }))

    # Add constraints
    for (i in 1:n_samp_latent) {
      for (j in 1:n_samp_latent) {
        if (i != j) {
          # Check if samples are from same trap and occasion
          if (all(idx[i, 1:2] == idx[j, 1:2])) {
            constraints[i, j] <- 0 # prevent combination
            constraints[j, i] <- 0 # ensure symmetry
            binconstraints <- TRUE
          }
        }
      }
    }
  }

  return(list(
    constraints = constraints,
    binconstraints = binconstraints
  ))

}

#' Generate K-Constraints Matrix for Capture-Recapture Data
#'
#' @description
#' Generates a constraint matrix that defines valid associations between marked individuals
#' and capture events based on marking status and sighting occasions. This function
#' is particularly useful in spatial capture-recapture studies where individuals may
#' have different marking statuses across sampling occasions.
#'
#' @param data List. Contains study data, optionally including 'markedS' matrix for
#'             marking status information.
#' @param y_sight_latent Array. 3D array of latent capture histories with dimensions
#'                       n_samp_latent × n_occasions × n_traps.
#' @param n_marked Integer. Number of marked individuals in the study.
#' @param M Integer. Total number of individuals (marked + unmarked) in the population.
#' @param n_samp_latent Integer. Number of latent samples.
#' @param K2 Integer. Number of secondary sampling occasions.
#'
#' @return Matrix. An M × n_samp_latent constraint matrix where:
#' \itemize{
#'   \item 0: Individual-sample combination not allowed
#'   \item 1: Individual-sample combination allowed (standard marking)
#'   \item 2: Individual-sample combination allowed (alternative marking status)
#' }
#'
#' @details
#' The function performs the following operations:
#'
#' 1. Input Validation:
#'    - Checks data structure and types
#'    - Validates dimensions and positive values
#'    - Ensures compatibility of markedS matrix if provided
#'
#' 2. Constraint Generation:
#'    - Initializes M × n_samp_latent constraint matrix
#'    - Processes each marked individual against all samples
#'    - Uses first sighting occasion to determine constraints
#'    - Incorporates marking status from markedS matrix
#'
#' If markedS is not provided in the data list, all marked individuals are assumed
#' to have standard marking status (1) for all occasions.
#'
generate_Kconstraints <- function(data, y_sight_latent, n_marked, M, n_samp_latent, K2) {
  # Input validation
  if (!is.list(data)) {
    stop("'data' must be a list")
  }
  if (!is.array(y_sight_latent) || length(dim(y_sight_latent)) != 3) {
    stop("'y_sight_latent' must be a 3D array")
  }
  if (!all(sapply(list(n_marked, M, n_samp_latent, K2), is.numeric))) {
    stop("n_marked, M, n_samp_latent, and K2 must be numeric")
  }
  if (!all(sapply(list(n_marked, M, n_samp_latent, K2), function(x) x > 0))) {
    stop("n_marked, M, n_samp_latent, and K2 must be positive")
  }

  # Initialize constraint matrix
  Kconstraints <- matrix(0, nrow = M, ncol = n_samp_latent)

  # Set up marked status matrix
  if ("markedS" %in% names(data)) {
    markedS <- data$markedS
    if (!is.matrix(markedS) || nrow(markedS) != n_marked || ncol(markedS) != K2) {
      stop("markedS must be a matrix with dimensions n_marked x K2")
    }
  } else {
    markedS <- matrix(1, nrow = n_marked, ncol = K2)
  }

  # Create constraints
  for (i in 1:n_marked) {
    for (j in 1:n_samp_latent) {
      # Find occasions where sightings occurred
      occ <- which(apply(y_sight_latent[j, , ], 2, sum) > 0)

      # Skip if no sightings
      if (length(occ) == 0) next

      # Set constraints based on marking status
      if (markedS[i, occ[1]] == 1) {  # Using first sighting occasion
        Kconstraints[i, j] <- 1
      } else if (markedS[i, occ[1]] == 2) {
        Kconstraints[i, j] <- 2
      }
    }
  }

  return(Kconstraints)
}

#' Combine Genetic Data with Different Observation Types
#'
#' @description
#' Combines genetic data from different observation types (marked, unmarked, unknown)
#' in capture-recapture studies. The function handles various combinations of data
#' types based on the specified parameters and validates data consistency across
#' different observation types.
#'
#' @param data List. A data structure containing the following components:
#' \itemize{
#'   \item G_unmarked: Matrix of genetic data for unmarked individuals
#'   \item G_marked: Matrix of genetic data for marked individuals
#'   \item G_unk: Optional matrix of genetic data for unknown status individuals
#'   \item G_marked_noID: Optional matrix of genetic data for marked individuals without IDs
#'   \item y_sight_unmarked: Capture histories for unmarked individuals
#'   \item y_sight_unk: Optional capture histories for unknown status individuals
#'   \item y_sight_marked_noID: Optional capture histories for marked individuals without IDs
#' }
#' @param ncat Integer. Number of genetic categories/markers. Default is 0.
#' @param nlevels Vector. Number of levels for each genetic marker. Default is NULL.
#' @param useUnk Logical. Whether to include unknown status individuals. Default is FALSE.
#' @param useMarkednoID Logical. Whether to include marked individuals without IDs. Default is FALSE.
#'
#' @return A list containing:
#' \itemize{
#'   \item G_use: Combined genetic data matrix
#'   \item G_marked: Modified genetic data for marked individuals
#'   \item y_sight_latent: Combined capture histories
#'   \item ncat: Updated number of genetic categories
#'   \item nlevels: Updated number of levels per category
#'   \item status: Vector indicating status of each individual (NULL if neither useUnk nor useMarkednoID)
#' }
#'
#' @details
#' The function performs the following operations:
#'
#' 1. Input Validation:
#'    - Checks for required data components
#'    - Validates consistency between genetic data and capture histories
#'    - Ensures matching dimensions across matrices
#'
#' 2. Data Combination:
#'    - Combines data based on useUnk and useMarkednoID parameters
#'    - Handles three scenarios:
#'      * Unknown status only (useUnk = TRUE, useMarkednoID = FALSE)
#'      * Marked without ID only (useUnk = FALSE, useMarkednoID = TRUE)
#'      * Both types (useUnk = TRUE, useMarkednoID = TRUE)
#'    - Adds status column to genetic data when appropriate
#'
#' Status codes in the output:
#' - 0: Unknown status
#' - 1: Marked
#' - 2: Unmarked
#'
#' @examples
#' \dontrun{
#' # Create example data
#' data <- list(
#'   G_unmarked = matrix(c(1,2,1, 2,1,2), nrow=2, ncol=3),
#'   G_marked = matrix(c(1,1,1, 2,2,2), nrow=2, ncol=3),
#'   G_unk = matrix(c(1,1,2), nrow=1, ncol=3),
#'   y_sight_unmarked = array(0, dim=c(2,3,4)),  # 2 unmarked, 3 traps, 4 occasions
#'   y_sight_unk = array(0, dim=c(1,3,4))        # 1 unknown, 3 traps, 4 occasions
#' )
#'
#' # Combine data including unknown status individuals
#' result <- combine_genetic_data(
#'   data = data,
#'   ncat = 3,
#'   nlevels = c(2,2,2),
#'   useUnk = TRUE,
#'   useMarkednoID = FALSE
#' )
#'}
combine_genetic_data <- function(
    data,
    ncat = 0,
    nlevels = NULL,
    useUnk = FALSE,
    useMarkednoID = FALSE
) {
  # Input validation
  if (is.null(data$G_unmarked) || is.null(data$G_marked) || is.null(data$y_sight_unmarked)) {
    stop("G_unmarked, G_marked, and y_sight_unmarked are required")
  }

  if (!is.null(data$G_unk) && is.null(data$y_sight_unk)) {
    stop("y_sight_unk must be provided when G_unk is present")
  }

  if (!is.null(data$G_marked_noID) && is.null(data$y_sight_marked_noID)) {
    stop("y_sight_marked_noID must be provided when G_marked_noID is present")
  }

  # Check matrix dimensions
  if (nrow(data$G_unmarked) != nrow(data$y_sight_unmarked)) {
    stop("Dimensions mismatch between G_unmarked and y_sight_unmarked")
  }

  G_marked <- data$G_marked

  if (useUnk & !useMarkednoID) {
    G_use <- rbind(data$G_unmarked, data$G_unk)
    status <- c(rep(2, nrow(data$G_unmarked)), rep(0, nrow(data$G_unk)))
    G_use <- cbind(G_use, status)
    G_marked <- cbind(G_marked, rep(1, nrow(G_marked)))
    ncat <- ncat + 1
    y_sight_latent <- combine_matrices(data$y_sight_unmarked, data$y_sight_unk)
  } else if (!useUnk & useMarkednoID) {
    G_use <- rbind(data$G_unmarked, data$G_marked_noID)
    status <- c(rep(2, nrow(data$G_unmarked)), rep(1, nrow(data$G_marked_noID)))
    G_use <- cbind(G_use, status)
    G_marked <- cbind(G_marked, rep(1, nrow(G_marked)))
    ncat <- ncat + 1
    y_sight_latent <- combine_matrices(data$y_sight_unmarked, data$y_sight_marked_noID)
  } else if (useUnk & useMarkednoID) {
    G_use <- rbind(data$G_unmarked, data$G_unk, data$G_marked_noID)
    status <- c(rep(2, nrow(data$G_unmarked)), rep(0, nrow(data$G_unk)), rep(1, nrow(data$G_marked_noID)))
    G_use <- cbind(G_use, status)
    G_marked <- cbind(G_marked, rep(1, nrow(G_marked)))
    ncat <- ncat + 1
    nlevels <- c(nlevels, 2)
    y_sight_latent <- combine_matrices(data$y_sight_unmarked, data$y_sight_unk, data$y_sight_marked_noID)
  } else {
    G_use <- data$G_unmarked
    y_sight_latent <- data$y_sight_unmarked
    status <- NULL
  }

  # Return results as a list
  return(list(
    G_use = G_use,
    G_marked = G_marked,
    y_sight_latent = y_sight_latent,
    ncat = ncat,
    nlevels = nlevels,
    status = status
  ))
}

#' Combine Matrices Along First Dimension
#'
#' @description
#' Combines multiple arrays or matrices along their first dimension (typically rows),
#' while preserving all other dimensions. This function handles arrays of arbitrary
#' dimensionality (2D matrices, 3D arrays, etc.) as long as all dimensions except
#' the first match across inputs.
#'
#' @param ... Arrays or matrices to be combined. All inputs must have the same number
#'           of dimensions and matching dimensions except for the first.
#'
#' @return An array or matrix containing all input arrays stacked along the first
#'         dimension, with other dimensions preserved.
#'
#' @details
#' The function performs the following operations:
#'
#' 1. Input Validation:
#'    - Checks that at least one matrix is provided
#'    - Verifies all inputs have same number of dimensions
#'    - Confirms all dimensions except the first match across inputs
#'
#' 2. Array Combination:
#'    - Calculates total size needed for first dimension
#'    - Creates zero-filled result array
#'    - Copies each input array into the appropriate position
#'    - Handles special cases for 2D, 3D, and 4D arrays
#'    - Uses dynamic indexing for higher dimensions
#'
#' The function supports arrays of any dimension, though explicit handling is
#' optimized for 2D-4D cases. Higher dimensional arrays are handled using
#' more general but potentially slower methods.
#'
#' @examples
#' \dontrun{
#' # 2D matrices
#' mat1 <- matrix(1:6, nrow=2, ncol=3)
#' mat2 <- matrix(7:12, nrow=2, ncol=3)
#' combined_2d <- combine_matrices(mat1, mat2)
#'
#' # 3D arrays
#' arr1 <- array(1:12, dim=c(2,3,2))
#' arr2 <- array(13:24, dim=c(2,3,2))
#' combined_3d <- combine_matrices(arr1, arr2)
#'
#' # Mixed number of rows
#' mat3 <- matrix(1:9, nrow=3, ncol=3)
#' mat4 <- matrix(10:15, nrow=2, ncol=3)
#' combined_mixed <- combine_matrices(mat3, mat4)
#'}
#' @examples
#' \dontrun{
#' # Will raise an error - different number of columns
#' mat1 <- matrix(1:6, nrow=2, ncol=3)
#' mat2 <- matrix(7:10, nrow=2, ncol=2)
#' combined <- combine_matrices(mat1, mat2)
#'}
combine_matrices <- function(...) {
  # Get list of input matrices
  matrices <- list(...)

  if (length(matrices) == 0) {
    stop("At least one matrix must be provided")
  }

  # Get dimensions of all input matrices
  dims_list <- lapply(matrices, dim)

  # Check if all matrices have same number of dimensions
  if (length(unique(lapply(dims_list, length))) != 1) {
    stop("All input matrices must have the same number of dimensions")
  }

  # Check if all matrices have matching dimensions (except first)
  ref_dims <- dims_list[[1]]
  for (i in seq_along(dims_list)) {
    if (any(dims_list[[i]][-1] != ref_dims[-1])) {
      stop("All matrices must have matching dimensions (except the first dimension)")
    }
  }

  # Calculate total required size for first dimension
  total_rows <- sum(sapply(dims_list, `[`, 1))

  # Create target dimensions
  target_dims <- ref_dims
  target_dims[1] <- total_rows

  # Create zero-filled matrix with target dimensions
  result <- array(0, dim = target_dims)

  # Fill the result matrix with input matrices
  current_row <- 1
  for (i in seq_along(matrices)) {
    current_matrix <- matrices[[i]]
    rows_to_add <- dims_list[[i]][1]
    index_range <- current_row:(current_row + rows_to_add - 1)

    # Handle different dimensionality cases
    if (length(target_dims) == 2) {
      result[index_range, ] <- current_matrix
    } else if (length(target_dims) == 3) {
      result[index_range, , ] <- current_matrix
    } else if (length(target_dims) == 4) {
      result[index_range, , , ] <- current_matrix
    } else {
      # For higher dimensions, use do.call with array indexing
      indices <- rep(list(quote(expr = )), length(target_dims))
      indices[[1]] <- index_range
      result <- do.call(`[<-`, c(list(result), indices, list(current_matrix)))
    }

    current_row <- current_row + rows_to_add
  }

  return(result)
}

#' Define State Space for Spatial Analysis
#'
#' @description
#' Defines the spatial state space for capture-recapture analysis using either
#' explicitly provided vertices or a buffer around sampling locations. The function
#' supports two methods of state space definition: polygon vertices or buffer-based
#' boundaries.
#'
#' @param data List. Must contain either:
#' \itemize{
#'   \item vertices: Matrix of coordinates defining the state space polygon
#'   \item buff: Numeric buffer distance to extend around sampling locations
#' }
#' @param Xall Matrix. Optional. Matrix of sampling location coordinates x, y.
#'             Required when using the buffer method.
#'
#' @return A list containing:
#' \itemize{
#'   \item vertices: Matrix of coordinates defining state space boundaries
#'   \item xlim: Vector of min, max x-coordinates
#'   \item ylim: Vector of min, max y-coordinates
#'   \item useverts: Logical indicating whether vertices method was used
#' }
#'
#' @details
#' The function supports two methods of defining the state space:
#'
#' 1. Vertices Method:
#'    - Uses explicitly provided polygon vertices
#'    - Calculates bounds directly from vertex coordinates
#'    - Suitable for irregular study areas
#'
#' 2. Buffer Method:
#'    - Creates rectangular state space around sampling locations
#'    - Extends boundaries by specified buffer distance
#'    - Requires sampling locations (Xall) to be provided
#'    - Suitable for regular study areas
#'
#' @examples
#' \dontrun{
#' # Using vertices method
#' data <- list(
#'   vertices = matrix(c(
#'     0,0,    # Bottom-left
#'     10,0,   # Bottom-right
#'     10,10,  # Top-right
#'     0,10    # Top-left
#'   ), ncol=2, byrow=TRUE)
#' )
#' state_space_vertices <- define_state_space(data)
#'
#' # Using buffer method
#' data <- list(buff = 500)  # 500 unit buffer
#' Xall <- matrix(c(
#'   100,100,  # Trap 1
#'   200,150,  # Trap 2
#'   150,200   # Trap 3
#' ), ncol=2, byrow=TRUE)
#' state_space_buffer <- define_state_space(data, Xall)
#'}
define_state_space <- function(data, Xall = NULL) {
  # Validate input
  if (!is.list(data)) {
    stop("data must be a list")
  }

  # Create state space based on vertices or buffer
  if ("vertices" %in% names(data)) {
    vertices <- data$vertices
    xlim <- c(min(vertices[, 1]), max(vertices[, 1]))
    ylim <- c(min(vertices[, 2]), max(vertices[, 2]))
    useverts <- TRUE
  } else if ("buff" %in% names(data)) {
    # Ensure Xall is provided when using buffer method
    if (is.null(Xall)) {
      stop("Xall must be provided when using buffer method")
    }

    buff <- data$buff
    xlim <- c(min(Xall[, 1]), max(Xall[, 1])) + c(-buff, buff)
    ylim <- c(min(Xall[, 2]), max(Xall[, 2])) + c(-buff, buff)
    vertices <- cbind(xlim, ylim)
    useverts <- FALSE
  } else {
    stop("user must supply either 'buff' or 'vertices' in data object")
  }

  # Return results as a list
  return(list(
    vertices = vertices,
    xlim = xlim,
    ylim = ylim,
    useverts = useverts
  ))
}

#' Validate Input Parameters for Data Analysis
#'
#' @description
#' Validates input parameters for data analysis, specifically checking the compatibility
#' between ID update methods and observation types in capture-recapture studies.
#'
#' @param input List. Must contain:
#' \itemize{
#'   \item IDup: Character. ID update method, must be either "MH" (Metropolis-Hastings) or "Gibbs"
#'   \item obstype: Character vector. Observation types, where the second element indicates
#'                  the observation model (e.g., "bernoulli")
#' }
#'
#' @return No return value. Throws an error if validation fails.
#'
#' @details
#' The function performs two validation checks:
#' 1. Verifies that IDup is either "MH" or "Gibbs"
#' 2. Ensures that Gibbs updates are not used with bernoulli observation types
#'
#' @examples
#' \dontrun{
#' # Valid input with MH updates and bernoulli observations
#' input <- list(
#'   IDup = "MH",
#'   obstype = c("standard", "bernoulli")
#' )
#' input_check(input)  # Passes validation
#'}
input_check <- function(input) {

  if (!input$IDup %in% c("MH", "Gibbs")) {
    stop("IDup must be MH or Gibbs")
  }
  if (input$obstype[2] == "bernoulli" & input$IDup == "Gibbs") {
    stop("Must use MH IDup for bernoulli data")
  }
}

#' Validate Data Structure for Capture-Recapture Analysis
#'
#' @description
#' Validates and standardizes data structures for capture-recapture analysis,
#' checking genetic data, capture histories, and associated metadata. The function
#' ensures data consistency and proper formatting while correcting matrix formats
#' where possible.
#'
#' @param data List. Must contain:
#' \itemize{
#'   \item G_marked: Matrix or vector of genetic data for marked individuals
#'   \item G_unmarked: Matrix or vector of genetic data for unmarked individuals
#'   \item IDlist: List containing:
#'     - IDcovs: List of genetic covariates
#'     - ncat: Integer number of genetic categories
#'   \item y_mark: Optional 3D array of marking capture histories
#'   \item y_sight_marked: 3D array of sighting histories for marked individuals
#'   \item y_sight_unmarked: 3D array of sighting histories for unmarked individuals
#' }
#'
#' @return The validated and potentially reformatted data list. Matrix formats are
#' corrected where necessary while maintaining the original data structure.
#'
#' @examples
#' \dontrun{
#' # Create valid data structure
#' data <- list(
#'   G_marked = matrix(c(1,2,1, 2,1,2), nrow=2, ncol=3),
#'   G_unmarked = matrix(c(1,1,1, 2,2,2), nrow=2, ncol=3),
#'   IDlist = list(
#'     IDcovs = list(c(1,2), c(1,2), c(1,2)),
#'     ncat = 3
#'   ),
#'   y_sight_marked = array(0, dim=c(2,4,3)),    # 2 individuals, 4 occasions, 3 traps
#'   y_sight_unmarked = array(0, dim=c(2,4,3)),  # 2 individuals, 4 occasions, 3 traps
#'   y_mark = array(0, dim=c(2,4,3))             # 2 individuals, 4 occasions, 3 traps
#' )
#'
#' # Validate data
#' validated_data <- data_check(data)
#' }
data_check <- function(data) {

  if (any(is.na(data$G_marked)) | any(is.na(data$G_unmarked))) {
    stop("Code missing IDcovs with a 0")
  }
  if (!is.matrix(data$G_marked) & !is.null(data$G_marked)) {
    data$G_marked <- matrix(data$G_marked)
  }
  if (!is.matrix(data$G_unmarked)) {
    data$G_marked <- matrix(data$G_unmarked)
  }
  if (!is.list(data$IDlist$IDcovs)) {
    stop("IDcovs must be a list")
  }
  if (ncol(data$G_marked) != data$IDlist$ncat) {
    stop("G_marked needs ncat number of columns")
  }
  if (ncol(data$G_unmarked) != data$IDlist$ncat) {
    stop("G_unmarked needs ncat number of columns")
  }
  if (length(dim(data$y_mark)) != 3 & !is.null(data$y_mark)) {
    stop("dim(y_mark) must be 3. Reduced to 2 during initialization")
  }
  if (length(dim(data$y_sight_marked)) != 3) {
    stop("dim(y_sight_marked) must be 3. Reduced to 2 during initialization")
  }
  if (length(dim(data$y_sight_unmarked)) != 3) {
    stop("dim(y_sight_unmarked) must be 3. Reduced to 2 during initialization")
  }

  return(data)
}

#' Calculate Truncated Normal Log Likelihood
#'
#' @description
#' Calculates the log likelihood of a value under a truncated normal distribution,
#' given the center, bounds, and standard deviation parameters.
#'
#' @param value Numeric. The value at which to calculate the log likelihood.
#' @param center Numeric. The center (mean) of the normal distribution before truncation.
#' @param bounds Numeric vector. A vector of length 2 specifying the
#'              bounds of the truncation.
#' @param sigma Numeric. The standard deviation of the normal distribution.
#'
#' @return Numeric. The log likelihood value of the truncated normal distribution
#'         at the specified point.
#'
#' @examples
#' \dontrun{
#' # Calculate log likelihood for a value within bounds
#' value <- 1.5
#' center <- 1.0
#' bounds <- c(0, 2)
#' sigma <- 0.5
#' log_lik <- calculate_log_likelihood(value, center, bounds, sigma)
#'
#' # Calculate for multiple values
#' values <- seq(0.1, 1.9, by=0.2)
#' log_liks <- sapply(values, calculate_log_likelihood,
#'                    center=1, bounds=c(0,2), sigma=0.5)
#'}
calculate_log_likelihood <- function(value, center, bounds, sigma) {
  log(
    dnorm(value, center, sigma) /
      (pnorm(bounds[2], center, sigma) - pnorm(bounds[1], center, sigma))
  )
}

#' Calculate Log Likelihood for Bernoulli or Poisson Observations
#'
#' @description
#' Calculates the log likelihood for capture-recapture data under either a Bernoulli
#' or Poisson observation model, accounting for detection probability and availability.
#'
#' @param obstype Character. Observation type, either "bernoulli" or "poisson".
#' @param y Numeric. Number of observed captures/detections.
#' @param K Numeric. Number of sampling occasions or effort.
#' @param lambda Numeric. Baseline detection rate.
#' @param z Numeric. Availability indicator (typically 0 or 1).
#'
#' @return Numeric. The log likelihood value under the specified observation model.
#'
#' @details
#' The function computes log likelihoods using:
#'
#' For Bernoulli model:
#' - Detection probability p = 1 - exp(-lambda)
#' - Log likelihood = log(Binomial(y | K, pz))
#'
#' For Poisson model:
#' - Log likelihood = log(Poisson(y | Klambda*z))
#'
#' The z parameter typically indicates whether an individual is available
#' for detection (z=1) or not (z=0).
#'
#' @examples
#' \dontrun{
#' # Poisson example
#' ll_pois <- calculate_ll_bern_pois(
#'   obstype = "poisson",
#'   y = 5,        # 5 detections
#'   K = 10,       # 10 occasions
#'   lambda = 0.5, # baseline rate
#'   z = 1         # individual available
#' )
#'}
calculate_ll_bern_pois <- function(obstype, y, K, lambda, z) {
  if (obstype == "bernoulli") {
    pd <- 1 - exp(-lambda)
    ll <- dbinom(y, K, pd * z, log = TRUE)
  } else {
    ll <- dpois(y, K * lambda * z, log = TRUE)
  }
  return(ll)
}

#' Calculate Pairwise Euclidean Distances Between Points
#'
#' @description
#' Calculates the Euclidean distances between two sets of 2D points. When only one
#' set is provided, calculates distances between all pairs of points within that set.
#'
#' @param x Matrix, data frame, or vector. First set of coordinates. Must be a 2-column
#'          matrix/data frame of coordinates, or a length-2 vector representing a single point.
#' @param y Matrix, data frame, or vector, or NULL. Optional second set of coordinates.
#'          Must be a 2-column matrix/data frame, or a length-2 vector. If NULL,
#'          distances are calculated between all pairs of points in x.
#'
#' @return Matrix of pairwise distances where:
#' \itemize{
#'   \item Rows correspond to points in x
#'   \item Columns correspond to points in y (or x if y is NULL)
#'   \item Each element i,j is the Euclidean distance between point i in x and point j in y
#' }
#'
#' @details
#' The function performs the following operations:
#'
#' 1. Input Validation:
#'    - Ensures x and y (if provided) contain 2D coordinates
#'    - Converts single points to 1-row matrices
#'
#' 2. Distance Calculation:
#'    - Uses Euclidean distance formula: sqrt((x₁-x₂)² + (y₁-y₂)²)
#'    - Efficiently computes all pairwise distances without loops
#'    - Returns result as a matrix
#'
#' When y is NULL, the function produces a symmetric matrix of
#' distances between all pairs of points in x.
#'
#' @examples
#' \dontrun{
#' # Single point to multiple points
#' point <- c(0, 0)
#' points <- matrix(c(1,0, 0,1, 1,1), ncol=2, byrow=TRUE)
#' distances <- pairwise_distances(point, points)
#'
#' # Multiple points to multiple points
#' points1 <- matrix(c(0,0, 1,0), ncol=2, byrow=TRUE)
#' points2 <- matrix(c(0,1, 1,1, 2,1), ncol=2, byrow=TRUE)
#' distances <- pairwise_distances(points1, points2)
#'
#' # All pairwise distances within one set
#' points <- matrix(c(0,0, 1,0, 0,1, 1,1), ncol=2, byrow=TRUE)
#' distances <- pairwise_distances(points)
#'}
pairwise_distances <- function(x, y = NULL) {
  # Ensure x is a 2-column matrix
  if (is.null(dim(x)) && length(x) == 2) {
    x <- matrix(x, nrow = 1)
  }
  if (ncol(x) != 2) {
    stop("Argument 'x' must be a 2-column matrix or data frame, or a length-2 vector.", call. = FALSE)
  }

  # If y is NULL, set y to x
  if (is.null(y)) {
    y <- x
  } else {
    # Ensure y is a 2-column matrix
    if (is.null(dim(y)) && length(y) == 2) {
      y <- matrix(y, nrow = 1)
    }
    if (ncol(y) != 2) {
      stop("Argument 'y' must be a 2-column matrix or data frame, or a length-2 vector.", call. = FALSE)
    }
  }

  # Compute pairwise distances
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = FALSE)
}

#' Check if Point is Within State Space
#'
#' @description
#' Determines whether a 2D point lies within a defined state space, which can be
#' either a rectangular boundary or an arbitrary polygon defined by vertices.
#'
#' @param point Numeric vector. A vector of length 2 containing x, y coordinates
#'             of the point to check.
#' @param state_space List. A list defining the state space with components:
#' \itemize{
#'   \item useverts: Logical. If TRUE, uses polygon vertices; if FALSE, uses rectangular bounds
#'   \item xlim: Numeric vector of length 2. min, max x-coordinates for rectangular bounds
#'   \item ylim: Numeric vector of length 2. min, max y-coordinates for rectangular bounds
#'   \item vertices: Matrix. Required if useverts=TRUE. Matrix of polygon vertex coordinates
#' }
#'
#' @return Logical. TRUE if the point lies within the state space, FALSE otherwise.
#'
#' @details
#' The function supports two methods of checking point inclusion:
#'
#' 1. Rectangular Bounds (useverts = FALSE):
#'    - Uses simple coordinate comparison
#'    - Checks if point lies within xlim and ylim ranges
#'    - Faster but limited to rectangular areas
#'
#' 2. Polygon Bounds (useverts = TRUE):
#'    - Uses point-in-polygon algorithm
#'    - Supports arbitrary polygon shapes
#'    - More flexible but computationally more intensive
#'
#' @examples
#' \dontrun{
#' # Check point in rectangular state space
#' rect_space <- list(
#'   useverts = FALSE,
#'   xlim = c(0, 10),
#'   ylim = c(0, 10)
#' )
#' point_in_area(c(5, 5), rect_space)  # TRUE
#' point_in_area(c(-1, 5), rect_space) # FALSE
#'
#' # Check point in polygon state space
#' poly_space <- list(
#'   useverts = TRUE,
#'   vertices = matrix(c(
#'     0,0,
#'     10,0,
#'     10,10,
#'     5,15,
#'     0,10
#'   ), ncol=2, byrow=TRUE)
#' )
#' point_in_area(c(5, 5), poly_space)   # TRUE
#' point_in_area(c(5, 12), poly_space)  # TRUE
#' point_in_area(c(-1, 5), poly_space)  # FALSE
#' }
point_in_area <- function(point,
                          state_space) {
  # Validate inputs
  if (!is.numeric(point) || length(point) != 2) {
    stop("point must be a numeric vector of length 2")
  }

  # Check using rectangular boundaries
  if (!state_space$useverts) {
    # Require both xlim and ylim for rectangular check
    if (is.null(state_space$xlim) || is.null(state_space$ylim) ||
        length(state_space$xlim) != 2 || length(state_space$ylim) != 2) {
      stop("For rectangular check, provide valid xlim and ylim")
    }

    return(
      point[1] < state_space$xlim[2] &
        point[1] > state_space$xlim[1] &
        point[2] < state_space$ylim[2] &
        point[2] > state_space$ylim[1]
    )
  }

  # Check using polygon vertices
  if (state_space$useverts) {
    # Require vertices
    if (is.null(state_space$vertices)) {
      stop("When useverts = TRUE, vertices must be provided")
    }

    # Determine if point is in polygon
    return(in_poly(point, state_space$vertices))
  }
}

#' Check if Point is Inside Polygon Using Ray-Casting Algorithm
#'
#' @description
#' Determines whether a 2D point lies within a polygon using the ray-casting
#' algorithm (also known as the even-odd rule or crossing number algorithm). The
#' function automatically closes the polygon by connecting the last vertex to the first.
#'
#' @param point Numeric vector. A vector of length 2 containing x, y coordinates
#'             of the point to check.
#' @param vertices Matrix or data frame. A two-column matrix or data frame containing
#'                the x, y coordinates of the polygon vertices in order. The polygon
#'                does not need to be closed (the function handles this automatically).
#'
#' @return Logical. TRUE if the point lies within the polygon, FALSE otherwise.
#'
#' @details
#' The ray-casting algorithm works by:
#'
#' 1. Casting a horizontal ray from the test point to the right
#' 2. Counting the number of times this ray intersects the polygon edges
#' 3. Using the even-odd rule to determine inclusion:
#'    - Odd number of intersections: point is inside
#'    - Even number of intersections: point is outside
#'
#' Implementation details:
#' 1. Input Validation:
#'    - Ensures point is a length-2 numeric vector
#'    - Verifies vertices is a 2-column matrix/data frame
#'
#' 2. Polygon Preparation:
#'    - Automatically closes the polygon by appending first vertex
#'    - Processes edges sequentially
#'
#' 3. Intersection Testing:
#'    - Tests if horizontal ray intersects each edge
#'    - Handles special cases like vertices and horizontal edges
#'
#' @examples
#' \dontrun{
#' # Create a triangular polygon
#' vertices <- matrix(c(
#'   0,0,   # Bottom-left
#'   2,0,   # Bottom-right
#'   1,2    # Top
#' ), ncol=2, byrow=TRUE)
#'
#' # Test points
#' in_poly(c(1,1), vertices)    # TRUE - point inside triangle
#' in_poly(c(0,2), vertices)    # FALSE - point outside triangle
#' in_poly(c(1,0), vertices)    # TRUE - point on edge
#'
#' # Create a complex polygon
#' vertices <- matrix(c(
#'   0,0,    # Start
#'   2,0,
#'   2,2,
#'   1,1,
#'   0,2     # End
#' ), ncol=2, byrow=TRUE)
#'
#' # Test points in complex polygon
#' in_poly(c(0.5,0.5), vertices)  # TRUE
#' in_poly(c(1.5,1.5), vertices)  # TRUE
#' in_poly(c(0.5,1.5), vertices)  # FALSE
#'}
in_poly <- function(point, vertices) {
  # Input validation
  if (!is.numeric(point) || length(point) != 2) {
    stop("point must be a numeric vector of length 2")
  }

  if (!is.matrix(vertices) && !is.data.frame(vertices)) {
    stop("vertices must be a matrix or data frame")
  }

  if (ncol(vertices) != 2) {
    stop("vertices must have exactly two columns (x and y)")
  }

  # Ensure vertices form a closed polygon by appending the first vertex to the end
  vertices <- rbind(vertices, vertices[1,])

  # Ray-casting algorithm
  intersect_count <- 0
  x <- point[1]
  y <- point[2]

  for (i in 1:(nrow(vertices) - 1)) {
    # Get current edge vertices
    x1 <- vertices[i, 1]
    y1 <- vertices[i, 2]
    x2 <- vertices[i+1, 1]
    y2 <- vertices[i+1, 2]

    # Check if the ray intersects the edge
    # Condition 1: y-coordinate is within the edge's y range
    # Condition 2: horizontal ray from point crosses the edge
    if (((y1 > y) != (y2 > y)) &&
        (x < (x2 - x1) * (y - y1) / (y2 - y1) + x1)) {
      intersect_count <- intersect_count + 1
    }
  }

  # Odd number of intersections means point is inside
  return(intersect_count %% 2 == 1)
}


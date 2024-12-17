# NEED TO INCORPORATE CONVENTIONAL MODEL (OMIT ALL MARKING PROCESSES)


#' Fit the generalized categorical spatial mark resight model. SMR_move fits the activity centers to cameras and snares separately whereas SMR fits them concurrently.
#'
#' Allowing for activity center relocation (act_center = "move") means the model
#' incorporates a process for updating activity centers between the marking and
#' sighting phases. If the data are sparse, there is no telemetry data, or if
#' the study period is short, it may better to not allow for activity center
#' relocation (act_center = "no move").
#'
#' @param data a data list as formatted by sim.genCatSMR(). See description for more details.
#' @param input Model definitions and inputs
#' @description This function fits the generalized categorical spatial mark resight model.
#' Modelling the marking process relaxes the assumption that the distribution of marked
#' individuals across the landscape is spatially uniform.
#'
#' the data list should be formatted to match the list outputted by
#' sim.genCatSMR(), but not all elements of that object are necessary.
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
    niter = 2400, nburn = 1200, nthin = 5, M = 200, act_center = "no move",
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
  act_center <- input$act_center

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

  # Initializations
  y_mark <- data$y_mark
  y_sight_marked <- data$y_sight_marked
  y_sight_unmarked <- data$y_sight_unmarked
  X1 <- as.matrix(data$X1)
  X2 <- as.matrix(data$X2)
  J1 <- nrow(X1)
  J2 <- nrow(X2)
  K1 <- data$K1
  K2 <- data$K2
  ncat <- data$IDlist$ncat
  nallele <- data$IDlist$nallele
  IDcovs <- data$IDlist$IDcovs
  buff <- data$buff
  Xall <- rbind(X1, X2)
  n_marked <- data$n_marked
  G_marked <- data$G_marked
  tf1 <- data$tf1
  tf2 <- data$tf2

  tune_parms <- c("lam0_mark", "lam0_sight", "sigma_p", "sigma_d", "s1", "s2")

  accept_rates <- tibble::tibble(
    acc_iter = rep(1:niter, each = length(proppars)),
    label = rep(names(proppars), niter),
    accept = dplyr::case_when(label %in% tune_parms ~ 0,
                              .default = 0.45)
  )

  tune.check <- 50
  batch_n <- 10

  if ("G_unmarked" %in% names(data)) {
    G_unmarked <- data$G_unmarked
    if (length(dim(y_sight_unmarked)) != 3) {
      stop("dim(y_sight_unmarked) must be 3. Reduced to 2 during initialization")
    }
    useUM <- TRUE
  } else {
    G_unmarked <- matrix(0, nrow = 0, ncol = ncat)
    useUM <- FALSE
  }

  if (any(is.na(G_marked)) | any(is.na(G_unmarked))) {
    stop("Code missing IDcovs with a 0")
  }
  if (!is.matrix(G_marked)) {
    G_marked <- matrix(G_marked)
  }
  if (!is.matrix(G_unmarked)) {
    G_marked <- matrix(G_unmarked)
  }
  if (!is.list(IDcovs)) {
    stop("IDcovs must be a list")
  }
  nlevels <- unlist(lapply(IDcovs, length))
  if (ncol(G_marked) != ncat) {
    stop("G_marked needs ncat number of columns")
  }
  if (ncol(G_unmarked) != ncat) {
    stop("G_unmarked needs ncat number of columns")
  }
  # Are there unknown marked status guys?
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
  # Are there marked no ID guys?
  useMarkednoID <- FALSE
  if ("G_marked.noID" %in% names(data)) {
    if (!is.na(data$G_marked.noID[1])) {
      G_marked.noID <- data$G_marked.noID
      if (!is.matrix(G_marked.noID)) {
        G_marked.noID <- matrix(G_marked.noID)
      }
      if (ncol(G_marked.noID) != ncat) {
        stop("G_marked.noID needs ncat number of columns")
      }
      y_sight_marked_noID <- data$y_sight_marked_noID
      useMarkednoID <- TRUE
    }
  }

  # data checks
  if (length(dim(y_mark)) != 3 & !is.null(y_mark)) {
    stop("dim(y_mark) must be 3. Reduced to 2 during initialization")
  }
  if (length(dim(y_sight_marked)) != 3) {
    stop("dim(y_sight_marked) must be 3. Reduced to 2 during initialization")
  }
  if (length(dim(y_sight_unmarked)) != 3) {
    stop("dim(y_sight_unmarked) must be 3. Reduced to 2 during initialization")
  }
  if (useUnk) {
    if (length(dim(y_sight_unk)) != 3) {
      stop("dim(y_sight_unk) must be 3. Reduced to 2 during initialization")
    }
  }
  if (useMarkednoID) {
    if (length(dim(y_sight_marked_noID)) != 3) {
      stop("dim(y_sight_marked_noID) must be 3. Reduced to 2 during initialization")
    }
  }
  if (!IDup %in% c("MH", "Gibbs")) {
    stop("IDup must be MH or Gibbs")
  }
  if (obstype[2] == "bernoulli" & IDup == "Gibbs") {
    stop("Must use MH IDup for bernoulli data")
  }

  # If using polygon state space
  if ("vertices" %in% names(data)) {
    vertices <- data$vertices
    useverts <- TRUE
    xlim <- c(min(vertices[, 1]), max(vertices[, 1]))
    ylim <- c(min(vertices[, 2]), max(vertices[, 2]))
  } else if ("buff" %in% names(data)) {
    buff <- data$buff
    xlim <- c(min(Xall[, 1]), max(Xall[, 1])) + c(-buff, buff)
    ylim <- c(min(Xall[, 2]), max(Xall[, 2])) + c(-buff, buff)
    vertices <- cbind(xlim, ylim)
    useverts <- FALSE
  } else {
    stop("user must supply either 'buff' or 'vertices' in data object")
  }

  ## pull out initial values
  psi <- inits$psi
  lam0_mark <- inits$lam0_mark
  lam0_sight <- inits$lam0_sight
  sigma_d <- inits$sigma_d
  sigma_p <- inits$sigma_p
  gamma <- inits$gamma

  if (useUnk & !useMarkednoID) {
    G_use <- rbind(G_unmarked, G_unk)
    status <- c(rep(2, nrow(G_unmarked)), rep(0, nrow(G_unk)))
    G_use <- cbind(G_use, status)
    G_marked <- cbind(G_marked, rep(1, nrow(G_marked)))
    ncat <- ncat + 1
    y_sight_latent <- abind::abind(y_sight_unmarked, y_sight_unk, along = 1)
  } else if (!useUnk & useMarkednoID) {
    G_use <- rbind(G_unmarked, G_marked.noID)
    status <- c(rep(2, nrow(G_unmarked)), rep(1, nrow(G_marked.noID)))
    G_use <- cbind(G_use, status)
    G_marked <- cbind(G_marked, rep(1, nrow(G_marked)))
    ncat <- ncat + 1
    y_sight_latent <- abind::abind(y_sight_unmarked, y_sight_marked_noID, along = 1)
  } else if (useUnk & useMarkednoID) {
    G_use <- rbind(G_unmarked, G_unk, G_marked.noID)
    status <- c(rep(2, nrow(G_unmarked)), rep(0, nrow(G_unk)), rep(1, nrow(G_marked.noID)))
    G_use <- cbind(G_use, status)
    G_marked <- cbind(G_marked, rep(1, nrow(G_marked)))
    ncat <- ncat + 1
    nlevels <- c(nlevels, 2)
    y_sight_latent <- abind::abind(y_sight_unmarked, y_sight_unk, y_sight_marked_noID, along = 1)
  } else {
    G_use <- G_unmarked
    y_sight_latent <- y_sight_unmarked
  }
  n_samp_latent <- nrow(y_sight_latent) # number of unmarked individuals + uncertainty
  if (is.na(nswap)) {
    nswap <- round(n_samp_latent / 2)
    # warning("nswap not specified, using round(n_samp_latent/2)")
  }

  # make constraints for data initialization
  constraints <- matrix(1, nrow = n_samp_latent, ncol = n_samp_latent)
  for (i in 1:n_samp_latent) {
    for (j in 1:n_samp_latent) {
      guys1 <- which(G_use[i, ] != 0)
      guys2 <- which(G_use[j, ] != 0)
      comp <- guys1[which(guys1 %in% guys2)]
      if (any(G_use[i, comp] != G_use[j, comp])) {
        constraints[i, j] <- 0
      }
    }
  }
  # If bernoulli data, add constraints that prevent y.true[i,j,k]>1
  binconstraints <- FALSE
  if (obstype[2] == "bernoulli") {
    idx <- t(apply(y_sight_latent, 1, function(x) {
      which(x > 0, arr.ind = TRUE)
    }))
    for (i in 1:n_samp_latent) {
      for (j in 1:n_samp_latent) {
        if (i != j) {
          if (all(idx[i, 1:2] == idx[j, 1:2])) {
            constraints[i, j] <- 0 # can't combine samples from same trap and occasion in binomial model
            constraints[j, i] <- 0
            binconstraints <- TRUE
          }
        }
      }
    }
  }

  #######################
  # not in concat
  # marking occasion order constraints
  Kconstraints <- matrix(0, nrow = M, ncol = n_samp_latent)
  if ("markedS" %in% names(data)) {
    markedS <- data$markedS
  } else {
    markedS <- matrix(1, nrow = n_marked, ncol = K2)
  }
  for (i in 1:n_marked) {
    for (j in 1:n_samp_latent) {
      occ <- which(apply(y_sight_latent[j, , ], 2, sum) > 0)
      if (markedS[i, occ] == 1) {
        Kconstraints[i, j] <- 1
      } else if (markedS[i, occ] == 2) {
        Kconstraints[i, j] <- 2
      }
    }
  }
  #######################

  # Build y_sight_true
  y_sight_true <- array(0, dim = c(M, J2, K2))
  y_sight_true[1:n_marked, , ] <- y_sight_marked
  ID <- rep(NA, n_samp_latent)
  idx <- n_marked + 1

  # assign all unmarked and unk samples to unmarked guys first
  for (i in 1:n_samp_latent) {
    if (useMarkednoID) {
      if (status[i] == 1) next
    }
    if (idx > M) {
      stop("Need to raise M to initialize y.true")
    }
    traps <- which(rowSums(y_sight_latent[i, , ]) > 0)
    y_sight_true2D <- apply(y_sight_true, c(1, 2), sum)
    if (length(traps) == 1) {
      cand <- which(y_sight_true2D[, traps] > 0) # guys caught at same traps
    } else {
      cand <- which(rowSums(y_sight_true2D[, traps]) > 0) # guys caught at same traps
    }
    cand <- cand[cand > n_marked]
    if (length(cand) > 0) {
      if (length(cand) > 1) { # if more than 1 ID to match to, choose first one
        cand <- cand[1]
      }
      # Check constraint matrix
      cands <- which(ID %in% cand) # everyone assigned this ID
      if (all(constraints[i, cands] == 1)) {
        # focal consistent with all partials already assigned and consistent with marked guy capture occasion
        y_sight_true[cand, , ] <- y_sight_true[cand, , ] + y_sight_latent[i, , ]
        ID[i] <- cand
      } else { # focal not consistent
        y_sight_true[idx, , ] <- y_sight_latent[i, , ]
        ID[i] <- idx
        idx <- idx + 1
      }
    } else { # no assigned samples at this trap
      y_sight_true[idx, , ] <- y_sight_latent[i, , ]
      ID[i] <- idx
      idx <- idx + 1
    }
  }

  ###################################
  # assign marked unknown ID guys
  ###################################
  if (useMarkednoID) { # Need to initialize these guys to marked guys
    fix <- which(status == 1)
    meanloc <- matrix(NA, nrow = n_marked, ncol = 2)
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
    for (i in 1:nrow(G_marked.noID)) {
      trap <- which(rowSums(y_sight_latent[i, , ]) > 0)
      compatible <- rep(FALSE, n_marked)
      for (j in 1:n_marked) {
        nonzero1 <- G_marked[j, 1:(ncat - 1)] != 0
        nonzero2 <- G_marked.noID[i, ] != 0
        nonzero <- which(nonzero1 & nonzero2)
        if (all(G_marked[j, nonzero] == G_marked.noID[i, nonzero])) {
          compatible[j] <- TRUE
        }
      }
      if (all(compatible == FALSE)) {
        stop(paste("No G_marked compatible with G_marked.noID "), i)
      }
      dists <- sqrt((X2[trap, 1] - meanloc[, 1])^2 + (X2[trap, 2] - meanloc[, 2])^2)
      dists[!compatible] <- Inf
      dists[which(Kconstraints[1:n_marked, fix[i]] == 0)] <- Inf # Exclude guys not marked yet
      if (all(is.finite(dists) == FALSE)) {
        stop(paste("No G_marked compatible with G_marked.noID "), i)
      }
      ID[fix[i]] <- which(dists == min(dists, na.rm = TRUE))[1]
      y_sight_true[ID[fix[i]], , ] <- y_sight_true[ID[fix[i]], , ] + y_sight_latent[fix[i], , ]
    }
  }
  if (binconstraints) {
    if (any(y_sight_true > 1)) stop("bernoulli data not initialized correctly")
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

  #######################
  # not in concat
  # Check marked guy k constraints
  if (any(data$markedS == 0 | data$markedS == 2)) {
    check <- which(ID <= n_marked)
    for (i in check) {
      if (Kconstraints[ID[i], i] == 0) {
        stop("ID initialized improperly for marked no ID samples")
      }
    }
  }
  #######################

  y_sight_true <- apply(y_sight_true, c(1, 2), sum)
  #######################
  # not in concat
  y_mark <- pad_matrix(y_mark, dims = c(M, J1, K1))
  y_mark2D <- apply(y_mark, c(1, 2), sum)
  #######################
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
  s1 <- cbind(runif(M, xlim[1], xlim[2]), runif(M, ylim[1], ylim[2])) # assign random locations
  y.all2D <- cbind(y_mark2D, y_sight_true)
  idx <- which(rowSums(y.all2D) > 0) # switch for those actually caught
  for (i in idx) {
    trps <- matrix(Xall[y.all2D[i, ] > 0, 1:2], ncol = 2, byrow = FALSE)
    if (nrow(trps) > 1) {
      s1[i, ] <- c(mean(trps[, 1]), mean(trps[, 2]))
    } else {
      s1[i, ] <- trps
    }
  }
  if (useverts == TRUE) {
    inside <- rep(NA, nrow(s1))
    for (i in 1:nrow(s1)) {
      inside[i] <- point_in_area(s1[i, ], vertices)
    }
    idx <- which(inside == FALSE)
    if (length(idx) > 0) {
      for (i in 1:length(idx)) {
        while (inside[idx[i]] == FALSE) {
          s1[idx[i], ] <- c(runif(1, xlim[1], xlim[2]), runif(1, ylim[1], ylim[2]))
          inside[idx[i]] <- point_in_area(s1[idx[i], ], vertices)
        }
      }
    }
  }
  s2 <- s1

  # collapse unmarked data to 2D
  y_sight_latent <- apply(y_sight_latent, c(1, 2), sum)

  # Initialize G_true
  G_true <- matrix(0, nrow = M, ncol = ncat)
  G_true[1:n_marked, ] <- G_marked
  for (i in unique(ID)) {
    idx <- which(ID == i)
    if (length(idx) == 1) {
      G_true[i, ] <- G_use[idx, ]
    } else {
      if (ncol(G_use) > 1) {
        G_true[i, ] <- apply(G_use[idx, ], 2, max) # consensus
      } else {
        G_true[i, ] <- max(G_use[idx, ])
      }
    }
  }
  if (useUnk | useMarkednoID) { # augmented guys are unmarked.
    if (max(ID) < M) {
      G_true[(max(ID) + 1):M, ncol(G_true)] <- 2
    }
    unkguys <- which(G_use[, ncol(G_use)] == 0)
  }

  G_latent <- G_true == 0 # Which genos can be updated?
  if (!(useUnk | useMarkednoID)) {
    for (j in 1:(ncat)) {
      fix <- G_true[, j] == 0
      G_true[fix, j] <- sample(IDcovs[[j]], sum(fix), replace = TRUE, prob = gamma[[j]])
    }
  } else {
    for (j in 1:(ncat - 1)) {
      fix <- G_true[, j] == 0
      G_true[fix, j] <- sample(IDcovs[[j]], sum(fix), replace = TRUE, prob = gamma[[j]])
    }
    # Split marked status back off
    Mark.obs <- G_use[, ncat]
    # Mark.status=G_true[,ncat]
    ncat <- ncat - 1
    G_use <- G_use[, 1:ncat]
    G_true <- G_true[, 1:ncat]
  }
  if (!is.matrix(G_use)) {
    G_use <- matrix(G_use, ncol = 1)
  }
  if (!is.matrix(G_true)) {
    G_true <- matrix(G_true, ncol = 1)
  }
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
  if (!is.null(y_mark)) {
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
  if (act_center == "move" & !is.null(y_mark)){
    ll_s2 <- log(dnorm(s2[, 1], s1[, 1], sigma_p) /
                   (pnorm(xlim[2], s1[, 1], sigma_p) - pnorm(xlim[1], s1[, 1], sigma_p)))
    ll_s2 <- ll_s2 + log(dnorm(s2[, 2], s1[, 2], sigma_p) /
                           (pnorm(ylim[2], s1[, 2], sigma_p) - pnorm(ylim[1], s1[, 2], sigma_p)))
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
          # llytrapsum <- sum(ll_y_mark_cand)
          accept_rates <- accept_rates %>%
            dplyr::mutate(
              accept = dplyr::case_when(
                acc_iter == iter & label == "lam0_mark" ~ 1,
                .default = accept
              )
            )
        }
      }
    } else {
      lam0_mark <- 0
      lamd_trap <- 0
      ll_y_mark <- 0
      # llytrapsum <- 0
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
        accept_rates <- accept_rates %>%
          dplyr::mutate(accept = dplyr::case_when(acc_iter == iter &
                                                    label == "lam0_sight" ~
                                                    1,
                                                  .default = accept))
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
        lamd_trap_cand <- 0
        ll_y_mark_cand <- 0
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
        accept_rates <- accept_rates %>%
          dplyr::mutate(accept = dplyr::case_when(acc_iter == iter &
                                                    label == "sigma_d" ~
                                                    1,
                                                  .default = accept))
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
          if (Mark.obs[l] == 2) { # This is an unmarked sample
            if (any(data$markedS == 0 | data$markedS == 2)) {
              possible <- possible[which(Kconstraints[possible, l] == 0)] # k marked status constraints
            } else {
              possible <- possible[possible > n_marked] # Can't swap to a marked guy
            }
          }
          if (Mark.obs[l] == 1) { # This is a marked sample
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
          if (Mark.obs[l] == 2) { # This is an unmarked sample
            if (any(data$markedS == 0 | data$markedS == 2)) {
              possible <- possible[which(Kconstraints[possible, l] == 0)] # k marked status constraints
            } else {
              possible <- possible[possible > n_marked] # Can't swap to a marked guy
            }
          }
          if (Mark.obs[l] == 1) { # This is a marked sample
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
    pbar.trap <- (1 - pd_trap)^K2D1
    pbar.sight <- (1 - pd_sight)^K2D2
    prob0.trap <- exp(rowSums(log(pbar.trap)))
    prob0.sight <- exp(rowSums(log(pbar.sight)))
    prob0 <- prob0.trap * prob0.sight


    fc <- prob0 * psi / (prob0 * psi + 1 - psi)
    z[known_vector == 0] <- rbinom(sum(known_vector == 0), 1, fc[known_vector == 0])
    if (obstype[1] == "bernoulli") {
      ll_y_mark <- dbinom(y_mark2D, K2D1, pd_trap * z, log = TRUE)
    } else {
      ll_y_mark <- dpois(y_mark2D, K2D1 * lamd_trap * z, log = TRUE)
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
    if (act_center == "move") {
      s_2_accept <- 0
      s_1_accept <- 0
      # s1 - activity centers for marked individuals
      for (i in 1:M) {
        Scand <- c(rnorm(1, s1[i, 1], proppars$s1), rnorm(1, s1[i, 2], proppars$s1))

        inbox <- point_in_area(Scand, xlim, ylim, vertices, useverts)
        if (inbox) {
          d1tmp <- sqrt((Scand[1] - X1[, 1])^2 + (Scand[2] - X1[, 2])^2)
          lamd_trap_cand[i, ] <- lam0_mark *
            exp(-d1tmp * d1tmp / (2 * sigma_d * sigma_d))

          # update log likelihood for s2
          ll_s2_cand[i] <-
            calculate_log_likelihood(s2[i, 1], Scand[1], xlim, sigma_p) +
            calculate_log_likelihood(s2[i, 2], Scand[2], ylim, sigma_p)

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

        inbox <- point_in_area(Scand, xlim, ylim, vertices, useverts)
        if (inbox) {
          d2tmp <- sqrt((Scand[1] - X2[, 1])^2 + (Scand[2] - X2[, 2])^2)
          lamd_sight_cand[i, ] <- lam0_sight *
            exp(-d2tmp * d2tmp / (2 * sigma_d * sigma_d))

          # movement likelihood
          ll_s2_cand[i] <-
            calculate_log_likelihood(Scand[1], s1[i, 1], xlim, sigma_p) +
            calculate_log_likelihood(Scand[2], s1[i, 2], ylim, sigma_p)

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

        inbox <- point_in_area(Scand, xlim, ylim, vertices, useverts)
        if (inbox) {
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
            # ll_s2[i] <- ll_s2_cand[i]
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

    accept_rates <- accept_rates %>%
      dplyr::mutate(accept = dplyr::case_when(acc_iter == iter &
                                                label == "s1" ~
                                                sum(s_1_accept) / M,
                                              .default = accept))
    accept_rates <- accept_rates %>%
      dplyr::mutate(accept = dplyr::case_when(acc_iter == iter &
                                                label == "s2" ~
                                                sum(s_2_accept) / M,
                                              .default = accept))



    # update sigma_p
    if (act_center == "move" & !is.null(y_mark)){
      sigma_p_cand <- rnorm(1, sigma_p, proppars$sigma_p)
      if (sigma_p_cand > 0) {
        ll_s2_cand <- log(dnorm(s2[, 1], s1[, 1], sigma_p_cand) /
                            (pnorm(xlim[2], s1[, 1], sigma_p_cand) - pnorm(xlim[1], s1[, 1], sigma_p_cand)))
        ll_s2_cand <- ll_s2_cand + log(dnorm(s2[, 2], s1[, 2], sigma_p_cand) /
                                         (pnorm(ylim[2], s1[, 2], sigma_p_cand) - pnorm(ylim[1], s1[, 2], sigma_p_cand)))


        if (runif(1) < exp(sum(ll_s2_cand) - sum(ll_s2))) {
          sigma_p <- sigma_p_cand
          ll_s2 <- ll_s2_cand
          accept_rates <- accept_rates %>%
            dplyr::mutate(accept = dplyr::case_when(acc_iter == iter &
                                                      label == "sigma_p" ~
                                                      1,
                                                    .default = accept))
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
    if(iter%%tune.check == 0){
      batch_n <- batch_n+1
      delta_n <- batch_n^-1
      # delta_n <- min(0.001,batch_n^-1)
      mean_accept <- accept_rates %>%
        dplyr::group_by(label) %>%
        dplyr::filter(acc_iter > (iter - tune.check + 1) & acc_iter <= iter) %>%
        dplyr::summarise(val = sum(accept)/tune.check) #Need to add/subtract 1 to tune check?

      # # Account for augmented data
      # mean_accept$val[mean_accept$label == "s1"] <- mean_accept$val[mean_accept$label == "s1"]/M
      # mean_accept$val[mean_accept$label == "s2"] <- mean_accept$val[mean_accept$label == "s2"]/M

      # # Add missing proppars
      # noname_index <- which(!(names(proppars) %in% mean_accept$label))
      # if (length(noname_index) > 0) {
      #   mean_accept <- mean_accept %>%
      #     dplyr::add_row(label = names(proppars)[noname_index],
      #             val = 0.45)
      # }
      mean_accept <- mean_accept %>%
        dplyr::mutate(label = factor(label, names(proppars))) %>%
        dplyr::arrange(label)

      # Adjust tuning parms with nonpropotional values
      # proppars[mean_accept$val > 0.45] <- unlist(proppars[mean_accept$val > 0.45]) + delta_n
      # proppars[mean_accept$val < 0.45] <- unlist(proppars[mean_accept$val < 0.45]) - delta_n

      # Adjust tuning parms proportional to acceptance rates
      proppars[mean_accept$val > 0.45] <- unlist(proppars[mean_accept$val > 0.45])/(1-abs(mean_accept$val[mean_accept$val > 0.45] - 0.45))
      proppars[mean_accept$val < 0.45] <- unlist(proppars[mean_accept$val < 0.45])*(1-abs(mean_accept$val[mean_accept$val < 0.45] - 0.45))
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
        Mark.obs <- rep(2, n_samp_latent)
      }
      for (l in 1:length(ID)) {
        cons <- Kconstraints[IDout[, l], l]
        if (any(cons == 1) & Mark.obs[l] == 2) {
          stop("Unmarked samples incorrectly assigned in MCMC!")
        }
        if (any(cons == 0) & Mark.obs[l] == 1) {
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

# Notes:
# Marked individuals are captured with snares and may be used with telemetry data. Relevant data parms are X1.in (snare locations), snares.in (individual tags and tag ids, capture history for each snare). Effort is recorded in effort.in
# Marked ("tagged") and unmarked ("unmarked") individuals are captured at cameras at camera locations, X2.in. Capture histories are recorded in cams.in with number captured (n.det), sex, marked (also tag.class), and gps (unused). Effort is extrapolated from the number of unique dates in cams.in

devtools::load_all()

# We can probably reorganize this...
library(readr)
library(tidyverse)
# library(AHMbook)
library(doParallel)
library(spdgt.sim)

# Download Google Drive files to temp folder
file_names <- c("capture_xy_all.csv",
                "BF_resight_xy.csv",
                "capture_obs_all.csv",
                "BF_resight_obs.csv",
                "BF_gps.csv",
                "capture_effort_all.csv"
)

# Create and download temp files for each csv file
tmpfl <- c()
for (ff in 1:length(file_names)) {
  tmpfl[ff] <- tempfile(fileext = ".csv")

  drive_folder <- googledrive::drive_find(
    file_names[ff],
    shared_drive = "NMDGF"
  ) %>%
    googledrive::drive_download(
      path = tmpfl[ff],
      overwrite = TRUE
    )
}

# Store all relevant outputs
all_out <- tibble(
  model = character(),
  parm = character(),
  means = double(),
  sds = double()
)

reps <- 1:3 # run 3 models

# Initializations

# MCMC Initializations
storeLatent <- TRUE
storeGamma <- TRUE
niter <- 30000
nburn <- 20000
nthin <- 1
# niter <- 5000
# nburn <- 1000
# nthin <- 1

# Collect data if on Google Drive
# Trap locations for marking (X1; snares) and sighting (X2; cameras) processes (required: x, y)
X1.in <- readr::read_csv(tmpfl[1], show_col_types = FALSE)
X2.in <- readr::read_csv(tmpfl[2], show_col_types = FALSE)
# Capture history of marked individuals classified by animal IDs (required: tag.id, IDs for all traps)
snares.in <- readr::read_csv(tmpfl[3], show_col_types = FALSE)
# Camera captures (required: week.id, tag.class, tag.id, location.id, n.det, sex.id)
cams.in <- readr::read_csv(tmpfl[4], show_col_types = FALSE)
# Telemetry data (required: tag.id, week.id, x, y)
locs.in <- readr::read_csv(tmpfl[5], show_col_types = FALSE)
# # Reorganize tag_id in location data to align with snare tag_id
locs.in$tag_id[locs.in$tag_id %in% 11:14] <- locs.in$tag_id[locs.in$tag_id %in% 11:14] - 1

# Effort (required: days)
effort.in <- readr::read_csv(tmpfl[6], show_col_types = FALSE)

# Trap locations for the marking process (X1; snares) and sighting process (X2; cameras)
x.scale <- 1000
y.scale <- 1000
X1 <- X1.in %>%
  summarise(x = x/x.scale,
            y = y/y.scale)
X2 <- X2.in %>%
  summarise(x = x/x.scale,
            y = y/y.scale)

# # Buffer is 1.5 times the average home range "diameter"
buff <- locs.in %>%
  group_by(tag) %>%
  summarise(mx = mean(c((max(x) - min(x))/x.scale,
                        (max(y) - min(y))/y.scale))/2) %>%
  ungroup() %>%
  summarise(means = mean(mx)*1.5) %>%
  pull(means)
# buff <- 100 # Distance to buffer to create the state space

# Number of traps
J1 <- nrow(X1)
J2 <- nrow(X2)

# Redefine data variables
snares <- as.data.frame(snares.in)
cams <- as.data.frame(cams.in)

K1 <- max(effort.in$days) # number of marking occasions
K2 <- length(unique(cams$week.id)) #number of sighting occasions

# Define tagged, untagged cams
# Delete hound traps (INCLUDE ALL TRAPS FOR NOW)
# X1=X1[1:31,]
cams.tagged <- cams[cams$tag.class == "tagged", ]
cams.untagged <- cams[cams$tag.class == "unmarked", ]
n.marked <- nrow(snares) # number of individuals captured and marked
n.unmarked <- nrow(cams.untagged)

# Augmented data. Note: large M could result in long run times
# M <- n.marked*10
M <- 160

# Complete sighting history for all captured individuals
y.mark <- array(0, dim = c(n.marked, J1, K1))
y.mark[, , 1] <- as.matrix(snares[, -c(1, 2)])

# Sighting history of the marked, observed marked status, and individually identified samples
y.sight.marked = array(0, dim = c(n.marked, J2, K2))

cams.tagged.marked <- cams.tagged %>%
  filter(n.det>0) %>%
  group_by(tag.id, location.id, week.id) %>%
  summarise(sum.n.det = sum(n.det),
            .groups = 'drop')

marked.indices <- cbind(cams.tagged.marked$tag.id,
                        cams.tagged.marked$location.id,
                        cams.tagged.marked$week.id)

y.sight.marked[marked.indices] <- cams.tagged.marked$sum.n.det

# sex for marked (1F, 2M)
G.marked <- unique(cams.tagged[ , c("tag.id", "sex.id")]) %>% pull(sex.id)

# sighting history of the observed marked status unmarked individual samples
y.sight.unmarked <- array(0, dim = c(n.unmarked, J2, K2))
unmarked.indices <- cbind(1:n.unmarked,
                          cams.untagged$location.id,
                          cams.untagged$week.id)
y.sight.unmarked[unmarked.indices] <- cams.untagged$n.det

# sex for unmarked
G.unmarked <- matrix(cams.untagged$sex.id, ncol = 1)
G.unmarked[is.na(G.unmarked)] <- 0

# Location matrix for telemetry data
locs <- array(NA, dim = c(n.marked, K2, 2))

# # Remove captures for ind 14 since not captured in snares
# locs.in <- filter(locs.in, tag_id != 14)

locs.indices <- cbind(locs.in$tag_id, # tag_id uses underscore in this dataset...
                      locs.in$week.id)
locs[cbind(locs.indices, 1)] <- as.numeric(locs.in$x) / x.scale
locs[cbind(locs.indices, 2)] <- as.numeric(locs.in$y) / y.scale

# Number of covariates (genotypes)
ncat <- 1
gamma <- IDcovs <- vector("list", ncat)

nlevels <- rep(2, ncat) # number of IDcovs per cat
for (i in 1:ncat) {
  gamma[[i]] <- rep(1 / nlevels[i], nlevels[i])
  IDcovs[[i]] <- 1:nlevels[i]
}
IDlist <- list(ncat = ncat, IDcovs = IDcovs)

# Trap operation vector or matrix for the marking process
# Exposure to capture does not vary by indiviudal or by trap and individual
# tf1 <- NA
tf1 <- effort.in$days
# # Exposure to capture varies by individual or by trap and individual
# tf1 <- matrix(rep(effort.in$days, n.marked),
#               ncol = J1,
#               nrow = n.marked,
#               byrow = TRUE)
# # Append augmented data
# tf1 <- rbind(tf1, matrix(rep(tf1[1, ], M - n.marked),
#                          nrow = M - n.marked,
#                          ncol = J1,
#                          byrow = TRUE))

# Trap operation vector or matrix for the sighting process
# # Exposure to capture does not vary by indiviudal or by trap and individual
# tf2 <- NA
tf2 <- rep(K2, J2)
# Exposure to capture varies by individual or by trap and individual
# make individual x trap effort to account for dead individual
# tf2 <- matrix(K2, ncol = J2, nrow = n.marked)
#
# tf2 <- rbind(tf2, matrix(K2,
#                          ncol = J2,
#                          nrow = M - n.marked))

# Make mark order matrix
markedS <- matrix(1, ncol = K2, nrow = n.marked) # Assume all individuals are marked for all occasions
# markedS[7,1:4]=0 #unmarked samples can be matched to marked individuals before they are marked (0=um, 1=marked)
# markedS[14,1:14]=0
# markedS[15,]=2#dead individuals can't be allocated any latent id samples (2=dead)

data <- vector("list")
data$y_mark <- y.mark
data$y_sight_marked <- y.sight.marked
data$y_sight_unmarked <- y.sight.unmarked
data$G_marked <- as.integer(G.marked)
data$G_unmarked <- G.unmarked
data$IDlist <- IDlist
data$locs <- locs
data$X1 <- as.matrix(X1)
data$X2 <- as.matrix(X2)
data$K1 <- K1
data$K2 <- K2
data$buff <- buff
data$markedS <- markedS
data$n_marked <- n.marked
data$tf1 <- tf1
data$tf2 <- tf2
xlim <- c(min(c(X1$x, X2$x)), max(c(X1$x, X2$x))) + c(-buff, buff)
ylim <- c(min(c(X1$y, X2$y)), max(c(X1$y, X2$y))) + c(-buff, buff)
area <- (xlim[2] - xlim[1]) * (ylim[2] - ylim[1])

mean_N <- c()
SD_N <- c()

input <- list(
  niter = 2400, nburn = 1200, nthin = 5, M = 200, act_center = "no move",
  inits = list(lam0_mark = 0.05, lam0_sight = 0.009, sigma_d = 5,
               sigma_p = 5, s1 = 5, s2 = 5, psi = 0.5,
               gamma = vector("list", 1)),
  obstype = c("bernoulli", "poisson"), nswap = NA,
  proppars = list(lam0_mark = 0.05, lam0_sight = 0.1, sigma_d = 0.02,
                  sigma_p = 0.2, s1 = 0.5, s2 = 0.25, s2t = 0.1),
  max_proppars <- list(lam0_mark = 100, lam0_sight = 100, sigma_d = 100,
                       sigma_p = 100, s1 = 100, s2 = 100, s2t = 100),
  min_proppars <- list(lam0_mark = .001, lam0_sight = .001, sigma_d = .001,
                       sigma_p = .001, s1 = .001, s2 = .001, s2t = .001),
  storeLatent = TRUE, storeGamma = TRUE, IDup = "Gibbs"
)

# # Use vertices instead of buffer
# verts <- rbind(c(430, 4060), c(310, 4060), c(310, 3945), c(345, 3945), c(345, 3910), c(370, 3910), c(370, 3945), c(410, 3945), c(410, 4000), c(430, 4000), c(430, 4060))
# data$vertices <- verts
# area = 60*120+55*100+35*25 # area of vertexed shape

# Run models in parallel
# n_cores <- parallel::detectCores() # Check number of cores available
my_cluster <- makeCluster(6, type = "PSOCK")
#register cluster to be used by %dopar%
doParallel::registerDoParallel(cl = my_cluster)

system.time({
  all_out_i <- foreach::foreach(
    iter=1:3,
    .packages = "magrittr"
  ) %dopar% {

    devtools::load_all()

    # for (iter in reps) {
    if (iter == 1) {
      input$act_center <- "no move"

      out_1 <- mcmc_SMR(data, input)
      # MCMC_SPIM_out_1 <- coda::as.mcmc(do.call(cbind, list(out_1$out)))
      # summary(MCMC_SPIM_out_1)

      all_out_i <- tibble::as_tibble(out_1$out) %>%
        dplyr::mutate(model = as.character(iter))
    }

    ##########################################
    if (iter == 2) {
      input$act_center <- "move"

      out_2 <- mcmc_SMR(data, input)
      # MCMC_SPIM_out_2 <- coda::as.mcmc(do.call(cbind, list(out_2$out)))
      # summary(MCMC_SPIM_out_2)

      all_out_i <- tibble::as_tibble(out_2$out) %>%
        dplyr::mutate(model = as.character(iter))
    }
    ####################################
    if (iter == 3) {
      # sex and no telemetry mobile
      data_no_tele <- data
      data_no_tele$locs <- NA

      input$act_center <- "move"

      out_3 <- mcmc_SMR(data_no_tele, input)
      # MCMC_SPIM_out_3 <- coda::as.mcmc(do.call(cbind, list(out_3$out)))
      # summary(MCMC_SPIM_out_3)

      all_out_i <- tibble::as_tibble(out_3$out) %>%
        dplyr::mutate(model = as.character(iter))
    }
    all_out_i <- all_out_i
  }
}) # Time stop
# Stop cluster
stopCluster(my_cluster)

# Is there a better way to select each element in the list?
all_out <- tibble::tibble()
for (ii in 1:3) {
  all_out <- all_out %>%
    dplyr::bind_rows(
      all_out_i[[ii]] %>%
        dplyr::mutate(model = as.character(ii))
    )
    # dplyr::add_row(model = as.character(ii),
    #                parm = colnames(all_out_i[[ii]]),
    #                means = colMeans(all_out_i[[ii]]),
    #                sds = apply(all_out_i[[ii]], 2, sd)
    # )
}

# all_out_unlist <- all_out_i[[2]]
# plot(all_out_unlist$sigma_p)

all_out <- dplyr::bind_rows(all_out_i) %>%
  tidyr::pivot_longer(cols = -model, names_to = "Parameter", values_to = "Mean")

# sigma_d == sigma
# all_out$parm[all_out$parm == "sigma_d"] <- "sigma"

all_out %>%
  dplyr::filter(Parameter == "N") %>%
  ggplot2::ggplot(ggplot2::aes(x = model, y = Mean, fill = model, group = model)) +
  ggplot2::geom_boxplot() +
  ggplot2::labs(y = "Abundance Estimates", "Model") +
  ggplot2::theme_minimal()

all_out %>%
  dplyr::filter(Parameter == "N") %>%
  dplyr::mutate(Mean = Mean * 100 / area) %>%
  ggplot2::ggplot(ggplot2::aes(x = model, y = Mean, fill = model, group = model)) +
  ggplot2::geom_boxplot() +
  labs(y = "Density Estimates", "Model") +
  theme_minimal()

all_out %>%
  dplyr::filter(Parameter %in% c("lam0_mark", "lam0_sight")) %>%
  ggplot2::ggplot(ggplot2::aes(x = Parameter, y = Mean, fill = model)) +
  ggplot2::geom_boxplot(position = ggplot2::position_dodge(width = 1)) +
  labs(y = "Mean Parameter Outputs", x = "Parameter") +
  theme_minimal()

all_out %>%
  dplyr::filter(Parameter %in% c("sigma_d", "sigma_p")) %>%
  ggplot2::ggplot(ggplot2::aes(x = Parameter, y = Mean, fill = model)) +
  ggplot2::geom_boxplot(position = ggplot2::position_dodge(width = 1)) +
  labs(y = "Mean Parameter Outputs", x = "Parameter") +
  theme_minimal()


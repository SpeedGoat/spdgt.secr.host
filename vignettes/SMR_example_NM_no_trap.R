# Notes:
# Marked individuals may be captured with snares or hounds and may be used with telemetry data. Relevant data parms are X1.in (snare locations), snares.in (individual tags and tag ids, capture history for each snare). Effort is recorded in effort.in
# Marked ("tagged") and unmarked ("unmarked") individuals are captured at cameras at camera locations, X.in. Capture histories are recorded in cams.in with number captured (n.det), sex, marked (also tag.class), and gps (unused). Effort is extrapolated from the number of unique dates in cams.in

# No snares were used for this run. Conventional SMR is used (no trap markings)
devtools::load_all()

library(dplyr)
library(doParallel)

# Download Google Drive files to temp folder
file_names <- tibble::tibble(
  names = c("resight_xy.csv",
            "resight_obs.csv",
            "gps.csv"),
  ID = c("1xOVJzs5OUue5KVkJWcg_foJX6I_d72K8",
         "1-PFNi323MPBxcu8Bl_Le65ZtNfLAdxii",
         "15MXOs5ZHMxE5F9dPdGa6p-MF45Siglon")
)

# Create and download temp files for each csv file
tmpfl <- c()
for (ff in 1:dim(file_names)[1]) {
  tmpfl[ff] <- tempfile(fileext = ".csv")

  # Download with ID
  googledrive::drive_download(
    googledrive::as_id(as.character(file_names[ff, 2])),
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
activity_centers <-
  tibble(model = character(),
         zout = double(),
         sxout = double(),
         syout = double())

# MCMC Initializations
storeLatent <- TRUE
storeGamma <- TRUE
# niter <- 30000
# nburn <- 20000
# nthin <- 1
niter <- 2000
nburn <- 1000
nthin <- 1

# Collect data if on Google Drive
# Trap locations for marking (X1; snares) and sighting (X; cameras) processes (required: x, y)
X.in <- readr::read_csv(tmpfl[1], show_col_types = FALSE)
# Camera captures (required: week.id, tag.class, tag.id, location.id, n.det, sex.id)
cams.in <- readr::read_csv(tmpfl[2], show_col_types = FALSE)
# Telemetry data (required: tag.id, week.id, x, y)
locs.in <- readr::read_csv(tmpfl[3], show_col_types = FALSE) |>
  dplyr::mutate(
    x = as.numeric(x),
    y = as.numeric(y)
  ) |>
  dplyr::rename(tag.id = tag_id)

# Redefine IDs to be sequential
X.in <- X.in |>
  dplyr::mutate(loc.id.seq = 1:n())
week.id.seq <- cams.in |>
  dplyr::select(week.id) |>
  dplyr::distinct() |>
  dplyr::mutate(week.id.seq = 1:n())

locs.in <- locs.in |>
  dplyr::full_join(week.id.seq, by = "week.id")

cams <- as.data.frame(cams.in |>
                        dplyr::full_join(X.in |>
                                           dplyr::select(location.id, loc.id.seq),
                                         by = "location.id") |>
                        dplyr::full_join(week.id.seq, by = "week.id"))

# number of sighting occasions
K <- length(unique(cams$week.id))

# Define tagged, untagged cams
cams.tagged <- cams |>
  dplyr::filter(!is.na(tag))
cams.untagged <- cams |>
  dplyr::filter(is.na(tag))
n.marked <- length(unique(cams.tagged$tag)) # number of individuals captured and marked
n.unmarked <- nrow(cams.untagged)

# Trap locations for the marking process (X1; snares) and sighting process (X; cameras)
x.scale <- 1000
y.scale <- 1000

# Omit missing data
X <- X.in |>
  summarise(x = x/x.scale,
            y = y/y.scale)
# Number of cameras
J <- nrow(X)

# # Buffer is 1.5 times the average home range "diameter"
buff <- locs.in |>
  group_by(tag) |>
  summarise(mx = mean(c((max(x) - min(x))/x.scale,
                        (max(y) - min(y))/y.scale))/2) |>
  ungroup() |>
  summarise(means = mean(mx)*1.5) |>
  pull(means)
# buff <- 100 # Distance to buffer to create the state space

# Augmented data. Note: large M could result in long run times
M <- n.marked*10
# M <- 160

# Sighting history of the marked, observed marked status, and individually identified samples
y.sight.marked = array(0, dim = c(n.marked, J, K))

########################
# DO ALL IN ONE PIPE??
cams.tagged.marked <- cams.tagged |>
  filter(n.det>0) |>
  group_by(tag.id, loc.id.seq, week.id.seq) |>
  summarise(sum.n.det = sum(n.det),
            .groups = 'drop')

marked.indices <- cbind(cams.tagged.marked$tag.id,
                        cams.tagged.marked$loc.id.seq,
                        cams.tagged.marked$week.id.seq)

y.sight.marked[marked.indices] <- cams.tagged.marked$sum.n.det
########################

# sex for marked (1F, 2M)
G.marked <- unique(cams.tagged[ , c("tag.id", "sex.id")]) |>
  pull(sex.id)

# sighting history of the observed marked status unmarked individual samples
y.sight.unmarked <- array(0, dim = c(n.unmarked, J, K))
unmarked.indices <- cbind(1:n.unmarked,
                          cams.untagged$loc.id.seq,
                          cams.untagged$week.id.seq)
y.sight.unmarked[unmarked.indices] <- cams.untagged$n.det

# sex for unmarked
G.unmarked <- matrix(cams.untagged$sex.id, ncol = 1)
G.unmarked[is.na(G.unmarked)] <- 0

# Location matrix for telemetry data
locs <- array(NA, dim = c(n.marked, K, 2))

locs.indices <- cbind(locs.in$tag.id,
                      locs.in$week.id.seq)
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

# Trap operation vector or matrix for the sighting process
tf <- rep(K, J)

data <- vector("list")
data$n_marked <- n.marked
data$y_sight_marked <- y.sight.marked
data$y_sight_unmarked <- y.sight.unmarked
data$G_marked <- as.integer(G.marked)
data$G_unmarked <- G.unmarked
data$IDlist <- IDlist
data$locs <- locs
data$X2 <- as.matrix(X)
data$K2 <- K
data$buff <- buff
data$tf2 <- tf
xlim <- c(min(X$x), max(X$x)) + c(-buff, buff)
ylim <- c(min(X$y), max(X$y)) + c(-buff, buff)
area <- (xlim[2] - xlim[1]) * (ylim[2] - ylim[1])

mean_N <- c()
SD_N <- c()

input <- list(
  niter = niter, nburn = nburn, nthin = 1, M = M, act_center = F,
  inits = list(lam0_mark = NA, lam0_sight = 0.05, sigma_d = 5,
               sigma_p = NA, s1 = NA, s2 = 5, psi = 0.5,
               gamma = gamma),
  obstype = c("bernoulli", "poisson"), nswap = NA,
  proppars = list(lam0_mark = 0.05, lam0_sight = 0.1, sigma_d = 0.02,
                  sigma_p = 0.2, s1 = 0.5, s2 = 0.25, s2t = 0.1),
  max_proppars <- list(lam0_mark = 100, lam0_sight = 100, sigma_d = 100,
                       sigma_p = 100, s1 = 100, s2 = 100, s2t = 100),
  min_proppars <- list(lam0_mark = .001, lam0_sight = .001, sigma_d = .001,
                       sigma_p = .001, s1 = .001, s2 = .001, s2t = .001),
  storeLatent = TRUE, storeGamma = TRUE, IDup = "Gibbs",
  model_choices = list(1:4)
)

# # Use vertices instead of buffer
# verts <- rbind(c(430, 4060), c(310, 4060), c(310, 3945), c(345, 3945), c(345, 3910), c(370, 3910), c(370, 3945), c(410, 3945), c(410, 4000), c(430, 4000), c(430, 4060))
# data$vertices <- verts
# area = 60*120+55*100+35*25 # area of vertexed shape

# Run models in parallel
final_out <- SMR_wrapper(data, input)

final_out %>%
  dplyr::filter(Parameter == "N") %>%
  ggplot2::ggplot(ggplot2::aes(x = model, y = Mean, fill = model, group = model)) +
  ggplot2::geom_boxplot() +
  ggplot2::labs(y = "Abundance Estimates", "Model") +
  ggplot2::theme_minimal()

final_out %>%
  dplyr::filter(Parameter == "N") %>%
  dplyr::mutate(Mean = Mean * 100 / area) %>%
  ggplot2::ggplot(ggplot2::aes(x = model, y = Mean, fill = model, group = model)) +
  ggplot2::geom_boxplot() +
  ggplot2::labs(y = "Density Estimates", "Model") +
  ggplot2::theme_minimal()

final_out %>%
  dplyr::filter(Parameter %in% c("lam0_mark", "lam0_sight")) %>%
  ggplot2::ggplot(ggplot2::aes(x = Parameter, y = Mean, fill = model)) +
  ggplot2::geom_boxplot(position = ggplot2::position_dodge(width = 1)) +
  ggplot2::labs(y = "Mean Parameter Outputs", x = "Parameter") +
  ggplot2::theme_minimal()

final_out %>%
  dplyr::filter(Parameter %in% c("sigma_d", "sigma_p")) %>%
  ggplot2::ggplot(ggplot2::aes(x = Parameter, y = Mean, fill = model)) +
  ggplot2::geom_boxplot(position = ggplot2::position_dodge(width = 1)) +
  ggplot2::labs(y = "Mean Parameter Outputs", x = "Parameter") +
  ggplot2::theme_minimal()


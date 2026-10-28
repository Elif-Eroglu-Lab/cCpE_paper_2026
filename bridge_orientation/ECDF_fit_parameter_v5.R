# Christoph Jüschke
# IMBA - Institute of Molecular Biotechnology. All rights reserved.
# Edited for bridge orientation use case by Sabrina Syamsul
#
# Analysis of spindle orientation angles
#
# This script:
# - Reads Excel input
# - Fits Vehicle and cCPE separately
# - Plots experimental ECDFs
# - Plots fitted model curves
# - Adds 95% bootstrap confidence bands
# - Uses 1-degree resolution for model curves / CI
# - Plots theoretical isotropic orientation distribution
# - Tests Vehicle vs theoretical distribution
# - Tests cCPE vs theoretical distribution
# - Saves lambda values and p-values to CSV
# - Saves main graph without legend/statistical text
# - Saves legend as a separate PNG and PDF
#
# Excel columns:
# treatment | individual | bridge | alpha | alphaSD | normSD
#
# Run from Terminal:
#
# R --slave --args "vers450angles_Raw_angles.xlsx" < fit_parameter_v5.R
#

set.seed(1)

# ============================================================
# PACKAGE
# ============================================================

library(readxl)

# ============================================================
# READ COMMAND-LINE ARGUMENT
# ============================================================

cmd_args <- commandArgs(T)

if (length(cmd_args) < 1) {
  stop(
    "Please provide the Excel file as a command-line argument."
  )
}

file <- cmd_args[1]

# ============================================================
# LOAD EXCEL DATA
# ============================================================

res <- read_excel(file)

# ============================================================
# CONVERT RELEVANT COLUMNS
# ============================================================

res$treatment <- as.character(res$treatment)

res$alpha <- as.numeric(
  as.character(res$alpha)
)

# ============================================================
# CHECK DATA
# ============================================================

cat("\nData structure:\n")
print(str(res))

cat("\nTreatments found:\n")
print(unique(res$treatment))

cat("\nNumber of observations per treatment:\n")
print(table(res$treatment))

cat("\nFirst 10 alpha values:\n")
print(head(res$alpha, 10))

# ============================================================
# PARAMETERS
# ============================================================

# Initial values for sigma_v and sigma_h
sigma_initial <- c(1, 1)

# Number of simulations for the main fitted model
n <- 50000

# Number of bootstrap repetitions
n_boot <- 500

# Number of simulations used during bootstrap fitting
n_boot_sim <- 10000

# Number of simulations used when creating
# the final bootstrap confidence curves
n_ci_sim <- 20000

# ============================================================
# COLOURS
# ============================================================

vehicle_color <- "#3373CC"

ccpe_color <- "#F2BF1A"

reference_color <- "#666666"

# ============================================================
# SIMULATE BASE NORMAL VECTORS
# ============================================================
#
# These vectors are used repeatedly during model fitting.
#
# Using fixed random vectors makes the optimization more stable
# because the objective function does not change due to new
# random numbers at every iteration.
#

C_vb <- array(
  data = c(
    rnorm(n),
    rnorm(n),
    rnorm(n)
  ),
  dim = c(
    n,
    3
  ),
  dimnames = list(
    cell = 1:n,
    coord = c(
      "x",
      "y",
      "z"
    )
  )
)

# ============================================================
# REFERENCE VECTOR POINTING UP
# ============================================================

P_n <- array(
  rep(
    c(0, 0, 1),
    n
  ),
  dim = c(
    3,
    n
  ),
  dimnames = list(
    coord = c(
      "X",
      "Y",
      "Z"
    ),
    cell = 1:n
  )
)

# ============================================================
# FUNCTION TO SIMULATE MODEL ANGLES
# ============================================================

simulate_angles <- function(
    sigma,
    n_sim
) {
  
  C_v <- array(
    data = c(
      rnorm(
        n_sim,
        0,
        sigma[1]
      ),
      rnorm(
        n_sim,
        0,
        sigma[1]
      ),
      rnorm(
        n_sim,
        0,
        sigma[2]
      )
    ),
    dim = c(
      n_sim,
      3
    )
  )
  
  # ----------------------------------------------------------
  # Normalize vectors
  # ----------------------------------------------------------
  
  C_n <- C_v / sqrt(
    apply(
      C_v,
      1,
      crossprod
    )
  )
  
  # ----------------------------------------------------------
  # Reference vector pointing along z-axis
  # ----------------------------------------------------------
  
  P <- array(
    rep(
      c(0, 0, 1),
      n_sim
    ),
    dim = c(
      3,
      n_sim
    )
  )
  
  # ----------------------------------------------------------
  # Calculate angle relative to vertical
  # ----------------------------------------------------------
  
  alp <- abs(
    asin(
      colSums(
        P * t(C_n)
      )
    )
  )
  
  return(alp)
}

# ============================================================
# FUNCTION TO FIT ONE TREATMENT
# ============================================================

fit_treatment <- function(
    treatment,
    res,
    C_vb,
    P_n,
    n
) {
  
  cat("\n")
  cat("============================================================\n")
  cat(
    "Fitting treatment:",
    treatment,
    "\n"
  )
  cat("============================================================\n")
  
  # ----------------------------------------------------------
  # Extract alpha values
  # ----------------------------------------------------------
  
  alpha_data <- res$alpha[
    res$treatment == treatment
  ]
  
  alpha_data <- as.numeric(
    alpha_data
  )
  
  alpha_data <- alpha_data[
    !is.na(alpha_data)
  ]
  
  if (
    length(alpha_data) == 0
  ) {
    
    stop(
      paste(
        "No valid alpha values found for treatment:",
        treatment
      )
    )
  }
  
  cat(
    "Number of observations:",
    length(alpha_data),
    "\n"
  )
  
  # ----------------------------------------------------------
  # Convert degrees to radians
  # ----------------------------------------------------------
  
  data <- alpha_data / 180 * pi
  
  # ----------------------------------------------------------
  # Experimental ECDF
  # ----------------------------------------------------------
  
  edat <- ecdf(
    data
  )
  
  dat <- edat(
    data
  )
  
  # ==========================================================
  # FIT FUNCTION
  # ==========================================================
  
  fitfun <- function(sigma) {
    
    C_v <- C_vb
    
    C_v[, "x"] <- C_vb[, "x"] * sigma[1]
    
    C_v[, "y"] <- C_vb[, "y"] * sigma[1]
    
    C_v[, "z"] <- C_vb[, "z"] * sigma[2]
    
    # --------------------------------------------------------
    # Normalize vectors
    # --------------------------------------------------------
    
    C_n <- C_v / sqrt(
      apply(
        C_v,
        1,
        crossprod
      )
    )
    
    # --------------------------------------------------------
    # Model angular distribution
    # --------------------------------------------------------
    
    alp5 <- ecdf(
      abs(
        asin(
          colSums(
            P_n * t(C_n)
          )
        )
      )
    )
    
    # --------------------------------------------------------
    # Difference between model and experimental ECDF
    # --------------------------------------------------------
    
    sum(
      (
        alp5(data) - dat
      )^2
    )
  }
  
  # ==========================================================
  # OPTIMIZE SIGMA
  # ==========================================================
  
  fit <- optim(
    sigma_initial,
    fitfun,
    method = "L-BFGS-B",
    lower = c(
      0,
      0
    ),
    upper = c(
      1,
      1
    ),
    control = list(
      trace = 0
    )
  )
  
  sigma <- fit$par
  
  # ============================================================
  # RESIMULATE USING OPTIMIZED SIGMAS
  # ============================================================
  
  C_v <- array(
    data = c(
      rnorm(
        n,
        0,
        sigma[1]
      ),
      rnorm(
        n,
        0,
        sigma[1]
      ),
      rnorm(
        n,
        0,
        sigma[2]
      )
    ),
    dim = c(
      n,
      3
    ),
    dimnames = list(
      cell = 1:n,
      coord = c(
        "x",
        "y",
        "z"
      )
    )
  )
  
  # ----------------------------------------------------------
  # Normalize
  # ----------------------------------------------------------
  
  C_n <- C_v / sqrt(
    apply(
      C_v,
      1,
      crossprod
    )
  )
  
  # ----------------------------------------------------------
  # Calculate fitted angles
  # ----------------------------------------------------------
  
  alp5 <- abs(
    asin(
      colSums(
        P_n * t(C_n)
      )
    )
  )
  
  # ============================================================
  # KS TEST:
  # DATA VS FITTED MODEL
  # ============================================================
  
  ks_fit <- ks.test(
    alp5,
    data
  )
  
  # ============================================================
  # LAMBDA VALUES
  # ============================================================
  
  lambda_v <- -log(
    sigma[1]
  )
  
  lambda_h <- -log(
    sigma[2]
  )
  
  # ============================================================
  # PRINT RESULTS
  # ============================================================
  
  cat("\n")
  cat(
    "Results for:",
    treatment,
    "\n"
  )
  cat("--------------------------------------------\n")
  
  cat(
    "Number of observations:",
    length(data),
    "\n"
  )
  
  cat(
    "sigma_v:",
    sigma[1],
    "\n"
  )
  
  cat(
    "sigma_h:",
    sigma[2],
    "\n"
  )
  
  cat(
    "lambda_v:",
    round(
      lambda_v,
      4
    ),
    "\n"
  )
  
  cat(
    "lambda_h:",
    round(
      lambda_h,
      4
    ),
    "\n"
  )
  
  cat(
    "KS D (data vs fitted model):",
    round(
      ks_fit$statistic,
      4
    ),
    "\n"
  )
  
  cat(
    "KS p-value (data vs fitted model):",
    ks_fit$p.value,
    "\n"
  )
  
  # ============================================================
  # RETURN RESULTS
  # ============================================================
  
  return(
    list(
      treatment = treatment,
      alpha_data = alpha_data,
      data = data,
      sigma = sigma,
      lambda_v = lambda_v,
      lambda_h = lambda_h,
      alp5 = alp5,
      ks_fit = ks_fit
    )
  )
}

# ============================================================
# FIT VEHICLE
# ============================================================

vehicle_fit <- fit_treatment(
  "Vehicle",
  res,
  C_vb,
  P_n,
  n
)

# ============================================================
# FIT cCPE
# ============================================================

ccpe_fit <- fit_treatment(
  "cCPE",
  res,
  C_vb,
  P_n,
  n
)

# ============================================================
# THEORETICAL DISTRIBUTION KS TESTS
#
# For an unbiased orientation in 3D:
#
# F(alpha) = sin(alpha)
#
# where alpha is in radians.
# ============================================================

vehicle_angles <- (
  vehicle_fit$alpha_data /
    180 *
    pi
)

ccpe_angles <- (
  ccpe_fit$alpha_data /
    180 *
    pi
)

# ------------------------------------------------------------
# Vehicle vs theoretical distribution
# ------------------------------------------------------------

ks_vehicle_random <- ks.test(
  vehicle_angles,
  function(x) sin(x)
)

# ------------------------------------------------------------
# cCPE vs theoretical distribution
# ------------------------------------------------------------

ks_ccpe_random <- ks.test(
  ccpe_angles,
  function(x) sin(x)
)

# ============================================================
# PRINT THEORETICAL DISTRIBUTION RESULTS
# ============================================================

cat("\n")
cat("============================================================\n")
cat("THEORETICAL DISTRIBUTION KS TESTS\n")
cat("============================================================\n")

cat(
  "\nVehicle vs theoretical orientation distribution:\n"
)

cat(
  "KS D =",
  round(
    ks_vehicle_random$statistic,
    4
  ),
  "\n"
)

cat(
  "KS p-value =",
  ks_vehicle_random$p.value,
  "\n"
)

cat(
  "\ncCPE vs theoretical orientation distribution:\n"
)

cat(
  "KS D =",
  round(
    ks_ccpe_random$statistic,
    4
  ),
  "\n"
)

cat(
  "KS p-value =",
  ks_ccpe_random$p.value,
  "\n"
)

# ============================================================
# BOOTSTRAP FUNCTION
#
# Resamples the experimental observations,
# refits the model,
# and stores sigma values.
#
# The same simulated random vectors are used for every
# bootstrap fit to reduce Monte-Carlo variation during
# optimization.
# ============================================================

bootstrap_fit <- function(
    alpha_data,
    n_boot,
    n_sim
) {
  
  boot_sigma <- matrix(
    NA,
    nrow = n_boot,
    ncol = 2
  )
  
  colnames(boot_sigma) <- c(
    "sigma_v",
    "sigma_h"
  )
  
  # ----------------------------------------------------------
  # Generate fixed random vectors for bootstrap fitting
  # ----------------------------------------------------------
  
  C_boot <- array(
    data = c(
      rnorm(n_sim),
      rnorm(n_sim),
      rnorm(n_sim)
    ),
    dim = c(
      n_sim,
      3
    )
  )
  
  P_boot <- array(
    rep(
      c(0, 0, 1),
      n_sim
    ),
    dim = c(
      3,
      n_sim
    )
  )
  
  # ----------------------------------------------------------
  # Bootstrap loop
  # ----------------------------------------------------------
  
  for (b in 1:n_boot) {
    
    if (
      b %% 25 == 0
    ) {
      
      cat(
        "Bootstrap:",
        b,
        "/",
        n_boot,
        "\n"
      )
    }
    
    # --------------------------------------------------------
    # Resample observed angles
    # --------------------------------------------------------
    
    boot_alpha <- sample(
      alpha_data,
      size = length(alpha_data),
      replace = TRUE
    )
    
    boot_data <- (
      boot_alpha /
        180 *
        pi
    )
    
    # --------------------------------------------------------
    # Experimental ECDF
    # --------------------------------------------------------
    
    boot_ecdf <- ecdf(
      boot_data
    )
    
    boot_dat <- boot_ecdf(
      boot_data
    )
    
    # ========================================================
    # Bootstrap objective function
    # ========================================================
    
    boot_fitfun <- function(sigma) {
      
      C_v <- C_boot
      
      C_v[, 1] <- (
        C_boot[, 1] *
          sigma[1]
      )
      
      C_v[, 2] <- (
        C_boot[, 2] *
          sigma[1]
      )
      
      C_v[, 3] <- (
        C_boot[, 3] *
          sigma[2]
      )
      
      # ------------------------------------------------------
      # Normalize vectors
      # ------------------------------------------------------
      
      C_n <- C_v / sqrt(
        apply(
          C_v,
          1,
          crossprod
        )
      )
      
      # ------------------------------------------------------
      # Model ECDF
      # ------------------------------------------------------
      
      model_ecdf <- ecdf(
        abs(
          asin(
            colSums(
              P_boot * t(C_n)
            )
          )
        )
      )
      
      # ------------------------------------------------------
      # Difference between bootstrap data and model
      # ------------------------------------------------------
      
      sum(
        (
          model_ecdf(
            boot_data
          ) -
            boot_dat
        )^2
      )
    }
    
    # ========================================================
    # Optimize bootstrap dataset
    # ========================================================
    
    boot_fit <- try(
      optim(
        sigma_initial,
        boot_fitfun,
        method = "L-BFGS-B",
        lower = c(
          0,
          0
        ),
        upper = c(
          1,
          1
        )
      ),
      silent = TRUE
    )
    
    # --------------------------------------------------------
    # Save successful result
    # --------------------------------------------------------
    
    if (
      !inherits(
        boot_fit,
        "try-error"
      )
    ) {
      
      if (
        all(
          is.finite(
            boot_fit$par
          )
        )
      ) {
        
        boot_sigma[b, ] <- boot_fit$par
      }
    }
  }
  
  return(
    boot_sigma
  )
}

# ============================================================
# BOOTSTRAP VEHICLE
# ============================================================

cat("\n")
cat("============================================================\n")
cat("BOOTSTRAPPING VEHICLE MODEL\n")
cat("============================================================\n")

vehicle_boot <- bootstrap_fit(
  vehicle_fit$alpha_data,
  n_boot,
  n_boot_sim
)

# ============================================================
# BOOTSTRAP cCPE
# ============================================================

cat("\n")
cat("============================================================\n")
cat("BOOTSTRAPPING cCPE MODEL\n")
cat("============================================================\n")

ccpe_boot <- bootstrap_fit(
  ccpe_fit$alpha_data,
  n_boot,
  n_boot_sim
)

# ============================================================
# CREATE BOOTSTRAP MODEL CURVES
#
# The confidence interval is calculated around the fitted
# MODEL distribution.
#
# Each bootstrap:
# 1. Takes a bootstrap-fitted sigma
# 2. Generates a model distribution
# 3. Calculates its ECDF
# 4. Evaluates the ECDF at every degree
#
# The 2.5th and 97.5th percentiles give a pointwise
# 95% confidence interval.
# ============================================================

create_boot_curves <- function(
    boot_sigma,
    x,
    n_sim = 20000
) {
  
  # ----------------------------------------------------------
  # Remove failed bootstrap fits
  # ----------------------------------------------------------
  
  original_n <- nrow(
    boot_sigma
  )
  
  boot_sigma <- boot_sigma[
    complete.cases(
      boot_sigma
    ),
    ,
    drop = FALSE
  ]
  
  cat(
    "Successful bootstrap fits:",
    nrow(boot_sigma),
    "of",
    original_n,
    "\n"
  )
  
  if (
    nrow(boot_sigma) < 10
  ) {
    
    stop(
      "Too few successful bootstrap fits to calculate a confidence band."
    )
  }
  
  # ----------------------------------------------------------
  # Matrix to store bootstrap curves
  # ----------------------------------------------------------
  
  boot_curves <- matrix(
    NA,
    nrow = nrow(boot_sigma),
    ncol = length(x)
  )
  
  # ==========================================================
  # Generate ONE fixed set of random vectors
  #
  # This means differences between bootstrap curves mainly
  # reflect differences in the fitted sigma values rather
  # than new Monte-Carlo noise.
  # ==========================================================
  
  C_boot_curve <- array(
    data = c(
      rnorm(n_sim),
      rnorm(n_sim),
      rnorm(n_sim)
    ),
    dim = c(
      n_sim,
      3
    )
  )
  
  P_boot_curve <- array(
    rep(
      c(0, 0, 1),
      n_sim
    ),
    dim = c(
      3,
      n_sim
    )
  )
  
  # ==========================================================
  # Calculate model curve for every bootstrap fit
  # ==========================================================
  
  for (
    i in 1:nrow(boot_sigma)
  ) {
    
    sigma <- boot_sigma[i, ]
    
    # --------------------------------------------------------
    # Apply bootstrap sigma values
    # --------------------------------------------------------
    
    C_v <- C_boot_curve
    
    C_v[, 1] <- (
      C_boot_curve[, 1] *
        sigma[1]
    )
    
    C_v[, 2] <- (
      C_boot_curve[, 2] *
        sigma[1]
    )
    
    C_v[, 3] <- (
      C_boot_curve[, 3] *
        sigma[2]
    )
    
    # --------------------------------------------------------
    # Normalize vectors
    # --------------------------------------------------------
    
    C_n <- C_v / sqrt(
      apply(
        C_v,
        1,
        crossprod
      )
    )
    
    # --------------------------------------------------------
    # Calculate model angles
    # --------------------------------------------------------
    
    sim <- abs(
      asin(
        colSums(
          P_boot_curve * t(C_n)
        )
      )
    )
    
    sim_deg <- (
      sim *
        180 /
        pi
    )
    
    # --------------------------------------------------------
    # Model ECDF
    # --------------------------------------------------------
    
    model_ecdf <- ecdf(
      sim_deg
    )
    
    # --------------------------------------------------------
    # Evaluate model ECDF at every degree
    # --------------------------------------------------------
    
    boot_curves[i, ] <- model_ecdf(
      x
    )
  }
  
  # ==========================================================
  # Pointwise 95% confidence interval
  # ==========================================================
  
  lower <- apply(
    boot_curves,
    2,
    quantile,
    probs = 0.025,
    na.rm = TRUE
  )
  
  upper <- apply(
    boot_curves,
    2,
    quantile,
    probs = 0.975,
    na.rm = TRUE
  )
  
  # ==========================================================
  # Return results
  # ==========================================================
  
  return(
    list(
      lower = as.numeric(
        lower
      ),
      upper = as.numeric(
        upper
      )
    )
  )
}

# ============================================================
# ANGLE GRID
#
# One point for every degree from 0 to 90.
# ============================================================

x <- 0:90

# ============================================================
# CREATE VEHICLE CONFIDENCE BAND
# ============================================================

cat("\n")
cat(
  "Creating Vehicle 95% confidence band...\n"
)

vehicle_ci <- create_boot_curves(
  vehicle_boot,
  x,
  n_sim = n_ci_sim
)

# ============================================================
# CREATE cCPE CONFIDENCE BAND
# ============================================================

cat(
  "Creating cCPE 95% confidence band...\n"
)

ccpe_ci <- create_boot_curves(
  ccpe_boot,
  x,
  n_sim = n_ci_sim
)

# ============================================================
# FITTED MODEL ECDFS
# ============================================================

vehicle_model_ecdf <- ecdf(
  vehicle_fit$alp5 *
    180 /
    pi
)

ccpe_model_ecdf <- ecdf(
  ccpe_fit$alp5 *
    180 /
    pi
)

# ============================================================
# EXPERIMENTAL ECDFS
# ============================================================

vehicle_ecdf <- ecdf(
  vehicle_fit$alpha_data
)

ccpe_ecdf <- ecdf(
  ccpe_fit$alpha_data
)

x_vehicle <- sort(
  unique(
    c(
      0,
      vehicle_fit$alpha_data,
      90
    )
  )
)

x_ccpe <- sort(
  unique(
    c(
      0,
      ccpe_fit$alpha_data,
      90
    )
  )
)

# ============================================================
# THEORETICAL ISOTROPIC ORIENTATION
#
# F(alpha) = sin(alpha)
# ============================================================

random_curve <- sin(
  x /
    180 *
    pi
)

# ============================================================
# PLOT SETTINGS
# ============================================================
#
# These settings are based on v4.
#
# cex.axis = 2
#   Axis tick labels are larger.
#
# cex.lab = 2.5
#   Axis titles are larger.
#
# font.axis = 1
#   Tick labels are regular rather than bold.
#
# font.lab = 1
#   Axis titles are regular rather than bold.
#
# mgp = c(5.5, 1.8, 0)
#   Gives substantially more space between axis titles,
#   tick labels and axes.
#
# tcl = -1.0
#   Longer tick marks.
#
# ============================================================

plot_margins <- c(
  8,
  8,
  2,
  2
)

plot_mgp <- c(
  5.5,
  1.8,
  0
)

# ============================================================
# FUNCTION TO DRAW THE MAIN PLOT
#
# Using one plotting function guarantees that the PNG and
# PDF versions are visually identical.
# ============================================================

draw_main_plot <- function() {
  
  # ==========================================================
  # PLOT SETTINGS
  # ==========================================================
  
  par(
    mar = plot_margins,
    las = 1,
    mgp = plot_mgp,
    cex.axis = 2,
    cex.lab = 2.5,
    font.axis = 1,
    font.lab = 1,
    tcl = -1.0
  )
  
  # ==========================================================
  # EMPTY PLOT
  #
  # plotmath is used so that alpha and degree render correctly
  # in both PNG and PDF.
  # ==========================================================
  
  plot(
    NA,
    NA,
    xlim = c(
      0,
      90
    ),
    ylim = c(
      0,
      1
    ),
    xlab = expression(
      "Angle " *
        alpha *
        " (" *
        degree *
        ")"
    ),
    ylab = "Fraction (%)",
    main = "",
    xaxt = "n",
    yaxt = "n",
    font.axis = 1,
    font.lab = 1,
    font.main = 1
  )
  
  # ==========================================================
  # AXES
  # ==========================================================
  
  axis(
    1,
    at = seq(
      0,
      90,
      by = 15
    ),
    lwd = 5,
    lwd.ticks = 5
  )
  
  axis(
    2,
    at = seq(
      0,
      1,
      by = 0.2
    ),
    labels = seq(
      0,
      100,
      by = 20
    ),
    lwd = 5,
    lwd.ticks = 5
  )
  
  box(
    lwd = 5
  )
  
  # ==========================================================
  # THEORETICAL ISOTROPIC ORIENTATION
  #
  # F(alpha) = sin(alpha)
  # ==========================================================
  
  lines(
    x,
    random_curve,
    col = reference_color,
    lwd = 5,
    lty = 2
  )
  
  # ==========================================================
  # VEHICLE 95% BOOTSTRAP CONFIDENCE BAND
  #
  # The band is deliberately drawn before the experimental
  # ECDF and fitted curve so those remain visually prominent.
  # ==========================================================
  
  polygon(
    c(
      x,
      rev(x)
    ),
    c(
      vehicle_ci$lower,
      rev(
        vehicle_ci$upper
      )
    ),
    col = adjustcolor(
      vehicle_color,
      alpha.f = 0.18
    ),
    border = NA
  )
  
  # ==========================================================
  # cCPE 95% BOOTSTRAP CONFIDENCE BAND
  # ==========================================================
  
  polygon(
    c(
      x,
      rev(x)
    ),
    c(
      ccpe_ci$lower,
      rev(
        ccpe_ci$upper
      )
    ),
    col = adjustcolor(
      ccpe_color,
      alpha.f = 0.18
    ),
    border = NA
  )
  
  # ==========================================================
  # VEHICLE EXPERIMENTAL ECDF
  #
  # Kept at lwd = 2, as in v4.
  # ==========================================================
  
  lines(
    x_vehicle,
    vehicle_ecdf(
      x_vehicle
    ),
    type = "s",
    col = "#000000",
    lwd = 2
  )
  
  points(
    vehicle_fit$alpha_data,
    vehicle_ecdf(
      vehicle_fit$alpha_data
    ),
    pch = 16,
    col = vehicle_color,
    cex = 1.5
  )
  
  # ==========================================================
  # cCPE EXPERIMENTAL ECDF
  #
  # Kept at lwd = 2, as in v4.
  # ==========================================================
  
  lines(
    x_ccpe,
    ccpe_ecdf(
      x_ccpe
    ),
    type = "s",
    col = "#000000",
    lwd = 2
  )
  
  points(
    ccpe_fit$alpha_data,
    ccpe_ecdf(
      ccpe_fit$alpha_data
    ),
    pch = 16,
    col = ccpe_color,
    cex = 1.5
  )
  
  # ==========================================================
  # VEHICLE FITTED MODEL CURVE
  # ==========================================================
  
  lines(
    x,
    vehicle_model_ecdf(
      x
    ),
    col = vehicle_color,
    lwd = 5
  )
  
  # ==========================================================
  # cCPE FITTED MODEL CURVE
  # ==========================================================
  
  lines(
    x,
    ccpe_model_ecdf(
      x
    ),
    col = ccpe_color,
    lwd = 5
  )
}

# ============================================================
# CREATE MAIN PNG
# ============================================================

png(
  filename = "Vehicle_cCPE_orientation_v5.png",
  width = 4000,
  height = 4000,
  res = 600
)

draw_main_plot()

dev.off()

# ============================================================
# CREATE MAIN PDF
#
# Vector output suitable for Illustrator, Inkscape,
# PowerPoint, publications, etc.
# ============================================================

pdf(
  file = "Vehicle_cCPE_orientation_v5.pdf",
  width = 6.67,
  height = 6.67,
  useDingbats = FALSE
)

draw_main_plot()

dev.off()

# ============================================================
# CREATE SEPARATE LEGEND PNG
# ============================================================

png(
  filename = "Vehicle_cCPE_legend_v5.png",
  width = 4000,
  height = 4000,
  res = 600
)

par(
  mar = c(
    0,
    0,
    0,
    0
  )
)

plot(
  NA,
  NA,
  xlim = c(
    0,
    1
  ),
  ylim = c(
    0,
    1
  ),
  axes = FALSE,
  xlab = "",
  ylab = "",
  bty = "n"
)

legend(
  "center",
  legend = c(
    "Theoretical isotropic orientation",
    "Vehicle",
    "cCpE"
  ),
  col = c(
    reference_color,
    vehicle_color,
    ccpe_color
  ),
  lwd = c(
    5,
    5,
    5
  ),
  lty = c(
    2,
    1,
    1
  ),
  bty = "n",
  cex = 1.5,
  seg.len = 2.5
)

dev.off()

# ============================================================
# CREATE SEPARATE LEGEND PDF
# ============================================================

pdf(
  file = "Vehicle_cCPE_legend_v5.pdf",
  width = 6.67,
  height = 6.67,
  useDingbats = FALSE
)

par(
  mar = c(
    0,
    0,
    0,
    0
  )
)

plot(
  NA,
  NA,
  xlim = c(
    0,
    1
  ),
  ylim = c(
    0,
    1
  ),
  axes = FALSE,
  xlab = "",
  ylab = "",
  bty = "n"
)

legend(
  "center",
  legend = c(
    "Theoretical isotropic orientation",
    "Vehicle",
    "cCpE"
  ),
  col = c(
    reference_color,
    vehicle_color,
    ccpe_color
  ),
  lwd = c(
    5,
    5,
    5
  ),
  lty = c(
    2,
    1,
    1
  ),
  bty = "n",
  cex = 1.5,
  seg.len = 2.5
)

dev.off()

# ============================================================
# SAVE RESULTS TO CSV
# ============================================================
#
# One row per treatment.
#
# Columns:
#
# treatment
# n
# sigma_v
# sigma_h
# lambda_v
# lambda_h
# KS_D_fitted_model
# KS_p_fitted_model
# KS_D_theoretical
# KS_p_theoretical
#
# The two KS tests remain separate so it is clear what each
# p-value represents.
#
# ============================================================

results_csv <- data.frame(
  
  treatment = c(
    "Vehicle",
    "cCpE"
  ),
  
  n = c(
    length(
      vehicle_fit$data
    ),
    length(
      ccpe_fit$data
    )
  ),
  
  sigma_v = c(
    vehicle_fit$sigma[1],
    ccpe_fit$sigma[1]
  ),
  
  sigma_h = c(
    vehicle_fit$sigma[2],
    ccpe_fit$sigma[2]
  ),
  
  lambda_v = c(
    vehicle_fit$lambda_v,
    ccpe_fit$lambda_v
  ),
  
  lambda_h = c(
    vehicle_fit$lambda_h,
    ccpe_fit$lambda_h
  ),
  
  KS_D_fitted_model = c(
    as.numeric(
      vehicle_fit$ks_fit$statistic
    ),
    as.numeric(
      ccpe_fit$ks_fit$statistic
    )
  ),
  
  KS_p_fitted_model = c(
    vehicle_fit$ks_fit$p.value,
    ccpe_fit$ks_fit$p.value
  ),
  
  KS_D_theoretical = c(
    as.numeric(
      ks_vehicle_random$statistic
    ),
    as.numeric(
      ks_ccpe_random$statistic
    )
  ),
  
  KS_p_theoretical = c(
    ks_vehicle_random$p.value,
    ks_ccpe_random$p.value
  ),
  
  stringsAsFactors = FALSE
)

# ============================================================
# WRITE CSV
# ============================================================

write.csv(
  results_csv,
  file = "Vehicle_cCPE_results_v5.csv",
  row.names = FALSE
)

# ============================================================
# FINAL SUMMARY
# ============================================================

cat("\n")
cat("============================================================\n")
cat("FINAL SUMMARY\n")
cat("============================================================\n")

# ------------------------------------------------------------
# Vehicle
# ------------------------------------------------------------

cat("\nVehicle:\n")

cat(
  "n =",
  length(
    vehicle_fit$data
  ),
  "\n"
)

cat(
  "sigma_v =",
  round(
    vehicle_fit$sigma[1],
    4
  ),
  "\n"
)

cat(
  "sigma_h =",
  round(
    vehicle_fit$sigma[2],
    4
  ),
  "\n"
)

cat(
  "lambda_v =",
  round(
    vehicle_fit$lambda_v,
    3
  ),
  "\n"
)

cat(
  "lambda_h =",
  round(
    vehicle_fit$lambda_h,
    3
  ),
  "\n"
)

cat(
  "KS D vs fitted model =",
  round(
    vehicle_fit$ks_fit$statistic,
    3
  ),
  "\n"
)

cat(
  "KS p vs fitted model =",
  vehicle_fit$ks_fit$p.value,
  "\n"
)

cat(
  "KS D vs theoretical distribution =",
  round(
    ks_vehicle_random$statistic,
    3
  ),
  "\n"
)

cat(
  "KS p vs theoretical distribution =",
  ks_vehicle_random$p.value,
  "\n"
)

# ------------------------------------------------------------
# cCPE
# ------------------------------------------------------------

cat("\ncCPE:\n")

cat(
  "n =",
  length(
    ccpe_fit$data
  ),
  "\n"
)

cat(
  "sigma_v =",
  round(
    ccpe_fit$sigma[1],
    4
  ),
  "\n"
)

cat(
  "sigma_h =",
  round(
    ccpe_fit$sigma[2],
    4
  ),
  "\n"
)

cat(
  "lambda_v =",
  round(
    ccpe_fit$lambda_v,
    3
  ),
  "\n"
)

cat(
  "lambda_h =",
  round(
    ccpe_fit$lambda_h,
    3
  ),
  "\n"
)

cat(
  "KS D vs fitted model =",
  round(
    ccpe_fit$ks_fit$statistic,
    3
  ),
  "\n"
)

cat(
  "KS p vs fitted model =",
  ccpe_fit$ks_fit$p.value,
  "\n"
)

cat(
  "KS D vs theoretical distribution =",
  round(
    ks_ccpe_random$statistic,
    3
  ),
  "\n"
)

cat(
  "KS p vs theoretical distribution =",
  ks_ccpe_random$p.value,
  "\n"
)

# ------------------------------------------------------------
# Bootstrap information
# ------------------------------------------------------------

cat("\nBootstrap information:\n")

cat(
  "Requested bootstrap fits =",
  n_boot,
  "\n"
)

cat(
  "Successful Vehicle fits =",
  sum(
    complete.cases(
      vehicle_boot
    )
  ),
  "\n"
)

cat(
  "Successful cCPE fits =",
  sum(
    complete.cases(
      ccpe_boot
    )
  ),
  "\n"
)

# ------------------------------------------------------------
# Output files
# ------------------------------------------------------------

cat("\n")
cat("Output files saved:\n")
cat("--------------------------------------------\n")

cat(
  "Vehicle_cCPE_orientation_v5.png\n"
)

cat(
  "Vehicle_cCPE_orientation_v5.pdf\n"
)

cat(
  "Vehicle_cCPE_legend_v5.png\n"
)

cat(
  "Vehicle_cCPE_legend_v5.pdf\n"
)

cat(
  "Vehicle_cCPE_results_v5.csv\n"
)

cat("============================================================\n")


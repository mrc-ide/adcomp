library(threemc)
library(dplyr)
library(ggplot2)
library(sf)
library(TMB)

compile("threemc.cpp")
dyn.load(dynlib("threemc"))
# remove circumcisions with missing type?
rm_missing_type <- FALSE

#### Metadata to run the models ####
# set country
cntry <- "SWZ"

k_dt <- 5 # Age knot spacing
# start_year <-  2002
# start_year <- 1985 # earlier start year to capture generational change in TMC
if (cntry == "LBR") cens_age <- 29 else cens_age <- 59
N <- 1000
forecast_year <- 2021
# paed_age_cutoff <- 10
paed_age_cutoff <- NULL
rw_order <- NULL # use RW or AR temporal prior?
# inc_time_tmc <- TRUE # include time random effect for TMC
inc_time_tmc <- TRUE


#### Reading in data ####
# Revert to using planar rather than spherical geometry in `sf`
sf::sf_use_s2(FALSE)

# Updated to redefine space here. Previously space variable reset within
# area_level. However, as we will be estimating the rates/coverage across
# multiple admin boundaries (area_level) internally in the model, this needs
# to be reset so that the hazard/integration matrices can distinguish between
# areas/area_levels

root_dir <- here::here("../../circumcision-coverage")
# read in data, filter for specific country and male surveys only
filters <- c("iso3" = cntry, sex = "male")
areas <- read_circ_data(file.path(root_dir, "data/areas.geojson"), filters) %>%
  dplyr::mutate(space = 1:dplyr::n()) # add space column to areas
survey_circumcision <- read_circ_data(
  file.path(root_dir, "data/survey_circumcision.csv.gz"), filters
) %>%
  mutate(survey_year = as.numeric(substr(survey_id, 4, 7)))
populations <- read_circ_data(
  file.path(root_dir,"data/population_singleage_aggr.csv.gz"),
  filters
)

# find survey years, rm max year (and also second largest year, if appropriate)

# survey_years <- sort(unique(survey_circumcision$survey_year))
# max_years <- max(survey_years)
# if (length(survey_years) > 1) {
#     if ((max_years - survey_years[length(survey_years) - 1]) == 1) {
#         max_years <- c(max_years, max_years - 1)
#     }
#     # record removed surveys
#     removed_surveys <- survey_circumcision %>%
#         filter(survey_year %in% max_years) %>%
#         distinct(survey_id) %>%
#         pull()
#     survey_circumcision <- survey_circumcision %>%
#         filter(!survey_year %in% max_years)
# } else {
#     # don't continue if there is only one survey available for a country
#     stop(paste0(
#         "Only one survey available for ",
#         cntry,
#         ", so OOS Validation not possible")
#     )
# }

# re-calculate start as earliest year - 50 (only where TMC needs to vary)
survey_years <- unique(survey_circumcision$survey_year)
if (inc_time_tmc == TRUE) {
  start_year <- first(survey_years) - 50
} else start_year <- 2002

# Fill in any NAs in populations with
min_pop_year <- min(populations$year)
if (start_year < min_pop_year) {
  missing_years <- start_year:(min_pop_year - 1)
  missing_rows <- tidyr::crossing(
    select(populations, -c(year, population)),
    "year"       = missing_years,
    "population" = NA
  )
  populations <- bind_rows(populations, missing_rows) %>%
    arrange(iso3, area_id, area_level, age, year) %>%
    group_by(iso3, area_id, area_level, age) %>%
    tidyr::fill(population, .direction = "downup") %>%
    ungroup()
}

# pull recommended area hierarchy for target country
area_lev <- threemc::datapack_psnu_area_level %>%
  filter(iso3 == cntry) %>%
  pull(psnu_area_level)

# don't model at the country level
if (length(area_lev) > 0 && area_lev == 0) area_lev <- NULL

# if area_level is missing, assume most common area lev in surveys
if (length(area_lev) == 0) {
  area_lev <- table(as.numeric(substr(survey_circumcision$area_id, 5, 5)))
  area_lev <- as.numeric(names(area_lev)[area_lev == max(area_lev)])
}

# area_lev <- 0 # run national level
# area_lev <- 1

# save location
# save_loc <- paste0("Runs/", cntry, "/")
# save_loc <- paste0("Runs/", cntry, "_", start_year, "/")
save_loc <- paste0("Runs/", cntry, "_", start_year, "/")
# save_loc <- paste0("Runs/", cntry, "_", start_year, "_rm_missing_type/")
# save_loc <- paste0("Runs/", cntry, "_", start_year, "_only_tmc_pre_2000/")
# save_loc <- paste0("Runs/", cntry, "_", start_year, "_oos/")
threemc::create_dirs_r(save_loc)


#### Preparing circumcision data ####

# pull latest and first censoring year from survey_id
cens_year <- max(survey_years)
start_year <- min(c(survey_years - 2, start_year)) # have lower bound on start

# Prepare circ data, and normalise survey weights and apply Kish coefficients.
survey_circ_preprocess <- prepare_survey_data(
  areas               = areas,
  survey_circumcision = survey_circumcision,
  area_lev            = area_lev,
  start_year          = start_year,
  cens_year           = cens_year,
  cens_age            = cens_age,
  rm_missing_type     = rm_missing_type,
  norm_kisk_weights   = TRUE
)

if (nrow(survey_circ_preprocess) == 0) {
  message("no valid surveys at this level") # move inside function!
}

# include indicator to determine whether there is any type distinction for cntry
if (all(is.na(survey_circ_preprocess$circ_who) &
        is.na(survey_circ_preprocess$circ_where))) {
  print("No type distinction made in valid surveys for this country")
  is_type <- FALSE
  paed_age_cutoff <- NULL # don't make design matrices piecewise on age
} else is_type <- TRUE


#### Shell dataset to estimate empirical rate ####

# Skeleton dataset

# Shell dataset creation changed in the following ways:
#    1) Single age population counts are now added to the output data set.
#       This is needed early on, as we will be aggregating the estimated
#       probabilities and cumulative incidence from the model on the district
#       level in order to include survey data not on the district level (or
#       administrative level of interest). These will be weighted by population.
#    2) Now produced for multiple levels of area_level. Function sets
#       up the shell dataset for all admin boundaries between national (admin 0)
#       and the district level (or administrative level of interest) rather
#       than letting survey_circumcision dictate one level of interest
#
# The internal aggregating to get the obs_mmc etc. still works as the functions
# now uses the new "space" variable defined above. These functions treat each
# "space" as a stratification variable and therefore self-contained. This has
# implications later where we have to specify the administrative boundaries
# we are primarily modelling on.


out <- create_shell_dataset(
  survey_circumcision = survey_circ_preprocess,
  populations         = populations,
  areas               = areas,
  area_lev            = area_lev,
  start_year          = start_year,
  end_year            = forecast_year,
  time1               = "time1",
  time2               = "time2",
  strat               = "space",
  age                 = "age",
  circ                = "indweight_st"
)


# Temp: fill in NAs in population with next known value
# (not very correct but will do for now!!)
if (!all(!is.na(out$population))) {
  message("Filling in missing populations with earliest known value")
  out <- out %>%
    group_by(area_id, area_name, area_level, space, circ_age, age) %>%
    tidyr::fill(population, .direction = "downup") %>%
    ungroup()
}


#### Dataset for modelling ####

dat_tmb <- threemc_prepare_model_data(
  out               = out,
  areas             = areas,
  area_lev          = area_lev,
  aggregated        = TRUE,
  weight            = "population",
  k_dt_time         = NULL,
  k_dt_age          = 5,
  rw_order          = rw_order,
  paed_age_cutoff   = paed_age_cutoff,
  inc_time_tmc      = inc_time_tmc
)


#### Modelling circumcision probabilities ####

# specify TMB model
if (is_type == TRUE & is.null(paed_age_cutoff)) {
  mod <- "Surv_SpaceAgeTime_ByType_withUnknownType"
} else if (is_type == TRUE) {
  mod <- "Surv_SpaceAgeTime_ByType_withUnknownType_Const_Paed_MMC"
} else {
  mod <- "Surv_SpaceAgeTime"
}

# remove mmc time correlation parameters, if fitting with RW precision matrix
if (!is.null(dat_tmb$Q_time)) {
  mod <- paste0(mod, "_RW")
}

# dummy paediatric MMC matrices
if (is.null(paed_age_cutoff)) {
  X_fixed_mmc_paed <- X_age_mmc_paed <- X_space_mmc_paed <- data.frame(0)
}

# dummy time TMC matrices
if (inc_time_tmc == FALSE) {
  X_time_tmc <- data.frame(0)
}

# Initial values
parameters <- with(
  dat_tmb,
  list(
    # intercept
    "u_fixed_mmc"            = rep(-5, ncol(X_fixed_mmc)),
    "u_fixed_mmc_paed"       = rep(-5, ncol(X_fixed_mmc_paed)),
    "u_fixed_tmc"            = rep(-5, ncol(X_fixed_tmc)),
    # age random effect
    "u_age_mmc"              = rep(0, ncol(X_age_mmc)),
    "u_age_mmc_paed"         = rep(0, ncol(X_age_mmc_paed)),
    "u_age_tmc"              = rep(0, ncol(X_age_tmc)),
    # time random effect for (non-paed) MMC
    "u_time_mmc"             = rep(0, ncol(X_time_mmc)),
    # time random effect for TMC
    "u_time_tmc"             = rep(0, ncol(X_time_tmc)),
    # Space random effect (district)
    "u_space_mmc"            = rep(0, ncol(X_space_mmc)),
    "u_space_mmc_paed"       = rep(0, ncol(X_space_mmc_paed)),
    "u_space_tmc"            = rep(0, ncol(X_space_tmc)),
    # Interactions for MMC
    "u_agetime_mmc"          = matrix(0, ncol(X_age_mmc), ncol(X_time_mmc)),
    "u_agespace_mmc"         = matrix(0, ncol(X_age_mmc), ncol(X_space_mmc)),
    "u_spacetime_mmc"        = matrix(0, ncol(X_time_mmc), ncol(X_space_mmc)),
    "u_agespace_mmc_paed"    = matrix(0, ncol(X_age_mmc_paed), ncol(X_space_mmc_paed)),
    # Interactions for TMC
    "u_agespace_tmc"         = matrix(0, ncol(X_age_tmc), ncol(X_space_tmc)),
    # Autocorrelation parameters for priors
    # Variance
    "logsigma_age_mmc"            = 0,
    "logsigma_age_mmc_paed"       = 0,
    "logsigma_time_mmc"           = 0,
    "logsigma_space_mmc"          = 0,
    "logsigma_space_mmc_paed"     = 0,
    "logsigma_agetime_mmc"        = 0,
    "logsigma_agespace_mmc"       = 0,
    "logsigma_agespace_mmc_paed"  = 0,
    "logsigma_spacetime_mmc"      = 0,
    "logsigma_age_tmc"            = 0,
    "logsigma_time_tmc"           = 0,
    "logsigma_space_tmc"          = 0,
    "logsigma_agespace_tmc"       = 0,
    # Mean
    "logitrho_mmc_time1"          = 2,
    "logitrho_mmc_time2"          = 2,
    "logitrho_mmc_time3"          = 2,
    "logitrho_mmc_age1"           = 2,
    "logitrho_mmc_paed_age1"      = 2,
    "logitrho_mmc_age2"           = 2,
    "logitrho_mmc_paed_age2"      = 2,
    "logitrho_mmc_age3"           = 2,
    "logitrho_tmc_time1"          = 2,
    "logitrho_tmc_age1"           = 2,
    "logitrho_tmc_age2"           = 2
  )
)

# remove paed-related parameters if not desired
if (is.null(paed_age_cutoff)) {
  # remove paed-related parameters
  parameters <- parameters[!grepl("paed", names(parameters))]
}

# remove mmc time correlation parameters, if fitting with RW precision matrix
if ("Q_time" %in% names(dat_tmb)) {
  parameters <- parameters[!grepl("logitrho_mmc_time", names(parameters))]
}

# remove time tmc terms, if not fitting model with non-constant tmc over time
if (inc_time_tmc == FALSE) {
  parameters <- parameters[
    # !grepl("time", names(parameters)) & !grepl("tmc", names(parameters))
    !names(parameters) %in% c("u_time_tmc", "logsigma_time_tmc", "logitrho_tmc_time1")
  ]
} else {
  mod <- paste0(mod, "2")
}

library(data.table)
fit_model <- function(fit = NULL,
                      dat_tmb = NULL,
                      mod = NULL,
                      parameters = NULL,
                      maps = NULL,
                      randoms = c(
                        "u_time_mmc", "u_age_mmc", "u_space_mmc",
                        "u_agetime_mmc", "u_agespace_mmc",
                        "u_spacetime_mmc", "u_age_tmc",
                        "u_space_tmc", "u_agespace_tmc"
                      ),
                      sample = TRUE,
                      smaller_fit_obj = FALSE,
                      sdreport = FALSE,
                      N = 1000,
                      verbose = TRUE,
                      ...) {


  # If model is not specified, allow function to choose based on dat_tmb
  # This also abstracts esoteric model specification from the user
  if (is.null(mod)) { # should be it's own function! Can test more easily then

    if (!is.null(parameters)) {
      param_names <- names(parameters)
    } else if (!is.null(fit)) {
      param_names <- names(fit$par)
      # add mapped parameters which won't be in fit$par, if appropriate
      if (!is.null(maps)) param_names <- c(param_names, names(maps))
    } else {
      stop("Please provide one of `parameters` or `fit`")
    }

    # Start with model with no type information
    mod <- "Surv_SpaceAgeTime"

    # determine whether there was information on circumcision type in `out`
    cond <- "type_info" %in% names(dat_tmb) && dat_tmb$type_info == TRUE ||
      # if we don't define dat_tmb, and instead have fit to be re-sampled from
      (!is.null(fit) && any(grepl("mmc", param_names)))

    if (cond) {
      mod <- paste0(mod, "_ByType_withUnknownType")
    }

    if (any(grepl("paed", param_names))) {
      mod <- paste0(mod, "_Const_Paed_MMC")
    }

    # if there are no (MMC for type mod) time autocorr hyperpars, use RW model
    if (!grepl("_ByType_withUnknownType", mod)) {
      cond <- grepl("logitrho", param_names) & grepl("time", param_names)
    } else {
      cond <- grepl("logitrho_mmc", param_names) & grepl("time", param_names)
    }
    # if there is a TMC autocorr hyperpar, use RW model for MMC and AR for TMC
    tmc_cond <- "logitrho_tmc_time1" %in% param_names


    if (all(cond == FALSE)) {
      if (tmc_cond) {
        mod <- paste0(mod, "_RW_MMC")
      } else {
        mod <- paste0(mod, "_RW")
      }
    }

    # if there is a time term for TMC, use the model with non-constant TMC
    # if (dat_tmb$type_info == TRUE && "u_time_tmc" %in% param_names) {
    cond <- dat_tmb$type_info == TRUE
    if (length(cond) == 0) cond <- "u_time_tmc" %in% param_names
    cond <- cond && "u_time_tmc" %in% param_names
    if (cond) {
      mod <- paste0(mod, "2")
    }

    # "RW_MMC" mod only valid where time TMC effect is used
    if (grepl("_RW_MMC", mod) && !grepl("2", mod)) {
      stop(paste(
        "Model with RW for MMC temporal parameters only available when ",
        "including time TMC effect"
      ))
    }
    message("mod not supplied, mod used = ", mod)
  }

  if (is.null(mod)) stop("Please provide one of `mod`, `parameters` or `fit`")

  # for specified "smaller fit" object (i.e. fit which requires resampling)
  if (!is.null(fit)) {
    if (!is.null(fit$sample)) {
      message("Sample already present in fit object, returning `fit`")
      return(fit)
    }
    if (!is.null(dat_tmb)) {
      message(paste0(
        "No need to specify dat_tmb or parameters for non-null fit, as they",
        " are replaced by those stored in fit"
      ))
    }
    # pull dat_tmb and parameters from small fit
    dat_tmb <- fit$tmb_data
    parameters <- split(fit$par.full, names(fit$par.full))
    init_params <- fit$par_init
    # pull different pars depending on whether the model has mmc/tmc split
    if (!mod %in% c("Surv_SpaceAgeTime", "Surv_SpaceAgeTime_RW")) {
      fit$par_init <- fit$par_init[names(fit$par_init) %in% names(parameters)]
      parameters <- parameters[names(fit$par_init)]

      if (any(is.na(names(parameters)))) {
        message("Removing NA parameters, may want to check specifications...")
        parameters <- parameters[!is.na(names(parameters))]
      }
    } else {
      # only need names and lengths, not values
      names(init_params) <- stringr::str_remove_all(
        names(init_params), "_mmc|_tmc"
      )
      init_params <- init_params[!duplicated(names(init_params))]
      parameters <- parameters[names(init_params)]
      # remove duplicate parameters
      parameters <- parameters[!duplicated(names(parameters))]
    }

    if (!is.null(maps)) {
      # ensure mapped parameters are in the same order as parameters for model
      mapped_pars <- is.na(names(parameters))
      param_order <- names(init_params)[mapped_pars]
      maps <- maps[match(names(maps), param_order)]

      # replace NAs in parameters with mapped parameters in par_init
      parameters[mapped_pars] <- init_params[
        names(init_params) %chin% names(maps)
      ]
      names(parameters)[mapped_pars] <- names(maps)
    }

    is_matrix <- vapply(init_params, is.matrix, logical(1))
    parameters[is_matrix] <- Map(matrix,
                                 parameters[is_matrix],
                                 nrow = lapply(init_params[is_matrix], nrow),
                                 ncol = lapply(init_params[is_matrix], ncol)
    )
    # if no fit == NULL, must have non-null dat_tmb & parameters
  } else {
    if (is.null(dat_tmb) || is.null(parameters)) {
      stop("Please specify non-null dat_tmb and parameters")
    }
  }

  # remove "mmc" from parameter & matrix names if required
  if (mod %in% c("Surv_SpaceAgeTime", "Surv_SpaceAgeTime_RW")) {
    remove_type_distinction <- function(x) {
      names(x) <- stringr::str_remove(names(x), "_mmc")
      x <- x[!grepl("_tmc", names(x))]
    }

    dat_tmb <- remove_type_distinction(
      dat_tmb[!names(dat_tmb) %chin% c("A_mmc", "A_tmc")]
    )
    names(dat_tmb)[names(dat_tmb) == "A_mc"] <- "A"

    parameters <- remove_type_distinction(parameters)
    randoms <- unique(stringr::str_remove(randoms, "_tmc|_mmc"))
  }

  # Only have named random parameters
  randoms <- randoms[randoms %chin% names(parameters)]
  if (length(randoms) == 0) randoms <- NULL

  # remove null parameters
  null_pars <- vapply(parameters, is.null, FUN.VALUE = logical(1))
  if (any(null_pars)) {
    message("Removing NULL parameters, check specification...")
    parameters <- parameters[!null_pars]
  }

  if (verbose) message("Creating TMB object with `TMB::MakeADFun`...")
  # Create TMB object
  obj <- TMB::MakeADFun(
    dat_tmb,
    parameters,
    random = randoms,
    map = maps,
    method = "BFGS",
    hessian = TRUE,
    DLL = mod,
    ...
  )
  # for specified fit, simply resample and return
  if (!is.null(fit)) {

    if (verbose) message("Resampling from `fit`...")
    fit$obj <- obj
    fit$obj$fn()
    fit <- circ_sample_tmb(
      fit = fit, obj = obj, nsample = N, sdreport = sdreport
    )
    fit$tmb_data <- fit$par_init <- NULL # make fit object smaller for saving
    return(fit)
  }

  # Run optimiser (use optim if all pars are fixed, nlminb otherwise)
  if (verbose) message("Optimising...")
  if (length(obj$par) == 0) {
    message("using optim")
    opt <- do.call(stats::optim, obj, ...)
  } else {
    opt <- stats::nlminb(
      start   = obj$par,
      obj     = obj$fn,
      gr      = obj$gr,
      control = list(trace = 1),
      ...
    )
  }

  # sample from TMB fit
  if (sample == TRUE) {
    if (verbose) message("Sampling...")
    fit <- circ_sample_tmb(
      obj = obj, opt = opt, nsample = N, sdreport = sdreport
    )
    # return smaller fit object
    if (smaller_fit_obj == TRUE) {
      if (verbose) message("Minimising fit object...")
      fit <- minimise_fit_obj(fit, dat_tmb, parameters)
    }
    return(fit)
  } else {
    return(opt)
  }
}

# fit model with TMB
fit <- fit_model(
  dat_tmb    = dat_tmb,
  mod        = "threemc",
  parameters = parameters,
  randoms    = c(
    "u_time_mmc", "u_age_mmc", "u_age_mmc_paed", "u_space_mmc",
    "u_agetime_mmc", "u_agespace_mmc", "u_agespace_mmc_paed",
    "u_spacetime_mmc",
    "u_time_tmc", "u_age_tmc", "u_space_tmc", "u_agespace_tmc"
  ),
  N          = N
)

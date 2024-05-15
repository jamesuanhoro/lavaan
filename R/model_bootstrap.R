do_model_based_para_boot <- function(
    fit,
    boot_model = "rmsea",
    replications = 1000L,
    parallel = "no",
    ncpus = max(c(1, parallel::detectCores()) - 1),
    iseed = NULL) {
  tmp_fit <- update(
    fit,
    se = "bootstrap",
    bootstrap = replications,
    parallel = parallel,
    ncpus = ncpus,
    iseed = iseed,
    do.fit = FALSE
  )
  tmp_fit@Options$do.fit <- TRUE

  boot_param <- 0
  err_var_vec <- NULL
  if (boot_model == "rmsea") {
    boot_param <- fitmeasures(fit, "rmsea")["rmsea"]
  } else if (boot_model == "crmr") {
    boot_param <- lavResiduals(fit, type = "crmr")$summary["ucrmr", ]
    err_var_vec <- diag(lavInspect(fit, "est")$theta)
  } else {
    warning("performing the standard parametric boostrap")
    boot_model <- NULL
  }

  boot_vcov <- lav_model_vcov(
    lavmodel = fit@Model,
    lavsamplestats = fit@SampleStats,
    lavoptions = tmp_fit@Options,
    lavdata = fit@Data,
    lavimplied = fit@implied,
    lavh1 = fit@h1,
    lavcache = fit@Cache,
    lavpartable = fit@ParTable,
    boot_type = "parametric",
    boot_model = boot_model,
    boot_param = boot_param,
    err_var_vec = err_var_vec
  )

  obj_sub <- create_lav_sub_objects_mod_boot(
    se = tmp_fit@Options$se,
    information = tmp_fit@Options$information,
    par_table = fit@ParTable,
    model = fit@Model,
    boot_vcov = boot_vcov
  )

  result <- create_lav_object_fit_mod_boot(
    fit = fit,
    tmp_fit = tmp_fit,
    obj_sub = obj_sub
  )
  return(result)
}

create_lav_sub_objects_mod_boot <- function(
    se, information, par_table, model, boot_vcov) {
  lavvcov <- list(
    se = se,
    information = information,
    vcov = boot_vcov
  )
  lavboot <- list()
  lavboot$coef <- attr(boot_vcov, "BOOT.COEF")
  par_tab <- par_table
  par_tab$se <- lav_model_vcov_se(
    lavmodel = model,
    lavpartable = par_tab,
    VCOV = boot_vcov,
    BOOT = lavboot$coef
  )
  ret <- list(
    lavvcov = lavvcov, lavboot = lavboot, par_tab = par_tab
  )
  return(ret)
}

create_lav_object_fit_mod_boot <- function(fit, tmp_fit, obj_sub) {
  boot_fit <- new(
    "lavaan",
    version      = as.character(packageVersion("lavaan")),
    call         = tmp_fit@call, # match.call
    timing       = tmp_fit@timing, # list
    Options      = tmp_fit@Options, # list
    ParTable     = obj_sub$par_tab, # list
    pta          = tmp_fit@pta, # list
    Data         = tmp_fit@Data, # S4 class
    SampleStats  = tmp_fit@SampleStats, # S4 class
    Model        = fit@Model, # S4 class
    Cache        = tmp_fit@Cache, # list
    Fit          = fit@Fit, # S4 class
    boot         = obj_sub$lavboot, # list
    optim        = fit@optim, # list
    implied      = fit@implied, # list
    loglik       = fit@loglik, # list
    vcov         = obj_sub$lavvcov, # list
    test         = fit@test, # list
    h1           = fit@h1, # list
    baseline     = fit@baseline, # list
    internal     = list(), # empty list
    external     = list() # empty list
  )
  return(boot_fit)
}

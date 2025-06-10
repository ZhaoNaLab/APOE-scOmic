setup_future <- function(object, plan=c('multicore', 'multisession'), workers=4) {
  plan = match.arg(plan)
  library(future)
  options(future.globals.maxSize = object.size(object)*10)
  options(future.resolve.recursive = Inf)
  options(future.wait.timeout = 10)
  plan(strategy=plan, workers=workers)
  rlog::log_info(paste('util/future.r setup_future: Running', plan, 'future with', workers, 'workers.'))
}

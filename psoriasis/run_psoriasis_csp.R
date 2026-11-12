# Fresh execution of the supplied, unchanged CSP code. Run from the package folder:
# Rscript run_psoriasis_csp.R [output_directory]
args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args)) args[1] else "reproduced/psoriasis_R"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
source("csp_functions.R")
raw <- read.csv("pso-psa-data.csv", stringsAsFactors = FALSE)
raw <- subset(raw, indication == "psoriasis")
dat <- data.frame(study = paste0("Study_", raw$study),
                  id = as.numeric(as.factor(raw$study)),
                  t = as.character(raw$treatment), r = raw$r, n = raw$n)
stopifnot(nrow(dat) == 56, length(unique(dat$id)) == 23,
          length(unique(dat$t)) == 10, !anyDuplicated(dat[,c("id","t")]))
for (model in c("fixed", "random")) {
  name <- if (model == "fixed") "common" else "random"
  result <- fit_csp_nma(dat, model = model, verbose = FALSE)
  write.csv(data.frame(target = names(result$theta_hat), theta = result$theta_hat,
                        tau2 = result$tau2), file.path(out,paste0("psoriasis_fit_",name,".csv")), row.names=FALSE)
  write.csv(result$overall,file.path(out,paste0("psoriasis_overall_",name,".csv")),row.names=FALSE)
  write.csv(data.frame(target=rownames(result$P),result$P,check.names=FALSE),
            file.path(out,paste0("psoriasis_P_",name,".csv")),row.names=FALSE)
  ee <- list()
  for (target in result$pair_names) {
    for (blk in result$blocks) {
      pk <- result$P[target,blk$row_names];names(pk) <- blk$local_pair_names
      e <- .canonical_endpoint_edges(pk,blk$arms,target,blk$study,blk$id,result$treatments)
      if (nrow(e)) {e$target <- target;ee[[length(ee)+1]] <- e}
    }
  }
  write.csv(do.call(rbind,ee),file.path(out,paste0("psoriasis_endpoint_edges_",name,".csv")),row.names=FALSE)
  write.csv(q_decomposition(result),file.path(out,paste0("psoriasis_Q_",name,".csv")),row.names=FALSE)
  saveRDS(result,file.path(out,paste0("psoriasis_fit_",name,".rds")))
  print(result$overall);cat("Model:",name,"tau2:",result$tau2,"\n")
}
# Compatibility check for the supplied illustrative data (same arm counts as COVID).
# Evaluate only the data assignment; the user's plotting demonstration remains unchanged.
for (ex in parse("csp_example.R")) {
  if (is.call(ex) && identical(ex[[1]],as.name("<-")) && identical(ex[[2]],as.name("example_data"))) eval(ex)
}
stopifnot(exists("example_data"))
for (name in c("common","random")) {
  result <- if (name=="common") fit_csp_nma(example_data,verbose=FALSE) else
    fit_csp_nma(example_data,model="random",tau2=0.002103984,verbose=FALSE)
  write.csv(data.frame(target=names(result$theta_hat),theta=result$theta_hat),
            file.path(out,paste0("covid_CSP_fit_",name,".csv")),row.names=FALSE)
  write.csv(data.frame(target=rownames(result$P),result$P,check.names=FALSE),
            file.path(out,paste0("covid_CSP_P_",name,".csv")),row.names=FALSE)
  ee <- list()
  for (target in result$pair_names) for (blk in result$blocks) {
    pk <- result$P[target,blk$row_names];names(pk) <- blk$local_pair_names
    e <- .canonical_endpoint_edges(pk,blk$arms,target,blk$study,blk$id,result$treatments)
    if (nrow(e)) {e$target <- target;ee[[length(ee)+1]] <- e}
  }
  write.csv(do.call(rbind,ee),file.path(out,paste0("covid_CSP_edges_",name,".csv")),row.names=FALSE)
}
writeLines(capture.output(sessionInfo()),file.path(out,"sessionInfo.txt"))
write.csv(installed.packages()[,c("Package","Version")],file.path(out,"installed_packages.csv"),row.names=FALSE)

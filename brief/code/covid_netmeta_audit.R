# Rerun of the netmeta contribution audit on the five-study COVID-19 network.
# Common-effect log odds ratio; corrected logits (r+1/2)/(n-r+1/2), arm variances 1/(r+1/2)+1/(n-r+1/2).
suppressPackageStartupMessages(library(netmeta))
cat("R", R.version.string, "| netmeta", as.character(packageVersion("netmeta")),
    "| cccp", as.character(packageVersion("cccp")), "\n")
arms <- read.csv(commandArgs(TRUE)[1])
arms$lo <- log((arms$r + .5) / (arms$n - arms$r + .5))
arms$v  <- 1/(arms$r + .5) + 1/(arms$n - arms$r + .5)
# all within-study pairs, TE = treat1 vs treat2 = logit(treat1) - logit(treat2)
pw <- do.call(rbind, lapply(split(arms, arms$study), function(s) {
  s <- s[order(s$t), ]; cmb <- combn(nrow(s), 2)
  data.frame(studlab = s$study[1], treat1 = s$t[cmb[1, ]], treat2 = s$t[cmb[2, ]],
             TE = s$lo[cmb[1, ]] - s$lo[cmb[2, ]], seTE = sqrt(s$v[cmb[1, ]] + s$v[cmb[2, ]]))
}))
net <- netmeta(TE, seTE, treat1, treat2, studlab, data = pw, sm = "OR", common = TRUE, random = FALSE)
trts <- net$trts
# pooled direct estimates with the multi-arm adjusted (Kron) weights w_{k,jl} = g_j g_l / sum g
g <- 1/arms$v
direct <- list()
for (s in split(arms, arms$study)) {
  gs <- 1/s$v; dk <- sum(gs); cmb <- combn(nrow(s), 2)
  for (c in seq_len(ncol(cmb))) {
    i <- cmb[1, c]; j <- cmb[2, c]; key <- paste(sort(c(s$t[i], s$t[j])), collapse = ":")
    a <- which(sort(c(s$t[i], s$t[j]))[1] == s$t[c(i, j)]); b <- 3 - a  # a = index of alphabetically first
    te <- s$lo[c(i, j)][a] - s$lo[c(i, j)][b]                            # first vs second
    w <- gs[i] * gs[j] / dk
    direct[[key]] <- rbind(direct[[key]], c(w = w, te = te))
  }
}
ybar <- sapply(direct, function(m) sum(m[, "w"] * m[, "te"]) / sum(m[, "w"]))
# compare with netmeta's own direct estimates
nd <- sapply(names(ybar), function(k) { ab <- strsplit(k, ":")[[1]]; net$TE.direct.common[ab[1], ab[2]] })
cat("max |ybar(Kron) - netmeta TE.direct.common| =", max(abs(ybar - nd)), "\n")
# aggregate (short) hat matrix: rows = all comparisons, columns = direct comparisons
H <- hatmatrix(net, method = "Davies", type = "full")$common
cat("hat matrix dims", dim(H), "\n")
colkeys <- colnames(H)
# map column labels "A:B" to ybar orientation (first vs second as labelled)
ycol <- sapply(colkeys, function(k) { ab <- strsplit(k, ":")[[1]]; key <- paste(sort(ab), collapse = ":")
  if (!key %in% names(ybar)) return(0)
  if (ab[1] == sort(ab)[1]) ybar[[key]] else -ybar[[key]] })
nodirect <- colkeys[!sapply(colkeys, function(k) paste(sort(strsplit(k, ":")[[1]]), collapse = ":") %in% names(ybar))]
cat("columns without direct evidence:", nodirect, "| max |H| in those columns:", max(abs(H[, nodirect, drop = FALSE])), "\n")
targets <- t(combn(trts, 2)); rows <- paste(targets[, 1], targets[, 2], sep = ":")
fitted <- sapply(seq_len(nrow(targets)), function(i) net$TE.common[targets[i, 1], targets[i, 2]])
hat_recon <- sapply(rows, function(r) sum(H[r, ] * ycol))
cat("max |hat row * direct(Kron) - fitted| =", max(abs(hat_recon - fitted)), "\n")
ndcol <- sapply(colkeys, function(k) { ab <- strsplit(k, ":")[[1]]; v <- net$TE.direct.common[ab[1], ab[2]]; if (is.na(v)) 0 else v })
cat("max |hat row * direct(netmeta TE.direct.common) - fitted| =", max(abs(sapply(rows, function(r) sum(H[r, ] * ndcol)) - fitted)), "\n")
methods <- c("shortestpath", "randomwalk", "cccp", "pseudoinverse")
res <- data.frame(Target = rows, Fitted = fitted)
for (m in methods) {
  cm <- netcontrib(net, method = m, hatmatrix.F1000 = FALSE)$common
  cc <- colnames(cm); stopifnot(all(cc %in% colkeys))
  sub <- sapply(rows, function(r) sum(sign(H[r, cc]) * cm[r, cc] * ycol[cc]))
  res[[m]] <- sub - fitted
  cat(m, ": row sums of percentages range", range(rowSums(cm)), "\n")
}
print(res, digits = 6)
cat("median |error|:", sapply(methods, function(m) median(abs(res[[m]]))), "\n")
cat("max |error|:   ", sapply(methods, function(m) max(abs(res[[m]]))), "\n")
write.csv(res, "covid_netmeta_audit_results.csv", row.names = FALSE)

# ---- label-orientation (unsigned) substitution, as in the archived October table
res_u <- data.frame(Target = rows, Fitted = fitted)
for (m in methods) {
  cm <- netcontrib(net, method = m, hatmatrix.F1000 = FALSE)$common; cc <- colnames(cm)
  res_u[[m]] <- sapply(rows, function(r) sum(cm[r, cc] * ycol[cc])) - fitted
}
print(res_u, digits = 6)
cat("unsigned median |error|:", sapply(methods, function(m) median(abs(res_u[[m]]))), "\n")
cat("unsigned max |error|:   ", sapply(methods, function(m) max(abs(res_u[[m]]))), "\n")
write.csv(res_u, "covid_netmeta_audit_results_labelorientation.csv", row.names = FALSE)
fmt <- function(x) ifelse(x < 0, sprintf("$%.6f$", x), sprintf("%.6f", x))
con <- file("covid_netmeta_rows_signed.tex", "w")
for (i in seq_len(nrow(res))) writeLines(paste(res$Target[i], "&", fmt(res$Fitted[i]), "&", paste(fmt(unlist(res[i, methods])), collapse = " & "), "\\\\"), con)
close(con)
con <- file("covid_netmeta_rows_unsigned.tex", "w")
for (i in seq_len(nrow(res_u))) writeLines(paste(res_u$Target[i], "&", fmt(res_u$Fitted[i]), "&", paste(fmt(unlist(res_u[i, methods])), collapse = " & "), "\\\\"), con)
close(con)

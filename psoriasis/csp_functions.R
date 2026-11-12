suppressPackageStartupMessages({
  library(Matrix)
  library(MASS)
  library(dplyr)
  library(ggplot2)
  library(plotly)
})

SAFE_STUDY_SEP <- " ||| "

# ============================================================
# Internal helpers
# ============================================================
.parse_pair <- function(pair_name) {
  x <- strsplit(pair_name, ":", fixed = TRUE)[[1]]
  list(a = x[1], b = x[2], u = x[1], v = x[2])
}

.get_block_map <- function(res) {
  mp <- res$blocks
  names(mp) <- sapply(mp, function(b) as.character(b$study))
  mp
}

.clean_study_label <- function(x, max_chars = 35) {
  x <- as.character(x)
  x <- gsub("\r", "", x, fixed = TRUE)
  x <- gsub("\n", " ", x, fixed = TRUE)
  x <- trimws(x)
  x <- sub("\\s*\\(.*?\\)", "", x)
  x <- gsub("\\s+", " ", x)
  ifelse(nchar(x) > max_chars, paste0(substr(x, 1, max_chars - 3), "..."), x)
}

.make_display_labels_unique <- function(x) {
  x <- as.character(x)
  out <- character(length(x))
  seen <- new.env(parent = emptyenv())
  for (i in seq_along(x)) {
    lab <- x[i]
    if (!exists(lab, envir = seen, inherits = FALSE)) {
      assign(lab, 1L, envir = seen)
      out[i] <- lab
    } else {
      k <- get(lab, envir = seen, inherits = FALSE) + 1L
      assign(lab, k, envir = seen)
      out[i] <- paste0(lab, " <", k, ">")
    }
  }
  out
}

.strip_unique_suffix <- function(x) sub(" <[0-9]+>$", "", x)

.split_path_studies <- function(studies_str, study_sep = SAFE_STUDY_SEP) {
  studies_str <- as.character(studies_str)
  if (length(studies_str) == 0 || is.na(studies_str) || !nzchar(studies_str)) return(character(0))
  if (grepl(study_sep, studies_str, fixed = TRUE)) {
    strsplit(studies_str, study_sep, fixed = TRUE)[[1]]
  } else {
    studies_str
  }
}

.canonical_endpoint_edges <- function(p_k, arms, target_pair, study_name, study_id,
                                      global_treatments, tol = 1e-12) {
  arms <- sort(unique(as.character(arms)))
  b <- setNames(rep(0, length(arms)), arms)

  for (nm in names(p_k)) {
    w <- as.numeric(p_k[[nm]])
    if (abs(w) <= tol) next
    uv <- strsplit(nm, ":", fixed = TRUE)[[1]]
    u <- uv[1]
    v <- uv[2]
    if (!(u %in% arms) || !(v %in% arms)) {
      stop("Local contrast contains a treatment not present in the study: ", nm)
    }
    b[[u]] <- b[[u]] + w
    b[[v]] <- b[[v]] - w
  }

  s <- pmax(b, 0)
  d <- pmax(-b, 0)
  names(s) <- names(d) <- arms

  edges <- list()
  frag <- 0L

  add_edge <- function(u, v, w, type) {
    if (w <= tol) return(invisible(NULL))
    frag <<- frag + 1L
    edges[[frag]] <<- data.frame(
      frag_id = frag,
      id = study_id,
      study = study_name,
      from = u,
      to = v,
      weight = as.numeric(w),
      type = type,
      stringsAsFactors = FALSE
    )
    invisible(NULL)
  }

  pp <- .parse_pair(target_pair)
  a <- pp$u
  b_target <- pp$v

  if (a %in% arms && b_target %in% arms) {
    w_direct <- min(s[[a]], d[[b_target]])
    if (w_direct > tol) {
      add_edge(a, b_target, w_direct, "Dir")
      s[[a]] <- s[[a]] - w_direct
      d[[b_target]] <- d[[b_target]] - w_direct
      if (s[[a]] <= tol) s[[a]] <- 0
      if (d[[b_target]] <= tol) d[[b_target]] <- 0
    }
  }

  max_iter <- 2L * length(arms) + 10L
  iter <- 0L

  repeat {
    srcs <- names(s)[s > tol]
    sinks <- names(d)[d > tol]
    if (length(srcs) == 0L && length(sinks) == 0L) break
    if (length(srcs) == 0L || length(sinks) == 0L) {
      stop("Canonical balance reduction failed: unmatched supply or demand in study ", study_name)
    }

    cand <- expand.grid(u = srcs, v = sinks, stringsAsFactors = FALSE)
    cand$lambda <- mapply(
      function(u, v) min(s[[u]], d[[v]]),
      cand$u,
      cand$v
    )

    lambda_max <- max(cand$lambda)
    tie_tol <- tol * max(1, abs(lambda_max))
    cand <- cand[abs(cand$lambda - lambda_max) <= tie_tol, , drop = FALSE]

    u_rank <- match(cand$u, global_treatments)
    v_rank <- match(cand$v, global_treatments)
    if (anyNA(u_rank) || anyNA(v_rank)) {
      stop("Treatment ordering is inconsistent with the fitted network.")
    }
    cand$rho <- (u_rank - 1L) * length(global_treatments) + v_rank
    cand <- cand[order(cand$rho, cand$u, cand$v), , drop = FALSE]

    u <- cand$u[1]
    v <- cand$v[1]
    w <- cand$lambda[1]
    add_edge(u, v, w, "Ind")

    s[[u]] <- s[[u]] - w
    d[[v]] <- d[[v]] - w
    if (s[[u]] <= tol) s[[u]] <- 0
    if (d[[v]] <= tol) d[[v]] <- 0

    iter <- iter + 1L
    if (iter > max_iter) {
      stop("Canonical balance reduction exceeded its finite-allocation bound in study ", study_name)
    }
  }

  if (length(edges) == 0L) {
    out <- data.frame(
      frag_id = integer(0), id = numeric(0), study = character(0),
      from = character(0), to = character(0), weight = numeric(0),
      type = character(0), stringsAsFactors = FALSE
    )
  } else {
    out <- do.call(rbind, edges)
    rownames(out) <- NULL
  }

  b_check <- setNames(rep(0, length(arms)), arms)
  if (nrow(out) > 0L) {
    for (i in seq_len(nrow(out))) {
      b_check[[out$from[i]]] <- b_check[[out$from[i]]] + out$weight[i]
      b_check[[out$to[i]]] <- b_check[[out$to[i]]] - out$weight[i]
    }
  }

  balance_error <- max(abs(b_check - b))
  if (balance_error > tol * max(1, max(abs(b)))) {
    stop(sprintf(
      "Canonical treatment-balance reconstruction failed in study %s: error = %.3e",
      study_name, balance_error
    ))
  }

  out
}

.collapse_path_segments <- function(seg_df) {
  if (nrow(seg_df) == 0) return(seg_df)
  seg_df <- seg_df[order(seg_df$segment_order), , drop = FALSE]
  rownames(seg_df) <- NULL
  groups <- integer(nrow(seg_df))
  g <- 1L
  groups[1] <- g
  if (nrow(seg_df) >= 2) {
    for (i in 2:nrow(seg_df)) {
      same_study <- seg_df$study[i] == seg_df$study[i - 1]
      connected <- seg_df$from[i] == seg_df$to[i - 1]
      if (same_study && connected) {
        groups[i] <- g
      } else {
        g <- g + 1L
        groups[i] <- g
      }
    }
  }
  out <- lapply(split(seq_len(nrow(seg_df)), groups), function(idx) {
    dd <- seg_df[idx, , drop = FALSE]
    data.frame(study = dd$study[1], from = dd$from[1], to = dd$to[nrow(dd)], stringsAsFactors = FALSE)
  })
  do.call(rbind, out)
}

.find_widest_reporting_path <- function(edge_df, source, sink, tol = 1e-12) {
  if (nrow(edge_df) == 0L) return(NULL)

  nodes <- sort(unique(c(source, sink, edge_df$from, edge_df$to)))
  cap <- setNames(rep(0, length(nodes)), nodes)
  pred_edge <- setNames(rep(NA_integer_, length(nodes)), nodes)
  visited <- setNames(rep(FALSE, length(nodes)), nodes)
  cap[[source]] <- Inf

  repeat {
    available <- nodes[!visited & (is.infinite(cap) | cap > tol)]
    if (length(available) == 0L) break

    vals <- cap[available]
    max_cap <- max(vals)
    if (is.infinite(max_cap)) {
      tied <- available[is.infinite(vals)]
    } else {
      tie_tol <- tol * max(1, abs(max_cap))
      tied <- available[abs(vals - max_cap) <= tie_tol]
    }
    u <- sort(tied)[1]
    visited[[u]] <- TRUE
    if (u == sink) break

    out_idx <- which(edge_df$from == u & edge_df$weight > tol)
    if (length(out_idx) == 0L) next
    out_idx <- out_idx[order(
      edge_df$to[out_idx],
      edge_df$study[out_idx],
      edge_df$id[out_idx],
      edge_df$frag_id[out_idx]
    )]

    for (ix in out_idx) {
      v <- edge_df$to[ix]
      if (visited[[v]]) next
      cand_cap <- min(cap[[u]], edge_df$weight[ix])
      improve_tol <- tol * max(1, abs(cand_cap), abs(cap[[v]]))
      better <- cand_cap > cap[[v]] + improve_tol
      tied_cap <- abs(cand_cap - cap[[v]]) <= improve_tol
      deterministic_tie <- tied_cap && (is.na(pred_edge[[v]]) || ix < pred_edge[[v]])
      if (better || deterministic_tie) {
        cap[[v]] <- cand_cap
        pred_edge[[v]] <- ix
      }
    }
  }

  if (!is.finite(cap[[sink]]) && source != sink) return(NULL)
  if (cap[[sink]] <= tol) return(NULL)

  idx_rev <- integer(0)
  current <- sink
  seen <- character(0)
  while (current != source) {
    if (current %in% seen) stop("Cycle encountered while reconstructing widest path.")
    seen <- c(seen, current)
    ix <- pred_edge[[current]]
    if (is.na(ix)) return(NULL)
    idx_rev <- c(idx_rev, ix)
    current <- edge_df$from[ix]
  }

  rev(idx_rev)
}

.compute_path_stats <- function(res, studies_str, path_str, weight, conf_level = 0.95) {
  zcrit <- qnorm(1 - (1 - conf_level) / 2)
  block_map <- .get_block_map(res)
  verts <- strsplit(path_str, "->", fixed = TRUE)[[1]]
  studies <- .split_path_studies(studies_str, study_sep = res$safe_study_sep)
  if ((length(verts) - 1L) != length(studies)) stop("Mismatch between path and studies for ", path_str)
  segs <- data.frame(study = studies, from = verts[-length(verts)], to = verts[-1], stringsAsFactors = FALSE)
  est_total <- 0
  var_total <- 0
  by_study <- split(segs, segs$study)
  for (st in names(by_study)) {
    blk <- block_map[[st]]
    dd <- by_study[[st]]
    idx <- integer(nrow(dd))
    sgn <- numeric(nrow(dd))
    est <- numeric(nrow(dd))
    for (i in seq_len(nrow(dd))) {
      nm1 <- paste(dd$from[i], dd$to[i], sep = ":")
      nm2 <- paste(dd$to[i], dd$from[i], sep = ":")
      if (nm1 %in% blk$local_pair_names) {
        idx[i] <- match(nm1, blk$local_pair_names); sgn[i] <- 1; est[i] <- blk$y[idx[i]]
      } else if (nm2 %in% blk$local_pair_names) {
        idx[i] <- match(nm2, blk$local_pair_names); sgn[i] <- -1; est[i] <- -blk$y[idx[i]]
      } else {
        stop("Edge not found in study block: ", st)
      }
    }
    est_total <- est_total + sum(est)
    Vsub <- blk$V[idx, idx, drop = FALSE]
    var_total <- var_total + as.numeric(t(sgn) %*% Vsub %*% sgn)
  }
  se_total <- sqrt(max(var_total, 0))
  data.frame(
    estimate = est_total,
    se = se_total,
    lower = est_total - zcrit * se_total,
    upper = est_total + zcrit * se_total,
    weight = weight,
    stringsAsFactors = FALSE
  )
}

.compute_direct_study_stats <- function(res, target, study_name, weight, conf_level = 0.95) {
  zcrit <- qnorm(1 - (1 - conf_level) / 2)
  blk <- .get_block_map(res)[[study_name]]
  pp <- .parse_pair(target)
  a <- pp$a; b <- pp$b
  nm1 <- paste(a, b, sep = ":")
  nm2 <- paste(b, a, sep = ":")
  if (nm1 %in% blk$local_pair_names) {
    idx <- match(nm1, blk$local_pair_names)
    est <- blk$y[idx]; var <- blk$V[idx, idx]
  } else if (nm2 %in% blk$local_pair_names) {
    idx <- match(nm2, blk$local_pair_names)
    est <- -blk$y[idx]; var <- blk$V[idx, idx]
  } else {
    stop("Direct contrast not found in study ", study_name)
  }
  se <- sqrt(max(var, 0))
  data.frame(
    estimate = est,
    se = se,
    lower = est - zcrit * se,
    upper = est + zcrit * se,
    weight = weight,
    stringsAsFactors = FALSE
  )
}

.build_direct_coef_vector <- function(res, target) {
  p_all <- res$P[target, , drop = TRUE]
  q_dir <- setNames(rep(0, length(p_all)), names(p_all))
  dir_df <- res$decomposition[[target]]$direct_paths
  if (is.null(dir_df) || nrow(dir_df) == 0) return(q_dir)
  block_map <- .get_block_map(res)
  pp <- .parse_pair(target)
  a <- pp$a; b <- pp$b
  for (i in seq_len(nrow(dir_df))) {
    st <- dir_df$studies[i]
    w <- dir_df$weight[i]
    blk <- block_map[[st]]
    nm1 <- paste(a, b, sep = ":")
    nm2 <- paste(b, a, sep = ":")
    if (nm1 %in% blk$local_pair_names) {
      idx <- match(nm1, blk$local_pair_names)
      q_dir[blk$row_names[idx]] <- q_dir[blk$row_names[idx]] + w
    } else if (nm2 %in% blk$local_pair_names) {
      idx <- match(nm2, blk$local_pair_names)
      q_dir[blk$row_names[idx]] <- q_dir[blk$row_names[idx]] - w
    } else {
      stop("Direct contrast not found in study ", st)
    }
  }
  q_dir
}

.compute_overall_rows <- function(res, target, conf_level = 0.95, tol = 1e-12) {
  zcrit <- qnorm(1 - (1 - conf_level) / 2)
  y_all <- unlist(lapply(res$blocks, `[[`, "y"))
  row_names <- unlist(lapply(res$blocks, `[[`, "row_names"))
  names(y_all) <- row_names
  V_all <- as.matrix(Matrix::bdiag(lapply(res$blocks, `[[`, "V")))
  rownames(V_all) <- row_names
  colnames(V_all) <- row_names
  p_all <- res$P[target, , drop = TRUE]
  q_dir <- .build_direct_coef_vector(res, target)
  q_ind <- p_all - q_dir
  obj <- res$decomposition[[target]]
  w_dir <- obj$overall$direct_weight
  w_ind <- obj$overall$indirect_weight
  C_dir <- sum(q_dir * y_all)
  C_ind <- sum(q_ind * y_all)
  theta_net <- unname(res$theta_hat[target])
  var_dir_C <- as.numeric(t(q_dir) %*% V_all %*% q_dir)
  var_ind_C <- as.numeric(t(q_ind) %*% V_all %*% q_ind)
  var_net <- as.numeric(t(p_all) %*% V_all %*% p_all)
  theta_dir <- if (abs(w_dir) > tol) C_dir / w_dir else NA_real_
  theta_ind <- if (abs(w_ind) > tol) C_ind / w_ind else NA_real_
  se_dir <- if (abs(w_dir) > tol) sqrt(max(var_dir_C, 0)) / w_dir else NA_real_
  se_ind <- if (abs(w_ind) > tol) sqrt(max(var_ind_C, 0)) / w_ind else NA_real_
  se_net <- sqrt(max(var_net, 0))
  bind_rows(
    data.frame(label = "Overall NMA", component = "Summary", estimate = theta_net, se = se_net,
               lower = theta_net - zcrit * se_net, upper = theta_net + zcrit * se_net, weight = 1,
               stringsAsFactors = FALSE),
    data.frame(label = "Overall direct", component = "Summary", estimate = theta_dir, se = se_dir,
               lower = theta_dir - zcrit * se_dir, upper = theta_dir + zcrit * se_dir, weight = w_dir,
               stringsAsFactors = FALSE),
    data.frame(label = "Overall indirect", component = "Summary", estimate = theta_ind, se = se_ind,
               lower = theta_ind - zcrit * se_ind, upper = theta_ind + zcrit * se_ind, weight = w_ind,
               stringsAsFactors = FALSE)
  )
}


.log_pdet <- function(M, tol = 1e-10) {
  ev <- eigen(M, symmetric = TRUE, only.values = TRUE)$values
  ev <- ev[ev > tol]
  if (length(ev) == 0) return(-Inf)
  sum(log(ev))
}

.make_tau_block <- function(M, tau2) {
  0.5 * tau2 * (M %*% t(M))
}

.estimate_tau2_reml <- function(tilde_y, tilde_X, B, V_within_blocks, M_blocks,
                                lower = 0, upper = NULL, tol = 1e-10,
                                max_upper_expansions = 20L) {
  XB <- tilde_X %*% B

  make_V <- function(tau2) {
    mats <- lapply(seq_along(V_within_blocks), function(i) {
      Vk <- V_within_blocks[[i]]
      Mk <- M_blocks[[i]]
      Vk + .make_tau_block(Mk, tau2)
    })
    as.matrix(Matrix::bdiag(mats))
  }

  neg_reml <- function(tau2) {
    V_tau <- make_V(tau2)
    V_plus <- MASS::ginv(V_tau)
    I_beta <- t(XB) %*% V_plus %*% XB
    beta_hat <- MASS::ginv(I_beta) %*% t(XB) %*% V_plus %*% tilde_y
    resid <- tilde_y - XB %*% beta_hat

    val <- 0.5 * (
      .log_pdet(V_tau, tol = tol) +
      .log_pdet(I_beta, tol = tol) +
      as.numeric(t(resid) %*% V_plus %*% resid)
    )
    if (!is.finite(val)) val <- .Machine$double.xmax / 100
    val
  }

  if (is.null(upper)) {
    yvar <- stats::var(as.numeric(tilde_y), na.rm = TRUE)
    if (!is.finite(yvar)) yvar <- 0
    upper <- max(1, yvar)
  }

  lower <- max(as.numeric(lower)[1], 0)
  upper <- max(as.numeric(upper)[1], lower + 1e-6)
  max_upper_expansions <- max(0L, as.integer(max_upper_expansions)[1])

  opt <- NULL
  for (j in 0:max_upper_expansions) {
    opt <- tryCatch(
      optimize(neg_reml, interval = c(lower, upper)),
      error = function(e) NULL
    )

    if (is.null(opt) || !is.finite(opt$objective) || !is.finite(opt$minimum)) {
      warning("REML optimization failed; tau^2 was set to 0.")
      return(0)
    }

    if (opt$minimum < 0.95 * upper || j == max_upper_expansions) break
    upper <- 2 * upper
  }

  if (opt$minimum >= 0.95 * upper) {
    warning("The REML estimate is near the upper search boundary; consider supplying tau2 directly if a wider search is needed.")
  }

  max(as.numeric(opt$minimum), 0)
}


# ============================================================
# 1. Fit CSP NMA and return P matrix
# ============================================================

fit_csp_nma <- function(dat, cc = 0.5, tol = 1e-12, verbose = TRUE,
                        model = c("fixed", "random"),
                        tau2 = NULL,
                        tau_method = c("REML", "fixed")) {
  model <- match.arg(model)
  tau_method <- match.arg(tau_method)

  dat <- dat[, c("study", "id", "t", "r", "n")]
  dat$study <- as.character(dat$study)
  dat$id <- as.numeric(dat$id)
  dat$t <- as.character(dat$t)
  dat$r <- as.numeric(dat$r)
  dat$n <- as.numeric(dat$n)

  dat <- dat[
    !is.na(dat$study) & !is.na(dat$id) & !is.na(dat$t) &
      !is.na(dat$r) & !is.na(dat$n),
    ,
    drop = FALSE
  ]

  trts <- sort(unique(dat$t))
  ref_trt <- trts[1]
  global_pairs <- t(combn(trts, 2))
  pair_names <- apply(global_pairs, 1, paste, collapse = ":")

  build_full_study_block <- function(dat_study, global_pair_names, cc = 0.5) {
    dat_study <- dat_study[order(dat_study$t), , drop = FALSE]
    arms <- dat_study$t
    r <- dat_study$r
    n <- dat_study$n
    study_name <- unique(dat_study$study)
    study_id <- unique(dat_study$id)
    k <- length(arms)

    r_adj <- r + cc
    n_adj <- n + 2 * cc
    z <- qlogis(r_adj / n_adj)
    s2 <- 1 / r_adj + 1 / (n_adj - r_adj)

    S <- diag(s2, nrow = k, ncol = k)
    cmb <- t(combn(seq_len(k), 2))
    n_pairs <- nrow(cmb)

    M <- matrix(0, nrow = n_pairs, ncol = k)
    y <- numeric(n_pairs)
    local_pair_names <- character(n_pairs)
    row_names <- character(n_pairs)

    for (i in seq_len(n_pairs)) {
      a <- cmb[i, 1]
      b <- cmb[i, 2]
      M[i, a] <- -1
      M[i, b] <-  1
      y[i] <- z[b] - z[a]
      local_pair_names[i] <- paste(arms[a], arms[b], sep = ":")
      row_names[i] <- paste0(study_id, ".", study_name, "::", local_pair_names[i])
    }

    V_within <- M %*% S %*% t(M)

    X <- matrix(0, nrow = n_pairs, ncol = length(global_pair_names))
    colnames(X) <- global_pair_names
    rownames(X) <- row_names
    for (i in seq_len(n_pairs)) {
      X[i, match(local_pair_names[i], global_pair_names)] <- 1
    }

    list(
      study = study_name,
      id = study_id,
      arms = arms,
      y = y,
      M=M,
      V_within = V_within,
      V = V_within,
      X = X,
      local_pair_names = local_pair_names,
      row_names = row_names,
      edge_names_in_tildeP = row_names
    )
  }

  study_list <- split(dat, dat$id)
  blocks <- lapply(study_list, build_full_study_block, global_pair_names = pair_names, cc = cc)

  tilde_y <- unlist(lapply(blocks, `[[`, "y"))
  tilde_X <- do.call(rbind, lapply(blocks, `[[`, "X"))
  proper_edge_names <- rownames(tilde_X)

  B <- matrix(0, nrow = nrow(global_pairs), ncol = length(trts) - 1)
  for (i in seq_len(nrow(global_pairs))) {
    u <- global_pairs[i, 1]
    v <- global_pairs[i, 2]
    if (u == ref_trt) {
      B[i, match(v, setdiff(trts, ref_trt))] <- 1
    } else {
      B[i, match(u, setdiff(trts, ref_trt))] <- -1
      B[i, match(v, setdiff(trts, ref_trt))] <- 1
    }
  }

  if (model == "random") {
    if (is.null(tau2)) {
      if (tau_method != "REML") {
        stop("For model = \"random\", either supply tau2 or use tau_method = \"REML\".")
      }
      tau2 <- .estimate_tau2_reml(
        tilde_y = tilde_y,
        tilde_X = tilde_X,
        B = B,
        V_within_blocks = lapply(blocks, `[[`, "V_within"),
        M_blocks = lapply(blocks, `[[`, "M"),
        tol = tol
      )
    }
    tau2 <- max(as.numeric(tau2)[1], 0)
  } else {
    tau2 <- 0
  }

  blocks <- lapply(blocks, function(blk) {
    blk$V <- blk$V_within + .make_tau_block(blk$M, tau2)
    blk
  })

  tilde_V <- as.matrix(bdiag(lapply(blocks, `[[`, "V")))
  rownames(tilde_V) <- proper_edge_names
  colnames(tilde_V) <- proper_edge_names

  V_plus <- ginv(tilde_V)
  XB <- tilde_X %*% B
  I_beta <- t(XB) %*% V_plus %*% XB
  P <- B %*% ginv(I_beta) %*% t(XB) %*% V_plus
  colnames(P) <- proper_edge_names
  rownames(P) <- pair_names

  theta_hat <- as.vector(P %*% tilde_y)
  names(theta_hat) <- pair_names

  empty_path_df <- function() {
    data.frame(target = character(0), type = character(0), studies = character(0),
               path = character(0), weight = numeric(0), stringsAsFactors = FALSE)
  }
  empty_edge_df <- function() {
    data.frame(frag_id = integer(0), id = integer(0), study = character(0),
               from = character(0), to = character(0), weight = numeric(0),
               stringsAsFactors = FALSE)
  }
  
  get_study_pk <- function(target_pair, blk) {
    edge_names <- blk$edge_names_in_tildeP
    p_k <- as.numeric(P[target_pair, edge_names, drop = TRUE])
    names(p_k) <- blk$local_pair_names
    p_k
  }
  
  extract_study_collapsed_direct_indirect <- function(target_pair, blk, tol = 1e-12) {
    p_k <- get_study_pk(target_pair, blk)
    endpoint_edges <- .canonical_endpoint_edges(
      p_k = p_k,
      arms = blk$arms,
      target_pair = target_pair,
      study_name = blk$study,
      study_id = blk$id,
      global_treatments = trts,
      tol = tol
    )

    if (nrow(endpoint_edges) == 0L) {
      return(list(
        direct_paths = empty_path_df(),
        direct_edges = empty_edge_df(),
        indirect_edges = empty_edge_df()
      ))
    }

    direct_edges <- endpoint_edges[endpoint_edges$type == "Dir", , drop = FALSE]
    indirect_edges <- endpoint_edges[endpoint_edges$type == "Ind", , drop = FALSE]

    direct_df <- if (nrow(direct_edges) > 0L) {
      data.frame(
        target = target_pair,
        type = "Dir",
        studies = direct_edges$study,
        path = paste(direct_edges$from, direct_edges$to, sep = "->"),
        weight = direct_edges$weight,
        stringsAsFactors = FALSE
      )
    } else {
      empty_path_df()
    }

    direct_edges$type <- NULL
    indirect_edges$type <- NULL

    list(
      direct_paths = direct_df,
      direct_edges = direct_edges,
      indirect_edges = indirect_edges
    )
  }
  
  decompose_indirect_from_residuals <- function(target_pair, residual_edges, tol = 1e-12) {
    if (nrow(residual_edges) == 0) {
      return(list(indirect_paths = empty_path_df(), residual_edges = residual_edges))
    }
    
    pp <- .parse_pair(target_pair)
    source <- pp$u
    sink <- pp$v
    work <- residual_edges[residual_edges$weight > tol, , drop = FALSE]
    path_rows <- list()
    pid <- 0L
    
    repeat {
      idx_path <- .find_widest_reporting_path(work, source = source, sink = sink, tol = tol)
      if (is.null(idx_path) || length(idx_path) == 0) break
      
      seg_df <- work[idx_path, , drop = FALSE]
      w_path <- min(seg_df$weight)
      
      seg_table <- data.frame(
        segment_order = seq_len(nrow(seg_df)),
        study = seg_df$study,
        from = seg_df$from,
        to = seg_df$to,
        stringsAsFactors = FALSE
      )
      
      collapsed <- .collapse_path_segments(seg_table)
      study_str <- paste(collapsed$study, collapse = SAFE_STUDY_SEP)
      verts <- c(collapsed$from[1], collapsed$to)
      path_str <- paste(verts, collapse = "->")
      
      pid <- pid + 1L
      path_rows[[pid]] <- data.frame(
        target = target_pair,
        type = "Ind",
        studies = study_str,
        path = path_str,
        weight = w_path,
        stringsAsFactors = FALSE
      )
      
      work$weight[idx_path] <- work$weight[idx_path] - w_path
      work$weight[abs(work$weight) < tol] <- 0
      work <- work[work$weight > tol, , drop = FALSE]
    }

    if (nrow(work) > 0L) {
      stop(sprintf(
        "Indirect path extraction for %s left %d positive residual edge(s).",
        target_pair, nrow(work)
      ))
    }
    
    ind_df <- if (length(path_rows)) do.call(rbind, path_rows) else empty_path_df()
    list(indirect_paths = ind_df, residual_edges = work)
  }
  
  combine_identical_paths <- function(path_df, tol = 1e-12) {
    if (nrow(path_df) == 0) return(path_df)
    agg <- aggregate(weight ~ target + type + studies + path, data = path_df, FUN = sum)
    agg <- agg[, c("target", "type", "studies", "path", "weight")]
    agg <- agg[order(agg$type, -agg$weight, agg$path, agg$studies), , drop = FALSE]
    rownames(agg) <- NULL
    agg
  }
  
  decompose_target <- function(target_pair, tol = 1e-12) {
    direct_parts <- list()
    pooled_indirect_edges_parts <- list()
    
    for (k in seq_along(blocks)) {
      blk <- blocks[[k]]
      tmp <- extract_study_collapsed_direct_indirect(target_pair, blk, tol = tol)
      if (nrow(tmp$direct_paths) > 0) direct_parts[[length(direct_parts) + 1]] <- tmp$direct_paths
      if (nrow(tmp$indirect_edges) > 0) pooled_indirect_edges_parts[[length(pooled_indirect_edges_parts) + 1]] <- tmp$indirect_edges
    }
    
    direct_df <- if (length(direct_parts)) do.call(rbind, direct_parts) else empty_path_df()
    pooled_indirect_edges <- if (length(pooled_indirect_edges_parts)) do.call(rbind, pooled_indirect_edges_parts) else empty_edge_df()
    ind_obj <- decompose_indirect_from_residuals(target_pair, pooled_indirect_edges, tol = tol)
    indirect_df <- ind_obj$indirect_paths
    all_paths <- rbind(direct_df, indirect_df)
    all_paths <- combine_identical_paths(all_paths, tol = tol)
    
    overall <- data.frame(
      target = target_pair,
      direct_weight = if (nrow(direct_df)) sum(direct_df$weight) else 0,
      indirect_weight = if (nrow(indirect_df)) sum(indirect_df$weight) else 0,
      total_path_weight = if (nrow(all_paths)) sum(all_paths$weight) else 0,
      stringsAsFactors = FALSE
    )
    
    list(
      target = target_pair,
      direct_paths = combine_identical_paths(direct_df, tol = tol),
      indirect_paths = combine_identical_paths(indirect_df, tol = tol),
      all_paths = all_paths,
      overall = overall
    )
  }
  
  decomposition <- setNames(vector("list", length(pair_names)), pair_names)
  for (i in seq_along(pair_names)) {
    if (verbose) cat(sprintf("Processing %d / %d : %s\n", i, length(pair_names), pair_names[i]))
    decomposition[[i]] <- decompose_target(pair_names[i], tol = tol)
  }
  
  overall <- do.call(rbind, lapply(decomposition, function(x) x$overall))
  
  list(
    data = dat,
    treatments = trts,
    pair_names = pair_names,
    blocks = blocks,
    P = P,
    theta_hat = theta_hat,
    safe_study_sep = SAFE_STUDY_SEP,
    decomposition = decomposition,
    overall = overall,
    model = model,
    tau2 = tau2,
    tau_method = if (model == "random" && is.null(match.call()$tau2)) tau_method else "fixed"
  )
}

# ============================================================
# 2. Generalized Cochran Q decomposition
# ============================================================
q_decomposition <- function(fit, tol = 1e-8) {
  required <- c("treatments", "blocks", "pair_names", "theta_hat")
  missing <- required[!vapply(required, function(x) !is.null(fit[[x]]), logical(1))]
  if (length(missing) > 0L) {
    stop("The fitted object is missing required component(s): ", paste(missing, collapse = ", "))
  }

  tilde_y <- unlist(lapply(fit$blocks, `[[`, "y"), use.names = FALSE)
  tilde_X <- do.call(rbind, lapply(fit$blocks, `[[`, "X"))
  tilde_V <- as.matrix(Matrix::bdiag(lapply(fit$blocks, `[[`, "V")))
  row_names <- unlist(lapply(fit$blocks, `[[`, "row_names"), use.names = FALSE)
  rownames(tilde_X) <- row_names
  rownames(tilde_V) <- row_names
  colnames(tilde_V) <- row_names

  V_plus <- MASS::ginv(tilde_V)

  theta_hat_vec <- as.numeric(fit$theta_hat[fit$pair_names])
  yhat_C <- as.vector(tilde_X %*% theta_hat_vec)

  design_key <- vapply(
    fit$blocks,
    function(blk) paste(sort(as.character(blk$arms)), collapse = "|||"),
    character(1)
  )
  design_levels <- unique(design_key)
  first_block <- match(design_levels, design_key)
  design_arms <- lapply(
    first_block,
    function(k) sort(as.character(fit$blocks[[k]]$arms))
  )
  p_by_design <- vapply(design_arms, function(x) length(x) - 1L, integer(1))

  col_end <- cumsum(p_by_design)
  col_start <- c(1L, head(col_end, -1L) + 1L)
  Z <- matrix(0, nrow = length(tilde_y), ncol = sum(p_by_design))

  row0 <- 0L
  for (k in seq_along(fit$blocks)) {
    blk <- fit$blocks[[k]]
    arms <- sort(as.character(blk$arms))
    ref <- arms[1]
    nonref <- arms[-1]
    local_pairs <- blk$local_pair_names

    Zk <- matrix(0, nrow = length(local_pairs), ncol = length(nonref))
    for (i in seq_along(local_pairs)) {
      uv <- strsplit(local_pairs[i], ":", fixed = TRUE)[[1]]
      u <- uv[1]
      v <- uv[2]
      if (u != ref) Zk[i, match(u, nonref)] <- -1
      if (v != ref) Zk[i, match(v, nonref)] <- 1
    }

    d <- match(design_key[k], design_levels)
    nr <- nrow(Zk)
    rows <- row0 + seq_len(nr)
    cols <- col_start[d]:col_end[d]
    Z[rows, cols] <- Zk
    row0 <- row0 + nr
  }

  p_D <- as.integer(Matrix::rankMatrix(Z)[1])
  p_D_expected <- sum(p_by_design)
  if (p_D != p_D_expected) {
    warning(sprintf(
      "The design-specific model has rank %d, whereas sum_d(t_d-1) = %d; degrees of freedom use rank(Z).",
      p_D, p_D_expected
    ))
  }

  I_D <- t(Z) %*% V_plus %*% Z
  gamma_D <- MASS::ginv(I_D) %*% t(Z) %*% V_plus %*% tilde_y
  yhat_D <- as.vector(Z %*% gamma_D)

  quad <- function(x) as.numeric(t(x) %*% V_plus %*% x)

  Q_net <- quad(tilde_y - yhat_C)
  Q_within <- quad(tilde_y - yhat_D)
  Q_inc_geom <- quad(yhat_D - yhat_C)

  scale_Q <- max(1, abs(Q_net), abs(Q_within), abs(Q_inc_geom))
  geom_error <- Q_net - Q_within - Q_inc_geom
  if (abs(geom_error) > tol * scale_Q) {
    stop(sprintf(
      "Q decomposition failed numerically: Q_net - Q_within - Q_inc = %.3e",
      geom_error
    ))
  }

  Q_inc <- Q_net - Q_within
  if (Q_inc < 0 && abs(Q_inc) <= tol * scale_Q) Q_inc <- 0
  if (Q_inc < 0) stop("Computed inconsistency component is negative beyond numerical tolerance.")

  rank_V <- as.integer(Matrix::rankMatrix(tilde_V)[1])
  T_n <- length(fit$treatments)
  df_net <- rank_V - (T_n - 1L)
  df_within <- rank_V - p_D
  df_inc <- p_D - (T_n - 1L)

  if (any(c(df_net, df_within, df_inc) < 0L)) {
    stop(sprintf(
      "Negative Q degrees of freedom: net=%d, within=%d, inc=%d",
      df_net, df_within, df_inc
    ))
  }
  if (df_net != df_within + df_inc) {
    stop("Q degrees of freedom do not satisfy df_net = df_within + df_inc.")
  }

  is_random <- !is.null(fit$model) && identical(fit$model, "random")

  if (is_random) {
    return(data.frame(
      Treatments = T_n,
      Designs = length(design_levels),
      rank_V = rank_V,
      p_D = p_D,
      tau2 = if (!is.null(fit$tau2)) as.numeric(fit$tau2)[1] else NA_real_,
      Q_net_RE = Q_net,
      Q_error = Q_within,
      Q_inc_RE = Q_inc,
      stringsAsFactors = FALSE
    ))
  }

  q_pvalue <- function(Q, df) {
    if (df > 0L) stats::pchisq(Q, df = df, lower.tail = FALSE) else NA_real_
  }

  data.frame(
    Treatments = T_n,
    Designs = length(design_levels),
    rank_V = rank_V,
    p_D = p_D,
    Q_net = Q_net,
    df_net = df_net,
    p_net = q_pvalue(Q_net, df_net),
    Q_het = Q_within,
    df_het = df_within,
    p_het = q_pvalue(Q_within, df_within),
    Q_inc = Q_inc,
    df_inc = df_inc,
    p_inc = q_pvalue(Q_inc, df_inc),
    stringsAsFactors = FALSE
  )
}

# ============================================================
# 3. Contrast decomposition table with estimates and weights
# ============================================================
contrast_decomposition_table <- function(fit, target, conf_level = 0.95) {
  obj <- fit$decomposition[[target]]
  if (is.null(obj)) stop("Target not found: ", target)
  
  overall_rows <- .compute_overall_rows(fit, target, conf_level = conf_level)
  
  direct_rows <- NULL
  if (!is.null(obj$direct_paths) && nrow(obj$direct_paths) > 0) {
    direct_rows <- lapply(seq_len(nrow(obj$direct_paths)), function(i) {
      st <- .compute_direct_study_stats(
        res = fit,
        target = target,
        study_name = obj$direct_paths$studies[i],
        weight = obj$direct_paths$weight[i],
        conf_level = conf_level
      )
      data.frame(
        section = "Direct",
        label = obj$direct_paths$studies[i],
        path = obj$direct_paths$path[i],
        weight = st$weight,
        estimate = st$estimate,
        se = st$se,
        lower = st$lower,
        upper = st$upper,
        stringsAsFactors = FALSE
      )
    })
    direct_rows <- bind_rows(direct_rows)
  }
  
  indirect_rows <- NULL
  if (!is.null(obj$indirect_paths) && nrow(obj$indirect_paths) > 0) {
    indirect_rows <- lapply(seq_len(nrow(obj$indirect_paths)), function(i) {
      st <- .compute_path_stats(
        res = fit,
        studies_str = obj$indirect_paths$studies[i],
        path_str = obj$indirect_paths$path[i],
        weight = obj$indirect_paths$weight[i],
        conf_level = conf_level
      )
      data.frame(
        section = "Indirect",
        label = obj$indirect_paths$studies[i],
        path = obj$indirect_paths$path[i],
        weight = st$weight,
        estimate = st$estimate,
        se = st$se,
        lower = st$lower,
        upper = st$upper,
        stringsAsFactors = FALSE
      )
    })
    indirect_rows <- bind_rows(indirect_rows)
  }
  
  summary_rows <- data.frame(
    section = "Overall",
    label = overall_rows$label,
    path = NA_character_,
    weight = overall_rows$weight,
    estimate = overall_rows$estimate,
    se = overall_rows$se,
    lower = overall_rows$lower,
    upper = overall_rows$upper,
    stringsAsFactors = FALSE
  )
  
  bind_rows(summary_rows, direct_rows, indirect_rows)
}

# 4. Forest plot
# ============================================================
plot_csp_forest <- function(fit, target, show_indirect_paths = FALSE,
                            point_size_range = c(2, 6),
                            title = NULL,
                            study_label_max_chars = 35) {
  obj <- fit$decomposition[[target]]
  if (is.null(obj)) stop("Target not found: ", target)
  
  # 1. Prepare Summary Rows
  summary_df <- .compute_overall_rows(fit, target)
  summary_df$component <- "Summary"
  
  # 2. Prepare Direct Study Rows
  direct_rows <- NULL
  if (!is.null(obj$direct_paths) && nrow(obj$direct_paths) > 0) {
    direct_rows <- lapply(seq_len(nrow(obj$direct_paths)), function(i) {
      st <- .compute_direct_study_stats(fit, target, obj$direct_paths$studies[i], obj$direct_paths$weight[i])
      data.frame(
        label = .clean_study_label(obj$direct_paths$studies[i], study_label_max_chars),
        component = "Direct study",
        estimate = st$estimate,
        se = st$se,
        lower = st$lower,
        upper = st$upper,
        weight = st$weight,
        stringsAsFactors = FALSE
      )
    })
    direct_rows <- bind_rows(direct_rows) %>% arrange(desc(weight), label)
  }
  
  # 3. Prepare Indirect Path Rows (Truncated 3-letter names with "/")
  path_rows <- NULL
  if (show_indirect_paths && !is.null(obj$indirect_paths) && nrow(obj$indirect_paths) > 0) {
    path_rows <- lapply(seq_len(nrow(obj$indirect_paths)), function(i) {
      st <- .compute_path_stats(fit, obj$indirect_paths$studies[i], obj$indirect_paths$path[i], obj$indirect_paths$weight[i])
      
      studies_vec <- .split_path_studies(obj$indirect_paths$studies[i], study_sep = fit$safe_study_sep)
      studies_short <- substr(trimws(studies_vec), 1, 3)
      studies_str <- paste(studies_short, collapse = "/")
      
      # Combined label for Y-axis (Indented studies on second line)
      combined_label <- sprintf("%s\n   %s", obj$indirect_paths$path[i], studies_str)
      
      data.frame(
        label = combined_label,
        component = "Indirect path",
        estimate = st$estimate,
        se = st$se,
        lower = st$lower,
        upper = st$upper,
        weight = st$weight,
        stringsAsFactors = FALSE
      )
    })
    path_rows <- bind_rows(path_rows) %>% arrange(desc(weight), label)
  }
  
  # 4. Combine all data
  df <- bind_rows(
    summary_df %>% filter(label == "Overall NMA"),
    summary_df %>% filter(label == "Overall direct"),
    direct_rows,
    summary_df %>% filter(label == "Overall indirect"),
    path_rows
  ) %>%
    mutate(
      label = as.character(label),
      proportion = weight,
      est_lab = ifelse(is.na(estimate), "", sprintf("%.3f", estimate)),
      wt_lab = ifelse(is.na(proportion), "", sprintf("%.3f", proportion))
    )
  
  # Create unique labels and add a dummy Header level for the top of the plot
  df$label_plot <- .make_display_labels_unique(df$label)
  unique_labels <- rev(df$label_plot)
  plot_levels <- c(unique_labels, "HEADER_ROW")
  df$label_plot <- factor(df$label_plot, levels = plot_levels)
  
  summary_part <- df %>% filter(component == "Summary")
  xmin <- min(df$lower, na.rm = TRUE)
  xmax <- max(df$upper, na.rm = TRUE)
  xrng <- xmax - xmin
  if (!is.finite(xrng) || xrng <= 0) xrng <- 1
  
  # Column X-axis positions
  x_left_text <- xmin - 0.10 * xrng
  x_right_text <- xmax + 0.18 * xrng
  
  if (is.null(title)) title <- paste("Forest plot for", target)
  
  # 5. Build Plot
  ggplot(df, aes(x = estimate, y = label_plot)) +
    geom_vline(data = summary_part %>% filter(label == "Overall NMA"),
               aes(xintercept = estimate), linetype = 2, color = "red", linewidth = 0.5) +
    geom_errorbar(aes(xmin = lower, xmax = upper), orientation = "y", width = 0.2, linewidth = 0.4, na.rm = TRUE) +
    geom_point(data = df %>% filter(component != "Summary"), aes(size = proportion),
               shape = 21, fill = "white", color = "black", stroke = 0.5, na.rm = TRUE) +
    geom_point(data = df %>% filter(component == "Summary"), aes(size = proportion),
               shape = 22, fill = "#b3d3df", color = "black", stroke = 0.6, na.rm = TRUE) +
    
    # Column Headers (Estimate and Weight)
    annotate("text", x = x_left_text, y = "HEADER_ROW", label = "Estimate", 
             hjust = 1, size = 3.5, fontface = "bold") +
    annotate("text", x = x_right_text, y = "HEADER_ROW", label = "Weight", 
             hjust = 0, size = 3.5, fontface = "bold") +
    
    # Numeric Data Labels
    geom_text(aes(x = x_left_text, label = est_lab), hjust = 1, size = 3) +
    geom_text(aes(x = x_right_text, label = wt_lab), hjust = 0, size = 3) +
    
    scale_size_continuous(range = point_size_range, limits = c(0, 1)) +
    # Suppress the dummy HEADER_ROW label on the actual Y-axis
    scale_y_discrete(labels = function(x) ifelse(x == "HEADER_ROW", "", .strip_unique_suffix(x))) +
    coord_cartesian(xlim = c(xmin - 0.25 * xrng, xmax + 0.35 * xrng), clip = "off") +
    theme_bw(base_size = 11) +
    theme(
      legend.position = "none",
      axis.text.y = element_text(size = 9, hjust = 0, lineheight = 0.9), 
      plot.margin = margin(10, 80, 10, 10), 
      panel.grid.minor = element_blank()
    ) +
    labs(
      title = title,
      subtitle = "95% CIs shown; red dashed line = NMA estimate",
      x = "Estimate",
      y = NULL
    )
}
# ============================================================
# 5. Tension plot
# ============================================================
plot_csp_tension <- function(fit, baseline = "A",
                             title = "Direct / indirect / network comparisons") {
  compute_tension_one <- function(target_name, tol = 1e-12) {
    y_all <- unlist(lapply(fit$blocks, `[[`, "y"))
    row_names <- unlist(lapply(fit$blocks, `[[`, "row_names"))
    names(y_all) <- row_names
    V_all <- as.matrix(Matrix::bdiag(lapply(fit$blocks, `[[`, "V")))
    rownames(V_all) <- row_names
    colnames(V_all) <- row_names
    
    p_all <- fit$P[target_name, , drop = TRUE]
    q_dir <- .build_direct_coef_vector(fit, target_name)
    q_ind <- p_all - q_dir
    obj <- fit$decomposition[[target_name]]
    w_dir <- obj$overall$direct_weight
    w_ind <- obj$overall$indirect_weight
    
    C_dir <- sum(q_dir * y_all)
    C_ind <- sum(q_ind * y_all)
    C_net <- sum(p_all * y_all)
    
    theta_dir <- if (abs(w_dir) > tol) C_dir / w_dir else NA_real_
    theta_ind <- if (abs(w_ind) > tol) C_ind / w_ind else NA_real_
    theta_net <- C_net
    
    var_dir <- if (abs(w_dir) > tol) as.numeric(t(q_dir) %*% V_all %*% q_dir) / (w_dir^2) else NA_real_
    var_ind <- if (abs(w_ind) > tol) as.numeric(t(q_ind) %*% V_all %*% q_ind) / (w_ind^2) else NA_real_
    var_net <- as.numeric(t(p_all) %*% V_all %*% p_all)
    
    out <- data.frame(
      target = target_name,
      component = c("Direct", "Indirect", "Network"),
      estimate = c(theta_dir, theta_ind, theta_net),
      se = sqrt(c(var_dir, var_ind, var_net)),
      weight = c(w_dir, w_ind, 1),
      stringsAsFactors = FALSE
    )
    out$lower <- out$estimate - 1.96 * out$se
    out$upper <- out$estimate + 1.96 * out$se
    out
  }
  
  keep <- grep(paste0("^", baseline, ":"), fit$pair_names, value = TRUE)
  tension_df <- bind_rows(lapply(keep, compute_tension_one)) %>%
    mutate(
      target = factor(target, levels = rev(unique(target))),
      component = factor(component, levels = c("Direct", "Indirect", "Network"))
    )
  
  dodge <- position_dodge(width = 0.55)
  
  ggplot(tension_df, aes(x = estimate, y = target, colour = component, shape = component)) +
    geom_vline(xintercept = 0, linetype = "dotted", linewidth = 0.5, colour = "black") +
    geom_errorbar(
      aes(xmin = lower, xmax = upper),
      orientation = "y",
      width = 0.18,
      linewidth = 0.5,
      position = dodge,
      na.rm = TRUE
    ) +
    geom_point(
      aes(size = weight),
      position = dodge,
      na.rm = TRUE
    ) +
    scale_shape_manual(values = c("Direct" = 16, "Indirect" = 17, "Network" = 15)) +
    scale_colour_manual(values = c("Direct" = "#F8766D", "Indirect" = "#7CAE00", "Network" = "#00BFC4")) +
    scale_size_continuous(range = c(2.5, 7), limits = c(0, 1), breaks = c(0.25, 0.5, 0.75, 1.0)) +
    labs(
      title = title,
      x = "Log-odds ratio",
      y = NULL,
      colour = NULL,
      shape = NULL,
      size = "Weight"
    ) +
    theme_bw(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "top",
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(colour = "grey85"),
      panel.grid.major.x = element_line(colour = "grey85")
    )
}

# ============================================================
# 6. 3D plot
# ============================================================
plot_csp_3d <- function(fit, title = "Canonical direct-study and overall indirect weights") {
  direct_tbl <- do.call(
    rbind,
    lapply(fit$decomposition, function(obj) {
      dd <- obj$direct_paths
      if (is.null(dd) || nrow(dd) == 0) return(NULL)
      data.frame(
        target = dd$target,
        source = dd$studies,
        type = "Direct",
        weight = as.numeric(dd$weight),
        stringsAsFactors = FALSE
      )
    })
  )
  
  indirect_tbl <- fit$overall %>%
    transmute(
      target = target,
      source = "Indirect",
      type = "Indirect",
      weight = as.numeric(indirect_weight)
    )
  
  fig_dat <- bind_rows(direct_tbl, indirect_tbl) %>%
    mutate(target = factor(target, levels = fit$pair_names))
  
  target_levels <- levels(fig_dat$target)
  source_levels <- c(sort(unique(fig_dat$source[fig_dat$source != "Indirect"])), "Indirect")
  
  fig_dat <- fig_dat %>%
    mutate(
      x = match(target, target_levels),
      y = match(source, source_levels)
    )
  
  p <- plot_ly()
  for (i in seq_len(nrow(fig_dat))) {
    clr <- if (fig_dat$type[i] == "Direct") "#1f77b4" else "#ff7f0e"
    p <- add_trace(
      p,
      type = "scatter3d",
      mode = "lines+markers",
      x = c(fig_dat$x[i], fig_dat$x[i]),
      y = c(fig_dat$y[i], fig_dat$y[i]),
      z = c(0, fig_dat$weight[i]),
      line = list(width = 16, color = clr),
      marker = list(size = 3, color = clr),
      hovertemplate = paste0(
        "Target: ", fig_dat$target[i], "<br>",
        "Source: ", fig_dat$source[i], "<br>",
        "Type: ", fig_dat$type[i], "<br>",
        "Weight: ", sprintf("%.4f", fig_dat$weight[i]),
        "<extra></extra>"
      ),
      showlegend = FALSE
    )
  }
  
  layout(
    p,
    scene = list(
      xaxis = list(title = "Target comparison", tickvals = seq_along(target_levels), ticktext = target_levels),
      yaxis = list(title = "Study / indirect", tickvals = seq_along(source_levels), ticktext = source_levels),
      zaxis = list(title = "Weight", range = c(0, max(fig_dat$weight) * 1.05)),
      camera = list(eye = list(x = 1.7, y = 1.6, z = 1.2))
    ),
    title = list(text = title)
  )
}

###############################################################################
##  SA-ASPC : complete replication pipeline  (v2, locked configuration)
##
##  Changes from v1, both required by the manuscript:
##    (1) Algorithm 1' uses the weight 2*beta*rho_bar, not beta*rho_bar.
##        Proposition 2.2 gives 2*rho_bar = d Pen / d W, so only this factor
##        makes the surrogate gradient equal the gradient of J_beta and the
##        descent argument of Corollary 2.3 valid.  This CHANGES the fits.
##    (2) The stability coordinate of criterion (12) uses the anchored map
##        z4 = (ST + 1)/2 on the adjusted-Rand functional, per Section 3.3,
##        not min-max.  The min-max variant is still computed and reported
##        as the sensitivity analysis of Section 5.5.
##
##  Run top to bottom.  Set RUN_* in BLOCK 0 to skip expensive blocks.
###############################################################################

## ===========================================================================
## BLOCK 0 : CONFIGURATION
## ===========================================================================

PATH_PANEL  <- "C:/Users/lenovo/Documents/data/panel_jateng.csv"
PATH_SHP    <- "C:/Users/lenovo/Downloads/3300.shp"
PATH_W <- "C:/Users/lenovo/Downloads/output/Wq_frozen.csv"
DIR_OUT     <- "output"
dir.create(DIR_OUT, showWarnings = FALSE, recursive = TRUE)

COL_REGION  <- "WADMKK"
COL_TIME    <- "TAHUN"
SHP_NAMECOL <- "KABKOT"
INDICATORS  <- c("APM_SD", "APM_SMP", "PRODPADI", "RLS", "SANITASI")

C_STAR      <- 8
M_FUZZ      <- 2.0
BETA_GRID   <- seq(0, 6, by = 0.25)
B_BANK      <- 100
SIGMA_PERT  <- 0.05
N_RESTART   <- 20
MAX_ITER    <- 200
TOL         <- 1e-7
T_MIN       <- 1e-6
DTW_NORM    <- TRUE
USE_ALG1P   <- TRUE

## criterion (12): weights on (FM, FS, XB, stability)
W_OBJ       <- c(FM = 0.35, FS = 0.25, XB = 0.20, ST = 0.20)
## "anchored" = z4 = (ST+1)/2 on the ARI functional      <- manuscript default
## "minmax"   = min-max on NDC                           <- Section 5.5 variant
STAB_MAP    <- "anchored"
R_ARI       <- 2            # intrinsic range of the adjusted Rand index

MASTER_SEED <- 20260905

RUN_CSELECT   <- TRUE
RUN_CALIB     <- TRUE
RUN_RULES     <- TRUE
RUN_WEIGHTS   <- TRUE
RUN_MCCERT    <- TRUE
RUN_BASELINE  <- TRUE
RUN_PARTITION <- TRUE
RUN_DRIFT     <- TRUE
RUN_CONVERGE  <- TRUE
RUN_RESTART   <- TRUE
RUN_RUNTIME   <- TRUE
RUN_DTWNORM   <- TRUE
RUN_REFINE    <- TRUE
RUN_FACTORIAL <- TRUE

set.seed(MASTER_SEED)
options(stringsAsFactors = FALSE)

for (p in c("cluster", "clue", "sf", "spdep", "parallel"))
  if (!requireNamespace(p, quietly = TRUE))
    stop(sprintf("install.packages('%s')", p))

###############################################################################
## BLOCK 1 : DATA AND SPATIAL WEIGHTS
###############################################################################

load_panel <- function(path, col_region, col_time, indicators) {
  raw <- utils::read.csv(path)
  stopifnot(all(c(col_region, col_time, indicators) %in% names(raw)))
  regions <- sort(unique(raw[[col_region]]))
  times   <- sort(unique(raw[[col_time]]))
  I <- length(regions); TT <- length(times); P <- length(indicators)
  if (any(table(raw[[col_region]]) != TT))
    stop("unbalanced panel: some regions do not have ", TT, " periods")
  X <- array(NA_real_, c(I, TT, P),
             dimnames = list(regions, as.character(times), indicators))
  for (i in seq_len(I)) {
    sub <- raw[raw[[col_region]] == regions[i], ]
    sub <- sub[order(sub[[col_time]]), ]
    X[i, , ] <- as.matrix(sub[, indicators, drop = FALSE])
  }
  if (anyNA(X)) stop("missing values in the panel")
  X
}

zscore_panel <- function(X) {
  for (p in seq_len(dim(X)[3])) {
    v <- X[, , p]; X[, , p] <- (v - mean(v)) / stats::sd(v)
  }
  X
}

build_W <- function(regions) {
  if (!is.na(PATH_W)) {
    W <- as.matrix(utils::read.csv(PATH_W, row.names = 1, check.names = FALSE))[regions, regions]
  } else {
    sf::sf_use_s2(FALSE)
    shp <- sf::st_read(PATH_SHP, quiet = TRUE)
    shp <- sf::st_make_valid(shp)              # NOT st_buffer(shp, 0)
    
    is_kota <- shp$KABKOTNO >= 71               # konvensi BPS: 71+ = kota, membedakan
    # dari kabupaten bernama sama (Semarang, dst.)
    nm  <- toupper(trimws(shp[[SHP_NAMECOL]]))
    nm[is_kota] <- paste("KOTA", nm[is_kota])
    
    rs  <- toupper(trimws(regions))
    rs  <- gsub("^KABUPATEN ", "", rs)          # sengaja TIDAK melucuti "^KOTA " lagi
    
    idx <- match(rs, nm)
    if (any(is.na(idx)))
      stop("not found in shapefile: ", paste(regions[is.na(idx)], collapse = ", "))
    if (any(duplicated(idx)))
      stop("duplicate shapefile match for: ",
           paste(regions[duplicated(idx) | duplicated(idx, fromLast = TRUE)], collapse = ", "))
    
    nb <- spdep::poly2nb(shp[idx, ], queen = TRUE)
    W  <- spdep::nb2mat(nb, style = "B", zero.policy = TRUE)
  }
  diag(W) <- 0; storage.mode(W) <- "double"
  dimnames(W) <- list(regions, regions)
  W
}

row_standardise <- function(W) {
  rs <- rowSums(W); z <- rs == 0; rs[z] <- 1
  P <- W / rs; P[z, ] <- 0; P
}

###############################################################################
## BLOCK 2 : MULTIVARIATE DTW, Eq. (1)
###############################################################################

dtw_pair <- function(A, B, normalise = TRUE) {
  n <- nrow(A); m <- nrow(B)
  Cst <- matrix(0, n, m)
  for (a in seq_len(n)) {
    dd <- sweep(B, 2, A[a, ], "-")
    Cst[a, ] <- sqrt(.rowSums(dd * dd, m, ncol(B)))
  }
  Acc <- matrix(Inf, n + 1, m + 1); Len <- matrix(Inf, n + 1, m + 1)
  Acc[1, 1] <- 0; Len[1, 1] <- 0
  for (a in seq_len(n)) for (b in seq_len(m)) {
    cand  <- c(Acc[a, b + 1], Acc[a + 1, b], Acc[a, b])
    lcand <- c(Len[a, b + 1], Len[a + 1, b], Len[a, b])
    k <- which.min(cand)
    Acc[a + 1, b + 1] <- Cst[a, b] + cand[k]
    Len[a + 1, b + 1] <- 1 + lcand[k]
  }
  if (normalise) Acc[n + 1, m + 1] / Len[n + 1, m + 1] else Acc[n + 1, m + 1]
}

dtw_matrix <- function(X, normalise = TRUE) {
  I <- dim(X)[1]
  traj <- lapply(seq_len(I), function(i) X[i, , , drop = TRUE])
  D <- matrix(0, I, I, dimnames = list(dimnames(X)[[1]], dimnames(X)[[1]]))
  for (i in seq_len(I - 1)) for (j in (i + 1):I) {
    v <- dtw_pair(traj[[i]], traj[[j]], normalise)
    D[i, j] <- v; D[j, i] <- v
  }
  D
}

###############################################################################
## BLOCK 3 : MODEL, Eqs. (2)-(5), Algorithms 1 and 1'
###############################################################################

## rho_ic(W) = sum_k p_ik (S_k - w_kc),  Eq. (3)
rho_mat <- function(Wm, Pm) {
  S <- rowSums(Wm)
  matrix(as.numeric(Pm %*% S), nrow(Wm), ncol(Wm)) - Pm %*% Wm
}

## Pen(W) = sum_i sum_k p_ik ( S_i S_k - <W_i, W_k> ),  Eq. (2a)
penalty_value <- function(Wm, Pm) {
  S <- rowSums(Wm)
  sum(Pm * (outer(S, S) - Wm %*% t(Wm)))
}

J_beta <- function(U, med, D, Pm, beta, m) {
  Wm <- U^m
  sum(Wm * D[, med, drop = FALSE]) + beta * penalty_value(Wm, Pm)
}

## Eq. (4)
membership_update <- function(dtil, m) {
  pow <- 1 / (m - 1)
  t(apply(dtil, 1, function(di) {
    z <- which(di <= .Machine$double.eps)
    if (length(z)) { u <- numeric(length(di)); u[z] <- 1 / length(z); u }
    else { r <- (1 / di)^pow; r / sum(r) }
  }))
}

## Eq. (5)
medoid_update <- function(U, D, m) {
  Wm <- U^m; chosen <- integer(0)
  for (cc in seq_len(ncol(U))) {
    cost <- as.numeric(crossprod(D, Wm[, cc]))
    cost[chosen] <- Inf
    chosen <- c(chosen, which.min(cost))
  }
  chosen
}

## one fit from one start.
##   alg1p = TRUE  -> Algorithm 1', symmetrised proximity AND the factor 2
##   alg1p = FALSE -> Algorithm 1 as published, weight rho with raw P
fcmd_stt_fit1 <- function(D, Pm, beta, C, m, U0, alg1p = USE_ALG1P,
                          max_iter = MAX_ITER, tol = TOL) {
  Pbar <- 0.5 * (Pm + t(Pm))
  U <- U0
  q <- medoid_update(U, D, m)
  Jold <- J_beta(U, q, D, Pm, beta, m)
  n_up <- 0L; n_it <- 0L; trace <- Jold
  for (it in seq_len(max_iter)) {
    n_it <- n_it + 1L
    if (alg1p) {
      ## factor 2 with Pbar: grad(surrogate) == grad(J_beta), Eq. (2c)
      rho  <- rho_mat(U^m, Pbar)
      dtil <- pmax(D[, q, drop = FALSE] + 2 * beta * rho, 0)
    } else {
      rho  <- rho_mat(U^m, Pm)
      dtil <- pmax(D[, q, drop = FALSE] + beta * rho, 0)
    }
    Ustar <- membership_update(dtil, m)
    if (alg1p) {
      tt <- 1
      repeat {
        Ucand <- (1 - tt) * U + tt * Ustar
        if (J_beta(Ucand, q, D, Pm, beta, m) <= Jold || tt <= T_MIN) break
        tt <- tt / 2
      }
      if (tt <= T_MIN) { trace <- c(trace, Jold); break }
      U <- Ucand
    } else U <- Ustar
    q <- medoid_update(U, D, m)
    Jnew <- J_beta(U, q, D, Pm, beta, m)
    if (Jnew > Jold + 1e-12) n_up <- n_up + 1L
    trace <- c(trace, Jnew)
    if (abs(Jold - Jnew) < tol * max(1, abs(Jold))) { Jold <- Jnew; break }
    Jold <- Jnew
  }
  list(U = U, medoids = q, J = Jold, iters = n_it,
       n_increase = n_up, trace = trace)
}

fcmd_stt <- function(D, Pm, beta, C, m = M_FUZZ, n_restart = N_RESTART,
                     alg1p = USE_ALG1P, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  I <- nrow(D); best <- NULL
  allJ <- numeric(n_restart); allInc <- integer(n_restart); allIt <- integer(n_restart)
  for (r in seq_len(n_restart)) {
    U0 <- matrix(stats::runif(I * C), I, C); U0 <- U0 / rowSums(U0)
    f <- fcmd_stt_fit1(D, Pm, beta, C, m, U0, alg1p = alg1p)
    allJ[r] <- f$J; allInc[r] <- f$n_increase; allIt[r] <- f$iters
    if (is.null(best) || f$J < best$J) best <- f
  }
  best$J_restarts <- allJ; best$inc_restarts <- allInc; best$it_restarts <- allIt
  best
}

###############################################################################
## BLOCK 4 : FUNCTIONALS, Eqs. (6)-(10)
###############################################################################

hard <- function(U) max.col(U, ties.method = "first")

fuzzy_moran <- function(U, Z, Pm, m) {
  Wm <- U^m; Zc <- scale(Z, TRUE, FALSE); dn <- rowSums(Zc * Zc)
  C <- ncol(U); FMc <- Mc <- numeric(C)
  for (cc in seq_len(C)) {
    w <- Wm[, cc]; S0 <- sum(outer(w, w) * Pm); Mc[cc] <- sum(w)
    if (S0 <= 0) { FMc[cc] <- 0; next }
    WZ <- Zc * w
    FMc[cc] <- (sum(w) / S0) * sum(Pm * (WZ %*% t(WZ))) / sum(w * dn)
  }
  sum(Mc * FMc) / sum(Mc)
}

spatial_coherence <- function(U, Pm) sum(Pm * (U %*% t(U))) / sum(Pm)

fuzzy_silhouette <- function(U, D, alpha = 1) {
  lab <- hard(U); I <- nrow(D); s <- numeric(I)
  for (i in seq_len(I)) {
    own <- setdiff(which(lab == lab[i]), i)
    a <- if (!length(own)) 0 else mean(D[i, own])
    oth <- setdiff(unique(lab), lab[i])
    b <- if (!length(oth)) 0 else min(vapply(oth, function(g) mean(D[i, lab == g]), 0))
    s[i] <- if (max(a, b) == 0) 0 else (b - a) / max(a, b)
  }
  Us <- t(apply(U, 1, sort, decreasing = TRUE))
  w <- (Us[, 1] - Us[, 2])^alpha
  sum(w * s) / sum(w)
}

xie_beni <- function(U, med, D, m) {
  dq <- D[med, med]; diag(dq) <- Inf
  sum((U^m) * (D[, med, drop = FALSE]^2)) / (nrow(U) * min(dq)^2)
}

ari <- function(x, y) {
  tb <- table(x, y); n <- sum(tb); cn2 <- function(z) sum(z * (z - 1) / 2)
  idx <- cn2(as.vector(tb)); ea <- cn2(rowSums(tb)); eb <- cn2(colSums(tb))
  ex <- ea * eb / cn2(n); mx <- (ea + eb) / 2
  if (mx == ex) 1 else (idx - ex) / (mx - ex)
}

E_overlap <- function(U) {
  I <- nrow(U); E <- matrix(0, I, I)
  for (i in seq_len(I - 1)) for (k in (i + 1):I) {
    v <- 1 - 0.5 * sum(abs(U[i, ] - U[k, ])); E[i, k] <- v; E[k, i] <- v
  }
  E
}
ndc <- function(U1, U2) {
  E1 <- E_overlap(U1); E2 <- E_overlap(U2); I <- nrow(U1)
  1 - sum(abs(E1[upper.tri(E1)] - E2[upper.tri(E2)])) / (I * (I - 1) / 2)
}

mean_entropy <- function(U) { u <- pmax(U, 1e-12); mean(-rowSums(u * log(u))) }

fuzzy_dunn <- function(U, med, D) {
  lab <- hard(U); dq <- D[med, med]; diag(dq) <- Inf
  diam <- vapply(unique(lab), function(g) {
    ii <- which(lab == g); if (length(ii) < 2) 0 else max(D[ii, ii]) }, 0)
  min(dq) / max(max(diam), .Machine$double.eps)
}

## number of units reassigned, after optimal label matching (Hungarian)
moved_matched <- function(U1, U2) {
  l1 <- hard(U1); l2 <- hard(U2)
  tb <- table(factor(l1, levels = seq_len(ncol(U1))),
              factor(l2, levels = seq_len(ncol(U2))))
  perm <- clue::solve_LSAP(tb, maximum = TRUE)
  inv <- integer(ncol(U2)); inv[as.integer(perm)] <- seq_along(perm)
  sum(l1 != inv[l2])
}

###############################################################################
## BLOCK 5 : PERTURBATION BANK, Proposition 3.4
###############################################################################

build_bank <- function(X, B = B_BANK, sigma = SIGMA_PERT,
                       normalise = DTW_NORM, seed = 12345) {
  set.seed(seed)
  lapply(seq_len(B), function(b) {
    Xb <- X + array(stats::rnorm(prod(dim(X)), 0, sigma), dim = dim(X))
    dtw_matrix(Xb, normalise)
  })
}

stability_at_beta <- function(bank, Pm, beta, C, m, Uref, n_restart = 3) {
  lab <- hard(Uref); B <- length(bank)
  a <- numeric(B); d <- numeric(B)
  for (b in seq_len(B)) {
    fb <- fcmd_stt(bank[[b]], Pm, beta, C, m, n_restart = n_restart)
    a[b] <- ari(lab, hard(fb$U)); d[b] <- ndc(Uref, fb$U)
  }
  list(ST = mean(a), ST_sd = stats::sd(a), NDC = mean(d), NDC_sd = stats::sd(d),
       ST_lo = mean(a) - 1.96 * stats::sd(a) / sqrt(B),
       ST_hi = mean(a) + 1.96 * stats::sd(a) / sqrt(B))
}

###############################################################################
## BLOCK 6 : RESCALING AND THE CRITERION, Eq. (12)
###############################################################################

mm01 <- function(x) { r <- range(x); if (diff(r) == 0) rep(1, length(x)) else (x - r[1]) / diff(r) }

## stab = "anchored" (Section 3.3 default) or "minmax" (Section 5.5 variant)
rescale_objectives <- function(tb, stab = STAB_MAP) {
  z4 <- switch(stab,
    anchored = (tb$ST + 1) / 2,          # anchored to the ARI range, R_A = 2
    minmax   = mm01(tb$NDC),             # min-max on the fuzzy concordance
    stop("unknown stab map"))
  cbind(FM = mm01(tb$FM), FS = mm01(tb$FS), XB = 1 - mm01(tb$XB), ST = z4)
}

d_ideal  <- function(Zn, w = W_OBJ) sqrt(as.numeric(((1 - Zn)^2) %*% w))
score_ws <- function(Zn, w = W_OBJ) as.numeric(Zn %*% w)

pareto_flag <- function(Zn) {
  vapply(seq_len(nrow(Zn)), function(i)
    !any(apply(Zn, 1, function(v) all(v >= Zn[i, ]) && any(v > Zn[i, ]))), TRUE)
}

knee_point <- function(beta, s) {
  x <- (beta - min(beta)) / diff(range(beta))
  y <- (s - min(s)) / diff(range(s))
  beta[which.max(y - x)]
}

###############################################################################
## BLOCK 7 : NUMBER OF CLUSTERS  (Table 5)
###############################################################################

select_C <- function(D, Pm, candidates = 2:10, m = M_FUZZ, n_boot = 50) {
  I <- nrow(D); res <- data.frame()
  for (Cc in candidates) {
    f <- fcmd_stt(D, Pm, 0, Cc, m, seed = 100 + Cc)
    lab <- hard(f$U); bs <- numeric(n_boot)
    for (b in seq_len(n_boot)) {
      idx <- unique(sort(sample(I, I, replace = TRUE)))
      if (length(idx) <= Cc) { bs[b] <- NA; next }
      Pb <- row_standardise((Pm[idx, idx, drop = FALSE] > 0) * 1)
      fb <- fcmd_stt(D[idx, idx, drop = FALSE], Pb, 0, Cc, m, n_restart = 3)
      bs[b] <- ari(lab[idx], hard(fb$U))
    }
    res <- rbind(res, data.frame(C = Cc, XB = xie_beni(f$U, f$medoids, D, m),
                                 FS = fuzzy_silhouette(f$U, D),
                                 Dunn = fuzzy_dunn(f$U, f$medoids, D),
                                 BootARI = mean(bs, na.rm = TRUE)))
  }
  res$Rank_XB <- rank(res$XB); res$Rank_FS <- rank(-res$FS)
  res$Rank_Dunn <- rank(-res$Dunn); res$Rank_Boot <- rank(-res$BootARI)
  res$TotalRank <- res$Rank_XB + res$Rank_FS + res$Rank_Dunn + res$Rank_Boot
  res[order(res$TotalRank), ]
}

###############################################################################
## BLOCK 8 : CALIBRATION GRID  (Table 6)
###############################################################################

calibrate <- function(D, Pm, Z, bank, grid = BETA_GRID, C = C_STAR, m = M_FUZZ,
                      verbose = TRUE) {
  out <- data.frame(); fits <- list()
  for (g in seq_along(grid)) {
    beta <- grid[g]; t0 <- proc.time()[3]
    f  <- fcmd_stt(D, Pm, beta, C, m, seed = 1000 + g)
    st <- stability_at_beta(bank, Pm, beta, C, m, f$U)
    out <- rbind(out, data.frame(
      beta = beta, FM = fuzzy_moran(f$U, Z, Pm, m),
      FS = fuzzy_silhouette(f$U, D), XB = xie_beni(f$U, f$medoids, D, m),
      ST = st$ST, ST_sd = st$ST_sd, NDC = st$NDC, NDC_sd = st$NDC_sd,
      ST_lo = st$ST_lo, ST_hi = st$ST_hi,
      SC = spatial_coherence(f$U, Pm), Entropy = mean_entropy(f$U),
      J = f$J, iters = f$iters, Runtime = proc.time()[3] - t0))
    fits[[as.character(beta)]] <- f
    if (verbose) cat(sprintf("beta=%.2f FM=%.4f ST=%.4f NDC=%.4f XB=%.4f\n",
                             beta, out$FM[g], out$ST[g], out$NDC[g], out$XB[g]))
  }
  list(table = out, fits = fits)
}

###############################################################################
## BLOCK 9 : SELECTION RULES  (Table 7)
###############################################################################

selection_rules <- function(tb, Zn, w = W_OBJ) {
  d <- d_ideal(Zn, w); s <- score_ws(Zn, w)
  Dm <- sqrt(as.numeric((Zn^2) %*% w)); cc <- Dm / (d + Dm)
  Hr <- apply(Zn, 2, function(z) { p <- pmax(z / sum(z), 1e-12)
                                   -sum(p * log(p)) / log(length(p)) })
  we <- (1 - Hr) / sum(1 - Hr)
  q0 <- Zn[1, "FS"] + Zn[1, "XB"]; qq <- Zn[, "FS"] + Zn[, "XB"]
  ok <- which(qq >= 0.9 * q0)
  data.frame(
    Rule = c("Distance to ideal point (SA-ASPC)", "Weighted sum", "TOPSIS",
             "Knee point", "Entropy weighting", "Equal weights",
             "ClustGeo-style 10% rule", "Fuzzy Moran maximisation"),
    beta = c(tb$beta[which.min(d)], tb$beta[which.max(s)], tb$beta[which.max(cc)],
             knee_point(tb$beta, s), tb$beta[which.min(d_ideal(Zn, we))],
             tb$beta[which.min(d_ideal(Zn, rep(0.25, 4)))],
             if (length(ok)) max(tb$beta[ok]) else tb$beta[1],
             tb$beta[which.max(tb$FM)]))
}

###############################################################################
## BLOCK 10 : WEIGHT SENSITIVITY  (Table 8)
###############################################################################

weight_sensitivity <- function(tb, Zn, n_draw = 5000, w = W_OBJ, seed = 2026) {
  set.seed(seed)
  Wd <- matrix(stats::rexp(n_draw * 4), n_draw, 4); Wd <- Wd / rowSums(Wd)
  sel <- tb$beta[apply((1 - Zn)^2 %*% t(Wd), 2, which.min)]
  fr <- as.data.frame(table(sel) / n_draw); names(fr) <- c("beta", "freq")
  fr$beta <- as.numeric(as.character(fr$beta)); fr <- fr[order(-fr$freq), ]
  fr$d_at_reported <- d_ideal(Zn, w)[match(fr$beta, tb$beta)]
  loo <- do.call(rbind, lapply(1:4, function(j) {
    ww <- w[-j] / sum(w[-j])
    data.frame(dropped = colnames(Zn)[j],
               beta = tb$beta[which.min(sqrt(as.numeric(((1 - Zn[, -j, drop = FALSE])^2) %*% ww)))])
  }))
  d <- d_ideal(Zn, w)
  list(freq = fr, loo = loo, gamma = sort(d)[2] - sort(d)[1],
       beta_hat = tb$beta[which.min(d)],
       equal_weights = tb$beta[which.min(d_ideal(Zn, rep(0.25, 4)))])
}

###############################################################################
## BLOCK 11 : MONTE CARLO CERTIFICATION  (Section 5.5)
###############################################################################
##  Section 3.5: population condition 2*kappa < gamma is implied by the
##  OBSERVABLE condition 4*kappa < gamma_hat, with kappa = 4 max(eta_r / R_r).

mc_certification <- function(tb, Zn, w = W_OBJ, B = B_BANK, stab = STAB_MAP) {
  if (stab == "anchored") { eta <- 1.96 * tb$ST_sd / sqrt(B); R <- R_ARI }
  else                    { eta <- 1.96 * tb$NDC_sd / sqrt(B); R <- diff(range(tb$NDC)) }
  kappa <- 4 * max(eta / R)
  d <- d_ideal(Zn, w); gam <- sort(d)[2] - sort(d)[1]
  list(map = stab, R = R, max_eta = max(eta), kappa = kappa,
       four_kappa = 4 * kappa, gamma_hat = gam,
       certified = (4 * kappa < gam), ratio = 4 * kappa / gam,
       conf_set = tb$beta[d <= min(d) + 2 * kappa],
       B_required = ceiling(B * (4 * kappa / gam)^2))
}

###############################################################################
## BLOCK 12 : BASELINES AND PARTITION COMPARISON  (Table 9)
###############################################################################

fcm_euclid <- function(Xf, C, m, n_start = 20, max_iter = 200) {
  I <- nrow(Xf); best <- NULL
  for (s in seq_len(n_start)) {
    U <- matrix(stats::runif(I * C), I, C); U <- U / rowSums(U); Dm <- NULL
    for (it in seq_len(max_iter)) {
      Wm <- U^m; V <- (t(Wm) %*% Xf) / colSums(Wm)
      Dm <- pmax(as.matrix(stats::dist(rbind(Xf, V)))[seq_len(I), I + seq_len(C)], 1e-12)
      Un <- membership_update(Dm^2, m)
      if (max(abs(Un - U)) < 1e-8) { U <- Un; break }
      U <- Un
    }
    obj <- sum((U^m) * Dm^2)
    if (is.null(best) || obj < best$obj) best <- list(U = U, obj = obj)
  }
  best
}

baseline_table <- function(D, Pm, Z, X, C, m, fits, beta_hat, beta_moran) {
  Xf <- matrix(X, nrow = dim(X)[1])
  crisp <- function(l, C) { U <- matrix(0, length(l), C); U[cbind(seq_along(l), l)] <- 1; U }
  add <- function(nm, U, med) data.frame(
    Method = nm, C = C,
    XB = if (is.null(med)) NA_real_ else xie_beni(U, med, D, m),
    FS = fuzzy_silhouette(U, D), FM = fuzzy_moran(U, Z, Pm, m),
    SC = spatial_coherence(U, Pm))
  km <- stats::kmeans(Xf, C, nstart = 50)
  pm <- cluster::pam(stats::as.dist(D), k = C)
  fc <- fcm_euclid(Xf, C, m)
  fd <- fcmd_stt(as.matrix(stats::dist(Xf)), Pm, 0, C, m, seed = 77)
  f0 <- fits[["0"]]; fa <- fits[[as.character(beta_hat)]]; fm <- fits[[as.character(beta_moran)]]
  tabout <- rbind(
    add("K-means", crisp(km$cluster, C), NULL),
    add("PAM (k-medoids)", crisp(pm$clustering, C), pm$id.med),
    add("FCM", fc$U, NULL),
    add("FCMd (Euclidean)", fd$U, fd$medoids),
    add("DTW-FCMd (beta = 0)", f0$U, f0$medoids),
    add(sprintf("DTW-FCMd-STT-FM (beta = %.2f)", beta_moran), fm$U, fm$medoids),
    add(sprintf("SA-ASPC (beta = %.2f)", beta_hat), fa$U, fa$medoids))
  agree <- data.frame(
    ARI_vs_moran = ari(hard(fa$U), hard(fm$U)),
    moved_vs_moran = moved_matched(fa$U, fm$U),
    ARI_vs_zero = ari(hard(fa$U), hard(f0$U)),
    moved_vs_zero = moved_matched(fa$U, f0$U))
  list(table = tabout, agreement = agree)
}

partition_profile <- function(fit, Xr) {
  U <- fit$U; lab <- hard(U); reg <- dimnames(Xr)[[1]]
  Xb <- apply(Xr, c(1, 3), mean)
  do.call(rbind, lapply(sort(unique(lab)), function(g) {
    ii <- which(lab == g)
    data.frame(Cluster = g, Medoid = reg[fit$medoids[g]], n = length(ii),
               t(round(colMeans(Xb[ii, , drop = FALSE]), 3)),
               MeanMaxMembership = round(mean(apply(U[ii, , drop = FALSE], 1, max)), 3))
  }))
}

partition_drift <- function(fits, tb) {
  f0 <- fits[["0"]]
  do.call(rbind, lapply(as.character(tb$beta), function(b) {
    f <- fits[[b]]
    data.frame(beta = as.numeric(b), moved_vs_zero = moved_matched(f0$U, f$U),
               ARI_vs_zero = ari(hard(f0$U), hard(f$U)),
               NDC_vs_zero = ndc(f0$U, f$U),
               mean_maxmemb = mean(apply(f$U, 1, max)))
  }))
}

###############################################################################
## BLOCK 13 : ALGORITHM 1 vs 1', AND RESTART BEHAVIOUR  (Tables 10, 11)
###############################################################################

convergence_study <- function(D, Pm, C = C_STAR, m = M_FUZZ,
                              betas = seq(0, 6, by = 0.5), n_restart = 20) {
  I <- nrow(D)
  do.call(rbind, lapply(betas, function(beta) {
    i1 <- t1 <- i2 <- t2 <- 0; ab <- numeric(n_restart)
    J1 <- J2 <- numeric(n_restart)
    for (r in seq_len(n_restart)) {
      set.seed(5000 + r)
      U0 <- matrix(stats::runif(I * C), I, C); U0 <- U0 / rowSums(U0)
      f1 <- fcmd_stt_fit1(D, Pm, beta, C, m, U0, alg1p = FALSE)
      f2 <- fcmd_stt_fit1(D, Pm, beta, C, m, U0, alg1p = TRUE)
      i1 <- i1 + f1$n_increase; t1 <- t1 + f1$iters
      i2 <- i2 + f2$n_increase; t2 <- t2 + f2$iters
      ab[r] <- ari(hard(f1$U), hard(f2$U)); J1[r] <- f1$J; J2[r] <- f2$J
    }
    data.frame(beta = beta, pct_increase_alg1 = 100 * i1 / t1,
               pct_increase_alg1p = 100 * i2 / t2,
               ARI_between = mean(ab), J_alg1 = min(J1), J_alg1p = min(J2))
  }))
}

fcmd_stt_pam <- function(D, Pm, beta, C, m = M_FUZZ) {
  med0 <- cluster::pam(stats::as.dist(D), k = C)$id.med
  U0 <- matrix(1 / C, nrow(D), C)
  rho0 <- rho_mat(U0^m, 0.5 * (Pm + t(Pm)))
  U0 <- membership_update(pmax(D[, med0, drop = FALSE] + 2 * beta * rho0, 0), m)
  fcmd_stt_fit1(D, Pm, beta, C, m, U0, alg1p = TRUE)
}

restart_study <- function(D, Pm, C = C_STAR, m = M_FUZZ,
                          betas = c(0, 0.5, 1, 2, 4, 6), n_restart = 100) {
  do.call(rbind, lapply(betas, function(b) {
    f <- fcmd_stt(D, Pm, b, C, m, n_restart = n_restart, seed = 42)
    p <- fcmd_stt_pam(D, Pm, b, C, m)
    data.frame(beta = b, J_pam = p$J, J_best_random = f$J,
               J_mean = mean(f$J_restarts), J_sd = stats::sd(f$J_restarts),
               n_distinct = length(unique(round(f$J_restarts, 6))),
               pam_rank = sum(f$J_restarts < p$J) + 1,
               n_steps_median = stats::median(f$it_restarts))
  }))
}

seed_stability <- function(D, Pm, Z, bank, C = C_STAR, seeds = 1:5) {
  do.call(rbind, lapply(seeds, function(s) {
    set.seed(s)
    cl <- calibrate(D, Pm, Z, bank, C = C, verbose = FALSE)
    Zs <- rescale_objectives(cl$table)
    data.frame(seed = s, beta_aspc = cl$table$beta[which.min(d_ideal(Zs))],
               beta_moran = cl$table$beta[which.max(cl$table$FM)],
               nondom = sum(pareto_flag(Zs)))
  }))
}

###############################################################################
## BLOCK 14 : RUNTIME, NORMALISATION, GRID REFINEMENT
###############################################################################

runtime_benchmark <- function(X, Pm, D, C = C_STAR, m = M_FUZZ,
                              grid = BETA_GRID, B = 20) {
  t0 <- proc.time()[3]
  bk <- build_bank(X, B = B)
  for (b in grid) {
    f <- fcmd_stt(D, Pm, b, C, m, n_restart = 3)
    invisible(stability_at_beta(bk, Pm, b, C, m, f$U, n_restart = 1))
  }
  t_re <- proc.time()[3] - t0
  t0 <- proc.time()[3]
  for (b in grid) {
    bk2 <- build_bank(X, B = B, seed = sample.int(1e6, 1))
    f <- fcmd_stt(D, Pm, b, C, m, n_restart = 3)
    invisible(stability_at_beta(bk2, Pm, b, C, m, f$U, n_restart = 1))
  }
  data.frame(G = length(grid), B = B, sec_reuse = t_re,
             sec_naive = proc.time()[3] - t0,
             speedup = (proc.time()[3] - t0) / t_re, predicted = length(grid))
}

normalisation_comparison <- function(X, Pm, Z, C = C_STAR, m = M_FUZZ) {
  do.call(rbind, lapply(c(TRUE, FALSE), function(nm) {
    Dn <- dtw_matrix(X, nm); bk <- build_bank(X, normalise = nm)
    cl <- calibrate(Dn, Pm, Z, bk, C = C, m = m, verbose = FALSE)
    Zn <- rescale_objectives(cl$table); d <- d_ideal(Zn)
    data.frame(normalised = nm, mean_dissim = mean(Dn[upper.tri(Dn)]),
               beta_aspc = cl$table$beta[which.min(d)],
               beta_moran = cl$table$beta[which.max(cl$table$FM)],
               nondom = sum(pareto_flag(Zn)),
               SC_lo = min(cl$table$SC), SC_hi = max(cl$table$SC))
  }))
}

grid_refinement <- function(D, Pm, Z, bank, C = C_STAR, m = M_FUZZ,
                            factors = c(1, 2, 4)) {
  do.call(rbind, lapply(factors, function(fc) {
    g <- seq(0, 6, by = 0.25 / fc)
    cl <- calibrate(D, Pm, Z, bank, grid = g, C = C, m = m, verbose = FALSE)
    Zn <- rescale_objectives(cl$table); d <- d_ideal(Zn)
    data.frame(refinement = fc, n_grid = length(g),
               beta_hat = g[which.min(d)], gamma = sort(d)[2] - sort(d)[1],
               nondom_fraction = mean(pareto_flag(Zn)))
  }))
}

###############################################################################
## BLOCK 15 : FACTORIAL STUDY  (Tables 12-17)
###############################################################################

FAST_B <- 30; FAST_RESTART <- 5; FAST_REFIT <- 1
CHECKPOINT_FILE <- file.path(DIR_OUT, "factorial_checkpoint_v2.rds")

one_cell <- function(D, Pm, Z, bank, sg, Cc, mm, rep, grid) {
  n <- length(grid)
  FM <- FS <- XB <- ST <- NDC <- numeric(n); Uh <- vector("list", n)
  for (g in seq_len(n)) {
    f  <- fcmd_stt(D, Pm, grid[g], Cc, mm, n_restart = FAST_RESTART)
    st <- stability_at_beta(bank, Pm, grid[g], Cc, mm, f$U, n_restart = FAST_REFIT)
    FM[g] <- fuzzy_moran(f$U, Z, Pm, mm); FS[g] <- fuzzy_silhouette(f$U, D)
    XB[g] <- xie_beni(f$U, f$medoids, D, mm)
    ST[g] <- st$ST; NDC[g] <- st$NDC; Uh[[g]] <- f$U
  }
  tb <- data.frame(beta = grid, FM = FM, FS = FS, XB = XB, ST = ST, NDC = NDC)
  Zn <- rescale_objectives(tb)          # honours STAB_MAP = "anchored"
  d <- d_ideal(Zn); s <- score_ws(Zn); i0 <- which.min(d)
  data.frame(sigma = sg, C = Cc, m = mm, rep = rep,
             beta_aspc = grid[i0], beta_moran = grid[which.max(FM)],
             beta_xb = grid[which.min(XB)], beta_knee = knee_point(grid, s),
             gamma = sort(d)[2] - sort(d)[1], nondom = mean(pareto_flag(Zn)),
             ST_at_hat = ST[i0], NDC_at_hat = NDC[i0],
             maxmemb = mean(apply(Uh[[i0]], 1, max)))
}

run_factorial <- function(X, Pm, Z, D,
                          sigmas = c(0.01, 0.03, 0.05, 0.07, 0.10),
                          Cs = c(4, 6, 8, 10), ms = c(1.5, 2.0, 2.5),
                          n_rep = 20, grid = BETA_GRID,
                          ncore = max(1, parallel::detectCores() - 1),
                          resume = TRUE) {
  cat("building", length(sigmas) * n_rep, "banks ...\n"); t0 <- proc.time()[3]
  banks <- setNames(vector("list", length(sigmas)), as.character(sigmas))
  for (si in seq_along(sigmas)) {
    banks[[si]] <- lapply(seq_len(n_rep), function(r) {
      set.seed(as.integer(1e5 * sigmas[si]) * 1000 + r)
      lapply(seq_len(FAST_B), function(b)
        dtw_matrix(X + array(stats::rnorm(prod(dim(X)), 0, sigmas[si]), dim(X)), DTW_NORM))
    })
    cat(sprintf("  sigma=%.2f (%.0f s)\n", sigmas[si], proc.time()[3] - t0))
  }
  cells <- expand.grid(rep = seq_len(n_rep), m = ms, C = Cs, sigma = sigmas,
                       KEEP.OUT.ATTRS = FALSE)[, c("sigma", "C", "m", "rep")]
  done <- NULL
  if (resume && file.exists(CHECKPOINT_FILE)) {
    done <- readRDS(CHECKPOINT_FILE)
    k <- paste(cells$sigma, cells$C, cells$m, cells$rep) %in%
         paste(done$sigma, done$C, done$m, done$rep)
    cells <- cells[!k, , drop = FALSE]
    cat("resuming:", nrow(done), "done,", nrow(cells), "left\n")
  }
  if (!nrow(cells)) return(done)
  job <- function(k) one_cell(D, Pm, Z, banks[[as.character(cells$sigma[k])]][[cells$rep[k]]],
                              cells$sigma[k], cells$C[k], cells$m[k], cells$rep[k], grid)
  cat("running", nrow(cells), "cells on", ncore, "cores ...\n"); t1 <- proc.time()[3]
  if (ncore > 1) {
    cl <- parallel::makeCluster(ncore); on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterExport(cl, ls(envir = globalenv()), envir = globalenv())
    parallel::clusterExport(cl, c("D", "Pm", "Z", "banks", "cells", "grid"), envir = environment())
    res <- parallel::parLapply(cl, seq_len(nrow(cells)), job)
  } else res <- lapply(seq_len(nrow(cells)), job)
  out <- do.call(rbind, res); if (!is.null(done)) out <- rbind(done, out)
  saveRDS(out, CHECKPOINT_FILE)
  cat(sprintf("done in %.1f min\n", (proc.time()[3] - t1) / 60))
  out
}

factorial_summaries <- function(fs, grid = BETA_GRID) {
  ends <- range(grid)
  pe <- function(v) c(100 * mean(v %in% ends), 100 * mean(!v %in% ends))
  endpoint <- data.frame(
    Rule = c("Fuzzy Moran maximisation", "XB minimisation", "SA-ASPC", "Knee point"),
    pct_endpoint = c(pe(fs$beta_moran)[1], pe(fs$beta_xb)[1],
                     pe(fs$beta_aspc)[1], pe(fs$beta_knee)[1]),
    pct_zero = c(100 * mean(fs$beta_moran == 0), 100 * mean(fs$beta_xb == 0),
                 100 * mean(fs$beta_aspc == 0), 100 * mean(fs$beta_knee == 0)),
    median_beta = c(stats::median(fs$beta_moran), stats::median(fs$beta_xb),
                    stats::median(fs$beta_aspc), stats::median(fs$beta_knee)))
  vd <- function(y) {
    a <- stats::anova(stats::lm(y ~ factor(sigma) * factor(C) * factor(m), data = fs))
    data.frame(term = rownames(a), pct = 100 * a[["Sum Sq"]] / sum(a[["Sum Sq"]]))
  }
  list(endpoint = endpoint, var_beta = vd(fs$beta_aspc),
       var_ST = vd(fs$ST_at_hat), var_NDC = vd(fs$NDC_at_hat),
       marginal = stats::aggregate(cbind(beta_aspc, beta_moran, ST_at_hat,
                                         NDC_at_hat, maxmemb, gamma) ~ sigma + C + m,
                                   fs, mean),
       distinct = stats::aggregate(beta_aspc ~ sigma + C + m, fs,
                                   function(z) length(unique(z))))
}

###############################################################################
## BLOCK 16 : MAIN
###############################################################################

wr <- function(x, nm) {
  utils::write.csv(x, file.path(DIR_OUT, paste0(nm, ".csv")), row.names = FALSE)
  cat("written:", nm, "\n")
}

X_raw   <- load_panel(PATH_PANEL, COL_REGION, COL_TIME, INDICATORS)
X       <- zscore_panel(X_raw)
regions <- dimnames(X)[[1]]
Wq      <- build_W(regions)
Pm      <- row_standardise(Wq)
D       <- dtw_matrix(X, DTW_NORM)
Zbar    <- apply(X, c(1, 3), mean)

cat(sprintf("I=%d T=%d p=%d | mean neighbours=%.2f | islands=%d | mean DTW=%.4f\n",
            dim(X)[1], dim(X)[2], dim(X)[3], mean(rowSums(Wq)),
            sum(rowSums(Wq) == 0), mean(D[upper.tri(D)])))
wr(data.frame(I = dim(X)[1], TT = dim(X)[2], p = dim(X)[3],
              mean_neighbours = mean(rowSums(Wq)), mean_dtw = mean(D[upper.tri(D)]),
              dtw_normalised = DTW_NORM, C_star = C_STAR, m = M_FUZZ,
              B = B_BANK, sigma = SIGMA_PERT, R = N_RESTART,
              stab_map = STAB_MAP, master_seed = MASTER_SEED), "T00_Setup")

if (RUN_CSELECT) { T5 <- select_C(D, Pm); wr(T5, "T05_ClusterNumber"); print(T5) }

if (RUN_CALIB) {
  bank <- build_bank(X)
  cal  <- calibrate(D, Pm, Zbar, bank)
  tb   <- cal$table
  Zn   <- rescale_objectives(tb, "anchored")
  Zm   <- rescale_objectives(tb, "minmax")
  tb$Score <- score_ws(Zn); tb$distance <- d_ideal(Zn)
  tb$distance_minmax <- d_ideal(Zm); tb$Pareto <- pareto_flag(Zn)
  beta_hat   <- tb$beta[which.min(tb$distance)]
  beta_mm    <- tb$beta[which.min(tb$distance_minmax)]
  beta_moran <- tb$beta[which.max(tb$FM)]
  wr(tb, "T06_CalibrationGrid")
  cat(sprintf("\nANCHORED: beta_hat=%.2f | MINMAX: beta_hat=%.2f | Moran=%.2f | nondom=%d/%d\n",
              beta_hat, beta_mm, beta_moran, sum(tb$Pareto), nrow(tb)))
  cat("entropy range:", round(range(tb$Entropy), 4), " log C =", round(log(C_STAR), 4), "\n")
  cat("Spearman ST~NDC:", round(stats::cor(tb$ST, tb$NDC, method = "spearman"), 3),
      "| FM~SC:", round(stats::cor(tb$FM, tb$SC, method = "spearman"), 3), "\n")
}

if (RUN_CALIB && RUN_RULES) {
  T7 <- selection_rules(tb, Zn)
  T7 <- cbind(T7, tb[match(T7$beta, tb$beta), c("FM", "FS", "XB", "ST", "NDC", "distance")])
  wr(T7, "T07_SelectionRules"); print(T7)
  T7m <- selection_rules(tb, Zm); wr(T7m, "T07b_SelectionRules_minmax"); print(T7m)
}

if (RUN_CALIB && RUN_WEIGHTS) {
  ws <- weight_sensitivity(tb, Zn)
  wr(ws$freq, "T08a_WeightSensitivity"); wr(ws$loo, "T08b_LeaveOneOut")
  cat(sprintf("gamma=%.6f | equal weights -> %.2f\n", ws$gamma, ws$equal_weights))
}

if (RUN_CALIB && RUN_MCCERT) {
  mca <- mc_certification(tb, Zn, stab = "anchored")
  mcm <- mc_certification(tb, Zm, stab = "minmax")
  cmp <- data.frame(map = c("anchored", "minmax"),
                    R = c(mca$R, mcm$R), kappa = c(mca$kappa, mcm$kappa),
                    four_kappa = c(mca$four_kappa, mcm$four_kappa),
                    gamma_hat = c(mca$gamma_hat, mcm$gamma_hat),
                    ratio = c(mca$ratio, mcm$ratio),
                    certified = c(mca$certified, mcm$certified),
                    B_required = c(mca$B_required, mcm$B_required),
                    conf_set_size = c(length(mca$conf_set), length(mcm$conf_set)))
  wr(cmp, "T_MonteCarloCertification"); print(cmp)
  cat("anchored confidence set:", mca$conf_set, "\n")
}

if (RUN_CALIB && RUN_BASELINE) {
  bl <- baseline_table(D, Pm, Zbar, X, C_STAR, M_FUZZ, cal$fits, beta_hat, beta_moran)
  wr(bl$table, "T09_Baselines"); wr(bl$agreement, "T09b_PartitionAgreement")
  print(bl$table); print(bl$agreement)
}

if (RUN_CALIB && RUN_PARTITION) {
  fh <- cal$fits[[as.character(beta_hat)]]
  wr(partition_profile(fh, X_raw), "T_CalibratedPartition")
  wr(data.frame(region = regions, cluster = hard(fh$U),
                max_membership = apply(fh$U, 1, max)), "T_RegionAssignments")
}

if (RUN_CALIB && RUN_DRIFT) {
  pd <- partition_drift(cal$fits, tb); wr(pd, "T_PartitionDrift"); print(pd, row.names = FALSE)
}

if (RUN_CONVERGE) { T10 <- convergence_study(D, Pm); wr(T10, "T10_Convergence"); print(T10) }
if (RUN_RESTART) {
  T11 <- restart_study(D, Pm); wr(T11, "T11_RestartDispersion"); print(T11)
  if (RUN_CALIB) { ss <- seed_stability(D, Pm, Zbar, bank); wr(ss, "T11b_SeedStability"); print(ss) }
}
if (RUN_RUNTIME) { rt <- runtime_benchmark(X, Pm, D); wr(rt, "T_Runtime"); print(rt) }
if (RUN_DTWNORM) { nc <- normalisation_comparison(X, Pm, Zbar); wr(nc, "T_DTWNormalisation"); print(nc) }
if (RUN_CALIB && RUN_REFINE) { gr <- grid_refinement(D, Pm, Zbar, bank); wr(gr, "T_GridRefinement"); print(gr) }

if (RUN_FACTORIAL) {
  fs <- run_factorial(X, Pm, Zbar, D, Cs = c(4, 6, 8, 10))
  wr(fs, "T12_FactorialRaw")
  fsm <- factorial_summaries(fs)
  wr(fsm$endpoint, "T13_EndpointRates");  wr(fsm$var_beta, "T14_VarianceBeta")
  wr(fsm$var_ST,   "T15_VarianceST");     wr(fsm$var_NDC,  "T16_VarianceNDC")
  wr(fsm$marginal, "T17_MarginalMeans");  wr(fsm$distinct, "T17b_DistinctBeta")
  print(fsm$endpoint); print(fsm$var_beta)
}

## ---- reproducibility record ---------------------------------------------
si <- sessionInfo()
repro <- data.frame(
  r_version = si$R.version$version.string,
  platform  = si$platform,
  cores     = parallel::detectCores(),
  sf        = as.character(utils::packageVersion("sf")),
  spdep     = as.character(utils::packageVersion("spdep")),
  cluster   = as.character(utils::packageVersion("cluster")),
  clue      = as.character(utils::packageVersion("clue")),
  master_seed = MASTER_SEED, stab_map = STAB_MAP,
  date      = format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
wr(repro, "T_Reproducibility"); print(t(repro))
saveRDS(list(cal = if (exists("cal")) cal else NULL, tb = if (exists("tb")) tb else NULL,
             D = D, Pm = Pm, Wq = Wq, Zbar = Zbar, X = X, X_raw = X_raw,
             C_STAR = C_STAR, STAB_MAP = STAB_MAP),
        file.path(DIR_OUT, "state_v2.rds"))
cat("\nAll blocks complete. Outputs in:", normalizePath(DIR_OUT), "\n")

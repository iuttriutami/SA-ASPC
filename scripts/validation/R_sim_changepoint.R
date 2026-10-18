## =============================================================================
## R_sim_changepoint.R — desain eksternal non-sinusoidal: trajektori dengan change point (Table S23)
##
## Menjawab reviewer (Section 6.5): satu keluarga arketipe yang berbeda dari sinusoid + drift,
## dan ukuran keberhasilan KEPUTUSAN selain ARI.
## Desain E3: lattice 6x6, queen, C0 = 3 pita kontigu, J = 3, T = 12, derau SAR (rho = 0.8, Gaussian).
##   Arketipe c: level konstan, lalu loncatan pada waktu tau_c = 4, 7, 10 (arah loncatan berbeda per klaster),
##   dikalikan faktor sinyal s. Proporsi pi unit diberi label acak (pi = 0 dan 0.15) seperti E1.
##   Sinyal "weak"/"strong" dipilih dengan pilot sehingga fit dangkal beta = 0 mencapai ARI ~0.30 / ~0.75.
## Protokol: deep (R = 100, enam salinan terperturbasi, dua sweep) atas grid {0, 0.5, ..., 6}, bank B = 30.
## Aturan: beta = 0, XB-only, Moran max, SA-ASPC (anchored), SA-ASPC tanpa ST, guard dominance, guard XB 2%,
##         dan oracle (grid point dengan ARI tertinggi).
## Ukuran keputusan: "keputusan benar" = (aturan memilih beta > 0) sama dengan (oracle memperoleh ARI
##   > ARI(beta=0) + 0.01). Dilaporkan juga proporsi lolos diagnostik empiris SA-ASPC.
##
## Input: SA_ASPC_pipeline_v2.R. Output: /kaggle/working/cp_out/tableS23_changepoint*.csv
## Checkpoint per replikasi (folder cp_ckpt). Jika sesi habis: pasang output sebagai Input lalu jalankan lagi.
## =============================================================================

N_SIDE <- 6; C0 <- 3; J <- 3; TT <- 12; RHO <- 0.8; TAU <- c(4, 7, 10)
PI_SET <- c(0, 0.15); N_REP <- 20
GRID <- seq(0, 6, by = 0.5); M_FUZZ <- 2
R_DEEP <- 100; N_PERT <- 6; SWEEPS <- 2; PERT_SD <- 0.30; R_SHALLOW <- 10
B_BANK <- 30; SIGMA <- 0.05; BANK_RESTART <- 3
W_OBJ <- c(FM = 0.35, FS = 0.25, XB = 0.20, ST = 0.20)
PILOT_S <- c(0.25, 0.5, 0.75, 1, 1.25, 1.5, 2, 2.5, 3); PILOT_REP <- 8; TARGET_ARI <- c(weak = 0.30, strong = 0.75)
CORES <- max(1, parallel::detectCores())

WORK <- "/kaggle/working"; CK <- file.path(WORK, "cp_ckpt"); OUT <- file.path(WORK, "cp_out")
dir.create(CK, showWarnings = FALSE); dir.create(OUT, showWarnings = FALSE)
inp <- function(pat) { f <- list.files("/kaggle/input", pat, recursive = TRUE, full.names = TRUE, ignore.case = TRUE); if (length(f)) f[1] else NA }
ck <- function(name, expr) {
  f <- file.path(CK, paste0(name, ".rds"))
  old <- list.files("/kaggle/input", paste0("^", name, "\\.rds$"), recursive = TRUE, full.names = TRUE)
  old <- if (length(old <- old[basename(dirname(old)) == basename(CK)])) old[1] else NA
  if (file.exists(f)) return(readRDS(f))
  if (!is.na(old)) { file.copy(old, f); return(readRDS(f)) }
  v <- expr; saveRDS(v, f); v
}

## ---------- pipeline ----------
pipe <- new.env()
for (e in parse(inp("SA_ASPC_pipeline_v2\\.R$"), keep.source = FALSE)) {
  if (is.call(e) && as.character(e[[1]]) %in% c("<-", "=") && is.name(e[[2]])) {
    nm <- as.character(e[[2]]); rhs <- e[[3]]
    if ((is.call(rhs) && identical(rhs[[1]], as.name("function"))) || grepl("^[A-Z][A-Z0-9_]*$", nm))
      if (!nm %in% c("X", "D")) tryCatch(eval(e, pipe), error = function(err) NULL)
  }
}
P_ <- function(nm) get(nm, envir = pipe)

## ---------- lattice, kebenaran, arketipe ----------
I <- N_SIDE^2
xy <- expand.grid(r = seq_len(N_SIDE), c = seq_len(N_SIDE))
Wq <- outer(seq_len(I), seq_len(I), function(a, b) pmax(abs(xy$r[a] - xy$r[b]), abs(xy$c[a] - xy$c[b])) == 1) * 1
Pm <- P_("row_standardise")(Wq)
truth0 <- as.integer(cut(xy$c, C0, labels = FALSE))                 # pita kontigu
JUMP <- rbind(c(1, -1, 0.5), c(-1, 0.5, 1), c(0.5, 1, -1))          # arah loncatan per klaster (C0 x J)
PRE  <- rbind(c(0.2, 0, -0.2), c(0, -0.2, 0.2), c(-0.2, 0.2, 0))     # level awal kecil per klaster
archetype <- function(c, s) sapply(seq_len(J), function(j) s * (PRE[c, j] + JUMP[c, j] * (seq_len(TT) >= TAU[c])))
Ainv <- solve(diag(I) - RHO * Pm)

gen <- function(s, pi_, seed) {
  set.seed(seed)
  truth <- truth0
  if (pi_ > 0) { sw <- sample.int(I, round(pi_ * I)); truth[sw] <- sample.int(C0, length(sw), replace = TRUE) }
  X <- array(0, c(I, TT, J), dimnames = list(sprintf("u%02d", seq_len(I)), seq_len(TT), paste0("V", seq_len(J))))
  for (i in seq_len(I)) X[i, , ] <- archetype(truth[i], s)
  for (t in seq_len(TT)) for (j in seq_len(J)) X[, t, j] <- X[, t, j] + as.vector(Ainv %*% rnorm(I))
  X <- tryCatch(P_("zscore_panel")(X), error = function(e) {   # cadangan: z-score per variabel
    for (j in seq_len(J)) X[, , j] <- (X[, , j] - mean(X[, , j])) / sd(X[, , j]); X })
  list(X = X, truth = truth)
}

## ---------- pilot sinyal ----------
pilot <- ck("pilot", do.call(rbind, lapply(PILOT_S, function(s) {
  a <- unlist(parallel::mclapply(seq_len(PILOT_REP), function(r) {
    d <- gen(s, 0, 5e6 + 100 * r + round(100 * s)); D <- P_("dtw_matrix")(d$X, TRUE)
    P_("ari")(d$truth, P_("hard")(P_("fcmd_stt")(D, Pm, 0, C0, M_FUZZ, n_restart = R_SHALLOW, seed = r)$U))
  }, mc.cores = CORES))
  data.frame(s = s, ari = mean(a))
})))
print(pilot)
S_LEVEL <- sapply(TARGET_ARI, function(a) pilot$s[which.min(abs(pilot$ari - a))])
cat("Sinyal terpilih:", paste(names(S_LEVEL), S_LEVEL, collapse = ", "), "\n")

## ---------- satu replikasi ----------
G <- length(GRID)
perturb <- function(U) { V <- U * exp(matrix(rnorm(length(U), 0, PERT_SD), nrow(U))); V / rowSums(V) }
deep_fit <- function(D, seed) {
  set.seed(seed)
  best <- lapply(seq_len(G), function(g) P_("fcmd_stt")(D, Pm, GRID[g], C0, M_FUZZ, n_restart = R_DEEP, seed = seed + g))
  for (sw in seq_len(SWEEPS)) for (ord in list(seq_len(G), rev(seq_len(G)))) for (k in seq_along(ord)[-1]) {
    g <- ord[k]; U0p <- best[[ord[k - 1]]]$U
    for (p in seq_len(N_PERT)) {
      f <- P_("fcmd_stt_fit1")(D, Pm, GRID[g], C0, M_FUZZ, if (p == 1) U0p else perturb(U0p))
      if (f$J < best[[g]]$J) best[[g]] <- f
    }
  }
  best
}
one_rep <- function(level, pi_, r) {
  seed <- 7e6 + 1e4 * match(level, names(S_LEVEL)) + 1e3 * match(pi_, PI_SET) + r
  d <- gen(S_LEVEL[[level]], pi_, seed); X <- d$X
  D <- P_("dtw_matrix")(X, TRUE); Zbar <- apply(X, c(1, 3), mean)
  fits <- deep_fit(D, seed)
  bank <- P_("build_bank")(X, B = B_BANK, sigma = SIGMA, normalise = TRUE, seed = seed + 99)
  st <- sapply(seq_len(G), function(g) {
    ref <- P_("hard")(fits[[g]]$U)
    sapply(bank, function(Db) P_("ari")(ref, P_("hard")(P_("fcmd_stt")(Db, Pm, GRID[g], C0, M_FUZZ, n_restart = BANK_RESTART)$U)))
  })
  tb <- data.frame(beta = GRID,
    FM = sapply(fits, function(f) P_("fuzzy_moran")(f$U, Zbar, Pm, M_FUZZ)),
    FS = sapply(fits, function(f) P_("fuzzy_silhouette")(f$U, D)),
    XB = sapply(fits, function(f) P_("xie_beni")(f$U, f$medoids, D, M_FUZZ)),
    ST = colMeans(st))
  Zn <- P_("rescale_objectives")(tb, "anchored")
  dd <- P_("d_ideal")(Zn, W_OBJ); o <- order(dd); k <- o[1]; gam <- dd[o[2]] - dd[o[1]]
  k_noST <- which.min(P_("d_ideal")(Zn, c(W_OBJ[1:3], 0) / sum(W_OBJ[1:3])))
  dom_ok <- !(tb$XB[k] > tb$XB[1] && tb$ST[k] < tb$ST[1])
  ari <- sapply(fits, function(f) P_("ari")(d$truth, P_("hard")(f$U)))
  pick <- c(beta0 = 1, XB_only = which.min(tb$XB), Moran_max = which.max(tb$FM), SA_ASPC = k, SA_ASPC_noST = k_noST,
            guard_dominance = if (dom_ok) k else 1, guard_XB2 = if (tb$XB[k] <= 1.02 * tb$XB[1]) k else 1,
            oracle = which.max(ari))
  emp_pass <- 4 * sqrt(W_OBJ[["ST"]]) * (1.96 * max(apply(st, 2, sd)) / sqrt(B_BANK)) / 2 < gam
  data.frame(level = level, pi = pi_, rep = r, rule = names(pick), beta = GRID[pick], ARI = ari[pick],
             ARI_beta0 = ari[1], penalty_pays = max(ari) > ari[1] + 0.01, gamma = gam, emp_pass = emp_pass)
}

cells <- expand.grid(level = names(S_LEVEL), pi = PI_SET, rep = seq_len(N_REP), stringsAsFactors = FALSE)
res <- do.call(rbind, parallel::mclapply(seq_len(nrow(cells)), function(i) with(cells[i, ],
  ck(sprintf("rep_%s_pi%02d_r%02d", level, round(100 * pi), rep), one_rep(level, pi, rep))),
  mc.cores = CORES, mc.preschedule = FALSE))
write.csv(res, file.path(OUT, "tableS23b_changepoint_replicates.csv"), row.names = FALSE)

## ---------- ringkasan: ARI dan keberhasilan keputusan ----------
summ <- do.call(rbind, lapply(split(res, list(res$level, res$pi, res$rule), drop = TRUE), function(s) {
  dif <- s$ARI - s$ARI_beta0
  h <- if (nrow(s) > 1) qt(0.975, nrow(s) - 1) * sd(dif) / sqrt(nrow(s)) else NA
  data.frame(level = s$level[1], pi = s$pi[1], rule = s$rule[1], mean_ARI = mean(s$ARI),
             diff_vs_beta0 = mean(dif), ci_half = h, worst_diff = min(dif), mean_beta = mean(s$beta),
             decision_correct = mean((s$beta > 0) == s$penalty_pays),
             share_penalty_pays = mean(s$penalty_pays))
}))
summ <- summ[order(summ$level, summ$pi, summ$rule), ]
num <- sapply(summ, is.numeric); summ[num] <- lapply(summ[num], function(v) round(v, 3))
print(summ, row.names = FALSE)
write.csv(summ, file.path(OUT, "tableS23_changepoint.csv"), row.names = FALSE)
sa <- res[res$rule == "SA_ASPC", ]
cat(sprintf("\nSA-ASPC: margin median %.4f; lolos diagnostik empiris %.0f%% replikasi\n",
            median(sa$gamma), 100 * mean(sa$emp_pass)))

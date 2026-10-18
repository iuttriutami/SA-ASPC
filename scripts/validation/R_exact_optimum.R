## =============================================================================
## R_exact_optimum.R — optimum global EKSAK pada beta = 0 sebagai tolok ukur protokol deep (Table S22)
##
## Menjawab reviewer 1.3. Pada beta = 0, untuk himpunan medoid q yang tetap, keanggotaan optimal
## punya bentuk tertutup (Proposition 1), sehingga untuk m = 2
##     J*(q) = sum_i 1 / sum_c 1/d(i, q_c)        (0 bila i sendiri medoid).
## Optimum global = min_q J*(q) atas semua himpunan C medoid. Kami mencarinya dengan
## branch-and-bound (Rcpp) dengan batas bawah yang sah:
##     J*(Q u R) >= sum_i 1 / (S_i(Q) + jumlah k invers-jarak terbesar milik i di kandidat sisa),
## sehingga hasilnya terbukti optimal bila pencarian selesai (status "proved").
## Dibandingkan dengan: protokol deep (N_MASTER master seed x R_DEEP restart) dan dangkal (best-of-20).
##
## Kasus: aplikasi 1 (C = 4, 8) dan aplikasi 2 (C = 4, 6).
## Input: Panel_JawaTengah_2022_2025.xlsx, Panel_App2_Jawa_2020_2024.csv, SA_ASPC_pipeline_v2.R
## Output: /kaggle/working/exact_out/tableS22_exact_optimum.csv
## Catatan: untuk beta > 0 langkah keanggotaan tidak konveks, jadi tidak ada sertifikat optimum;
## titik yang diuji di sini (beta = 0) adalah pilihan aplikasi 1 pada C = 8.
## =============================================================================

CASES <- list(list(app = 1, C = 4), list(app = 1, C = 8), list(app = 2, C = 4), list(app = 2, C = 6))
M_FUZZ <- 2; R_DEEP <- 200; N_MASTER <- 20; R_SHALLOW <- 20; N_SHALLOW <- 20
TIME_LIMIT <- 4 * 3600                     # detik per kasus untuk branch-and-bound
CORES <- max(1, parallel::detectCores())

WORK <- "/kaggle/working"; CK <- file.path(WORK, "exact_ckpt"); OUT <- file.path(WORK, "exact_out")
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

## ---------- data ----------
load_D <- function(app) ck(sprintf("D_app%d", app), {
  if (app == 1) {
    raw <- as.data.frame(readxl::read_excel(inp("Panel_JawaTengah.*\\.xlsx$")))
    reg <- "WADMKK"; ind <- c("APM_SD", "APM_SMP", "PRODPADI", "RLS", "SANITASI")
  } else {
    raw <- read.csv(inp("Panel_App2_Jawa.*\\.csv$"), check.names = FALSE)
    reg <- "KODE"; ind <- c("IPM", "MISKIN", "TPT", "SANITASI", "LOG_PDRB_KAP")
  }
  names(raw) <- toupper(trimws(names(raw)))
  tmp <- file.path(WORK, sprintf("panel%d_tmp.csv", app)); write.csv(raw, tmp, row.names = FALSE)
  X <- P_("zscore_panel")(P_("load_panel")(tmp, reg, "TAHUN", ind))
  P_("dtw_matrix")(X, TRUE)
})
Dl <- list(load_D(1), load_D(2))
cat(sprintf("Rerata DTW: app1 %.3f (naskah 2.958), app2 %.3f (naskah 2.854)\n",
            mean(Dl[[1]][upper.tri(Dl[[1]])]), mean(Dl[[2]][upper.tri(Dl[[2]])])))

## ---------- branch-and-bound (Rcpp) ----------
Rcpp::sourceCpp(code = '
#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <functional>
#include <ctime>
using namespace Rcpp;
// D: dissimilarities (n x n). Zero distances are handled exactly via counters (no 1/0).
// [[Rcpp::export]]
List bb_medoids(NumericMatrix D, int C, double ub, double time_limit) {
  const int n = D.nrow();
  std::vector<double> inv((size_t)n * n); std::vector<char> zero((size_t)n * n);
  for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) {
    zero[(size_t)i * n + j] = D(i, j) <= 0; inv[(size_t)i * n + j] = D(i, j) > 0 ? 1.0 / D(i, j) : 0.0; }
  std::vector<double> top((size_t)n * (n + 1) * (C + 1), 0.0); std::vector<char> zs((size_t)n * (n + 1), 0);
  auto T = [&](int i, int j, int k) -> double& { return top[((size_t)i * (n + 1) + j) * (C + 1) + k]; };
  for (int i = 0; i < n; i++) {
    std::vector<double> best;
    for (int j = n - 1; j >= 0; j--) {
      zs[(size_t)i * (n + 1) + j] = zs[(size_t)i * (n + 1) + j + 1] || zero[(size_t)i * n + j];
      best.push_back(inv[(size_t)i * n + j]); std::sort(best.begin(), best.end(), std::greater<double>());
      if ((int)best.size() > C) best.pop_back();
      double s = 0; for (int k = 1; k <= C; k++) { if (k <= (int)best.size()) s += best[k - 1]; T(i, j, k) = s; }
    }
  }
  std::vector<double> S(n, 0.0); std::vector<int> Z(n, 0), cur, bestset;
  double bestval = ub, nodes = 0; bool timeout = false; std::clock_t t0 = std::clock();
  std::function<void(int,int)> dfs = [&](int j, int depth) {
    if (timeout) return;
    nodes += 1;
    if (std::fmod(nodes, 1048576.0) == 0 && double(std::clock() - t0) / CLOCKS_PER_SEC > time_limit) { timeout = true; return; }
    int k = C - depth;
    if (k == 0) {
      double v = 0; for (int i = 0; i < n; i++) if (!Z[i]) v += 1.0 / S[i];
      if (v < bestval) { bestval = v; bestset = cur; }
      return;
    }
    double lb = 0;
    for (int i = 0; i < n; i++) if (!Z[i] && !zs[(size_t)i * (n + 1) + j]) lb += 1.0 / (S[i] + T(i, j, k));
    if (lb >= bestval) return;
    for (int l = j; l <= n - k; l++) {
      for (int i = 0; i < n; i++) { S[i] += inv[(size_t)i * n + l]; Z[i] += zero[(size_t)i * n + l]; }
      cur.push_back(l); dfs(l + 1, depth + 1); cur.pop_back();
      for (int i = 0; i < n; i++) { S[i] -= inv[(size_t)i * n + l]; Z[i] -= zero[(size_t)i * n + l]; }
      if (timeout) return;
    }
  };
  dfs(0, 0);
  IntegerVector bs(bestset.size()); for (size_t a = 0; a < bestset.size(); a++) bs[a] = bestset[a] + 1;
  return List::create(_["value"] = bestval, _["medoids"] = bs, _["nodes"] = nodes, _["proved"] = !timeout,
                      _["seconds"] = double(std::clock() - t0) / CLOCKS_PER_SEC);
}')

Jstar <- function(D, q) { Dq <- D[, q, drop = FALSE]; v <- 1 / rowSums(1 / Dq); v[apply(Dq <= 0, 1, any)] <- 0; sum(v) }
J_of_fit <- function(D, f) sum(f$U^M_FUZZ * D[, f$medoids])

res <- do.call(rbind, lapply(CASES, function(cs) {
  D <- Dl[[cs$app]]; C <- cs$C; tag <- sprintf("app%d_C%d", cs$app, C)
  cat("\n==", tag, "==\n")
  deep <- ck(paste0("deep_", tag), parallel::mclapply(seq_len(N_MASTER), function(s)
    P_("fcmd_stt")(D, diag(nrow(D)), 0, C, M_FUZZ, n_restart = R_DEEP, seed = 31e6 + 1000 * s), mc.cores = CORES))
  shal <- ck(paste0("shallow_", tag), parallel::mclapply(seq_len(N_SHALLOW), function(s)
    P_("fcmd_stt")(D, diag(nrow(D)), 0, C, M_FUZZ, n_restart = R_SHALLOW, seed = 41e6 + 1000 * s), mc.cores = CORES))
  Jd <- sapply(deep, `[[`, "J"); Js <- sapply(shal, `[[`, "J")
  bd <- deep[[which.min(Jd)]]
  chk <- c(J_reported = min(Jd), J_recomputed = J_of_fit(D, bd), J_closed_form = Jstar(D, bd$medoids))
  cat("Cek objektif (harus hampir sama):", round(chk, 6), "\n")
  ## urutkan kandidat: medoid deep dulu, lalu unit paling sentral (memperkuat pemangkasan)
  ord <- unique(c(bd$medoids, order(rowSums(D))))
  bb <- ck(paste0("bb_", tag), bb_medoids(D[ord, ord], C, ub = Jstar(D, bd$medoids) * (1 + 1e-9), time_limit = TIME_LIMIT))
  Jex <- if (length(bb$medoids)) bb$value else Jstar(D, bd$medoids)
  q_ex <- if (length(bb$medoids)) ord[bb$medoids] else bd$medoids
  hard_ex <- apply(D[, q_ex, drop = FALSE], 1, which.min)
  ari_deep <- P_("ari")(hard_ex, P_("hard")(bd$U))
  ari_shal <- sapply(shal, function(f) P_("ari")(hard_ex, P_("hard")(f$U)))
  tol <- 1e-6 * Jex
  data.frame(case = tag, n = nrow(D), C = C, combinations = choose(nrow(D), C),
             proved_optimal = bb$proved, bb_nodes = bb$nodes, bb_seconds = round(bb$seconds, 1),
             J_exact = Jex, J_deep_best = min(Jd),
             rel_gap_deep = (min(Jd) - Jex) / Jex,
             share_deep_seeds_at_optimum = mean(Jd <= Jex + tol),
             J_shallow_median = median(Js), rel_gap_shallow_median = (median(Js) - Jex) / Jex,
             rel_gap_shallow_max = (max(Js) - Jex) / Jex,
             share_shallow_at_optimum = mean(Js <= Jex + tol),
             same_medoids_deep = setequal(q_ex, bd$medoids),
             ARI_deep_vs_exact = ari_deep, ARI_shallow_vs_exact_min = min(ari_shal))
}))
num <- sapply(res, is.numeric); res[num] <- lapply(res[num], function(v) signif(v, 6))
print(res, row.names = FALSE)
write.csv(res, file.path(OUT, "tableS22_exact_optimum.csv"), row.names = FALSE)
cat("\nproved_optimal = TRUE berarti J_exact terbukti optimum global pada beta = 0.\n")
cat("rel_gap_deep = 0 dan same_medoids_deep = TRUE berarti protokol deep mencapai optimum global.\n")

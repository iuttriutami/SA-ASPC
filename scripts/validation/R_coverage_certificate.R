## =============================================================================
## R_coverage_certificate.R — validasi empiris sertifikat Monte Carlo (Section 3.4)
## Aplikasi 1 (35 kab/kota Jawa Tengah, C = 8, m = 2), fit kanonik deep naskah.
##
## Desain:
##   1. Bangun POOL anggota bank perturbasi (sigma = 0.05) dan hitung ARI tiap anggota
##      terhadap partisi kanonik di setiap beta (sama seperti ST di naskah).
##   2. Target beta* = pilihan SA-ASPC dengan ST dari seluruh POOL (pendekatan target Monte Carlo).
##   3. Untuk B = 25, 50, 100, 200, 400: ambil N_REP bank berukuran B dari POOL (tanpa pengembalian),
##      hitung beta-hat, margin, sertifikat eksak, diagnostik empiris, dan confidence set.
##   4. Laporkan: P(beta-hat = beta*), proporsi tersertifikasi, false certification
##      (tersertifikasi tetapi salah), cakupan confidence set, ukuran confidence set.
##
## Input: Panel_JawaTengah_2022_2025.xlsx, Wq_frozen.csv, E8_canonical_fits.rds, SA_ASPC_pipeline_v2.R
## Checkpoint per blok 100 anggota -> aman dari batas 12 jam (pasang output versi sebelumnya sebagai Input).
## =============================================================================

POOL      <- 2000                     # ukuran pool (target)
CHUNK     <- 100                      # anggota per checkpoint
B_SET     <- c(25, 50, 100, 200, 400)
N_REP     <- 200
SIGMA     <- 0.05
BANK_RESTART <- 3                     # sama dengan analisis utama
W_OBJ     <- c(FM = 0.35, FS = 0.25, XB = 0.20, ST = 0.20)
DELTA     <- 0.05
GRID      <- seq(0, 6, by = 0.25)
M_FUZZ    <- 2
CORES     <- max(1, parallel::detectCores())

WORK <- "/kaggle/working"; CK <- file.path(WORK, "cov_ckpt"); OUT <- file.path(WORK, "cov_out")
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

## ---------- data aplikasi 1 ----------
raw <- as.data.frame(readxl::read_excel(inp("Panel_JawaTengah.*\\.xlsx$")))
names(raw) <- toupper(trimws(names(raw)))
IND <- c("APM_SD", "APM_SMP", "PRODPADI", "RLS", "SANITASI")
tmp <- file.path(WORK, "panel1_tmp.csv"); write.csv(raw, tmp, row.names = FALSE)
X1 <- P_("zscore_panel")(P_("load_panel")(tmp, "WADMKK", "TAHUN", IND))
D1 <- ck("D_app1", P_("dtw_matrix")(X1, TRUE))
cat("Rerata DTW:", round(mean(D1[upper.tri(D1)]), 3), "(naskah 2.958)\n")
reg <- dimnames(X1)[[1]]
W1 <- as.matrix(read.csv(inp("Wq_frozen\\.csv$"), row.names = 1, check.names = FALSE))
nk <- function(v) { v <- toupper(trimws(sub("^X", "", v))); v <- sub("^KABUPATEN\\s+|^KAB\\.?\\s+", "", v); gsub("[^A-Z]", "", v) }
mr <- match(nk(reg), nk(rownames(W1))); mc <- match(nk(reg), nk(colnames(W1)))
if (anyNA(mr) || anyNA(mc)) stop("Nama unit panel tidak cocok dengan Wq_frozen.csv: ", paste(reg[is.na(mr)], collapse = ", "))
W1 <- W1[mr, mc]; dimnames(W1) <- list(reg, reg); Pm <- P_("row_standardise")(W1)
Zbar <- apply(X1, c(1, 3), mean)

canon <- readRDS(inp("E8_canonical_fits\\.rds$"))
stopifnot(length(canon) == length(GRID), nrow(canon[[1]]$U) == 35)
C1 <- ncol(canon[[1]]$U); cat("C =", C1, "| grid =", length(GRID), "titik\n")
ref_h <- sapply(canon, function(f) P_("hard")(f$U))

## objektif deterministik dari fit kanonik
FM <- sapply(canon, function(f) P_("fuzzy_moran")(f$U, Zbar, Pm, M_FUZZ))
FS <- sapply(canon, function(f) P_("fuzzy_silhouette")(f$U, D1))
XB <- sapply(canon, function(f) P_("xie_beni")(f$U, f$medoids, D1, M_FUZZ))

## ---------- 1. pool: ARI anggota bank terhadap partisi kanonik ----------
n_chunk <- POOL / CHUNK
st_pool <- do.call(rbind, lapply(seq_len(n_chunk), function(j) ck(sprintf("pool_%03d", j), {
  cat("blok", j, "/", n_chunk, "\n")
  bank <- P_("build_bank")(X1, B = CHUNK, sigma = SIGMA, normalise = TRUE, seed = 900000 + j)
  do.call(rbind, parallel::mclapply(seq_along(bank), function(b) {
    sapply(seq_along(GRID), function(g)
      P_("ari")(ref_h[, g], P_("hard")(P_("fcmd_stt")(bank[[b]], Pm, GRID[g], C1, M_FUZZ, n_restart = BANK_RESTART)$U)))
  }, mc.cores = CORES, mc.preschedule = FALSE))
})))
cat("Pool:", nrow(st_pool), "anggota x", ncol(st_pool), "titik grid\n")

## ---------- 2. seleksi dan sertifikat ----------
select_cert <- function(ST, B) {
  tb <- data.frame(beta = GRID, FM = FM, FS = FS, XB = XB, ST = ST)
  Zn <- P_("rescale_objectives")(tb, "anchored"); d <- P_("d_ideal")(Zn, W_OBJ)
  o <- order(d); k <- o[1]; gam <- d[o[2]] - d[o[1]]
  L <- log(2 * length(GRID) / DELTA); k_ex <- sqrt(W_OBJ[["ST"]]) * (2 * sqrt(L / (2 * B))) / 2
  list(k = k, gamma = gam, cert = 4 * k_ex < gam, set = which(d <= min(d) + 2 * k_ex), k_ex = k_ex)
}
target <- select_cert(colMeans(st_pool), nrow(st_pool))
k_star <- target$k
cat(sprintf("Target (POOL = %d): beta* = %.2f, margin = %.4f\n", nrow(st_pool), GRID[k_star], target$gamma))

## ---------- 3. bank kecil berulang ----------
set.seed(2026)
cov <- do.call(rbind, lapply(B_SET, function(B) {
  r <- t(replicate(N_REP, {
    idx <- sample(nrow(st_pool), B)
    s <- select_cert(colMeans(st_pool[idx, , drop = FALSE]), B)
    sd_max <- max(apply(st_pool[idx, , drop = FALSE], 2, sd))
    emp_ok <- 4 * sqrt(W_OBJ[["ST"]]) * (1.96 * sd_max / sqrt(B)) / 2 < s$gamma
    c(correct = s$k == k_star, certified = s$cert, false_cert = s$cert && s$k != k_star,
      covered = k_star %in% s$set, set_size = length(s$set), emp_pass = emp_ok, gamma = s$gamma)
  }))
  data.frame(B = B, P_correct = mean(r[, "correct"]), P_certified = mean(r[, "certified"]),
             P_false_cert = mean(r[, "false_cert"]), coverage = mean(r[, "covered"]),
             mean_set_size = mean(r[, "set_size"]), P_empirical_pass = mean(r[, "emp_pass"]),
             median_gamma = median(r[, "gamma"]))
}))
cov[, -1] <- round(cov[, -1], 3)
print(cov, row.names = FALSE)
write.csv(cov, file.path(OUT, "tableS20_certificate_coverage_app1.csv"), row.names = FALSE)
cat("\nNominal: coverage >= 0.95 dan false certification <= 0.05 (delta = 0.05).\n")

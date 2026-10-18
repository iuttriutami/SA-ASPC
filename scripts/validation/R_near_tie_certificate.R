## =============================================================================
## R_near_tie_certificate.R — validasi sertifikat pada keputusan hampir seri (Table S21)
##
## Menjawab reviewer 1.2: (a) target harus independen dari bank evaluasi, (b) uji harus
## mencakup margin kecil. Desain semi-sintetis:
##   * Populasi = distribusi empiris 2.000 anggota pool (ARI per titik grid) dari
##     R_coverage_certificate.R. Bank evaluasi diambil DENGAN pengembalian dari populasi ini,
##     sehingga target beta* dihitung eksak dari rata-rata populasi (bukan aproksimasi MC)
##     dan setiap bank independen dari target.
##   * Margin diatur dengan menggeser jarak ideal-point d(beta) secara deterministik
##     (offset konstan per titik grid). Offset tidak mengubah sebaran ST, sehingga asumsi
##     sertifikat (Section 3.4) tetap berlaku persis.
##   * Dua arah: "toward0" (beta* = 0, pesaing terdekat titik kedua) dan "reversed"
##     (pesaing dijadikan beta* dengan margin gamma, beta = 0 menjadi runner-up).
##
## Input: output versi R_coverage_certificate.R (folder cov_ckpt/ berisi pool_*.rds dan D_app1.rds),
##        Panel_JawaTengah_2022_2025.xlsx, Wq_frozen.csv, E8_canonical_fits.rds, SA_ASPC_pipeline_v2.R
## Output: /kaggle/working/neartie_out/tableS21_near_tie.csv  (+ tableS21b_detail.csv)
## Waktu: beberapa menit (tidak ada fitting ulang).
## =============================================================================

GAMMAS <- c(0.002, 0.005, 0.01, 0.02, 0.05, 0.10, 0.15, 0.20)
B_SET  <- c(25, 50, 100, 200, 400, 1000, 2000, 5000)
N_REP  <- 2000
W_OBJ  <- c(FM = 0.35, FS = 0.25, XB = 0.20, ST = 0.20)
DELTA  <- 0.05
GRID   <- seq(0, 6, by = 0.25)
M_FUZZ <- 2

WORK <- "/kaggle/working"; CK <- file.path(WORK, "cov_ckpt"); OUT <- file.path(WORK, "neartie_out")
dir.create(CK, showWarnings = FALSE); dir.create(OUT, showWarnings = FALSE)
inp <- function(pat) { f <- list.files("/kaggle/input", pat, recursive = TRUE, full.names = TRUE, ignore.case = TRUE); if (length(f)) f[1] else NA }
ck_load <- function(name) {                         # hanya MEMBACA checkpoint lama
  f <- list.files("/kaggle/input", paste0("^", name, "\\.rds$"), recursive = TRUE, full.names = TRUE)
  f <- f[basename(dirname(f)) == "cov_ckpt"]
  if (!length(f)) stop("Checkpoint ", name, " tidak ditemukan. Pasang output versi R_coverage_certificate.R sebagai Input.")
  readRDS(f[1])
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

## ---------- objektif deterministik (sama dengan R_coverage_certificate.R) ----------
raw <- as.data.frame(readxl::read_excel(inp("Panel_JawaTengah.*\\.xlsx$")))
names(raw) <- toupper(trimws(names(raw)))
IND <- c("APM_SD", "APM_SMP", "PRODPADI", "RLS", "SANITASI")
tmp <- file.path(WORK, "panel1_tmp.csv"); write.csv(raw, tmp, row.names = FALSE)
X1 <- P_("zscore_panel")(P_("load_panel")(tmp, "WADMKK", "TAHUN", IND))
D1 <- ck_load("D_app1"); reg <- dimnames(X1)[[1]]
W1 <- as.matrix(read.csv(inp("Wq_frozen\\.csv$"), row.names = 1, check.names = FALSE))
nk <- function(v) { v <- toupper(trimws(sub("^X", "", v))); v <- sub("^KABUPATEN\\s+|^KAB\\.?\\s+", "", v); gsub("[^A-Z]", "", v) }
W1 <- W1[match(nk(reg), nk(rownames(W1))), match(nk(reg), nk(colnames(W1)))]
dimnames(W1) <- list(reg, reg); Pm <- P_("row_standardise")(W1)
Zbar <- apply(X1, c(1, 3), mean)
canon <- readRDS(inp("E8_canonical_fits\\.rds$"))
FM <- sapply(canon, function(f) P_("fuzzy_moran")(f$U, Zbar, Pm, M_FUZZ))
FS <- sapply(canon, function(f) P_("fuzzy_silhouette")(f$U, D1))
XB <- sapply(canon, function(f) P_("xie_beni")(f$U, f$medoids, D1, M_FUZZ))

pool <- do.call(rbind, lapply(1:20, function(j) ck_load(sprintf("pool_%03d", j))))
cat("Pool:", nrow(pool), "x", ncol(pool), "\n")
G <- length(GRID); L <- log(2 * G / DELTA)
d_of <- function(ST) {
  tb <- data.frame(beta = GRID, FM = FM, FS = FS, XB = XB, ST = ST)
  P_("d_ideal")(P_("rescale_objectives")(tb, "anchored"), W_OBJ)
}
d_pop <- d_of(colMeans(pool))
o <- order(d_pop); k1 <- o[1]; k2 <- o[2]
cat(sprintf("Populasi: beta* = %.2f, runner-up = %.2f, margin asli = %.4f\n", GRID[k1], GRID[k2], d_pop[k2] - d_pop[k1]))

## offset deterministik yang menghasilkan margin populasi = gamma
make_offset <- function(gamma, dir) {
  off <- numeric(G)
  if (dir == "toward0") {                      # k1 tetap target; tarik k2 mendekat, titik lain ditahan >= gamma
    off[k2] <- (d_pop[k1] + gamma) - d_pop[k2]
    others <- setdiff(seq_len(G), c(k1, k2))
    off[others] <- pmax(0, (d_pop[k1] + gamma) - d_pop[others])
  } else {                                     # k2 jadi target dengan margin gamma di atas k1
    off[k2] <- (d_pop[k1] - gamma) - d_pop[k2]
    others <- setdiff(seq_len(G), c(k1, k2))
    off[others] <- pmax(0, d_pop[k1] - d_pop[others])
  }
  off
}

sim_cell <- function(gamma, dir, B) {
  off <- make_offset(gamma, dir); dp <- d_pop + off
  kstar <- which.min(dp); stopifnot(abs(sort(dp)[2] - sort(dp)[1] - gamma) < 1e-9)
  kex <- sqrt(W_OBJ[["ST"]]) * sqrt(L / (2 * B))             # = sqrt(w4) * eps / 2
  r <- t(replicate(N_REP, {
    idx <- sample.int(nrow(pool), B, replace = TRUE)
    stb <- pool[idx, , drop = FALSE]
    d <- d_of(colMeans(stb)) + off
    oo <- order(d); k <- oo[1]; gam <- d[oo[2]] - d[oo[1]]
    cert <- 4 * kex < gam
    emp <- 4 * sqrt(W_OBJ[["ST"]]) * (1.96 * max(apply(stb, 2, sd)) / sqrt(B)) / 2 < gam
    set <- which(d <= min(d) + 2 * kex)
    c(correct = k == kstar, cert = cert, false_cert = cert && k != kstar,
      emp = emp, false_emp = emp && k != kstar, covered = kstar %in% set, size = length(set))
  }))
  data.frame(direction = dir, gamma = gamma, B = B, B_required = ceiling(8 * W_OBJ[["ST"]] * L / gamma^2),
             P_correct = mean(r[, "correct"]), P_certified = mean(r[, "cert"]),
             P_false_cert = mean(r[, "false_cert"]),
             P_false_cert_given_cert = if (any(r[, "cert"] == 1)) sum(r[, "false_cert"]) / sum(r[, "cert"]) else NA,
             coverage = mean(r[, "covered"]), mean_set_size = mean(r[, "size"]),
             P_empirical_pass = mean(r[, "emp"]), P_false_empirical = mean(r[, "false_emp"]))
}

set.seed(20261018)
res <- do.call(rbind, lapply(c("toward0", "reversed"), function(dir)
  do.call(rbind, lapply(GAMMAS, function(g) do.call(rbind, lapply(B_SET, function(B) sim_cell(g, dir, B)))))))
num <- sapply(res, is.numeric); res[num] <- lapply(res[num], function(v) round(v, 4))
write.csv(res, file.path(OUT, "tableS21b_detail.csv"), row.names = FALSE)

## ringkasan per gamma (dipooled atas arah dan B)
summ <- do.call(rbind, lapply(split(res, list(res$direction, res$gamma), drop = TRUE), function(s) data.frame(
  direction = s$direction[1], gamma = s$gamma[1], B_required = s$B_required[1],
  min_coverage = min(s$coverage), max_false_cert = max(s$P_false_cert),
  max_false_empirical = max(s$P_false_empirical),
  P_certified_at_Bmax = s$P_certified[s$B == max(s$B)],
  P_correct_at_B100 = s$P_correct[s$B == 100])))
summ <- summ[order(summ$direction, summ$gamma), ]
print(summ, row.names = FALSE)
write.csv(summ, file.path(OUT, "tableS21_near_tie.csv"), row.names = FALSE)
cat("\nNominal (delta = 0.05): coverage >= 0.95 dan false certification <= 0.05 di SETIAP sel.\n")
cat("Sel terburuk: coverage min =", min(res$coverage), "| false cert maks =", max(res$P_false_cert),
    "| false empirical maks =", max(res$P_false_empirical), "\n")

## =====================================================================
## R_aplikasi_pertama_C.R
## Aplikasi pertama SA-ASPC (35 kab/kota Jawa Tengah, 2022-2025) dengan C ditetapkan,
## memakai pipeline yang sama dengan aplikasi kedua (fungsi ASLI SA_ASPC_pipeline_v2.R).
##
## Jalankan dua kali (Save Version per run, folder output terpisah otomatis):
##   C_FIXED <- 8  -> app1_ckpt_C8/, app1_out_C8/  (cek: harus mereproduksi beta_hat = 0 naskah)
##   C_FIXED <- 4  -> app1_ckpt_C4/, app1_out_C4/  (pembanding untuk Table C = 4 vs C = 8)
## Input yang dibutuhkan: Panel_JawaTengah_2022_2025.xlsx, Wq_frozen.csv, SA_ASPC_pipeline_v2.R
## =====================================================================

## ---------------- PENGATURAN -----------------------------------------
PANEL_PATTERN <- "Panel_JawaTengah.*\\.(csv|xlsx)$"
COL_REGION    <- "WADMKK"
COL_NAME      <- "WADMKK"
COL_TIME      <- "TAHUN"
INDICATORS    <- c("APM_SD", "APM_SMP", "PRODPADI", "RLS", "SANITASI")

W_MODE        <- "frozen"        # "frozen" (CSV 0/1, row.names = kunci) atau "shapefile"
W_FROZEN_PAT  <- "Wq_frozen\\.csv$"         # matriks queen beku aplikasi 1 (naskah)
SHP_PATTERN   <- "\\.shp$"
SHP_KEY       <- "KODE_KAB"

C_FIXED       <- 4               # jalankan dua kali: 8 (cek replikasi naskah) dan 4 (pembanding)
C_CAND        <- 2:10
M_FUZZ2       <- 2
GRID_MODE     <- "absolute"      # grid naskah {0, 0.25, ..., 6}
BASE_GRID     <- seq(0, 6, by = 0.25)
W_OBJ2        <- c(FM = 0.35, FS = 0.25, XB = 0.20, ST = 0.20)

R_DEEP2  <- 100; N_PERT2 <- 6; SWEEPS2 <- 2; N_MASTER2 <- 5; PERT_SD <- 0.30
R_SHALLOW <- 20; N_SHALLOW_SEEDS <- 20
B0 <- 100; B_MAX <- 400; SIGMA2 <- 0.05; BANK_RESTART <- 3; DELTA <- 0.05
CORES <- max(1, parallel::detectCores())
## ---------------------------------------------------------------------

WORK <- "/kaggle/working"
TAG  <- if (is.na(C_FIXED)) "" else sprintf("_C%d", C_FIXED)
CK   <- file.path(WORK, paste0("app1_ckpt", TAG)); OUT <- file.path(WORK, paste0("app1_out", TAG))
dir.create(CK, showWarnings = FALSE); dir.create(OUT, showWarnings = FALSE)
say  <- function(...) cat(sprintf(...), "\n")
gate <- function(lbl, ok, det = "") say("  [%s] %s %s", if (isTRUE(ok)) "PASS" else "CEK", lbl, det)
inp  <- function(pat) { f <- list.files("/kaggle/input", pat, recursive = TRUE, full.names = TRUE, ignore.case = TRUE); if (length(f)) f[1] else NA }
ck   <- function(name, expr) {                    # checkpoint: hitung sekali, simpan
  f <- file.path(CK, paste0(name, ".rds"))
  old <- list.files("/kaggle/input", paste0("^", name, "\\.rds$"), recursive = TRUE, full.names = TRUE)
  old <- if (length(old <- old[basename(dirname(old)) == basename(CK)])) old[1] else NA   # hanya folder dengan TAG yang sama
  if (file.exists(f)) return(readRDS(f))
  if (!is.na(old)) { file.copy(old, f); return(readRDS(f)) }
  v <- expr; saveRDS(v, f); v                     # 'expr' dievaluasi malas, hanya bila perlu
}
say("Run: %s (checkpoint %s, output %s)", if (TAG == "") "utama" else sub("_", "", TAG), basename(CK), basename(OUT))

## ---------------------------------------------------------------------
## 0. Fungsi pipeline
## ---------------------------------------------------------------------
pipe <- new.env()
for (e in parse(inp("SA_ASPC_pipeline_v2\\.R$"), keep.source = FALSE)) {
  if (is.call(e) && as.character(e[[1]]) %in% c("<-", "=") && is.name(e[[2]])) {
    nm <- as.character(e[[2]]); rhs <- e[[3]]
    if ((is.call(rhs) && identical(rhs[[1]], as.name("function"))) || grepl("^[A-Z][A-Z0-9_]*$", nm))
      if (!nm %in% c("X", "D")) tryCatch(eval(e, pipe), error = function(err) NULL)
  }
}
P_ <- function(nm) get(nm, envir = pipe)

## ---------------------------------------------------------------------
## S1. Data
## ---------------------------------------------------------------------
f_panel <- inp(PANEL_PATTERN); if (is.na(f_panel)) stop("Panel tidak ditemukan: ", PANEL_PATTERN)
raw <- if (grepl("xlsx$", f_panel, ignore.case = TRUE)) as.data.frame(readxl::read_excel(f_panel)) else read.csv(f_panel, check.names = FALSE)
names(raw) <- toupper(trimws(names(raw)))
say("File panel: %s", basename(f_panel)); say("Kolom panel: %s", paste(names(raw), collapse = ", "))
pick_col <- function(want, cands, label) {
  hit <- intersect(toupper(c(want, cands)), names(raw))
  if (!length(hit)) stop(sprintf("Kolom %s tidak ditemukan. Isi pengaturan dengan salah satu kolom di atas.", label))
  if (hit[1] != toupper(want)) say("  [INFO] kolom %s: '%s' tidak ada, dipakai '%s'", label, want, hit[1])
  hit[1]
}
COL_REGION <- pick_col(COL_REGION, c("KODE_KAB", "KODE", "KDKAB", "KD_KAB", "IDKAB", "KODE_BPS", "KAB_ID",
                                     "WADMKK", "KABKOT", "NAMA_KAB", "KABUPATEN", "NAMA"), "unit")
COL_TIME   <- pick_col(COL_TIME, c("TAHUN", "YEAR", "TH", "PERIODE"), "waktu")
COL_NAME   <- if (toupper(COL_NAME) %in% names(raw)) toupper(COL_NAME) else COL_REGION
INDICATORS <- toupper(INDICATORS)
miss_ind <- setdiff(INDICATORS, names(raw))
if (length(miss_ind)) stop("Indikator tidak ditemukan: ", paste(miss_ind, collapse = ", "),
                           ". Sesuaikan INDICATORS dengan nama kolom di atas.")
KEY_IS_NAME <- !grepl("KODE|KD|ID", COL_REGION)
norm_name <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x <- gsub("^KABUPATEN\\s+|^KAB\\.?\\s+", "", x)
  x <- gsub("^KOTA\\s+ADM(INISTRASI)?\\.?\\s+", "KOTA ", x)
  gsub("[^A-Z ]", "", gsub("\\s+", " ", x))
}
if (KEY_IS_NAME) {
  say("  [PERHATIAN] kunci unit berbasis nama ('%s'). Pastikan kabupaten dan kota bernama sama dibedakan.", COL_REGION)
  raw[[COL_REGION]] <- norm_name(raw[[COL_REGION]])
}
raw[[COL_REGION]] <- as.character(raw[[COL_REGION]])
tab <- table(raw[[COL_REGION]]); TT <- length(unique(raw[[COL_TIME]]))
say("Panel: %d unit, %d periode (%s), %d indikator", length(tab), TT,
    paste(range(raw[[COL_TIME]]), collapse = "-"), length(INDICATORS))
gate("panel seimbang", all(tab == TT)); gate("tanpa NA", !anyNA(raw[, INDICATORS]),
     sprintf("(%d NA)", sum(is.na(raw[, INDICATORS]))))
if (any(tab != TT) || anyNA(raw[, INDICATORS])) stop("Lengkapi/seimbangkan panel dulu (lihat daftar unit di atas).")
if (!(length(tab) == 35 && TT == 4)) stop("Ini bukan panel aplikasi 1 (harus 35 unit x 2022-2025).")
tmp <- file.path(WORK, "panel2_tmp.csv"); write.csv(raw, tmp, row.names = FALSE)
X_raw <- P_("load_panel")(tmp, COL_REGION, COL_TIME, INDICATORS)
X <- P_("zscore_panel")(X_raw); regions <- dimnames(X)[[1]]; I <- length(regions)
D <- ck("D", P_("dtw_matrix")(X, TRUE))
Zbar <- apply(X, c(1, 3), mean)
mD <- mean(D[upper.tri(D)]); say("I = %d, T = %d, rerata DTW off-diagonal = %.3f", I, TT, mD)
names_lookup <- if (COL_NAME %in% names(raw)) tapply(raw[[COL_NAME]], raw[[COL_REGION]], `[`, 1)[regions] else regions

## ---------------------------------------------------------------------
## S2. Matriks spasial (beku dari R_buat_W_jawa.R)
## ---------------------------------------------------------------------
Wq <- ck("Wq", {
  if (W_MODE == "frozen") {
    fW <- c(list.files(WORK, W_FROZEN_PAT, full.names = TRUE), inp(W_FROZEN_PAT)); fW <- fW[!is.na(fW)][1]
    if (is.na(fW)) stop("Matriks beku tidak ditemukan (", W_FROZEN_PAT, "); jalankan R_buat_W_jawa.R dulu.")
    W <- as.matrix(read.csv(fW, row.names = 1, check.names = FALSE))
    rownames(W) <- colnames(W) <- sub("^X", "", rownames(W))
    if (KEY_IS_NAME) rownames(W) <- colnames(W) <- norm_name(rownames(W))
    miss <- setdiff(regions, rownames(W))
    if (length(miss)) stop("Unit panel tidak ada di matriks beku: ", paste(miss, collapse = ", "))
    W[regions, regions]
  } else {
    library(sf); library(spdep)
    shp <- sf::st_read(inp(SHP_PATTERN), quiet = TRUE)
    names(shp)[names(shp) != attr(shp, "sf_column")] <- toupper(names(shp)[names(shp) != attr(shp, "sf_column")])
    say("Kolom shapefile: %s", paste(setdiff(names(shp), attr(shp, "sf_column")), collapse = ", "))
    SHP_KEY <- toupper(SHP_KEY)
    if (!SHP_KEY %in% names(shp)) {
      cand <- intersect(if (KEY_IS_NAME) c("WADMKK", "KABKOT", "NAMOBJ", "NAMA_KAB", "NAME_2", "KABUPATEN")
                        else c("KODE_KAB", "KDKAB", "KD_KAB", "KDPKAB", "KODE", "IDKAB", "KODE_BPS"), names(shp))
      if (!length(cand)) stop("Kolom kunci shapefile tidak ditemukan; isi SHP_KEY dengan salah satu kolom di atas.")
      say("  [INFO] SHP_KEY '%s' tidak ada, dipakai '%s'", SHP_KEY, cand[1]); SHP_KEY <- cand[1]
    }
    shp[[SHP_KEY]] <- as.character(shp[[SHP_KEY]])
    if (KEY_IS_NAME) {
      shp[[SHP_KEY]] <- norm_name(shp[[SHP_KEY]])
      kcol <- intersect(c("KABKOTNO", "KDKAB", "KODE_KAB", "KD_KAB"), names(shp))[1]
      if (!is.na(kcol)) {
        kk <- suppressWarnings(as.integer(substr(gsub("\\D", "", as.character(shp[[kcol]])), nchar(gsub("\\D", "", as.character(shp[[kcol]]))) - 1, 99)))
        is_kota <- !is.na(kk) & kk >= 71 & !grepl("^KOTA ", shp[[SHP_KEY]])
        shp[[SHP_KEY]][is_kota] <- paste("KOTA", shp[[SHP_KEY]][is_kota])
        say("  [INFO] %d unit kota dikenali dari kolom %s (kode >= 71)", sum(kk >= 71, na.rm = TRUE), kcol)
      }
    }
    miss <- setdiff(regions, shp[[SHP_KEY]])
    if (length(miss)) {
      extra <- setdiff(shp[[SHP_KEY]], regions)
      stop("Unit panel tidak ditemukan di shapefile: ", paste(miss, collapse = ", "),
           "\n  Kandidat di shapefile yang tidak terpakai: ", paste(head(extra, 30), collapse = ", "))
    }
    if (anyDuplicated(shp[[SHP_KEY]][shp[[SHP_KEY]] %in% regions]))
      stop("Kunci shapefile tidak unik untuk unit panel (mis. kab & kota bernama sama). Pakai kode BPS.")
    shp <- shp[match(regions, shp[[SHP_KEY]]), ]
    shp <- sf::st_make_valid(shp)
    nb <- spdep::poly2nb(shp, queen = TRUE, snap = 1e-6)
    W <- spdep::nb2mat(nb, style = "B", zero.policy = TRUE)
    iso <- which(rowSums(W) == 0)
    if (length(iso)) {
      cen <- sf::st_coordinates(sf::st_point_on_surface(sf::st_geometry(shp)))
      for (i in iso) { d <- sqrt(colSums((t(cen) - cen[i, ])^2)); d[i] <- Inf; j <- which.min(d); W[i, j] <- W[j, i] <- 1 }
      say("  %d unit terisolasi dihubungkan ke tetangga terdekat: %s", length(iso), paste(names_lookup[iso], collapse = ", "))
    }
    dimnames(W) <- list(regions, regions); W
  }
})
write.csv(Wq, file.path(OUT, "W_app1_frozen.csv"))
gate("matriks simetris tanpa unit terisolasi", isSymmetric(unname(Wq)) && all(rowSums(Wq) > 0),
     sprintf("(rerata derajat %.2f)", mean(rowSums(Wq))))
if (requireNamespace("spdep", quietly = TRUE)) {
  nc <- spdep::n.comp.nb(suppressWarnings(spdep::mat2listw(Wq, style = "B", zero.policy = TRUE))$neighbours)$nc
  gate("graf terhubung (1 komponen)", nc == 1, sprintf("(%d komponen)", nc))
}
Pm <- P_("row_standardise")(Wq)

## ---------------------------------------------------------------------
## S3. Jumlah klaster (rank sum kriteria internal di beta = 0, seperti naskah)
## ---------------------------------------------------------------------
csel <- ck("Cselect", P_("select_C")(D, Pm, candidates = C_CAND, m = M_FUZZ2))
print(csel)
write.csv(as.data.frame(csel), file.path(OUT, "tableS4_app1_Cselect.csv"), row.names = FALSE)
C2 <- if (!is.na(C_FIXED)) C_FIXED else {
  cs <- as.data.frame(csel)
  rc <- grep("rank.?sum|ranksum|totalrank|total.?rank", names(cs), value = TRUE, ignore.case = TRUE)[1]
  cc <- grep("^C$|^K$|candidate", names(cs), value = TRUE, ignore.case = TRUE)[1]
  if (is.na(rc) || is.na(cc)) stop("Struktur select_C tidak dikenali; isi C_FIXED secara manual.")
  tie_col <- intersect(c("BootARI", "BOOTARI"), names(cs))[1]          # pemecah seri: stabilitas tertinggi
  o <- if (!is.na(tie_col)) order(cs[[rc]], -cs[[tie_col]]) else order(cs[[rc]])
  as.integer(cs[[cc]][o[1]])
}
say("C = %d%s", C2, if (is.na(C_FIXED)) " (otomatis)" else " (ditetapkan)")

GRID <- if (GRID_MODE == "relative") BASE_GRID * mD / 2.958 else BASE_GRID
G <- length(GRID); say("Grid: %d titik, beta maks = %.3f (%s)", G, max(GRID), GRID_MODE)

## ---------------------------------------------------------------------
## S4. Protokol deep: restart + sweep kontinuasi, per master seed, lalu pool
## ---------------------------------------------------------------------
perturb <- function(U) { V <- U * exp(matrix(rnorm(length(U), 0, PERT_SD), nrow(U))); V / rowSums(V) }
deep_seed <- function(seed) {
  set.seed(seed)
  best <- lapply(seq_len(G), function(g) P_("fcmd_stt")(D, Pm, GRID[g], C2, M_FUZZ2, n_restart = R_DEEP2, seed = seed + g))
  for (sw in seq_len(SWEEPS2)) for (ord in list(seq_len(G), rev(seq_len(G)))) for (k in seq_along(ord)[-1]) {
    g <- ord[k]; U0p <- best[[ord[k - 1]]]$U
    for (p in seq_len(N_PERT2)) {
      f <- P_("fcmd_stt_fit1")(D, Pm, GRID[g], C2, M_FUZZ2, if (p == 1) U0p else perturb(U0p))
      if (f$J < best[[g]]$J) best[[g]] <- f
    }
  }
  best
}
seeds <- 20260905 + 1000 * seq_len(N_MASTER2)
say("S4: protokol deep (%d master seed, %d restart, %d sweep) ...", N_MASTER2, R_DEEP2, SWEEPS2)
deep_all <- parallel::mclapply(seq_along(seeds), function(k) ck(sprintf("deep_seed%02d", k), deep_seed(seeds[k])),
                               mc.cores = min(CORES, N_MASTER2), mc.preschedule = FALSE)
bad <- which(vapply(deep_all, function(x) inherits(x, "try-error") || is.null(x), logical(1)))
if (length(bad)) stop("Seed deep gagal: ", paste(bad, collapse = ", "), ". Jalankan ulang (checkpoint seed lain aman).")
Jmat <- sapply(deep_all, function(fs) sapply(fs, `[[`, "J"))           # G x seeds
canon <- lapply(seq_len(G), function(g) deep_all[[which.min(Jmat[g, ])]][[g]])
share_best <- rowMeans(abs(Jmat - apply(Jmat, 1, min)) <= 1e-8 * apply(Jmat, 1, min))
saveRDS(canon, file.path(OUT, "E8_app1_canonical_fits.rds"))

## ---------------------------------------------------------------------
## S5. Protokol dangkal (praktik umum) untuk dibandingkan
## ---------------------------------------------------------------------
say("S5: protokol dangkal (%d seed, %d restart) ...", N_SHALLOW_SEEDS, R_SHALLOW)
shallow <- ck("shallow", parallel::mclapply(seq_len(N_SHALLOW_SEEDS), function(s)
  lapply(seq_len(G), function(g) P_("fcmd_stt")(D, Pm, GRID[g], C2, M_FUZZ2, n_restart = R_SHALLOW, seed = 7e6 + 100 * s + g)),
  mc.cores = CORES))

## ---------------------------------------------------------------------
## S6. Bank, objektif, aturan seleksi
## ---------------------------------------------------------------------
bank_labels <- function(B, seed, tag) ck(tag, {
  bank <- P_("build_bank")(X, B = B, sigma = SIGMA2, normalise = TRUE, seed = seed)
  parallel::mclapply(seq_len(G), function(g) lapply(bank, function(Db)
    P_("hard")(P_("fcmd_stt")(Db, Pm, GRID[g], C2, M_FUZZ2, n_restart = BANK_RESTART)$U)), mc.cores = CORES)
})
st_from <- function(fits, labs) sapply(seq_len(G), function(g) {
  ref <- P_("hard")(fits[[g]]$U); sapply(labs[[g]], function(l) P_("ari")(ref, l))
})
objectives <- function(fits, st_rep) data.frame(
  beta = GRID,
  FM = sapply(fits, function(f) P_("fuzzy_moran")(f$U, Zbar, Pm, M_FUZZ2)),
  FS = sapply(fits, function(f) P_("fuzzy_silhouette")(f$U, D)),
  XB = sapply(fits, function(f) P_("xie_beni")(f$U, f$medoids, D, M_FUZZ2)),
  ST = colMeans(st_rep), ST_sd = apply(st_rep, 2, sd), NDC = 0,
  SC = sapply(fits, function(f) P_("spatial_coherence")(f$U, Pm)),
  J = sapply(fits, `[[`, "J"))
select_all <- function(tb, w = W_OBJ2) {
  Zn <- P_("rescale_objectives")(tb, "anchored"); d <- P_("d_ideal")(Zn, w)
  o <- order(d); k <- o[1]
  noST <- P_("d_ideal")(Zn, c(w[1:3], 0) / sum(w[1:3]))
  eq <- P_("d_ideal")(Zn, rep(0.25, 4))
  dom_ok <- !(tb$XB[k] > tb$XB[1] && tb$ST[k] < tb$ST[1])
  list(d = d, k = k, gamma = d[o[2]] - d[o[1]],
       rules = data.frame(rule = c("SA-ASPC", "SA-ASPC equal weights", "SA-ASPC without ST", "Moran max", "XB min",
                                   "Guard dominance", "Guard XB 2%"),
                          beta = c(GRID[k], GRID[which.min(eq)], GRID[which.min(noST)], GRID[which.max(tb$FM)],
                                   GRID[which.min(tb$XB)], if (dom_ok) GRID[k] else 0,
                                   if (tb$XB[k] <= 1.02 * tb$XB[1]) GRID[k] else 0)),
       nondominated = sum(P_("pareto_flag")(Zn)))
}

say("S6: bank B = %d ...", B0)
labs0 <- bank_labels(B0, 12345, "bank_main")
st0 <- st_from(canon, labs0)
tb  <- objectives(canon, st0); sel <- select_all(tb)
say("Deep: beta_hat = %.3f, margin = %.4f, nondominated %d/%d", GRID[sel$k], sel$gamma, sel$nondominated, G)

## ---------------------------------------------------------------------
## S7. Sertifikasi (konstanta anchored), bank diperbesar bila perlu
## ---------------------------------------------------------------------
certify2 <- function(tb, B, gam, d, w4 = W_OBJ2[["ST"]]) {
  L <- log(2 * G / DELTA); eps <- 2 * sqrt(L / (2 * B))
  k_ex <- sqrt(w4) * eps / 2
  k_em <- sqrt(w4) * (1.96 * max(tb$ST_sd) / sqrt(B)) / 2
  data.frame(B = B, gamma_hat = gam, fourk_emp = 4 * k_em, ratio_emp = 4 * k_em / gam,
             fourk_exact = 4 * k_ex, ratio_exact = 4 * k_ex / gam, certified_exact = 4 * k_ex < gam,
             conf_set_exact = paste(round(GRID[d <= min(d) + 2 * k_ex], 3), collapse = " "),
             B_required_exact = ceiling(8 * w4 * L / gam^2))
}
cert <- certify2(tb, B0, sel$gamma, sel$d); print(cert)
if (!cert$certified_exact && cert$B_required_exact <= B_MAX) {
  B_add <- cert$B_required_exact - B0
  say("Memperbesar bank: +%d anggota (total %d)", B_add, cert$B_required_exact)
  labs_add <- bank_labels(B_add, 54321, sprintf("bank_add_%d", B_add))
  st_add <- st_from(canon, labs_add)
  st_all <- rbind(st0, st_add)
  tb_big <- objectives(canon, st_all); sel_big <- select_all(tb_big)
  cert <- rbind(cert, certify2(tb_big, nrow(st_all), sel_big$gamma, sel_big$d)); print(cert)
} else if (!cert$certified_exact) say("B yang dibutuhkan (%d) > B_MAX (%d): laporkan sebagai tidak tersertifikasi.", cert$B_required_exact, B_MAX)

## ---------------------------------------------------------------------
## S8. Kedalaman optimisasi, bobot, korelasi, peta
## ---------------------------------------------------------------------
say("S8: perbandingan protokol, wilayah bobot, peta ...")
## (a) pilihan dangkal per seed (memakai bank yang sama -> beda hanya dari fit utama)
sh_tab <- ck("shallow_selection", do.call(rbind, lapply(seq_along(shallow), function(s) {
  st_s <- st_from(shallow[[s]], labs0)
  ss <- select_all(objectives(shallow[[s]], st_s))
  data.frame(seed = s, beta_hat = GRID[ss$k], gamma = ss$gamma, moran_beta = ss$rules$beta[4])
})))
deep_sel <- do.call(rbind, lapply(seq_along(deep_all), function(k) {
  ss <- select_all(objectives(deep_all[[k]], st_from(deep_all[[k]], labs0))); data.frame(seed = k, beta_hat = GRID[ss$k], gamma = ss$gamma)
}))

## (b) wilayah bobot dan radius (Theorem 14)
Zn <- P_("rescale_objectives")(tb, "anchored"); s2 <- (1 - Zn)^2
rdir <- function(n, a) { x <- matrix(rgamma(n * length(a), a), n, byrow = TRUE); x / rowSums(x) }
set.seed(2026)
reg_flat  <- table(factor(GRID[max.col(-(rdir(20000, rep(1, 4)) %*% t(s2)), "first")], levels = GRID)) / 20000
reg_local <- table(factor(GRID[max.col(-(rdir(20000, 50 * W_OBJ2) %*% t(s2)), "first")], levels = GRID)) / 20000
g0 <- sel$k
rad <- min(sapply(setdiff(seq_len(G), g0), function(h) { a <- s2[h, ] - s2[g0, ]; sum(W_OBJ2 * a) / sqrt(sum((a - mean(a))^2)) }))
spear <- round(cor(tb[, c("beta", "FM", "FS", "XB", "ST")], method = "spearman"), 2)

## (c) simpan semua tabel
write.csv(tb, file.path(OUT, "tableS5_app1_grid.csv"), row.names = FALSE)
write.csv(sel$rules, file.path(OUT, "tableS6_app1_rules.csv"), row.names = FALSE)
write.csv(cert, file.path(OUT, "tableS7_app1_certification.csv"), row.names = FALSE)
depth <- data.frame(protocol = c(rep("shallow", nrow(sh_tab)), rep("deep (single master seed)", nrow(deep_sel))),
                    seed = c(sh_tab$seed, deep_sel$seed), beta_hat = c(sh_tab$beta_hat, deep_sel$beta_hat),
                    gamma = c(sh_tab$gamma, deep_sel$gamma))
write.csv(depth, file.path(OUT, "tableS8_app1_depth.csv"), row.names = FALSE)
write.csv(data.frame(beta = GRID, share_seeds_best_J = share_best), file.path(OUT, "tableS9_app1_seed_agreement.csv"), row.names = FALSE)
write.csv(data.frame(beta = GRID, flat = as.numeric(reg_flat), local = as.numeric(reg_local)),
          file.path(OUT, "tableS10_app1_weight_regions.csv"), row.names = FALSE)
write.csv(spear, file.path(OUT, "tableS10b_app1_spearman.csv"))
write.csv(data.frame(region = regions, name = names_lookup,
                     cluster_beta0 = P_("hard")(canon[[1]]$U),
                     cluster_selected = P_("hard")(canon[[sel$k]]$U)),
          file.path(OUT, "app1_partitions.csv"), row.names = FALSE)

## (d) peta partisi
f_gpkg <- c(list.files(WORK, "batas_app1_jawa\\.gpkg$", full.names = TRUE), inp("batas_app1_jawa\\.gpkg$")); f_gpkg <- f_gpkg[!is.na(f_gpkg)][1]
if (requireNamespace("sf", quietly = TRUE) && (W_MODE == "shapefile" || !is.na(f_gpkg))) {
  if (W_MODE != "shapefile") { SHP_PATTERN <- "batas_app1_jawa\\.gpkg$"; SHP_KEY <- "KODE" }
  shp <- sf::st_read(if (W_MODE == "shapefile") inp(SHP_PATTERN) else f_gpkg, quiet = TRUE)
  names(shp)[names(shp) != attr(shp, "sf_column")] <- toupper(names(shp)[names(shp) != attr(shp, "sf_column")])
  shp[[SHP_KEY]] <- as.character(shp[[SHP_KEY]])
  if (KEY_IS_NAME) {
    shp[[SHP_KEY]] <- norm_name(shp[[SHP_KEY]])
    kcol <- intersect(c("KABKOTNO", "KDKAB", "KODE_KAB", "KD_KAB"), names(shp))[1]
    if (!is.na(kcol)) {
      dd <- gsub("\\D", "", as.character(shp[[kcol]])); kk <- suppressWarnings(as.integer(substr(dd, nchar(dd) - 1, 99)))
      ik <- !is.na(kk) & kk >= 71 & !grepl("^KOTA ", shp[[SHP_KEY]]); shp[[SHP_KEY]][ik] <- paste("KOTA", shp[[SHP_KEY]][ik])
    }
  }
  shp <- shp[match(regions, shp[[SHP_KEY]]), ]
  k_alt <- if (sel$k == 1) which.max(tb$FM) else 1
  png(file.path(OUT, "figS2_app1_map.png"), width = 2400, height = 1000, res = 200)
  op <- par(mfrow = c(1, 2), mar = c(0.5, 0.5, 2, 0.5)); pal <- hcl.colors(C2, "Set 2")
  for (k in c(sel$k, k_alt)) plot(sf::st_geometry(shp), col = pal[P_("hard")(canon[[k]]$U)], border = "grey40", lwd = 0.3,
                                  main = sprintf("beta = %.3f%s", GRID[k], if (k == sel$k) " (selected)" else ""))
  par(op); dev.off()
}

## ---------------------------------------------------------------------
## Ringkasan untuk naskah
## ---------------------------------------------------------------------
say("\n===== RINGKASAN APLIKASI PERTAMA%s =====", if (TAG == "") "" else paste0(" (", sub("_", "", TAG), ")"))
say("I = %d, T = %d, C = %d, rerata DTW = %.3f, derajat = %.2f", I, TT, C2, mD, mean(rowSums(Wq)))
say("Deep (pooled %d seed): beta_hat = %.3f, margin = %.4f", N_MASTER2, GRID[sel$k], sel$gamma)
say("Deep per seed: beta_hat = {%s}", paste(sort(unique(round(deep_sel$beta_hat, 3))), collapse = ", "))
say("Dangkal (%d seed): beta_hat = {%s}; proporsi beta_hat = deep: %.2f",
    nrow(sh_tab), paste(sort(unique(round(sh_tab$beta_hat, 3))), collapse = ", "), mean(abs(sh_tab$beta_hat - GRID[sel$k]) < 1e-9))
say("Moran max (deep) = %.3f", sel$rules$beta[4])
say("FS di beta = 0: %.3f; FS di beta_hat: %.3f", tb$FS[1], tb$FS[sel$k])
say("Bobot: beta_hat dipilih %.0f%% simplex, %.1f%% sekitar w0; radius %.3f vs tepi %.3f",
    100 * reg_flat[sel$k], 100 * reg_local[sel$k], rad, min(W_OBJ2) * sqrt(4 / 3))
print(cert[nrow(cert), ])
say("Output: %s", OUT)

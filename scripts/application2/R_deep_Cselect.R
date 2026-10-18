## =============================================================================
## Pemilihan C dengan protokol deep (beta = 0) — SA-ASPC
## Tempel di cell baru SETELAH R_aplikasi_kedua.R selesai dijalankan di sesi yang sama
## (butuh objek: P_, pipe, D, Pm, ck, perturb, M_FUZZ2, CORES, OUT).
##
## Desain: kriteria dihitung dengan implementasi yang SAMA untuk fit dangkal dan fit deep,
## sehingga perbedaan peringkat hanya berasal dari kedalaman optimisasi.
##   Bagian A: aplikasi 2 (Jawa, 114 unit)          -> langsung jalan
##   Bagian B: aplikasi 1 (Jawa Tengah 2022–2025)   -> jalan hanya bila panelnya ada di Input
## =============================================================================

## ---------- Bagian 0: siapkan objek bila skrip utama belum dijalankan di sesi ini ----------
## (memuat fungsi pipeline + D + matriks W dari checkpoint app2_ckpt; tidak menghitung ulang analisis utama)
if (!all(sapply(c("P_", "D", "Pm", "ck", "perturb"), exists))) {
  cat("Objek sesi belum ada -> memuat pipeline dan checkpoint ...\n")
  WORK <- "/kaggle/working"; CK <- file.path(WORK, "app2_ckpt"); OUT <- file.path(WORK, "app2_out")
  dir.create(CK, showWarnings = FALSE, recursive = TRUE); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
  M_FUZZ2 <- 2; PERT_SD <- 0.30; CORES <- max(1, parallel::detectCores())
  inp <- function(pat) { f <- list.files("/kaggle/input", pat, recursive = TRUE, full.names = TRUE, ignore.case = TRUE); if (length(f)) f[1] else NA }
  ck <- function(name, expr) {
    f <- file.path(CK, paste0(name, ".rds"))
    old <- list.files("/kaggle/input", paste0("^", name, "\\.rds$"), recursive = TRUE, full.names = TRUE)
    old <- if (length(old <- old[basename(dirname(old)) == basename(CK)])) old[1] else NA
    if (file.exists(f)) return(readRDS(f))
    if (!is.na(old)) { file.copy(old, f); return(readRDS(f)) }
    v <- expr; saveRDS(v, f); v
  }
  f_pipe <- inp("SA_ASPC_pipeline_v2\\.R$")
  if (is.na(f_pipe)) stop("SA_ASPC_pipeline_v2.R tidak ada di Input (dataset 'pipeline v2').")
  pipe <- new.env()
  for (e in parse(f_pipe, keep.source = FALSE)) {
    if (is.call(e) && as.character(e[[1]]) %in% c("<-", "=") && is.name(e[[2]])) {
      nm <- as.character(e[[2]]); rhs <- e[[3]]
      if ((is.call(rhs) && identical(rhs[[1]], as.name("function"))) || grepl("^[A-Z][A-Z0-9_]*$", nm))
        if (!nm %in% c("X", "D")) tryCatch(eval(e, pipe), error = function(err) NULL)
    }
  }
  P_ <- function(nm) get(nm, envir = pipe)
  perturb <- function(U) { V <- U * exp(matrix(rnorm(length(U), 0, PERT_SD), nrow(U))); V / rowSums(V) }

  D <- ck("D", {                                   # biasanya langsung dimuat dari checkpoint
    f_panel <- inp("Panel_App2_Jawa.*\\.csv$"); raw <- read.csv(f_panel, check.names = FALSE)
    names(raw) <- toupper(trimws(names(raw)))
    tmp <- file.path(WORK, "panel2_tmp.csv"); write.csv(raw, tmp, row.names = FALSE)
    P_("dtw_matrix")(P_("zscore_panel")(P_("load_panel")(tmp, "KODE", "TAHUN",
                      c("IPM", "MISKIN", "TPT", "SANITASI", "LOG_PDRB_KAP"))), TRUE)
  })
  Wq <- ck("Wq", {
    fW <- inp("W_app2_frozen\\.csv$")
    if (is.na(fW)) stop("Checkpoint Wq.rds dan W_app2_frozen.csv tidak ditemukan. Pasang output versi sebelumnya sebagai Input.")
    as.matrix(read.csv(fW, row.names = 1, check.names = FALSE))
  })
  Pm <- P_("row_standardise")(Wq)
  cat("Siap: D", paste(dim(as.matrix(D)), collapse = "x"), "| W", paste(dim(Wq), collapse = "x"), "\n")
}
stopifnot(nrow(as.matrix(D)) == 114, nrow(as.matrix(Pm)) == 114)
C_SET    <- 2:10        # kandidat C (sama dengan naskah)
N_MASTER <- 20          # master seed  (= Section 5.3)
R_DEEP   <- 200         # restart per master seed
N_PERT   <- 10          # warm start terperturbasi per master seed
R_SHAL   <- 20          # protokol dangkal: best-of-20
B_SUB    <- 30          # subsampel untuk kriteria stabilitas
SUB_FRAC <- 0.8         # proporsi unit per subsampel

## ---------- kriteria internal (implementasi sama untuk dangkal & deep) ----------
hard <- function(U) P_("hard")(U)
dunn_crisp <- function(h, Dm) {
  ks <- sort(unique(h)); if (length(ks) < 2) return(NA_real_)
  diam  <- max(sapply(ks, function(k) { i <- which(h == k); if (length(i) < 2) 0 else max(Dm[i, i]) }))
  sep   <- min(combn(ks, 2, function(p) min(Dm[h == p[1], h == p[2]])))
  sep / diam
}
sub_stability <- function(f, C, Dm, Pmat, tag) {
  h_full <- hard(f$U); set.seed(4242 + C)
  idx_list <- replicate(B_SUB, sort(sample(nrow(Dm), floor(SUB_FRAC * nrow(Dm)))), simplify = FALSE)
  aris <- parallel::mclapply(seq_len(B_SUB), function(b) {
    i <- idx_list[[b]]; Pi <- Pmat[i, i]; Pi <- Pi / pmax(rowSums(Pi), 1e-12)
    g <- P_("fcmd_stt")(Dm[i, i], Pi, 0, C, M_FUZZ2, n_restart = R_SHAL, seed = 5e6 + 100 * b + C)
    P_("ari")(hard(g$U), h_full[i])
  }, mc.cores = CORES, mc.preschedule = FALSE)
  mean(unlist(aris))
}
criteria <- function(f, C, Dm, Pmat, tag) c(
  XB      = P_("xie_beni")(f$U, f$medoids, Dm, M_FUZZ2),
  FS      = P_("fuzzy_silhouette")(f$U, Dm),
  Dunn    = dunn_crisp(hard(f$U), Dm),
  SubARI  = sub_stability(f, C, Dm, Pmat, tag))

## ---------- fit dangkal dan fit deep pada beta = 0 ----------
fit_shallow <- function(C, Dm, Pmat) P_("fcmd_stt")(Dm, Pmat, 0, C, M_FUZZ2, n_restart = R_SHAL, seed = 7e6 + C)
fit_deep <- function(C, Dm, Pmat) {
  fits <- parallel::mclapply(seq_len(N_MASTER), function(s) {
    f <- P_("fcmd_stt")(Dm, Pmat, 0, C, M_FUZZ2, n_restart = R_DEEP, seed = 31e6 + 1000 * s + C)
    for (p in seq_len(N_PERT)) {
      g <- P_("fcmd_stt_fit1")(Dm, Pmat, 0, C, M_FUZZ2, perturb(f$U)); if (g$J < f$J) f <- g
    }
    f
  }, mc.cores = CORES, mc.preschedule = FALSE)
  J <- sapply(fits, `[[`, "J"); best <- fits[[which.min(J)]]
  best$share_best <- mean(abs(J - min(J)) <= 1e-8 * abs(min(J)))
  best
}

rank_sum <- function(tab) {
  r <- cbind(rank(tab$XB), rank(-tab$FS), rank(-tab$Dunn), rank(-tab$SubARI))
  rowSums(r)
}

run_Cselect <- function(Dm, Pmat, tag) {
  rows <- lapply(C_SET, function(C) {
    cat(sprintf("[%s] C = %d ...\n", tag, C))
    fs <- ck(sprintf("Cshal_%s_C%02d", tag, C), fit_shallow(C, Dm, Pmat))
    fd <- ck(sprintf("Cdeep_%s_C%02d", tag, C), fit_deep(C, Dm, Pmat))
    cs <- ck(sprintf("Ccrit_shal_%s_C%02d", tag, C), criteria(fs, C, Dm, Pmat, tag))
    cd <- ck(sprintf("Ccrit_deep_%s_C%02d", tag, C), criteria(fd, C, Dm, Pmat, tag))
    data.frame(C = C, J_shallow = fs$J, J_deep = fd$J,
               J_drop_pct = 100 * (1 - fd$J / fs$J), share_seed_best = fd$share_best,
               ARI_shallow_deep = P_("ari")(hard(fs$U), hard(fd$U)),
               XB_s = cs["XB"], FS_s = cs["FS"], Dunn_s = cs["Dunn"], SubARI_s = cs["SubARI"],
               XB_d = cd["XB"], FS_d = cd["FS"], Dunn_d = cd["Dunn"], SubARI_d = cd["SubARI"])
  })
  tab <- do.call(rbind, rows); rownames(tab) <- NULL
  tab$Rank_shallow <- rank_sum(setNames(tab[, c("XB_s", "FS_s", "Dunn_s", "SubARI_s")], c("XB", "FS", "Dunn", "SubARI")))
  tab$Rank_deep    <- rank_sum(setNames(tab[, c("XB_d", "FS_d", "Dunn_d", "SubARI_d")], c("XB", "FS", "Dunn", "SubARI")))
  pick <- function(r, tie) tab$C[order(r, -tie)[1]]
  cat(sprintf("\n[%s] C terpilih  dangkal: %d | deep: %d\n", tag,
              pick(tab$Rank_shallow, tab$SubARI_s), pick(tab$Rank_deep, tab$SubARI_d)))
  cat(sprintf("[%s] FS maksimum  dangkal: C = %d | deep: C = %d\n", tag,
              tab$C[which.max(tab$FS_s)], tab$C[which.max(tab$FS_d)]))
  tab
}

## =============================================================================
## Bagian A — aplikasi 2 (Jawa, 114 unit)
## =============================================================================
t0 <- Sys.time()
tabA <- run_Cselect(as.matrix(D), as.matrix(Pm), "app2")
print(cbind(tabA[, 1:6], round(tabA[, -(1:6)], 4)), row.names = FALSE)
write.csv(tabA, file.path(OUT, "tableS19_app2_Cselect_deep.csv"), row.names = FALSE)
cat("Waktu bagian A:", round(difftime(Sys.time(), t0, units = "hours"), 2), "jam\n")

## =============================================================================
## Bagian B — aplikasi 1 (Jawa Tengah 2022–2025, 35 unit), bila panelnya ada di Input
## =============================================================================
find_app1_panel <- function() {
  cand <- list.files("/kaggle/input", "\\.(csv|xlsx)$", recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
  for (f in cand) {
    r <- tryCatch(if (grepl("xlsx$", f, ignore.case = TRUE)) as.data.frame(readxl::read_excel(f))
                  else read.csv(f, check.names = FALSE), error = function(e) NULL)
    if (is.null(r)) next
    names(r) <- toupper(trimws(names(r)))
    tc <- intersect(c("TAHUN", "YEAR", "TH", "PERIODE"), names(r))[1]; if (is.na(tc)) next
    yrs <- sort(unique(r[[tc]]))
    if (length(yrs) == 4 && all(yrs == 2022:2025) && nrow(r) == 140) return(list(file = f, raw = r, tc = tc))
  }
  NULL
}
p1 <- find_app1_panel()
if (is.null(p1)) {
  cat("\nPanel aplikasi 1 (35 unit x 2022–2025 = 140 baris) TIDAK ditemukan di /kaggle/input.\n",
      "Tambahkan lewat Add Input lalu jalankan ulang cell ini; bagian A tidak akan dihitung ulang.\n")
} else {
  raw1 <- p1$raw; tc1 <- p1$tc
  rc1  <- if ("WADMKK" %in% names(raw1)) "WADMKK" else
          names(raw1)[sapply(names(raw1), function(v) length(unique(raw1[[v]])) == 35)][1]
  ind1 <- c("APM_SD", "APM_SMP", "PRODPADI", "RLS", "SANITASI")   # 5 indikator naskah (AKSESKES tidak dipakai)
  stopifnot(all(ind1 %in% names(raw1)), !anyNA(raw1[, ind1]), length(unique(raw1[[rc1]])) == 35)
  cat("\nPanel aplikasi 1:", p1$file, "\nUnit:", rc1, "| Waktu:", tc1,
      "\nIndikator:", paste(ind1, collapse = ", "), "\n")
  tmp1 <- file.path("/kaggle/working", "panel1_tmp.csv"); write.csv(raw1, tmp1, row.names = FALSE)
  X1 <- P_("zscore_panel")(P_("load_panel")(tmp1, rc1, tc1, ind1))
  D1 <- ck("D_app1", P_("dtw_matrix")(X1, TRUE))
  fW <- list.files("/kaggle/input", "Wq_frozen\\.csv$", recursive = TRUE, full.names = TRUE)[1]
  if (is.na(fW)) stop("Wq_frozen.csv (W aplikasi 1) tidak ada di Input.")
  W1 <- as.matrix(read.csv(fW, row.names = 1, check.names = FALSE))
  reg1 <- dimnames(X1)[[1]]
  nk <- function(v) {                                                 # kunci nama dinormalkan
    v <- toupper(trimws(sub("^X", "", v)))
    v <- sub("^KABUPATEN\\s+|^KAB\\.?\\s+", "", v)
    gsub("[^A-Z]", "", v)
  }
  m_row <- match(nk(reg1), nk(rownames(W1))); m_col <- match(nk(reg1), nk(colnames(W1)))
  if (anyNA(m_row) || anyNA(m_col)) {
    cat("Unit panel yang tidak cocok dengan W:", paste(reg1[is.na(m_row)], collapse = ", "), "\n",
        "Contoh nama di W:", paste(head(rownames(W1), 8), collapse = ", "), "\n")
    stop("Kunci unit panel dan Wq_frozen.csv tidak cocok; samakan penamaannya.")
  }
  W1 <- W1[m_row, m_col]; dimnames(W1) <- list(reg1, reg1)          # urutan W = urutan D1
  stopifnot(nrow(W1) == 35, ncol(W1) == 35)
  P1 <- P_("row_standardise")(W1)
  cat("Rerata DTW off-diagonal aplikasi 1:", round(mean(D1[upper.tri(D1)]), 3), "(naskah: 2.958)\n")
  tabB <- run_Cselect(as.matrix(D1), as.matrix(P1), "app1")
  print(cbind(tabB[, 1:6], round(tabB[, -(1:6)], 4)), row.names = FALSE)
  write.csv(tabB, file.path(OUT, "tableS19b_app1_Cselect_deep.csv"), row.names = FALSE)
}

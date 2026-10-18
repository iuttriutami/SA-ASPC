## =====================================================================
## R_buat_W_jawa.R  (versi final)
## Matriks queen beku untuk 114 kab/kota Pulau Jawa di Panel_App2_Jawa_2020_2024.csv
## Unit panel mengikuti batas GADM 3.6.
##   - pemekaran (Kota Serang, Kota Tangerang Selatan, Bandung Barat) digabung ke induk
##   - waduk/danau/Kepulauan Seribu dibuang
##   - unit terisolasi -> tetangga terdekat
##   - Pulau Madura disambungkan via Suramadu (Bangkalan <-> Kota Surabaya)
## Output: /kaggle/working/W_app2_frozen.csv, /kaggle/working/batas_app2_jawa.gpkg
## =====================================================================
library(sf); library(spdep)
sf::sf_use_s2(FALSE)                     # GEOS planar: irisan batas stabil
WORK  <- "/kaggle/working"
say   <- function(...) cat(sprintf(...), "\n")
inp   <- function(pat) { f <- list.files("/kaggle/input", pat, recursive = TRUE, full.names = TRUE, ignore.case = TRUE); if (length(f)) f[1] else NA }
quiet <- function(expr) suppressMessages(suppressWarnings(expr))

pn <- read.csv(inp("Panel_App2_Jawa.*\\.csv$"))
units <- unique(pn[, c("KODE", "PROVINSI", "WADMKK", "TIPE")])
say("Unit panel: %d", nrow(units))

## ---- normalisasi kunci: provinsi | tipe | nama dasar
norm_prov <- function(x) {
  x <- toupper(gsub("[^A-Za-z ]", "", x))
  x[grepl("JAKARTA", x)] <- "JAKARTA"; x[grepl("YOGYA", x)] <- "YOGYAKARTA"
  gsub("\\s+", "", x)
}
norm_base <- function(x) gsub("[^A-Z]", "", sub("^KOTA\\s+(ADM(INISTRASI)?\\.?\\s+)?", "", toupper(trimws(x))))
units$key <- paste(norm_prov(units$PROVINSI), units$TIPE, norm_base(units$WADMKK), sep = "|")

## ---- batas wilayah
read_bnd <- function() {
  f <- inp("\\.(shp|gpkg|geojson|json)$")
  if (!is.na(f) && !grepl("3300|batas_app2", basename(f))) { say("Batas dari input: %s", basename(f)); return(sf::st_read(f, quiet = TRUE)) }
  urls <- c("https://geodata.ucdavis.edu/gadm/gadm3.6/gpkg/gadm36_IDN_gpkg.zip",
            "https://geodata.ucdavis.edu/gadm/gadm4.1/json/gadm41_IDN_2.json.zip")
  for (u in urls) {
    z <- file.path(WORK, basename(u))
    ok <- tryCatch({ download.file(u, z, mode = "wb", quiet = TRUE); TRUE }, error = function(e) FALSE)
    if (!ok) next
    ff <- unzip(z, exdir = file.path(WORK, "gadm")); say("Unduh: %s", basename(u))
    g <- ff[grepl("gpkg$", ff)]
    if (length(g)) return(sf::st_read(g[1], layer = "gadm36_IDN_2", quiet = TRUE))
    j <- ff[grepl("json$", ff)]; if (length(j)) return(sf::st_read(j[1], quiet = TRUE))
  }
  stop("Tidak ada file batas di input dan unduhan gagal. Nyalakan Internet di Kaggle atau unggah file GADM IDN level 2.")
}
b <- read_bnd()
names(b)[names(b) != attr(b, "sf_column")] <- toupper(names(b)[names(b) != attr(b, "sf_column")])
if (!all(c("NAME_1", "NAME_2") %in% names(b))) stop("File batas harus berkolom NAME_1 dan NAME_2 (format GADM). Kolom: ", paste(names(b), collapse = ", "))
b <- b[norm_prov(b$NAME_1) %in% unique(norm_prov(units$PROVINSI)), ]
tp <- toupper(paste(if ("TYPE_2" %in% names(b)) b$TYPE_2 else "", if ("ENGTYPE_2" %in% names(b)) b$ENGTYPE_2 else ""))
is_kota <- grepl("KOTA|CITY", tp) | grepl("^KOTA ", toupper(b$NAME_2)) | grepl("JAKARTA", toupper(b$NAME_2))
b$key <- paste(norm_prov(b$NAME_1), ifelse(is_kota, "KOTA", "KAB"), norm_base(b$NAME_2), sep = "|")
b <- b[!grepl("KEPULAUANSERIBU|WATERBODY|DANAU|WADUK", b$key), ]

miss <- setdiff(units$key, b$key)
if (length(miss)) {
  say("Belum cocok (%d): %s", length(miss), paste(miss, collapse = ", "))
  say("Kandidat di file batas: %s", paste(head(setdiff(b$key, units$key), 40), collapse = ", "))
  stop("Perbaiki ejaan di salah satu sisi (lihat dua daftar di atas), lalu jalankan ulang.")
}
extra <- setdiff(b$key, units$key)
if (length(extra)) say("[INFO] %d poligon tidak ada di panel, digabung ke induk/tetangga: %s", length(extra), paste(extra, collapse = ", "))

## ---- penggabungan poligon ekstra (pemekaran -> induk; lainnya -> batas bersama terpanjang)
b <- quiet(sf::st_make_valid(b))
if (length(extra)) {
  fixmap <- c("JAWABARAT|KAB|BANDUNGBARAT"   = "JAWABARAT|KAB|BANDUNG",
              "BANTEN|KOTA|SERANG"           = "BANTEN|KAB|SERANG",
              "BANTEN|KOTA|TANGERANGSELATAN" = "BANTEN|KAB|TANGERANG")
  for (k in extra) {
    gi   <- which(b$key == k)
    gg   <- quiet(sf::st_union(sf::st_geometry(b)[gi]))
    if (k %in% names(fixmap)) {
      tgt <- which(b$key == fixmap[[k]])
      if (length(tgt) != 1) stop("Target fixmap tidak ditemukan: ", fixmap[[k]])
    } else {
      cand <- which(b$key %in% units$key)
      cand <- cand[lengths(quiet(sf::st_intersects(sf::st_geometry(b)[cand], gg))) > 0]
      if (!length(cand)) { say("  %s: tidak bersentuhan, dibuang", k); next }
      sh <- vapply(cand, function(j) {
        x <- quiet(sf::st_intersection(sf::st_boundary(sf::st_geometry(b)[j]), sf::st_boundary(gg)))
        if (!length(x)) return(0)
        x <- quiet(sf::st_collection_extract(x, "LINESTRING"))
        if (!length(x)) 0 else as.numeric(sum(sf::st_length(x)))
      }, numeric(1))
      tgt <- cand[which.max(sh)]
    }
    sf::st_geometry(b)[tgt] <- quiet(sf::st_make_valid(sf::st_union(c(sf::st_geometry(b)[tgt], gg))))
    say("  %s -> %s", k, b$key[tgt])
  }
  b <- b[b$key %in% units$key, ]
}

## ---- matriks queen
b  <- b[match(units$key, b$key), ]
nb <- quiet(spdep::poly2nb(b, queen = TRUE, snap = 1e-6))
W  <- spdep::nb2mat(nb, style = "B", zero.policy = TRUE)

## (a) unit tanpa tetangga -> tetangga terdekat
iso <- which(rowSums(W) == 0)
if (length(iso)) {
  cen <- quiet(sf::st_coordinates(sf::st_point_on_surface(sf::st_geometry(b))))
  for (i in iso) { d <- sqrt(colSums((t(cen) - cen[i, ])^2)); d[i] <- Inf; j <- which.min(d); W[i, j] <- W[j, i] <- 1 }
  say("[INFO] unit terisolasi dihubungkan ke tetangga terdekat: %s", paste(units$WADMKK[iso], collapse = ", "))
}

## (b) komponen terpisah (Madura) -> Suramadu
link_pair <- function(a, b_) {
  i <- which(units$key == a); j <- which(units$key == b_)
  if (length(i) != 1 || length(j) != 1) {
    say("Kandidat: %s", paste(grep("BANGKALAN|SURABAYA", units$key, value = TRUE), collapse = ", "))
    stop("Kunci tidak ditemukan: ", a, " / ", b_)
  }
  W[i, j] <<- W[j, i] <<- 1
}
n_comp <- function(W) spdep::n.comp.nb(quiet(spdep::mat2listw(W, style = "B", zero.policy = TRUE))$neighbours)
comp <- n_comp(W)
say("Komponen sebelum link manual: %d", comp$nc)
if (comp$nc > 1) {
  for (g in seq_len(comp$nc)) say("  komponen %d: %s", g, paste(units$WADMKK[comp$comp.id == g], collapse = ", "))
  link_pair("JAWATIMUR|KAB|BANGKALAN", "JAWATIMUR|KOTA|SURABAYA")
  comp <- n_comp(W)
  say("Komponen setelah link Suramadu: %d", comp$nc)
  if (comp$nc > 1) stop("Masih ada komponen terpisah; tambahkan link_pair() sesuai daftar komponen di atas.")
}

## ---- simpan
dimnames(W) <- list(units$KODE, units$KODE)
stopifnot(isSymmetric(W), all(rowSums(W) > 0))
write.csv(W, file.path(WORK, "W_app2_frozen.csv"))
quiet(sf::st_write(sf::st_sf(units, geometry = sf::st_geometry(b)), file.path(WORK, "batas_app2_jawa.gpkg"), quiet = TRUE, delete_dsn = TRUE))
say("Unit berderajat 1: %s", paste(units$WADMKK[rowSums(W) == 1], collapse = ", "))
say("Selesai: %d unit, 1 komponen, rerata derajat %.2f, min %d, maks %d -> W_app2_frozen.csv, batas_app2_jawa.gpkg",
    nrow(W), mean(rowSums(W)), min(rowSums(W)), max(rowSums(W)))

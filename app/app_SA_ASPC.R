## =============================================================================
## SA-ASPC Explorer (R Shiny)
## Interactive companion to: "When does a spatial penalty pay? Calibrating the
## contiguity coefficient in fuzzy clustering of multivariate time trajectories"
## (Utami & Daniaty).
##
## The app recomputes the SA-ASPC selection and its Monte Carlo certificate from
## the calibration-grid table produced by the analysis scripts, so a reader can
## change the weights, the bank size and delta and see how the decision responds.
##
## Inputs (CSV files written by the analysis scripts):
##   * calibration grid (required): columns beta, FM, FS, XB, ST
##       e.g. tableS5_app2_grid.csv; T06_CalibrationGrid.csv for the first panel
##   * ARI matrix between grid points (optional): tableS11_app2_ARI_grid.csv
##   * fuzziness / coherence by beta (optional): tableS13_app2_fuzziness_join.csv
##   * partitions (optional): app2_partitions.csv (region, cluster_beta0, cluster_selected)
##   * boundaries (optional, needs package sf): batas_app2_jawa.gpkg
##
## Run:  install.packages(c("shiny", "ggplot2"))   # plus "sf" for the map tab
##       shiny::runApp("app_SA_ASPC.R")   or open in RStudio and click "Run App"
## =============================================================================

library(shiny)
library(ggplot2)

W_PAPER <- c(FM = 0.35, FS = 0.25, XB = 0.20, ST = 0.20)

## ---------- core computations (same formulas as Sections 3.3-3.4) -----------
minmax <- function(v) { r <- diff(range(v)); if (r <= 0) rep(1, length(v)) else (v - min(v)) / r }

rescale_anchored <- function(tb) cbind(
  FM = minmax(tb$FM),
  FS = minmax(tb$FS),
  XB = 1 - minmax(tb$XB),
  ST = (tb$ST + 1) / 2)                     # anchored stability map, R_A = 2

d_ideal <- function(Z, w) sqrt(as.vector(((1 - Z)^2) %*% (w / sum(w))))

select_beta <- function(tb, w) {
  Z <- rescale_anchored(tb); d <- d_ideal(Z, w); o <- order(d, tb$beta)
  list(d = d, k = o[1], gamma = d[o[2]] - d[o[1]], Z = Z)
}

certify <- function(beta, d, gamma, B, delta, w4) {
  G <- length(beta); L <- log(2 * G / delta)
  eps   <- 2 * sqrt(L / (2 * B))             # Hoeffding half-width for R_A = 2
  kappa <- sqrt(w4) * eps / 2                # anchored kappa_delta
  list(fourk = 4 * kappa, ratio = 4 * kappa / gamma, certified = 4 * kappa < gamma,
       B_required = ceiling(8 * w4 * L / gamma^2),
       conf_set = beta[d <= min(d) + 2 * kappa])
}

rules_table <- function(tb, w) {
  s  <- select_beta(tb, w); k <- s$k; Z <- s$Z
  eq <- d_ideal(Z, rep(0.25, 4)); noST <- d_ideal(Z, c(w[1:3], 0))
  dom_ok <- !(tb$XB[k] > tb$XB[1] && tb$ST[k] < tb$ST[1])
  data.frame(rule = c("SA-ASPC (current weights)", "SA-ASPC equal weights", "SA-ASPC without ST",
                      "Moran maximisation", "XB minimisation", "Guard: dominance", "Guard: XB tolerance 2%"),
             beta = round(c(tb$beta[k], tb$beta[which.min(eq)], tb$beta[which.min(noST)],
                            tb$beta[which.max(tb$FM)], tb$beta[which.min(tb$XB)],
                            if (dom_ok) tb$beta[k] else 0,
                            if (tb$XB[k] <= 1.02 * tb$XB[1]) tb$beta[k] else 0), 3))
}

rdirichlet <- function(n, a) { x <- matrix(rgamma(n * length(a), a), n, byrow = TRUE); x / rowSums(x) }

weight_regions <- function(tb, alpha, n = 20000, seed = 2026) {
  set.seed(seed)
  s2 <- (1 - rescale_anchored(tb))^2
  Wd <- rdirichlet(n, alpha)
  pick <- max.col(-(Wd %*% t(s2)), ties.method = "first")
  data.frame(beta = tb$beta, share = tabulate(pick, nrow(tb)) / n)
}

read_csv_safe <- function(f, ...) if (is.null(f)) NULL else read.csv(f$datapath, check.names = FALSE, ...)

## ---------------------------------- UI ---------------------------------------
ui <- fluidPage(
  titlePanel("SA-ASPC Explorer: spatial-penalty calibration and certification"),
  sidebarLayout(
    sidebarPanel(width = 3,
      h4("1. Data"),
      fileInput("f_grid", "Calibration grid (beta, FM, FS, XB, ST)", accept = ".csv"),
      fileInput("f_ari",  "ARI matrix between grid points (optional)", accept = ".csv"),
      fileInput("f_fuz",  "Fuzziness / coherence by beta (optional)", accept = ".csv"),
      fileInput("f_part", "Partitions (optional)", accept = ".csv"),
      fileInput("f_gpkg", "Boundaries .gpkg (optional, needs sf)", accept = ".gpkg"),
      hr(),
      h4("2. Weights (normalised to sum 1)"),
      sliderInput("wFM", "FM (fuzzy Moran)",       0, 1, W_PAPER[["FM"]], 0.01),
      sliderInput("wFS", "FS (fuzzy silhouette)",  0, 1, W_PAPER[["FS"]], 0.01),
      sliderInput("wXB", "XB (Xie-Beni)",          0, 1, W_PAPER[["XB"]], 0.01),
      sliderInput("wST", "ST (stability)",         0, 1, W_PAPER[["ST"]], 0.01),
      actionButton("reset_w", "Reset to paper weights"),
      hr(),
      h4("3. Certification"),
      numericInput("B", "Bank size B", 100, min = 10, step = 10),
      numericInput("delta", "delta", 0.05, min = 0.001, max = 0.5, step = 0.01),
      hr(),
      helpText("Grid points are re-scaled exactly as in Section 3.3 (min-max for FM, FS, XB;",
               "anchored (ST+1)/2 for stability). The certificate is the exact Hoeffding",
               "certificate of Section 3.4 for the anchored criterion.")
    ),
    mainPanel(width = 9,
      tabsetPanel(
        tabPanel("Calibration",
                 br(), uiOutput("summary_box"),
                 plotOutput("p_obj", height = 330), plotOutput("p_d", height = 260),
                 h4("Selection rules"), tableOutput("t_rules")),
        tabPanel("Certification",
                 br(), tableOutput("t_cert"),
                 plotOutput("p_ratio", height = 330),
                 helpText("The selection is certified when the ratio 4*kappa/gamma falls below 1.",
                          "The certificate concerns Monte Carlo error only, conditional on the data,",
                          "the grid, C, m and the weights.")),
        tabPanel("Weight sensitivity",
                 br(), sliderInput("conc", "Concentration of the local Dirichlet around the current weights",
                                   5, 200, 50, 5),
                 plotOutput("p_regions", height = 380),
                 helpText("Flat: weights drawn uniformly from the simplex. Local: Dirichlet centred at the",
                          "current weights. Bars give the share of draws that select each grid point (Theorem 14).")),
        tabPanel("Partition stability",
                 br(), plotOutput("p_ari", height = 520),
                 helpText("Pairwise ARI between the hard partitions at the grid points. Blocks of ARI = 1",
                          "are plateaus: grid points that give the same partition.")),
        tabPanel("Fuzziness and coherence",
                 br(), plotOutput("p_fuz", height = 420),
                 helpText("PC_norm: normalised partition coefficient (0 = uniform memberships).",
                          "join: share of neighbouring pairs in the same hard cluster.")),
        tabPanel("Map",
                 br(),
                 fluidRow(column(4, selectInput("map_col", "Partition column", choices = NULL)),
                          column(4, selectInput("gpkg_key", "Region key in boundaries", choices = NULL)),
                          column(4, selectInput("part_key", "Region key in partitions", choices = NULL))),
                 plotOutput("p_map", height = 560)),
        tabPanel("Data", br(), tableOutput("t_grid"))
      )
    )
  )
)

## -------------------------------- server -------------------------------------
server <- function(input, output, session) {

  observeEvent(input$reset_w, {
    for (nm in names(W_PAPER)) updateSliderInput(session, paste0("w", nm), value = W_PAPER[[nm]])
  })

  w <- reactive({
    v <- c(FM = input$wFM, FS = input$wFS, XB = input$wXB, ST = input$wST)
    validate(need(sum(v) > 0, "At least one weight must be positive."))
    v / sum(v)
  })

  grid <- reactive({
    validate(need(!is.null(input$f_grid), "Upload the calibration-grid CSV to start."))
    tb <- read_csv_safe(input$f_grid)
    names(tb)[tolower(names(tb)) == "beta"] <- "beta"
    miss <- setdiff(c("beta", "FM", "FS", "XB", "ST"), names(tb))
    validate(need(!length(miss), paste("Missing columns:", paste(miss, collapse = ", "))))
    tb[order(tb$beta), ]
  })

  sel  <- reactive(select_beta(grid(), w()))
  cert <- reactive({ s <- sel(); certify(grid()$beta, s$d, s$gamma, input$B, input$delta, w()[["ST"]]) })

  ## --- Calibration tab
  output$summary_box <- renderUI({
    s <- sel(); ce <- cert(); b <- grid()$beta
    wellPanel(HTML(sprintf(
      "<b>Selected penalty:</b> %.3f &nbsp;|&nbsp; <b>margin</b> &gamma;&#770; = %.5f &nbsp;|&nbsp; <b>certified at B = %d:</b> %s &nbsp;|&nbsp; <b>B required:</b> %s",
      b[s$k], s$gamma, input$B, if (ce$certified) "yes" else "no", format(ce$B_required, big.mark = ","))))
  })

  output$p_obj <- renderPlot({
    tb <- grid(); k <- sel()$k
    long <- rbind(data.frame(beta = tb$beta, value = tb$FM, obj = "FM (higher better)"),
                  data.frame(beta = tb$beta, value = tb$FS, obj = "FS (higher better)"),
                  data.frame(beta = tb$beta, value = tb$XB, obj = "XB (lower better)"),
                  data.frame(beta = tb$beta, value = tb$ST, obj = "ST (higher better)"))
    ggplot(long, aes(beta, value)) + geom_line() + geom_point(size = 1.4) +
      geom_vline(xintercept = tb$beta[k], linetype = 2, colour = "firebrick") +
      facet_wrap(~obj, scales = "free_y", nrow = 1) +
      labs(x = expression(beta), y = NULL, title = "Calibration objectives (dashed: selected penalty)") +
      theme_bw(base_size = 12)
  })

  output$p_d <- renderPlot({
    tb <- grid(); s <- sel(); ce <- cert()
    df <- data.frame(beta = tb$beta, d = s$d, in_set = tb$beta %in% ce$conf_set)
    ggplot(df, aes(beta, d)) + geom_line() +
      geom_point(aes(shape = in_set), size = 2.4) +
      scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 16), name = "In exact confidence set") +
      geom_vline(xintercept = tb$beta[s$k], linetype = 2, colour = "firebrick") +
      labs(x = expression(beta), y = "distance to ideal point",
           title = "Criterion (8) over the grid") + theme_bw(base_size = 12) + theme(legend.position = "bottom")
  })

  output$t_rules <- renderTable(rules_table(grid(), w()), digits = 3)

  ## --- Certification tab
  output$t_cert <- renderTable({
    ce <- cert(); s <- sel()
    data.frame(quantity = c("gamma-hat (margin)", "4 kappa_delta (exact)", "ratio 4 kappa / gamma",
                            "certified", "B required", "exact confidence set (grid points)"),
               value = c(signif(s$gamma, 4), signif(ce$fourk, 4), signif(ce$ratio, 4),
                         if (ce$certified) "yes" else "no", format(ce$B_required, big.mark = ","),
                         paste(round(ce$conf_set, 3), collapse = ", ")))
  })

  output$p_ratio <- renderPlot({
    s <- sel(); G <- nrow(grid())
    Bs <- unique(round(exp(seq(log(10), log(max(2000, 2 * min(cert()$B_required, 1e6))), length.out = 200))))
    r <- sapply(Bs, function(B) certify(grid()$beta, s$d, s$gamma, B, input$delta, w()[["ST"]])$ratio)
    ggplot(data.frame(B = Bs, ratio = r), aes(B, ratio)) + geom_line() +
      geom_hline(yintercept = 1, linetype = 2) + geom_vline(xintercept = input$B, colour = "firebrick") +
      scale_x_log10() + scale_y_log10() +
      labs(x = "bank size B (log scale)", y = "4 kappa / gamma (log scale)",
           title = sprintf("Exact certificate versus bank size (G = %d, delta = %.3f)", G, input$delta)) +
      theme_bw(base_size = 12)
  })

  ## --- Weight sensitivity tab
  output$p_regions <- renderPlot({
    tb <- grid()
    flat  <- weight_regions(tb, rep(1, 4)); flat$prior <- "flat (uniform on simplex)"
    local <- weight_regions(tb, input$conc * w()); local$prior <- sprintf("local (concentration %d)", input$conc)
    df <- rbind(flat, local)
    ggplot(df, aes(factor(round(beta, 2)), share)) + geom_col() + facet_wrap(~prior, ncol = 1) +
      labs(x = expression(beta), y = "share of weight draws selecting this grid point") +
      theme_bw(base_size = 12) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
  })

  ## --- Partition stability tab
  output$p_ari <- renderPlot({
    validate(need(!is.null(input$f_ari), "Upload the ARI matrix CSV (e.g. tableS11_app2_ARI_grid.csv)."))
    A <- as.matrix(read.csv(input$f_ari$datapath, row.names = 1, check.names = FALSE))
    lab <- rownames(A); n <- nrow(A)
    df <- data.frame(i = rep(seq_len(n), n), j = rep(seq_len(n), each = n), ari = as.vector(A))
    k <- which.min(abs(as.numeric(lab) - grid()$beta[sel()$k]))
    ggplot(df, aes(i, j, fill = ari)) + geom_tile() +
      scale_fill_gradientn(colours = hcl.colors(60, "YlOrRd", rev = TRUE), limits = c(min(0.6, min(df$ari)), 1), name = "ARI") +
      annotate("point", x = k, y = k, shape = 4, size = 5, stroke = 1.5) +
      scale_x_continuous(breaks = seq_len(n), labels = lab, expand = c(0, 0)) +
      scale_y_continuous(breaks = seq_len(n), labels = lab, expand = c(0, 0)) +
      labs(x = expression(beta), y = expression(beta), title = "Pairwise ARI between grid partitions (x: selected)") +
      coord_equal() + theme_minimal(base_size = 11) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
  })

  ## --- Fuzziness tab
  output$p_fuz <- renderPlot({
    validate(need(!is.null(input$f_fuz), "Upload the fuzziness CSV (e.g. tableS13_app2_fuzziness_join.csv)."))
    fz <- read.csv(input$f_fuz$datapath, check.names = FALSE)
    cols <- intersect(c("PC_norm", "gap12", "join", "join_exp"), names(fz))
    validate(need(length(cols) > 0 && "beta" %in% names(fz), "Expected columns beta and PC_norm / gap12 / join."))
    long <- do.call(rbind, lapply(cols, function(cl) data.frame(beta = fz$beta, value = fz[[cl]], measure = cl)))
    ggplot(long, aes(beta, value)) + geom_line() + geom_point(size = 1.3) +
      geom_vline(xintercept = grid()$beta[sel()$k], linetype = 2, colour = "firebrick") +
      facet_wrap(~measure, scales = "free_y") + labs(x = expression(beta), y = NULL) + theme_bw(base_size = 12)
  })

  ## --- Map tab
  part <- reactive(read_csv_safe(input$f_part))
  shp  <- reactive({
    if (is.null(input$f_gpkg)) return(NULL)
    validate(need(requireNamespace("sf", quietly = TRUE), "Install the sf package for the map."))
    sf::st_read(input$f_gpkg$datapath, quiet = TRUE)
  })
  observe({
    p <- part(); if (!is.null(p)) {
      updateSelectInput(session, "map_col", choices = setdiff(names(p), c("region", "name")),
                        selected = intersect(c("cluster_selected", names(p)), names(p))[1])
      updateSelectInput(session, "part_key", choices = names(p), selected = names(p)[1])
    }
    s <- shp(); if (!is.null(s)) {
      nm <- setdiff(names(s), attr(s, "sf_column"))
      updateSelectInput(session, "gpkg_key", choices = nm,
                        selected = intersect(c("KODE", "key", "KEY", nm), nm)[1])
    }
  })
  output$p_map <- renderPlot({
    p <- part(); s <- shp()
    validate(need(!is.null(p) && !is.null(s), "Upload both the partitions CSV and the boundaries .gpkg."),
             need(nzchar(input$map_col) && nzchar(input$gpkg_key) && nzchar(input$part_key), "Choose the columns."))
    key_s <- as.character(s[[input$gpkg_key]]); key_p <- as.character(p[[input$part_key]])
    s$cluster <- factor(p[[input$map_col]][match(key_s, key_p)])
    validate(need(any(!is.na(s$cluster)), "No region keys matched; check the two key columns."))
    ggplot(s) + geom_sf(aes(fill = cluster), colour = "grey30", linewidth = 0.1) +
      scale_fill_brewer(palette = "Set2", na.value = "grey90") +
      labs(title = sprintf("Hard partition: %s", input$map_col), fill = "Cluster") + theme_void(base_size = 12)
  })

  ## --- Data tab
  output$t_grid <- renderTable({
    tb <- grid(); s <- sel()
    cbind(tb[, intersect(c("beta", "FM", "FS", "XB", "ST", "ST_sd", "SC"), names(tb))],
          d = s$d, selected = seq_len(nrow(tb)) == s$k)
  }, digits = 4)
}

shinyApp(ui, server)

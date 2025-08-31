cat("\n======================================================================\n")
cat("FUNÇÃO: fatgraf.DBC\n")
cat("Descrição : Análises fatoriais para DBC automatizadas com gráficos e tabelas integradas.\n")
cat("Desenvolvimento: Walter Baida (Doutorando em Biometria e Estatística)\n")
cat("                 Wilhan Valasco (Doutorando em Biotecnologia)\n")
cat("Contato       : +55 (64) 99250-2925 | +55 (64) 99286-2039\n")
cat("======================================================================\n")

fatgraf.DBC <- function(data, variavel,
                    f1 = "Fator1", f2 = "Fator2", bloco = "Rep",
                    fat.nome = c("Fator 1", "Fator 2"),
                    Rank = FALSE, alpha = 0.05,
                    mcomp = c("sk", "tukey", "lsd", "snk"),
                    alfa.t = 0.05, alfa.f = 0.05,
                    y_label = NULL, x_label = NULL, titulo = NULL,
                    ylim_min = NULL, ylim_max = NULL,
                    cores = NULL,
                    ang = 0,
                    ang_rot = 0,
                    console = TRUE) {
  
  pb <- txtProgressBar(min = 0, max = 5, style = 3)
  
  pacotes <- c("agricolae", "AgroR", "ExpDes", "dplyr", "readr", "ggplot2")
  for (p in pacotes) {
    if (!require(p, character.only = TRUE, quietly = TRUE)) {
      install.packages(p, dependencies = TRUE)
      library(p, character.only = TRUE)
    }
  }
  setTxtProgressBar(pb, 1)
  
  erro_padrao <- function(x) {
    sd(x, na.rm = TRUE) / sqrt(length(na.omit(x)))
  }
  
  required_cols <- c(f1, f2, bloco, variavel)
  missing_cols <- required_cols[!required_cols %in% names(data)]
  if (length(missing_cols) > 0) {
    stop("As seguintes colunas estão ausentes no data frame: ", paste(missing_cols, collapse = ", "))
  }
  
  fator1 <- factor(data[[f1]])
  fator2 <- factor(data[[f2]])
  controle <- factor(data[[bloco]])
  y <- data[[variavel]]
  
  if (!(length(fator1) == length(fator2) && length(fator2) == length(controle) && length(controle) == length(y))) {
    stop("As variáveis fornecidas têm comprimentos diferentes. Verifique os dados.")
  }
  
  dat <- data.frame(fator1 = fator1, fator2 = fator2, controle = controle, y = y)
  
  if (Rank) {
    modelo <- with(dat, FAT2DBC(fator1, fator2, controle, rank(y),
                                norm = "sw", homog = "bt",
                                alpha.f = alfa.f, alpha.t = alfa.t,
                                quali = c(TRUE, TRUE),
                                names.fat = fat.nome, mcomp = mcomp,
                                plot.on = FALSE, print.on = F))
  } else {
    modelo <- with(dat, FAT2DBC(fator1, fator2, controle, y,
                                norm = "sw", homog = "bt",
                                alpha.f = alfa.f, alpha.t = alfa.t,
                                quali = c(TRUE, TRUE),
                                names.fat = fat.nome, mcomp = mcomp,
                                plot.on = FALSE, print.on = F))
  }
  
  rownames(modelo$anova) <- gsub("Fator1", fat.nome[1],
                                 gsub("Fator2", fat.nome[2],
                                      gsub("Fator1:Fator2", paste(fat.nome, collapse = ":"), rownames(modelo$anova))))
  
  setTxtProgressBar(pb, 2)
  
  tabela_medias <- modelo$plot$interaction$data
  mediana <- aggregate(y ~ fator1:fator2, data = dat, FUN = if (Rank) median else mean)
  colnames(mediana) <- c("f1", "f2", "mediana")
  erro <- aggregate(y ~ fator1:fator2, data = dat, FUN = erro_padrao)
  colnames(erro) <- c("f1", "f2", "erro")
  resultado_final <- merge(tabela_medias, mediana, by = c("f1", "f2"))
  resultado_final <- merge(resultado_final, erro, by = c("f1", "f2"))
  resultado_final$Letras <- paste0(resultado_final$letra, resultado_final$letra1)
  
  dados_plot <- resultado_final %>%
    dplyr::select(f2, f1, media_fito = mediana, erro, Letras) %>%
    mutate(
      label = paste0(round(media_fito, 2), " ", Letras),
      ymax = media_fito + erro,
      ymin = media_fito - erro
    )
  
  setTxtProgressBar(pb, 3)
  
  y_min <- ifelse(is.null(ylim_min), 0, ylim_min)
  y_max <- ifelse(is.null(ylim_max), max(dados_plot$ymax, na.rm = TRUE) * 1.5, ylim_max)
  
  if (is.null(cores)) {
    cores <- RColorBrewer::brewer.pal(n = length(unique(dados_plot$f2)), name = "Set2")
  }
  
  p <- ggplot(dados_plot, aes(x = f1, y = media_fito, fill = factor(f2))) +
    geom_col(position = position_dodge(width = 0.9), width = 0.8, color = NA) +
    geom_errorbar(aes(ymin = ymin, ymax = ymax),
                  position = position_dodge(width = 0.9), width = 0.2, color = "black") +
    geom_text(
      aes(label = label, y = ymax + 0.03 * max(ymax, na.rm = TRUE)),
      position = position_dodge(width = 0.9),
      angle = ang_rot,
      hjust = ifelse(ang_rot == 90, 0, 0.5),
      vjust = ifelse(ang_rot == 90, 0.5, 0),
      size = 3.5,
      fontface = "plain",
      family = "serif"
    ) +
    scale_y_continuous(limits = c(y_min, y_max), expand = c(0, 0)) +
    scale_fill_manual(values = cores) +
    labs(
      title = titulo %||% paste("Análise de", variavel),
      x = x_label %||% expression(italic("Tratamentos")),
      y = y_label %||% variavel,
      fill = fat.nome[2]
    ) +
    theme_classic(base_size = 8, base_family = "serif") +
    theme(
      axis.text.x = element_text(
        angle = ang,
        hjust = ifelse(ang == 90, 1, 0.5),
        vjust = ifelse(ang == 90, 0.5, 1),
        color = "black", size = 12
      ),
      axis.line = element_line(color = "black", size = 0.1),
      axis.text.y = element_text(color = "black"),
      axis.title.x = element_text(face = "plain", color = "black", size = 14),
      axis.title.y = element_text(face = "plain", color = "black", size = 14),
      panel.border = element_rect(color = "black", fill = NA, size = 0.4),
      panel.grid.major = element_line(color = "gray90", size = 0.1),
      panel.grid.minor = element_line(color = "gray90", size = 0.1),
      panel.background = element_rect(fill = "transparent"),
      plot.background = element_rect(fill = "transparent", color = NA),
      legend.position = c(0.05, 0.95),
      legend.justification = c("left", "top"),
      legend.direction = "horizontal",
      legend.background = element_rect(fill = "transparent", color = NA),
      legend.box.background = element_blank(),
      legend.text = element_text(size = rel(0.6))
    )
  
  
  setTxtProgressBar(pb, 4)
  
  qm_res <- modelo$anova[["QM"]][rownames(modelo$anova) == "Residuals"]
  media_geral <- mean(modelo$resp, na.rm = TRUE)
  cv <- 100 * sqrt(qm_res) / media_geral
  
  if (console) {
    print(p)
    cat(paste(rep("-", 80), collapse = ""), "\n")    
    cat("\nComparações múltiplas\n")
    cat(paste(rep("-", 80), collapse = ""), "\n")
    print(
      dados_plot %>%
        dplyr::select(f1, f2, media_fito, Letras, erro) %>%
        dplyr::arrange(f1, f2)
    )
    cat(paste(rep("-", 80), collapse = ""), "\n")   
    cat("\nAnálise de variância\n")
    cat(paste(rep("-", 80), collapse = ""), "\n")
    print(modelo$anova)
    cat(sprintf("\nCoeficiente de variação (CV): %.2f%%\n", cv))
    
    cat(paste(rep("-", 80), collapse = ""), "\n")
    cat("\nTeste de normalidade dos resíduos\n")
    cat(paste(rep("-", 80), collapse = ""), "\n")
    print(modelo$norm)
    
    cat(paste(rep("-", 80), collapse = ""), "\n")
    cat("\nHomogeneidade de variâncias (teste de homocedasticidade)\n")
    cat(paste(rep("-", 80), collapse = ""), "\n")
    print(modelo$homog)
  }
  
  setTxtProgressBar(pb, 5)
  close(pb)
  
  return(list(
    grafico = p,
    tabela = dados_plot %>%
      dplyr::select(f1, f2, media_fito, Letras, erro) %>%
      dplyr::arrange(f1, f2),
    anova = modelo$anova,
    CV = cv,
    normalidade = modelo$norm,
    homogeneidade = modelo$homog,
    modelo = modelo
  ))
}



cat("\nExemplo de uso com banco de dados fictício:\n")
cat(paste(rep("-", 80), collapse = ""), "\n")
cat(
  '# Gerar banco de dados fictício para o exemplo\n',
  'set.seed(123)\n',
  'Cultivar <- rep(c("C1", "C2", "C3", "C4"), each = 4 * 3)\n',
  'Mistura  <- rep(rep(c("M1", "M2", "M3", "M4"), each = 3), 4)\n',
  'Rep      <- rep(1:3, times = 16)\n',
  'Prod     <- rnorm(48, mean = 60, sd = 10)\n\n',
  'df <- data.frame(Cultivar, Mistura, Rep, Prod)\n\n',
  '# Exemplo de uso da função\n',
  'res <- fatgraf.DBC(\n',
  '  data = df,\n',
  '  variavel = "Prod",\n',
  '  f1 = "Cultivar",\n',
  '  f2 = "Mistura",\n',
  '  bloco = "Rep",\n',
  ' fat.nome = c("Cultivar", "Mistura"),\n',
  '  Rank = TRUE,\n',
  '  mcomp = "sk",\n',
  '  alfa.t = 0.1,\n',
  '  alfa.f = 0.99,\n',
  '  y_label = "Prod.",\n',
  '  x_label = "Cultivares",\n',
  '  titulo = "Cultivar x Mistura",\n',
  '  ylim_min = 0,\n',
  '  ylim_max = 100,\n',
  '  ang = 0,\n',
  '  ang_rot = 90,\n',
  '  cores = c("#1b9e77", "blue", "orange", "gray"),\n',
  '  console = TRUE\n',
  ')\n',
  sep = ""
)

#' @importFrom ggplot2 ggplot aes geom_point geom_errorbar geom_hline coord_flip
#' @importFrom ggplot2 scale_y_continuous labs theme theme_minimal element_text element_blank
NULL

utils::globalVariables(c("Object", "Ratio", "Ratio.lower", "Ratio.upper"))

#' Create unadjusted and adjusted OR/HR forest plot
#'
#' This function creates forest plots of OR/HR. It can make unadjusted if the entrance object is from compareGroups family, or adjusted if it is an add_adjusted_ratio object.
#'
#'
#' @param compare_object An object of class compareGroups family or add_adjusted_ratio with de OR/HR calculated for plotting.
#'
#' @param hline_forest_plot numeric value indicating the position of the horizontal reference line in the forest plot. Default is 1, representing the null effect. If you do not want to display it, you can set the value to NULL.
#'
#' @param log_scale Logical. If TRUE, the OR/HR will be displayed on a logarithmic scale in the forest plot.
#'
#'
#' @return A forest_plot with the OR/HR according to the input object.
#'
#'
#' @examples
#' \dontrun{
#'
#' library("compareGroups")
#' library("add2compareGroups")
#'
#' data(lung_patients)
#'
#' # Create inicial compareGroups objects
#' cg_death <- compareGroups(
#'   DEATH ~ GENDER + AGE + SMOKING + PACK_YEARS + RECURRENCE,
#'   data = lung_patients
#' ) # Create a basic compareGroups
#'
#' # Transform it to a createTable for better late interpretation
#' ct_death <- createTable(cg_death, show.ratio = TRUE)
#'
#' # Common adjusted variables
#' OR_common <- add_adjusted_ratio(ct_death, adjust_common = list("AGE", "GENDER"))
#' export2md(OR_common)
#'
#' # Individual adjusted variables
#' OR_individual <- add_adjusted_ratio(ct_death,
#'   adjust_common = list("AGE", "GENDER"),
#'   adjust = list("SMOKING" = list("PACK_YEARS", "PROGRESSION"), "RECURRENCE" = list("SURGERY")),
#'   show.n = TRUE
#' )
#' export2md(OR_individual)
#'
#' # Create the forest_plots
#' ratio_forest_plot(ct_death)
#' ratio_forest_plot(ct_death, log_scale = TRUE)
#' ratio_forest_plot(OR_individual, hline_forest_plot = NULL)
#' ratio_forest_plot(OR_individual, log_scale = TRUE)
#'
#' }
#' @export

# Inici de la funció ----

ratio_forest_plot <- function(compare_object, hline_forest_plot = 1, log_scale = FALSE){

  if (inherits(compare_object, "add_adjusted_ratio")){

    # és un objecte de la nostre class
    # extreiem la matriu_forest de l'objecte entrat

    is_adjust <- "Adjusted"

    data_util <- attr(compare_object, "data_util")
    var_dep <- attr(compare_object, "var_dep")
    var_indep_label <- attr(compare_object, "var_indep_label")

    matriu_forest <- attr(compare_object, "matriu_forest")

  } else  if(inherits(compare_object, c("createTable", "descrTable", "compareGroups"))){

    # En cas de ser d'quests objectes
    # és createTable o descrTable
    # treiem l'objecte comareGroups

    if(inherits(compare_object, c("createTable", "descrTable"))){
      compare_object <- attr(compare_object, "x")[[1]]
    }

    # Si es compareGropus hem de generar la matriu_forest
    # Fem un indicador que els ratios son crus per després posar-ho al gràfic
    # extreiem la matriu OR/HR de la mateixa manera que en la funció principal

    is_adjust <- "Raw"

    data_util <- attr(compare_object, "Xext")
    var_dep <- attr(compare_object, "yname.orig")
    var_indep_label <- attr(compare_object, "names")

    matriu_forest <- NULL

    for(i in (seq_along(var_indep_label))){

      if (is.factor(data_util[[var_dep]])){
        matriu_ratio <- attr(compare_object [[i]], "OR")
      } else if(is.Surv(data_util[[var_dep]])){
        matriu_ratio <- attr(compare_object [[i]], "HR")
      } else {stop("Dependent variable must be a binary factor or a survival object.")}

      fila_ref <- which(is.na(attr(compare_object[[i]], "p.ratio")))

      if (length(fila_ref) != 0){
        rnames <- rownames(matriu_ratio)[-fila_ref]
        matriu_ratio <- matrix(matriu_ratio[-fila_ref, ], ncol = 3)
        rownames(matriu_ratio) <- rnames
      } else {
        rownames(matriu_ratio) <- var_indep_label[i]
      }

      matriu_estim <- cbind(rownames(matriu_ratio), matriu_ratio)

      if (!(var_indep_label[i] %in% matriu_estim[, 1])){
        matriu_estim[, 1] <- paste(var_indep_label[i], "-", matriu_estim[, 1])
      }

      matriu_forest <- rbind(matriu_forest, matriu_estim)




    }
  } else {stop("The input object must be of family class 'compareGroups'. ")}


  df_forest <- as.data.frame(matriu_forest)
  df_forest[-1] <- lapply(df_forest[-1], as.numeric)

  colnames(df_forest) <- c("Object", "Ratio", "Ratio.lower", "Ratio.upper")


  ## Creació del forest plot ----

  # Creem el nom del eix hortizontal personalitzar per si és OR/HR i també si es escala log o lineal
  # També calculem els valors minim i maxim

  or_lower <- min(df_forest$Ratio.lower)
  or_upper <- max(df_forest$Ratio.upper)


  if (is.factor(data_util[[var_dep]])) {
    ratio <- "Odds Ratio"
  } else if (is.Surv(data_util[[var_dep]])) {
    ratio <- "Hazard Ratio"
  }

  if (log_scale) {
    trans_forest <- "log"
    h_label <- paste0(is_adjust, " ", ratio, " (log scale)")

    # Limits pel gràfic
    min_modif <- floor(or_lower*10)/10
    max_modif <- ceiling(or_upper*10)/10


    # Per calcular el by
    min_log <- floor(log(or_lower)*10)/10
    max_log <- ceiling(log(or_upper)*10)/10


    by_modif <- diff(range(min_log, max_log)) / 10


    by_modif <- ceiling(by_modif/0.05)*0.05

    seq_initial <- exp(seq(min_log, max_log, by = by_modif))

    dif_u <- seq_initial[which.min(abs(seq_initial-1))] - 1

    breaks <- round(seq_initial - dif_u, 2) # Transformació inversa: exp(x) perquè li fa el logaritme al fer el gràfic

  } else {
    trans_forest <- "identity"
    h_label <- paste0(is_adjust, " ", ratio, " (linear scale)")

    min_modif <- floor(or_lower*10)/10
    max_modif <- ceiling(or_upper*10)/10

    by_modif <- diff(range(min_modif, max_modif)) / 10


    by_modif <- ceiling(by_modif/0.05)*0.05

    seq_initial <- seq(min_modif, max_modif, by = by_modif)

    dif_u <- seq_initial[which.min(abs(seq_initial-1))] - 1


    breaks <- round(seq_initial - dif_u, 2)
  }

  df_forest$Object <- factor(df_forest$Object, levels = rev(unique(df_forest$Object)))

  # Creem el forest plot base
  plot_forest <- ggplot(df_forest, aes(x = Object, y = Ratio, ymin = Ratio.lower, ymax = Ratio.upper)) +
    geom_errorbar(
      width = 0.15
    ) +  # punts amb barres d'IC
    geom_point(
      size = 2.25
    ) +
    geom_hline(
      yintercept = hline_forest_plot,
      linetype = "dashed",
      color = "red",
      size = 0.5
    ) +  # línia de referència
    coord_flip() +  # gira els eixos per tenir variables verticals
    scale_y_continuous(
      limits = c(min_modif, max_modif),
      breaks = breaks,
      labels = if (log_scale) {  # Etiquetes personalitzades per escala log
        function(x) formatC(x, format = "f", digits = 2)
      } else {
        function(x) as.character(formatC(x, format = "f", digits = 2))
      },
      trans = trans_forest
    ) +
    labs(
      x = "Variables",
      y = h_label,
      title = "Forest Plot"
    ) +
    theme_minimal() +
    theme(
      axis.text.y = element_text(),
      axis.title.y = element_blank(),
      axis.text.x  = element_text(angle = 45, hjust = 1)
    )

  return(plot_forest)

}




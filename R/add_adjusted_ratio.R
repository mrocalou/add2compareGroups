#' @importFrom stats glm binomial coef confint drop1 update model.frame setNames as.formula
#' @importFrom survival is.Surv coxph
#' @importFrom compareGroups createTable
NULL

#' Add adjusted OR/HR columns to compareGroups objects
#'
#' This function adds new columns with adjusted Odds Ratios (OR) or Hazard Ratios (HR), along with their respective statistical parameters (such as p-values and p-overall), to an object from compareGroups family. Users can specify different adjustment variables for each descripted variable in the comparison or define a common set of covariables to be applied across all variables.
#'
#'
#' @param compare_object An object of class compareGroups, createTable or descrTable that contains the original comparisons.
#'
#' @param adjust_common A character list specifying a common set of covariables used to adjust all variables included in compare_object.
#'
#' @param adjust A names list where each element corresponds to a variable in compare_object, and its value is a character list containing the covariates used to adjust that specific variable.
#'
#' @param show.ratio Logical. If TRUE, adjusted OR/HR are shown in the table.
#'
#' @param show.p.ratio Logical. If TRUE, p-values corresponding to each adjusted OR/HR is diplayed in the table.
#'
#' @param show.p.overall Logical. If TRUE, the overall p-values for group significance ('p.overall' column) is displayed in the table.
#'
#' @param show.n Logical. If TRUE, the number of observations used to compute each OR/HR estimate is displayed in the table.
#'
#' #' @param cor_test Logical. If TRUE, checks for multicollinearity among adjustment variables using VIF and excludes those with VIF > 5.
#'
#'
#' @return An updated object of class add_adjusted_ratio and createTable, generated from the input compare_object, that now includes additional columns with the adjusted OR/HR values and their associated statistics.
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
#' OR_common <- add_adjusted_ratio(ct_death,
#'   adjust_common = list("AGE", "GENDER")
#' )
#' export2md(OR_common)
#'
#' # Individual adjusted variables
#' OR_individual <- add_adjusted_ratio(ct_death,
#'   adjust_common = list("AGE", "GENDER"),
#'    adjust = list("SMOKING" = list("PACK_YEARS", "PROGRESSION"), "RECURRENCE" = list("SURGERY")),
#'    show.n = TRUE
#' )
#' export2md(OR_individual)
#'
#' }
#' @export



# Inici de la funció ----

add_adjusted_ratio <- function(compare_object, adjust_common = NULL, adjust = NULL, show.ratio = TRUE, show.p.ratio = show.ratio, show.p.overall = TRUE, show.n = FALSE, cor_test = TRUE){

  #////////////////////////////////////////////

  ## Comprovacions inicials i estructura perque compare_object sigui un objecte compareGroups ----

  if(inherits(compare_object, c("cbind.createTable", "rbind.createTable"))){
    stop("This function is not applicable to objects of class 'cbind.createTable' or 'rbind.createTable'.")
  }


  if (!show.ratio && !show.p.ratio && !show.p.overall) {
    stop("At least one of 'show.ratio', 'show.p.ratio', or 'show.p.overall' must be TRUE for this function to run.")
  }

  ### 'base' -> objecte createTable o descrTable en cas que s'entri aquest.
  ### 'compare_object' -> objecte compareGroups

  if(inherits(compare_object, c("createTable", "descrTable"))){
    base <- compare_object
    compare_object <- attr(compare_object, "x")[[1]]

  } else if(inherits(compare_object, "compareGroups")){
    base <- createTable(compare_object)

  } else {stop("The input object must be of family class 'compareGroups'.")}



  #////////////////////////////////////////////

  ## Objectes que necesitarem ----

  # Secció on es creen objectes per retornar outputs o per treballar amb ells.

  ### Per fer el forest_plot ----
  matriu_forest <- NULL

  ### Per retornar els resultats ----
  results <- list()

  ### Pels missatges d'ajustament

  text_adj <- NULL

  ### Data base de l'objecte compareGroups entrat ----
  data_util <- attr(compare_object, "Xext")
  nrow_data <- nrow(data_util)

  ### Noms de les variables de tot el dataframe ----

  all_vars_data_util <- names(attr(compare_object, "Xext"))

  ### Les variables ----
  var_dep <- attr(compare_object, "yname.orig")
  #var_dep_label <- attr(compare_object, "yname")
  var_indep <- attr(compare_object, "varnames.orig")
  var_indep_label <- attr(compare_object, "names")



  # Mantener en una lista compacta.

  #////////////////////////////////////////////

  ## Comprovacions pre protocolos ----

  # A la hora de calcular los OR/HR para el la variable de interés, el cálculo predeterminado es sencillamenteun modelo crudo. Sin embargo, en ocasiones es necesario onbtenre unos OR/HR ajustado según un modelo  multivariante, que incluya la variable de interés  junto a otras convariables.
  # Para poder hacer esto: El parametro "adjust" es un ajuste en la fórmula para hacer los cálculos.
  # Reajustar la fórmula require de verificar varias cosas, de modo que hay que hacer checks con el parámetro "adjust".

  # creem una llista NULL per en cas que sigui necesari ordenar adjust:
  adjust_sort <- NULL

  # Si no hi ha ajustament, no segueix perquè no té sentit

  if (is.null(adjust) && is.null (adjust_common)){
    stop("At least one of 'adjust' or 'adjust_common' must be provided in order to execute the process.")
  }

  if (!is.null(adjust) && !is.null(adjust_common)){
    message("Since adjust and adjust_common exist, adjust will be used")
  }

  ## si "adjust" és NULL i existeix "adjust_common"
  ## adjust passa a ser una llista amb els noms de "var_indep" i amb la llista de "adjust_common" a cada variable

  if (is.null(adjust) && !is.null(adjust_common)) {

    adjust <- lapply(var_indep, function(x){
      adjust_list_comun <- adjust_common

      return(adjust_list_comun)
    })

    names(adjust) <- var_indep
  }


  # En cas que 'adjust' no sigui NULL, ha de passar unes comrpovacions:

  if(!is.null(adjust)){

    #////////////////////////////////////////////

    ### 'adjust' ha de ser una llista
    if (!is.list(adjust)) {
      stop("Parameter 'adjust' must be a list.")
    }

    ### guardem els noms del parametre 'adjust' en una llista

    names_adjust <- names(adjust)

    ### els noms de les llistes han de ser strings

    if (!all(sapply(names_adjust, is.character))) {
      stop("All names in the 'adjust' list must be strings.")
    }

    ### no hi pot haver cap nom de la llista repetit

    if (length(names_adjust) != length(unique(names_adjust))) {
      ## Arnau Lagarda: aqui noms_adjust hauria de ser names_adjust no?99
      duplicats <- names_adjust[duplicated(names_adjust)]
      stop(paste0("There are duplicated independent variables in 'adjust': ", paste(duplicats, collapse = ", ")))
    }

    ### per cada nom de 'adjust'

    for (var_name in names_adjust) {

      #### 'var_name' ha de ser una de les variables independents

      if(!var_name %in% var_indep){
        stop(paste0("Variable '", var_name, "' is not an independent variable."))
      }

      #### extreiem la llista de la variable per veure si compleix les condicions que volem

      var_name_list <- adjust[[var_name]]

      #### 'var_name_list' ha de ser una llista

      if (!is.list(var_name_list)) {
        stop(paste0("Element '", var_name, "' in 'adjust' is not a list."))
      }

      #### tots els elements de la llista han de ser strings

      if (!all(sapply(var_name_list, is.character))) {
        stop(paste0("All variables within '", var_name, "' must be provided as strings."))
      }

      #### la variable independent no es troba dins la seva propia llista

      if (var_name %in% var_name_list) {
        adjust[[var_name]] <- adjust[[var_name]][adjust[[var_name]] != var_name]
        message(paste0("Variable '", var_name, "' has been removed from its own adjustment list to avoid redundancy."))
      }

      #### la llista no pot tenir velors duplicats

      if (length(var_name_list) != length(unique(var_name_list))) {
        stop(paste0("There are duplicate variables in the adjustment list for '", var_name, "'."))
      }

      #### totes les variables de la llista es troben en en dataframe (objecte 'all_vars_data_util')

      if (!all(var_name_list %in% all_vars_data_util)) {
        faltants <- var_name_list[!var_name_list %in% all_vars_data_util]
        if(length(faltants) == 1){
          stop(paste0("Variable '", paste(faltants, collapse = ", "),
                      "' is not present in the original dataset."))
        } else{
          stop(paste0("Variables '", paste(faltants, collapse = ", "),
                      "' are not present in the original dataset."))
        }

      }

    }

  }

  ### En cas que totes aquestes condicions es compleixin, ordenem l'objecte adjust perquè tingui l'ordre de les variables independents:

  if(!identical(names(adjust), var_indep)){
    adjust_sort <- setNames(vector("list", length(var_indep)), var_indep)
    for (i in seq_along(var_indep)) {
      var <- var_indep[i]
      if (var %in% names(adjust)) {
        adjust_sort[[var]] <- adjust[[var]]
      } else {
        adjust_sort[[var]] <- NULL
        message(paste0("Variable '", var, "' does not have a specified adjustment. Set to NULL."))
      }
    }
  } else {adjust_sort <- adjust}  ### Arnau Lagarda, el parentesis del else no és correcte, haurien de ser claudators {} 99

  if (cor_test){
    # Transformem totes les vairabels a numeric perquè la correlació serà igual i així podem calcular el VIF
    data_VIF <- data.frame(suppressWarnings(lapply(data_util, as.numeric)))

    for (var in names(adjust_sort)){
      for (v_adj in unlist(adjust_sort[var])){

        adj_formula <- as.formula(paste0(var, " ~ ", v_adj))
        mod_lm <- lm(adj_formula, data = data_VIF)
        smod <- summary(mod_lm)
        r_mod <- (smod$r.squared)^2
        VIF <- 1/(1-r_mod)

        if (!is.na(VIF) && (VIF > 5 || VIF == "Inf")) {
          adjust_sort[[var]] <- adjust_sort[[var]][unlist(adjust_sort[[var]]) != v_adj]
          message(paste0("Removed '", v_adj, "' from '", var, "' adjustments due to high collinearity."))
        }
      }
    }
  }


  #////////////////////////////////////////////

  ## Càlcul de OR o HR ----

  # Seccion dedidca al protocolo del cálculo del OR/HR para la tabla del compareGroups

  ### Para variables categoricas ----
  #### Comprovació de que el factor depenent té només 2 ----
  # De moment, només es poden fer OR en cas que la variable dependent tingui només dos nivells.

  if (is.factor( data_util[[var_dep]]) && nlevels( data_util[[var_dep]] ) != 2) {
    stop("Dependent factor must have exactly two levels.")
  }

  #### Inici del protocol per calcular OR ----
  if (is.factor(data_util[[var_dep]]) && nlevels(data_util[[var_dep]]) == 2) {

    for(i in seq_along(var_indep)){

      ###### Check inicial ----

      # En cas que var_indep[i] no es trobi en els noms de adjust_sort
      # els OR que es veuran seran els crus.

      if(!(var_indep[i] %in% names(adjust_sort))){
        text_adj[i] <- paste0(var_indep[i], " is no adjusted.")

        fila_ref <- which(is.na(attr(compare_object[[i]], "p.ratio")))

        matriu_OR <- attr(compare_object [[i]], "OR")

        if(is.null(rownames(matriu_OR))){
          rownames(matriu_OR) <- var_indep_label[i]
        }

        if (length(fila_ref) != 0){
          rnames <- rownames(matriu_OR)[-fila_ref]
          matriu_OR <- matrix(matriu_OR[-fila_ref, ], ncol = 3)
          rownames(matriu_OR) <- rnames
        }

        matriu_estim <- cbind(rownames(matriu_OR), matriu_OR)

        if (!(var_indep_label[i] %in% matriu_estim[, 1])){
          matriu_estim[, 1] <- paste(var_indep_label[i], "-", matriu_estim[, 1])
        }
        matriu_forest <- rbind(matriu_forest, matriu_estim)



        # En cas contrari, es crea un vector "adjust_vars" amb les variables de dins.
        # I es procedeix al calcul dels OR ajustats

      } else if(var_indep[i] %in% names(adjust_sort)){

        adjust_vars <- unlist(adjust_sort[grep(var_indep[i], names(adjust_sort))])


        text_adj[i] <- paste0(var_indep[i], " is adjusted by ", paste(adjust_vars, collapse = ", "))


        ##### Crea la formula del modelo ajustado(que pude ser solo crudo también) y se ajusta el modelo correspondiente ----

        formula_mod <- as.formula(paste(var_dep, "~", paste(c(var_indep[i], adjust_vars), collapse = " + ")))
        formula_initial <- as.formula(paste(var_dep, "~", var_indep[i]))

        model_glm <- glm(formula_mod, data = data_util, family = binomial)



        model_glm_initial <- glm (formula_initial, data = data_util, family = binomial)

        ##### Numero de dades utilitzades ----

        N_mod <- nrow_data - length(model_glm$na.action)
        compare_object[[i]]$sam[1] <- N_mod

        ##### Extracció de coeficinetes y outputs necesarios para la tabla final ----

        # 'coef_limit' s'utiltiza per saber quants dels coeficients són part de la variable sola i així triar els coeficients que necesitem
        coef_limit <- length(coef(model_glm_initial))

        # guardem els objectes en una matriu per després afegir-los en el seu lloc

        fact_ratio <- attr(compare_object[[i]], "fact.ratio")

        matriu_estim <- cbind(
          exp(   coef(model_glm) * fact_ratio  )[2:coef_limit],
          exp( suppressMessages(confint(model_glm)) * fact_ratio )[2:coef_limit,1],
          exp( suppressMessages(confint(model_glm)) * fact_ratio )[2:coef_limit,2]
        )

        # extreiem p_valors

        p_values_mod <- summary(model_glm)$coefficients[2:coef_limit, "Pr(>|z|)"]

        # calculem el p-overall

        p_overall <- drop1(update(model_glm, data = model.frame(model_glm)), test= "Chisq")

        compare_object[[i]]$p.overall <- p_overall[var_indep[i],"Pr(>Chi)"]

        # busquem quina és la fila de referencia per tal de deixar-la com està (ho fem amb el p.ratio perquè només hem de buscar un NA)

        fila_ref <- which(is.na(attr(compare_object[[i]], "p.ratio")))

        #### Substituim els valors obtinguts a 'compare_object'----

        # Primer guardem els noms de les files i columnes per deixar els mateixos
        col_names_OR <- colnames(attr(compare_object [[i]], "OR"))
        row_names_OR <- rownames(attr(compare_object [[i]], "OR"))

        if(length(fila_ref) != 0){
          rnames <- rownames(attr(compare_object [[i]], "OR"))[-fila_ref ]
          attr(compare_object [[i]], "OR")[-fila_ref, ] <- matriu_estim
          attr(compare_object [[i]], "p.ratio")[-fila_ref] <- p_values_mod
          rownames(matriu_estim) <- rnames


        } else {
          attr(compare_object [[i]], "OR") <- matriu_estim
          attr(compare_object [[i]], "p.ratio") <- p_values_mod
          rownames(matriu_estim) <- var_indep_label[i]
        }

        colnames(attr(compare_object [[i]], "OR")) <- col_names_OR
        rownames(attr(compare_object [[i]], "OR")) <- row_names_OR

        #### Preparació  matriu de forest ----

        # Afegim la columna de noms a la 'matriu_estim' per així afegir-la a la 'matriu_forest' per després poder fer el forest_plot

        matriu_estim <- cbind(rownames(matriu_estim), matriu_estim)

        if (!(var_indep_label[i] %in% matriu_estim[, 1])){
          matriu_estim[, 1] <- paste(var_indep_label[i], "-", matriu_estim[, 1])
        }
        matriu_forest <- rbind(matriu_forest, matriu_estim)

      }

    }


    ### Supervivència ----

    ## Jo canviaria el is.Surv(data_util[[var_dep]]) per inherits(data_util[[var_dep]], "Surv")
  } else if  (is.Surv(data_util[[var_dep]])) {

    for(i in seq_along(var_indep)){

      ###### Check inicial ----

      # En cas que var_indep[i] no es trobi en els noms de adjust_sort
      # els OR que es veuran seran els crus.

      if(!(var_indep[i] %in% names(adjust_sort))){
        text_adj[i] <- paste0(var_indep[i], " is no adjusted.")

        fila_ref <- which(is.na(attr(compare_object[[i]], "p.ratio")))

        matriu_OR <- attr(compare_object [[i]], "HR")

        if(is.null(rownames(matriu_OR))){
          rownames(matriu_OR) <- var_indep_label[i]
        }

        if (length(fila_ref) != 0){
          rnames <- rownames(matriu_OR)[-fila_ref]
          matriu_OR <- matrix(matriu_OR[-fila_ref, ], ncol = 3)
          rownames(matriu_OR) <- rnames
        }

        matriu_estim <- cbind(rownames(matriu_OR), matriu_OR)

        if (!(var_indep_label[i] %in% matriu_estim[, 1])){
          matriu_estim[, 1] <- paste(var_indep_label[i], "-", matriu_estim[, 1])
        }
        matriu_forest <- rbind(matriu_forest, matriu_estim)

      }

      # En cas contrari, es crea un vector "adjust_vars" amb les variables de dins.
      # I es procedeix al calcul dels OR ajustats

      else if(var_indep[i] %in% names(adjust_sort)){
        adjust_vars <- unlist(adjust_sort[grep(var_indep[i], names(adjust_sort))])


        text_adj[i] <- paste0(var_indep[i], " is adjusted by ", paste(adjust_vars, collapse = ", "))


        ##### Crea la formula del modelo ajustado(que pude ser solo crudo también) y se ajusta el modelo correspondiente ----

        formula_mod <- as.formula(paste(var_dep, "~", paste(c(var_indep[i], adjust_vars), collapse = " + ")))
        formula_initial <- as.formula(paste(var_dep, "~", var_indep[i]))

        model_Surv <- coxph(formula_mod, data = data_util)
        model_Surv_initial <- coxph(formula_initial, data = data_util)

        ##### Numero de dades utilitzades ----

        N_mod <- nrow_data - length(model_Surv$na.action)
        compare_object[[i]]$sam[1] <- N_mod

        ##### Extracció de coeficinetes y outputs necesarios para la tabla final ----

        # 'coef_limit' s'utiltiza per saber quants dels coeficients són part de la variable sola i així triar els coeficients que necesitem

        coef_limit <- length(coef(model_Surv_initial))

        # guardem els objectes en una matriu per després afegir-los en el seu lloc, multipilcant pel fact.ratio

        fact_ratio <- attr(compare_object[[i]], "fact.ratio")

        matriu_estim <- cbind(
          exp(   coef(model_Surv) * fact_ratio  )[1:coef_limit],
          exp( suppressMessages(confint(model_Surv)) * fact_ratio )[1:coef_limit,1],
          exp( suppressMessages(confint(model_Surv)) * fact_ratio )[1:coef_limit,2]
        )

        # extreiem p_valors

        p_values_mod <- summary(model_Surv)$coefficients[1:coef_limit, "Pr(>|z|)"]

        # calculem el p-overall

        p_overall <- drop1(update(model_Surv, data = model.frame(model_Surv)), test= "Chisq")
        compare_object[[i]]$p.overall <- p_overall[var_indep[i],"Pr(>Chi)"]

        # busquem quina és la fila de referencia per tal de deixar-la com està

        fila_ref <- which(is.na(attr(compare_object[[i]], "p.ratio")))

        #### Substituim els valors obtinguts a 'compare_object'----

        # Primer guardem els noms de les files i columnes per deixar els mateixos
        col_names_HR <- colnames(attr(compare_object [[i]], "HR"))
        row_names_HR <- rownames(attr(compare_object [[i]], "HR"))


        if(length(fila_ref) != 0){
          rnames <- rownames(attr(compare_object [[i]], "HR"))[-fila_ref ]
          attr(compare_object [[i]], "HR")[-fila_ref, ] <- matriu_estim
          attr(compare_object [[i]], "p.ratio")[-fila_ref] <- p_values_mod
          rownames(matriu_estim) <- rnames

        } else {
          attr(compare_object [[i]], "HR") <- matriu_estim
          attr(compare_object [[i]], "p.ratio") <- p_values_mod
          rownames(matriu_estim) <- var_indep_label[i]
        }

        colnames(attr(compare_object [[i]], "HR")) <- col_names_HR
        rownames(attr(compare_object [[i]], "HR")) <- row_names_HR

        #### Preparació  matriu de forest ----

        # Afegim la columna de noms a la 'matriu_estim' per així afegir-la a la 'matriu_forest' per després poder fer el forest_plot

        matriu_estim <- cbind(rownames(matriu_estim), matriu_estim)

        if (!(var_indep_label[i] %in% matriu_estim[, 1])){
          matriu_estim[, 1] <- paste(var_indep_label[i], "-", matriu_estim[, 1])
        }
        matriu_forest <- rbind(matriu_forest, matriu_estim)

      }

    }

  } else {stop("Dependent variable must be a binary factor or a survival object.")}


  ## Creacció de la taula de resultat ----

  table_adjust <- cbind("Raw" = base,
                        "Adjusted" = createTable(compare_object,
                                                 show.ratio = show.ratio,
                                                 show.descr = FALSE,
                                                 hide = attr(base, "hide"),
                                                 hide.no = attr(base, "hide.no"),
                                                 show.n = show.n,
                                                 show.p.ratio = show.p.ratio,
                                                 show.p.overall = show.p.overall))


  # Guardar atributs ----

  results <- table_adjust

  attr(results, "matriu_forest") <- matriu_forest

  attr(results, "data_util") <- data_util

  attr(results, "var_dep") <- var_dep

  attr(results, "var_indep_label") <- var_indep_label

  attr(results, "text_adj") <- text_adj

  # Crear class propi de l'objecte ----

  class(results) <- c(class(results), "add_adjusted_ratio")

  # Retornar la taula ----

  return(results)

}



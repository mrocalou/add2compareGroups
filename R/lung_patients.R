#' Training data
#'
#' curated and semi-synthetic version of the original dataset from an oncology patient cohort, designed to support reproducible analysis in lung cancer reserch.
#'
#' @format A data frame with 19 columns
#' \describe{
#'   \item{\code{PACIENT_ID}}{Unique patient code}
#'   \item{\code{GENDER}}{Patient sex}
#'   \item{\code{AGE}}{Age at diagnosis}
#'   \item{\code{DEATH}}{Vital status at last follow-up}
#'   \item{\code{SURV_MONTHS}}{Months from diagnosis to last follow-up or death}
#'   \item{\code{SMOKING}}{Tabacco consumption status}
#'   \item{\code{PACK_YEARS}}{Estimated cumulative tabacco exposure}
#'   \item{\code{CHEMO}}{Chemotherapy administered}
#'   \item{\code{RADIO}}{UniqRadiotherapy administeredue}
#'   \item{\code{ECOG_C}}{Performance status}
#'   \item{\code{SURGERY}}{Whether the patient underwent surgery}
#'   \item{\code{DISEASE_FREE}}{Free of disease at last assessment}
#'   \item{\code{RECURRENCE}}{Evidence of recurrence}
#'   \item{\code{PROGRESSION}}{Evidence of disease progression}
#'   \item{\code{T_STAGE}}{Clinical T descriptor}
#'   \item{\code{N_STAGE}}{Clinical N descriptor}
#'   \item{\code{M_STAGE}}{Clinical M descriptor}
#'   \item{\code{STAGE}}{Clinical stage at diagnosis}
#'   \item{\code{SURG_INTENT}}{Surgical intention}
#' }
#'
"lung_patients"

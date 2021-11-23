#' import clinical data from Brennan 2016
#' @source \url{https://www.cell.com/cms/10.1016/j.cell.2013.09.034/attachment/9cefc2e8-caac-4225-bcdd-70f105ccf568/mmc7.xlsx}
#' @note CSV extracted from XLSX; names normalized by janitor::clean_names, empty column deleted,
#' `case_id` added as rownames.
#' @examples
#' head(get_brennan_clin)
#' @export
get_brennan_clin = function() {
 tmp = read.csv(system.file("extdata/clinsubtrimmed.csv", package="BiocAIML"))
 tmp = tmp[,-16L]
 names(tmp) = c("case_id", "secondary_or_recurrent", "age_at_procedure", "gender", 
   "path_dx", "mgmt_status", "x2012_methylation_class", "g_cimp_methylation", 
   "idh1_status", "expression_subclass", "therapy_class", "vital_status", 
   "os_days", "progression_status", "pfs_days")  # produced by janitor
 rownames(tmp) = tmp$case_id
 tmp
}

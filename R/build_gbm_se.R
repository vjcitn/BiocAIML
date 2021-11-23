#' use curatedTCGAData and Brennan metadata to build an annotated GBM MAE
#' @import curatedTCGAData
#' @importFrom S4Vectors DataFrame
#' @import BiocFileCache
#' @param cache if not NULL, location of BiocFileCache instance
#' @note Will cache the created SummarizedExperiment instance unless told not to.
#' `rname` is `"brenn_gbm_se"`.
#' @examples
#' gbmse = build_gbm_se()
#'
#' # check effect of MGMT methylation on survival
#'
#' xm = gbmse[, gbmse$mgmt_status !="" & gbmse$vital_status !=""]
#' library(survival)
#' ss = Surv(xm$os_days, 1*(xm$vital_status=="DECEASED"))
#' plot(survfit(ss~xm$mgmt_status), lty=1:2)
#' legend(1100, .95, lty=c(1,2), legend=c("MGMT methylated", "unmethylated"))
#' title("Time-on-study/vital status for 123 GBM patients analyzed in Brennan et al. PMID 24120142")
#' survdiff(ss~xm$mgmt_status)
#' @export
build_gbm_se = function(cache=BiocFileCache::BiocFileCache()) {
 if (!is.null(cache)) {
   chk = bfcquery(cache, "brenn_gbm_se")
   if (nrow(chk)>0) return(readRDS(chk$rpath))
   }
 mae = curatedTCGAData("GBM", assay="RNAseq2GeneNorm", dry.run=FALSE, version="2.0.1")
 se = experiments(mae)[[1]]
 colnames(se) = substr(colnames(se), 1, 12)
 br = get_brennan_clin()
 okids = intersect(colnames(se), rownames(br))
 br = br[okids,]
 se = se[,okids]
 colData(se) = DataFrame(br)
 assayNames(se) = "RNAseq2GeneNorm"
 if (is.null(cache)) return(se)
 tf = tempfile()
 saveRDS(se, tf)
 bfcadd(cache, rname="brenn_gbm_se", fpath=tf, action="copy")
 readRDS(bfcquery(cache, "brenn_gbm_se")$rpath)
}
 



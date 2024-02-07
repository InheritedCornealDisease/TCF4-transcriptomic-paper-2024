# VEP annotations instead of CoCoRV's ANNOVAR

case_snps = function(data, extraParamJason = NULL) {
  outIndex = logical(dim(data)[1])
  outIndex[ ] = F

  cadd <- data[, "CADD_PHRED"] > 15 & !is.na(data[, "CADD_PHRED"])
  eaf <- data[, "gnomADe_AF"] < 0.005 | !is.na(data[, "gnomADe_AF"])
  gaf <- data[, "gnomADg_AF"] < 0.005 | !is.na(data[, "gnomADg_AF"])
  outIndex <- (cadd & eaf & gaf) 

  return(outIndex)
}
context("Test prep")
source("../src/prep.R")
load.data <- function(){
  # Read the original gene set collection
  require("GSEABase") || stop("Package GSEABase is not available!")
  original.gsc.address <- "../data/GenesetCollectionGS1-5/OriginalGS1-GS5.gmt"
  original.gsc <- GSEABase::getGmt(original.gsc.address, 
                                   geneIdType=EntrezIdentifier(),
                                   collectionType=BroadCollection()
                                   )
  annotation <- "hu6800"
  expression.profile.address <- paste("../data/ArtificialExpressionData/",
                                      "ArtificialExpressionProfilewithEntrezID.csv",
                                      sep="")
  background <- rownames(read.csv(expression.profile.address))
  return(list(original.gsc=original.gsc, annotation=annotation, background=background))
}


test_that("Gene sets are filtered as they should.", {
  data <- load.data()
  original.gsc <- data$original.gsc
  annotation <- data$annotation
  background <- data$background
  # Filter Geneset1
  gsc1to5 <- prep.genesets(original.gsc, background,
                           min.size=11, max.size=Inf)
  expect_equal(names(gsc1to5), c("Geneset2","Geneset3","Geneset4","Geneset5"))
  # filter Geneset5
  gsc1to5 <- prep.genesets(original.gsc, background,
                           min.size=1, max.size=49)
  expect_equal(names(gsc1to5), c("Geneset1","Geneset2","Geneset3","Geneset4"))
  # filter Geneset1 and Geneset5
  gsc1to5 <- prep.genesets(original.gsc, background,
                           min.size=11, max.size=49)
  expect_equal(names(gsc1to5), c("Geneset2","Geneset3","Geneset4"))
  
})

test_that("Gene sets are loaded correctly.", {
  # Load gene set collection
  original.gsc.address <- "../data/GenesetCollectionGS1-5/OriginalGS1-GS5.gmt"
  gsc1to5 <- prep.loadGMT(original.gsc.address)
  genes <- geneIds(gsc1to5[["Geneset1"]])
  expect_equal(genes, c("6558", "2743", "2566", "1187", "2554", "2741", "2742",
                       "1811", "6569", "6559"))

  genes <- geneIds(gsc1to5[["Geneset2"]])
  expect_equal(genes, c("1120", "3931", "3930", "3612", "2222", "2224", "337",
                        "3422", "1119", "3158", "1584", "5336", "5724", "3157",
                        "9663", "2805", "47", "2171", "8644", "6307"))

  genes <- geneIds(gsc1to5[["Geneset3"]])
  expect_equal(genes, c("1033", "6714", "4221", "896", "5629", "4171", "1956",
                        "836", "4745", "595", "7023", "1032", "5728", "3087",
                        "2810", "11140", "8850", "890", "1030", "892", "1164",
                        "7026", "8452", "2965", "7128", "898", "1022", "1029",
                        "1050", "1647"))

  genes <- geneIds(gsc1to5[["Geneset4"]])
  expect_equal(genes, c("7520", "7514", "6051", "829", "5162", "4841", "3421",
                        "25814", "1964", "4582", "3475", "9341", "9686",
                        "10128", "4190", "10212", "5538", "3030", "2058",
                        "7417", "8887", "6721", "3066", "6601", "378", "648",
                        "23173", "1176", "4676", "821", "8655", "2958", "4012",
                        "9632", "10902", "10657", "9375", "23196", "2783",
                        "889"))

  genes <- geneIds(gsc1to5[["Geneset5"]])
  expect_equal(genes, c("3151", "3184", "10370", "10521", "3015", "10274",
                        "5292", "5422", "9937", "6223", "5981", "5833", "8554",
                        "1114", "23168", "2047", "9037", "688", "4091", "655",
                        "4753", "2146", "1616", "29902", "9824", "2176", "990",
                        "7353", "983", "7091", "5888", "7374", "23649", "2329",
                        "5122", "9230", "4678", "4500", "6567", "5511",
                        "23338", "10971", "22913", "6532", "9088", "23468",
                        "11073", "7465", "3178", "3925"))
  
})

test_that("Expression profile is loaded correctly.",{
  address <- "../data/ArtificialExpressionData/ArtificialExpressionProfile.csv"
  annotation <- "hu6800"
  contrast <- c("c", "c", "c", "c", "c", "d", "d", "d", "d", "d")
  eset <- prep.loadExpressionSet(address, contrast, annotation, sep=",")
  expect_equal(dim(exprs(eset)), c(150, 10))
  # Check the first row of expression profile
  expected <- c(0.833733170755138, 0.258213780476749, 1.15123169923349,
                1.28549945572469, 0.153504324579999, 3.67755590547646,
                1.06210247077319, 2.05952695737273, 2.43278278242153,
                1.47164238587154)
  difference <- as.matrix(exprs(eset)["U30246_at", ] - expected)
  expect_true(norm(difference) < 1E-7)
  # Check the last row of expression profile
  expected <- c(0.380864742480641, 2.20953097005084, 1.06034880325832,
                0.67559253538877, 0.54429758077997, 0.056446768842962,
                1.29547400508158, 0.478497082309476, 1.58601296107319,
                0.63542425470372)
  difference <- as.matrix(exprs(eset)["M31303_rna1_at", ] - expected)
  expect_true(norm(difference) < 1E-7)
  })

test_that("prep.makeFGeneSetCollection works as expected.", {
  genesets <- list()
  genesets[["Geneset1"]] <- c("U30246_at", "U33267_at", "X15376_at",
                            "Z30643_at", "X14766_at", "X52009_s_at",
                            "X52008_at", "L02785_at", "L13258_at",
                            "X91220_at")
  genesets[["Geneset2"]] <- c("U62317_rna3_at", "M12625_at",
                            "L25931_s_at", "X66922_at", "X69141_at",
                            "Z47055_s_at", "J02758_s_at",
                            "X17025_at", "D10704_at","X83618_at",
                            "D16154_at", "M32879_at", "M37238_s_at",
                            "D10202_at", "L25798_at", "D87436_at",
                            "M37400_at", "X64330_at", "M94856_at",
                            "D17793_at", "U60205_at")
  collection <- prep.makeFGeneSetCollection(genesets)
  expect_identical(names(collection), c("Geneset1", "Geneset2"))
  expect_identical(geneIds(collection[["Geneset1"]]),
                   genesets[["Geneset1"]])
  expect_identical(geneIds(collection[["Geneset2"]]),
                   genesets[["Geneset2"]])
  })
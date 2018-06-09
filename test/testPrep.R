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
                                      "ArtificialExpressionProfile.csv",
                                      sep="")
  background <- rownames(read.csv(expression.profile.address))
  return(list(original.gsc=original.gsc, annotation=annotation, background=background))
}

test_that("Test prep.annotation.pkg.name", {
  expect_identical(prep.annotation.pkg.name("hu6800"), "hu6800.db")
  expect_identical(prep.annotation.pkg.name("hu6800.db"), "hu6800.db")
  expect_identical(prep.annotation.pkg.name("hu6800.db"), "hu6800.db")
  expect_error(prep.annotation.pkg.name(""))
  expect_error(prep.annotation.pkg.name(NA))
  expect_error(prep.annotation.pkg.name(NULL))
  }
)

test_that("Gene Ids as translated as they should",{
  data <- load.data()
  original.gsc <- data$original.gsc
  annotation <- data$annotation
  background <- data$background
  gsc1to5 <- prep.genesets(original.gsc, annotation, background,
                           min.size=1, max.size=Inf)
  expect_identical(length(gsc1to5), length(original.gsc))
  expect_equal(names(gsc1to5), c("Geneset1","Geneset2","Geneset3","Geneset4","Geneset5"))
  expect_equal(gsc1to5[["Geneset1"]], c("U30246_at", "U33267_at", "X15376_at",
                                         "Z30643_at", "X14766_at",
                                         "X52009_s_at", "X52008_at",
                                         "L02785_at", "L13258_at",
                                         "X91220_at"))
  expect_true(all(gsc1to5[["Geneset2"]] %in% c("U62317_rna3_at", "M12625_at",
                                               "L25931_s_at", "X66922_at",
                                               "X69141_at", "Z47055_s_at",
                                               "J02758_s_at", "X17025_at",
                                               "D10704_at", "X83618_at",
                                               "D16154_at", "M32879_at",
                                               "M37238_s_at", "D10202_at",
                                               "L25798_at", "D87436_at",
                                               "M37400_at", "X64330_at",
                                               "M94856_at", "D17793_at",
                                               "U60205_at")))
  expect_true(all(gsc1to5[["Geneset3"]] %in% c("L25876_at", "K03218_at",
                                           "U93237_rna2_at", "M92287_at",
                                           "U44060_at", "D21063_at",
                                           "X00588_at", "U13737_at",
                                           "D83017_s_at", "X59798_at",
                                           "S73885_s_at", "U40343_at",
                                           "U92436_at", "X67235_s_at",
                                           "X57348_s_at", "U43077_at",
                                           "U57317_at", "X51688_at",
                                           "L36844_at", "M74091_at",
                                           "X54942_at", "M64497_at",
                                           "U58089_at", "M95809_at",
                                           "M59465_at", "M74093_at",
                                           "X95406_at", "L20320_at",
                                           "U26727_at", "U34070_s_at",
                                           "M60974_s_at")))
  # Mapping Ids is not one-to-one
  expect_true(all(gsc1to5[["Geneset4"]] %in% c("M30938_at", "Y08614_at",
                                               "U49278_at", "U56637_at",
                                               "D90086_at", "U02493_at",
                                               "Z68129_cds1_at", "Z93784_at",
                                               "L18960_at", "J05582_s_at",
                                               "M35093_s_at", "X52228_at",
                                               "Y10313_at", "U64520_at",
                                               "D50911_at", "M92439_at",
                                               "D55654_at", "U90426_at",
                                               "U44772_at", "D16480_at",
                                               "X54326_at", "L08666_at",
                                               "U33821_at", "U02031_at",
                                               "U31814_at", "U66616_at",
                                               "M36341_at", "L13689_at",
                                               "D42084_at", "U91932_at",
                                               "U77456_at", "L10284_at",
                                               "U32944_at", "U14193_at",
                                               "D50810_at", "D38555_at",
                                               "X87613_at", "M88108_at",
                                               "U81006_at", "D80005_at",
                                               "M36429_s_at", "U90268_at")))

  expect_true(all(gsc1to5[["Geneset5"]] %in% c("X13546_rna1_at", "M94630_at",
                                            "U65093_at", "U59321_at", "M37583_at",
                                            "Z75330_at", "M16750_s_at",
                                            "M54915_s_at", "X06745_at",
                                            "D42045_at", "M81757_at",
                                            "L14922_at", "L24783_at",
                                            "D84307_at", "U78524_at",
                                            "Y00064_at", "D87440_at",
                                            "L40636_at", "U52840_at",
                                            "D14520_at", "U59914_at",
                                            "X51801_at", "D83018_at",
                                            "U61145_at", "AF006041_at",
                                            "U79274_at", "D13638_s_at",
                                            "X66894_s_at", "U77949_at",
                                            "U64444_at", "X05360_at",
                                            "M99439_at", "D14134_at",
                                            "X89398_cds2_at", "L24559_at",
                                            "Z11737_at", "X64810_at",
                                            "X79780_at", "M97856_at",
                                            "X97261_at", "X97261_r_at",
                                            "U05321_at", "U14575_at",
                                            "D87076_at", "X56468_at",
                                            "L38696_at", "L05568_at",
                                            "U56816_at", "L07515_at",
                                            "D87448_at", "X62048_at",
                                            "X79536_at", "M31303_rna1_at")))
}
)

test_that("Gene sets are filtered as they should.", {
  data <- load.data()
  original.gsc <- data$original.gsc
  annotation <- data$annotation
  background <- data$background
  # Filter Geneset1
  gsc1to5 <- prep.genesets(original.gsc, annotation, background,
                           min.size=11, max.size=Inf)
  expect_equal(names(gsc1to5), c("Geneset2","Geneset3","Geneset4","Geneset5"))
  # filter Geneset5
  gsc1to5 <- prep.genesets(original.gsc, annotation, background,
                           min.size=1, max.size=49)
  expect_equal(names(gsc1to5), c("Geneset1","Geneset2","Geneset3","Geneset4"))
  # filter Geneset1 and Geneset5
  gsc1to5 <- prep.genesets(original.gsc, annotation, background,
                           min.size=11, max.size=49)
  expect_equal(names(gsc1to5), c("Geneset2","Geneset3","Geneset4"))
  
})
if(!require("testthat")){
  install.packages("testthat")
  library("testthat")
}
testthat::test_dir("test/")
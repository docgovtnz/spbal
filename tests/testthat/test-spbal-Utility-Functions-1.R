# Validate spbal functions, features and parameter validation.

testthat::test_that("1. Verify internal function functions correctly.", {
  sf_object <- sf::st_read(system.file("shape/nc.shp", package="sf"))
  expect_error(spbal::getPanel(sf_object, 1),
               "spbal(getPanel) Simple file object does not contain a feature named panel_id.", fixed=TRUE)
})

# Validate and compare Halton point generators.

#testthat::test_that("2. Compare cppRSHalton() and cppRSHalton_br() same n, seeds.", {
#  chp <- cppRSHalton(n = 1000, seeds = c(123, 456))
#  chpbr <- cppRSHalton_br(n = 1000, seeds = c(123, 456))
#  expect_equal(chp[,2:3], chpbr$pts)
#})

testthat::test_that("3. Compare cppBASpts() and cppRSHalton_br() same n, bases.", {
  chp <- cppBASpts(n = 1000, bases = c(2, 3))
  chpbr <- cppRSHalton_br(n = 1000, bases = c(2, 3), seeds = chp$seeds)
  expect_equal(chp$pts, chpbr$pts)
})

testthat::test_that("4. Ensure seeds are random for BAS - testing setBASSeed().", {
  set.seed(1)

  ## We will make up a long and skinny shape inside a really big bounding box.
  ## Then sample and make sure nothing weird happens.
  vals <- data.frame(X = runif(30, 0, 10), Y = runif(30, 0, 100))
  bb.df <- c("xmin" = -10.222, "ymin" = -20.1525, "xmax" = 133.54, "ymax" = 525.223)
  bb <- sf::st_as_sfc(sf::st_bbox(bb.df))
  shp <- sf::st_as_sf(vals, coords = c("X", "Y")) |> sf::st_buffer(4) |> sf::st_union()
  shp <- sf::st_cast(sf::st_combine(shp), "POLYGON")

  set.seed(234)
  res <- NULL
  for( i in 1:1000 ){
    res <- rbind(res, findBASSeed(shp, bb, n = 1))
  }
  duplicates <- which(duplicated(res))

  ## Throw an error if any updates in the BAS start to not create all unique points.
  expect_length(duplicates, 0)
})

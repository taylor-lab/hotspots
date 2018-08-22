library(hotspots)
library(testthat)

test_that("test that it works", {
    
    expect_equal(2+3, 5)
    expect_equal(0.0+0.0, 0)
    expect_equal(5 + -5, 0)
    ## expect errors if input not numeric
    expect_error('a' +2)
    expect_error(10 + 'b')
    expect_error('a' + 'b')
    
})

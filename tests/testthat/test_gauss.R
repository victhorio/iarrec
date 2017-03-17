library(iarrec)
context("Gaussian related functions")

test_that("lgauss", {
	expect_lt(abs(lgauss(0, 0, 1) - -0.9189385332046727805633), 1e-15)
	expect_lt(abs(lgauss(100, 0, 1) - -5000.918938533204709529), 1e-15)
	expect_lt(abs(lgauss(30, 30, 50) - -2.87495003591874587201), 1e-15)
	expect_lt(abs(lgauss(1e-5, 2, 1e-5) - -199993.1624808006745297), 1e-10)
	expect_equal(lgauss(1e-20, 2, 1e-20), -2e20)

	expect_lt(mean(abs(lgauss(seq(0.1, 100, 0.001), 5, 5))) - 287.8423224894216900793, 1e-15)
})

test_that("lbigauss", {
	expect_lt(abs(lbigauss(0, 0, 0, 0, 1, 1, 0) - -1.837877066409345339082), 1e-15)
	expect_lt(abs(lbigauss(0, 0, 0, 0, 1, 1, 0.5) - -1.694036030183454721865), 1e-15)
	expect_lt(abs(lbigauss(0, 0, 0, 0, 5000, 5000, -539) - -10.34922581299003141453), 1e-15)
	expect_lt(abs(lbigauss(32, 48, -25, -8, 60, 98, 30) - -37.72519254298495638977), 1e-14)
})

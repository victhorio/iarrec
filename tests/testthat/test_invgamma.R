library(iarrec)
context("Inverse gamma related functions")


test_that("dinvgamma", {
	expect_lt(abs(dinvgamma(3, 3, 3) - 0.06131324019524036), 1e-15)
	expect_lt(abs(dinvgamma(100000, 3, 200000) - 5.413411329464517562947e-06), 1e-15)
	expect_lt(abs(dinvgamma(1e-4, 0.1, 0.1)), 1e-15)
	expect_lt(abs(dinvgamma(1e-5, 2, 1e-5) - 36787.94411714417219628), 1e-10)
	expect_lt(abs(dinvgamma(1e-20, 2, 1e-20) - 36787944117143998464), 3e5)

	expect_lt(mean(abs(dinvgamma(seq(0.1, 100, 0.001), 5, 5))) - 0.01000990978570895709177, 1e-15)
})

test_that("linvgamma", {
	expect_lt(abs(linvgamma(3, 3, 3) - -2.791759469228055401402e+00), 1e-15)
	expect_lt(abs(linvgamma(100000, 3, 200000) - -12.12663110385033604643), 1e-15)
	expect_lt(abs(linvgamma(1e-4, 0.1, 0.1) - -992.3515967518598017705), 1e-15)
	expect_lt(abs(linvgamma(1e-5, 2, 1e-5) - 1.051292546497022684093e+01), 2e-15)
	expect_lt(abs(linvgamma(1e-20, 2, 1e-20) - 4.505170185988090736373e+01), 1e-14)

	expect_lt(mean(linvgamma(seq(0.1, 100, 0.001), 5, 5)) - -17.14920651169919096901, 1e-15)
})

# The generator route inside a confluence band (#18) takes the same dose the
# closed forms take: rates only into the depot or central compartment, and NA
# for a steady state with no valid dose where the closed form gives NA.  The
# 2-cmt IV parameters below have a rate below 1e-4 of the largest, so every
# call is in the band.

.bandCall <- function(
  rate,
  type = 0L,
  tau = 0,
  tinf = 0,
  amt = 0,
  bolusCmt = 0L,
  dt = 0.5
) {
  linCmtModelDouble(
    dt,
    1e-5,
    10,
    1e-4,
    30,
    0,
    0,
    0,
    rep(0, 2 + 2 * 4),
    rate,
    2L,
    0L,
    1L,
    TRUE,
    type,
    tau,
    tinf,
    amt,
    bolusCmt,
    30L,
    3L
  )
}

test_that("a banded step reads only the depot/central rate", {
  .a <- .bandCall(c(20, 0))
  .b <- .bandCall(c(20, 5))
  expect_true(is.finite(.a$val))
  expect_equal(.b$val, .a$val)
  expect_equal(.b$Jg, .a$Jg)
})

test_that("a banded steady state with no valid dose is NA, as the closed forms", {
  expect_true(is.na(.bandCall(c(0, 0), type = 1L)$val))
  expect_true(is.na(
    .bandCall(c(0, 0), type = 3L, tau = 24, amt = 100, bolusCmt = 1L)$val
  ))
  expect_true(is.finite(.bandCall(c(0, 0), type = 3L, tau = 24, amt = 100)$val))
})

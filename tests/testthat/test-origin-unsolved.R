# The per-origin reads (which1 = -9 dose time, -10 bioavailability) at a
# record the solve never advances: one at the same time as the record before
# it -- an observation at a dose time, the trough sample of most regimens.
# It has no cached decomposition of its own and used to read NA, which an
# estimator reading the -10 sensitivity at every observation (nlmixr2est's
# FOCEi with f(depot) on an eta) took for a failed subject.
.originPar <- c(tcl = log(4), tv = log(20), tka = log(1.1), eta_f = 0.05, eta_lag = 0.1)

.originFd <- function(m, e, par, h = 1e-6) {
  p1 <- .originPar
  p1[[par]] <- p1[[par]] + h
  p0 <- .originPar
  p0[[par]] <- p0[[par]] - h
  (rxode2::rxSolve(m, e, params = p1)$cp - rxode2::rxSolve(m, e, params = p0)$cp) / (2 * h)
}

.originRel <- function(a, b) max(abs(a - b)) / (max(abs(b)) + 1e-8)

test_that("the bioavailability read answers at an observation at a dose time", {
  m <- rxode2::rxode2({
    cl <- exp(tcl); v <- exp(tv); ka <- exp(tka)
    f(depot) <- exp(eta_f)
    cp <- linCmtB(rx__PTR__, t, 2, 1, 1, -1, -1, 1, cl, v, 0, 0, 0, 0, ka)
    dF <- linCmtB(rx__PTR__, t, 2, 1, 1, -10, 7, 1, cl, v, 0, 0, 0, 0, ka)
  })
  # observations at t = 0 and t = 12, both dose times
  e <- rxode2::et(amt = 100, cmt = "depot", ii = 12, addl = 1) |>
    rxode2::et(c(0, 0.5, 2, 6, 12, 12.5, 18, 24))
  s <- rxode2::rxSolve(m, e, params = .originPar)
  expect_false(anyNA(s$dF))
  expect_lt(.originRel(s$dF, .originFd(m, e, "eta_f")), 1e-6)
})

test_that("the dose-time read answers at an observation at another compartment's dose time", {
  # a lagged depot alongside an unlagged IV dose into central at t = 12: the
  # observation at t = 12 follows the IV dose at the same time
  m <- rxode2::rxode2({
    cl <- exp(tcl); v <- exp(tv); ka <- exp(tka)
    lag <- 2 * exp(eta_lag)
    alag(depot) <- lag
    cp <- linCmtB(rx__PTR__, t, 2, 1, 1, -1, -1, 1, cl, v, 0, 0, 0, 0, ka)
    d9 <- lag * linCmtB(rx__PTR__, t, 2, 1, 1, -9, 7, 1, cl, v, 0, 0, 0, 0, ka)
  })
  e <- rxode2::et(amt = 100, cmt = "depot") |>
    rxode2::et(amt = 40, cmt = "central", time = 12) |>
    rxode2::et(c(0.5, 4, 8, 12, 12.5, 18, 24))
  s <- rxode2::rxSolve(m, e, params = .originPar)
  expect_false(anyNA(s$d9))
  expect_lt(.originRel(s$d9, .originFd(m, e, "eta_lag")), 1e-6)
})

# The per-origin reads (which1 = -9 dose time, -10 bioavailability) at a
# record the solve never advances: one at the same time as the record before
# it -- an observation at a dose time.  It has no cached decomposition of its
# own and used to read NA, which nlmixr2est's FOCEi (reading -10 at every
# observation) took for a failed subject.

# which2 for the concentration through compartment q (0-based)
.ouW2 <- function(q) as.integer(q * 8L + 7L)

# central finite difference of the concentration in a scale (`scale`) or a
# shift (`shift`) of the doses into linCmt state `cmt`
.ouFd <- function(ev, p, ncmt, oral0, cmt, scale = FALSE, h = 1e-6) {
  .move <- function(s) {
    .e <- ev
    .d <- .e$cmt == cmt
    if (scale) {
      .e$amt[.d] <- .e$amt[.d] * (1 + s)
      .e$rate[.d] <- .e$rate[.d] * (1 + s)
    } else {
      .e$time[.d] <- .e$time[.d] + s
    }
    .hsRun(.e, p, ncmt, oral0)$lhs
  }
  (.move(h) - .move(-h)) / (2 * h)
}

.ouRel <- function(a, b) max(abs(a - b)) / (max(abs(b)) + 1e-8)

test_that("the bioavailability read answers at an observation at a dose time", {
  for (.ncmt in 1:3) {
    # depot boluses at 0 and 12 with observations at both dose times
    .ev <- .hsEv(1L, times = c(0, 0.5, 2, 6, 12, 12.5, 18, 24))
    .p <- .hsPar(.ncmt, 1L)
    .b <- .hsRun(.ev, .p, .ncmt, 1L, useB = TRUE, calls = cbind(-10L, .ouW2(0L)), originMask = 1L)
    expect_false(anyNA(.b$calls))
    expect_lt(.ouRel(.b$calls[, 1], .ouFd(.ev, .p, .ncmt, 1L, 0L, scale = TRUE)), 1e-6)
  }
})

test_that("the dose-time read answers at an observation at another compartment's dose time", {
  for (.ncmt in 1:3) {
    # a depot bolus at 0 and a central bolus at 12, observed at 12
    .ev <- data.frame(
      time = c(0, 0.5, 4, 8, 12, 12, 12.5, 18, 24),
      cmt = c(0L, -1L, -1L, -1L, 1L, -1L, -1L, -1L, -1L),
      amt = c(100, 0, 0, 0, 40, 0, 0, 0, 0),
      rate = 0
    )
    .p <- .hsPar(.ncmt, 1L)
    .b <- .hsRun(.ev, .p, .ncmt, 1L, useB = TRUE, calls = cbind(-9L, .ouW2(0L)), originMask = 1L)
    .ok <- .ev$time > 0
    expect_false(anyNA(.b$calls[.ok, 1]))
    expect_lt(.ouRel(.b$calls[.ok, 1], .ouFd(.ev, .p, .ncmt, 1L, 0L)[.ok]), 1e-6)
  }
})

test_that("the dose-time read at an infusion start uses the started rate", {
  # a depot bolus at 0 and a central infusion over 12 to 14, observed as it
  # starts: the read is the derivative after the start, a one-sided
  # difference moving the infusion earlier
  .ev <- data.frame(
    time = c(0, 4, 12, 12, 13, 14, 18),
    cmt = c(0L, -1L, 1L, -1L, -1L, 1L, -1L),
    amt = c(100, 0, 0, 0, 0, 0, 0),
    rate = c(0, 0, 50, 0, 0, -50, 0)
  )
  .p <- .hsPar(2L, 1L)
  .b <- .hsRun(.ev, .p, 2L, 1L, useB = TRUE, calls = cbind(-9L, .ouW2(1L)), originMask = 2L)
  .h <- 1e-6
  .e <- .ev
  .e$time[.e$cmt == 1L] <- .e$time[.e$cmt == 1L] - .h
  .fd <- (.hsRun(.ev, .p, 2L, 1L)$lhs - .hsRun(.e, .p, 2L, 1L)$lhs) / .h
  expect_equal(.b$calls[4, 1], .fd[4], tolerance = 1e-5)
  expect_equal(.b$calls[5, 1], .ouFd(.ev, .p, 2L, 1L, 1L)[5], tolerance = 1e-5)
})

test_that("the per-origin reads stay NA when the model tracks no origin", {
  .ev <- .hsEv(1L, times = c(0, 0.5, 2, 12, 12.5))
  .p <- .hsPar(1L, 1L)
  .b <- .hsRun(.ev, .p, 1L, 1L, useB = TRUE, calls = cbind(c(-9L, -10L), .ouW2(0L)))
  expect_true(all(is.na(.b$calls)))
})

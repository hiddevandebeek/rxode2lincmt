# The carry kept in a solve buffer (an rxode2 sending op->numLinCarry), against
# the fixed array and the R reference

test_that("the carry in a solve buffer matches the fixed array", {
  set.seed(45)
  .ev <- c(101L, 0L, .carryEvid(5L, 0L), 0L, 3L, 0L)
  .d <- c(100, 0, 2, 0, 0, 0)
  for (.fast in c(FALSE, TRUE)) {
    expect_identical(
      .carryRun1(.ev, .d, fast = .fast, inSolve = TRUE),
      .carryRun1(.ev, .d, fast = .fast)
    )
  }
  for (.sh in list(c(1L, 0L), c(2L, 1L), c(3L, 1L))) {
    .cs <- .carryCase(
      .carryEv(.carryEvid(4L, 0L), 40),
      sum(.sh),
      8L,
      thetaStep = 5L
    )
    for (.fast in c(FALSE, TRUE)) {
      .a <- .carryProbe(
        .cs$evid,
        .cs$dose,
        .cs$theta,
        .cs$M,
        .cs$add,
        .sh[1],
        .sh[2],
        fast = .fast
      )
      .b <- .carryProbe(
        .cs$evid,
        .cs$dose,
        .cs$theta,
        .cs$M,
        .cs$add,
        .sh[1],
        .sh[2],
        fast = .fast,
        inSolve = TRUE
      )
      expect_identical(.b, .a)
    }
  }
})

test_that("the carry in a solve buffer is not capped at 8 pairs", {
  set.seed(46)
  .cs <- .carryCase(.carryEv(.carryEvid(5L, 1L), 0.5), 4L, 40L, thetaStep = 4L)
  for (.fast in c(FALSE, TRUE)) {
    expect_lt(
      .carryErr(.cs, 3L, 1L, fast = .fast, inSolve = TRUE),
      1e-12,
      label = .fast
    )
  }
  expect_error(
    .carryProbe(.cs$evid, .cs$dose, .cs$theta, .cs$M, .cs$add, 3L, 1L),
    "1 to 8 pairs"
  )
})

test_that("the carry in a solve buffer sits after the other states", {
  set.seed(47)
  .cs <- .carryCase(.carryEv(.carryEvid(4L, 1L), 0), 3L, 3L, thetaStep = 3L)
  .a <- .carryProbe(.cs$evid, .cs$dose, .cs$theta, .cs$M, .cs$add, 2L, 1L)
  # the probe errors if a state before the carry is written
  .b <- .carryProbe(
    .cs$evid,
    .cs$dose,
    .cs$theta,
    .cs$M,
    .cs$add,
    2L,
    1L,
    inSolve = TRUE,
    pad = 5L
  )
  expect_identical(.b, .a)
})

test_that("a row with no carry makes the rest of the subject's carry NA", {
  .ev <- c(101L, 0L, 0L, 0L, 0L, 0L)
  .d <- c(100, 0, 0, 0, 0, 0)
  expect_identical(
    .carryRun1(.ev, .d, inSolve = TRUE, nAll = 3L),
    c(.carryRun1(.ev, .d)[1:3], rep(NA_real_, 3L))
  )
  # rows with no carry are not advanced
  .linCmtCarryFastStats(TRUE)
  .carryRun1(.ev, .d, inSolve = TRUE, nAll = 3L)
  expect_equal(.linCmtCarryFastStats(TRUE)[["advCalls"]], 3)
})

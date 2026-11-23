# The carry in this package's own per-subject storage (pointer slot 10,
# linCmtCarryOwn()), against the carry in the host's ind->linCmtCarryT.

.carryOwn <- function(own) {
  .Call("_rxode2lincmt_carryOwnSet", as.integer(own), PACKAGE = "rxode2lincmt")
}

test_that("owned storage gives the host array's carry to the bit", {
  .prev <- .carryOwn(0L)
  on.exit(.carryOwn(.prev))
  set.seed(7)
  for (.sh in list(c(1L, 0L), c(1L, 1L), c(2L, 1L), c(3L, 1L))) {
    .m <- sum(.sh)
    for (.np in c(1L, 5L, 32L)) {
      for (.c in seq_len(.m) - 1L) {
        .evs <- list(
          .carryEv(.carryEvid(4L, .c), 40),
          .carryEv(.carryEvid(5L, .c), 0.5),
          .carryEv(c(.carryEvid(4L, .c), 3L, .carryEvid(5L, .c)), c(30, 0, 2))
        )
        for (.ev in .evs) {
          for (.fast in c(FALSE, TRUE)) {
            .cs <- .carryCase(.ev, .m, .np, thetaStep = 6L)
            .run <- function() {
              .carryProbe(
                .cs$evid,
                .cs$dose,
                .cs$theta,
                .cs$M,
                .cs$add,
                .sh[1],
                .sh[2],
                fast = .fast
              )
            }
            .carryOwn(0L)
            .host <- .run()
            .carryOwn(1L)
            .own <- .run()
            expect_identical(.own, .host)
          }
        }
      }
    }
  }
})

test_that("owned storage grows past 32 columns to what the model writes", {
  .prev <- .carryOwn(1L)
  on.exit(.carryOwn(.prev))
  expect_identical(.linCmtCarryMaxPairs(), 1048576L)
  set.seed(11)
  for (.np in c(33L, 70L, 129L)) {
    .cs <- .carryCase(
      .carryEv(.carryEvid(5L, 1L), 0.5),
      3L,
      .np,
      thetaStep = 5L
    )
    expect_lt(.carryErr(.cs, 2L, 1L), 1e-12, label = .np)
    expect_lt(.carryErr(.cs, 2L, 1L, fast = TRUE), 1e-12, label = .np)
  }
})

test_that("an owned restart leaves the host array alone", {
  .prev <- .carryOwn(1L)
  on.exit(.carryOwn(.prev))
  .r <- .Call(`_rxode2lincmt_carryRestartProbe`, 8L, 3L)
  expect_true(all(.r[1:136] == 7))
  expect_true(is.nan(.r[137]))
})

test_that("the host array is the carry until rxode2 hands it over", {
  .prev <- .carryOwn(0L)
  on.exit(.carryOwn(.prev))
  .o <- .rxode2lincmtHostInfo()$offsets
  on.exit(.rxode2lincmtIniHost(list(.o)), add = TRUE)
  .w <- .o
  .w[2] <- 8L
  .rxode2lincmtIniHost(list(.w))
  expect_identical(.linCmtCarryMaxPairs(), 8L)
  expect_identical(.carryOwn(1L), 0L)
  expect_identical(.linCmtCarryMaxPairs(), 1048576L)
  expect_identical(.carryOwn(0L), 1L)
  expect_identical(.linCmtCarryMaxPairs(), 8L)
})

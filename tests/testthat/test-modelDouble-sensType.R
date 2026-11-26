# linCmtModelDouble() must accept the reverse-mode code 31 and the auto code
# 100 (which resolves to 31); an unknown code used to leave the Jacobian
# uninitialized.
.call3 <- function(sensType) {
  .a <- numeric(.lcNalast(3L, 1L, TRUE))
  .a[1:4] <- c(50, 20, 5, 2)
  linCmtModelDouble(
    0.7,
    1.0,
    20,
    2.0,
    40,
    0.5,
    60,
    1.1,
    .a,
    rep(0, 4),
    3L,
    1L,
    1L,
    TRUE,
    0L,
    0,
    0,
    0,
    0L,
    0L,
    as.integer(sensType),
    0.001
  )
}

test_that("sensType 31 and 100 give the forward-mode Jacobian", {
  .fwd <- .call3(30L)
  for (.st in c(3L, 31L, 100L)) {
    .r <- .call3(.st)
    expect_equal(.r$val, .fwd$val, tolerance = 1e-12)
    expect_equal(.r$J, .fwd$J, tolerance = 1e-10)
    expect_equal(.r$Jg, .fwd$Jg, tolerance = 1e-10)
  }
})

test_that("an unknown sensType errors instead of returning garbage", {
  expect_error(.call3(99L), "unsupported sensType")
})

.callS <- function(ncmt, oral0, sensType, ndiff = 0L, rate = c(0, 0)) {
  .a <- numeric(.lcNalast(ncmt, oral0, TRUE))
  .a[seq_len(ncmt + oral0)] <- c(50, 20, 5, 2)[seq_len(ncmt + oral0)]
  linCmtModelDouble(
    0.7,
    1.0,
    20,
    2.0,
    40,
    0.5,
    60,
    1.1,
    .a,
    rate[seq_len(1L + oral0)],
    as.integer(ncmt),
    as.integer(oral0),
    1L,
    TRUE,
    0L,
    0,
    0,
    0,
    0L,
    as.integer(ndiff),
    as.integer(sensType),
    1e-5
  )
}

for (.ncmt in 1:3) {
  for (.oral0 in 0:1) {
    test_that(sprintf("every sensType agrees with forward-mode AD (ncmt %d, oral %d)", .ncmt, .oral0), {
      skip_on_cran()
      for (.rate in list(c(0, 0), c(3, 4))) {
        .fwd <- .callS(.ncmt, .oral0, 3L, rate = .rate)
        expect_equal(.callS(.ncmt, .oral0, 32L, rate = .rate)$Jg, .fwd$Jg, tolerance = 1e-12)
        expect_equal(.callS(.ncmt, .oral0, 31L, rate = .rate)$Jg, .fwd$Jg, tolerance = 1e-10)
        for (.st in c(1L, 2L, 10L, 20L)) {
          .r <- .callS(.ncmt, .oral0, .st, rate = .rate)
          expect_equal(.r$val, .fwd$val, tolerance = 1e-12)
          expect_equal(.r$Jg, .fwd$Jg, tolerance = 1e-3, label = sprintf("sensType %d", .st))
        }
      }
    })
  }
}

test_that("a partial ndiff mask gives the requested columns of the full Jacobian", {
  # ka 1, p1 2, v1 4, p2 8, p3 16, p4 32, p5 64; columns are p1, v1, p2..p5, ka
  for (.ncmt in 1:3) {
    for (.oral0 in 0:1) {
      .full <- .callS(.ncmt, .oral0, 3L)
      for (.st in c(3L, 31L, 32L)) {
        .r <- .callS(.ncmt, .oral0, .st, ndiff = 2L + 4L + .oral0)
        expect_equal(.r$val, .full$val, tolerance = 1e-12)
        expect_equal(.r$J[, 1:2], .full$J[, 1:2], tolerance = 1e-10)
        if (.oral0 == 1L) {
          expect_equal(.r$J[, ncol(.r$J)], .full$J[, ncol(.full$J)], tolerance = 1e-10)
        }
      }
    }
  }
})

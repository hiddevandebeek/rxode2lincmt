test_that("clearance parameterization gives the expected micro-constants", {
  expect_equal(linCmtMicros(p1 = 2, v1 = 20, ncmt = 1), c(v = 20, k = 0.1))
  expect_equal(linCmtMicros(p1 = 2, v1 = 20, p2 = 3, p3 = 40, ncmt = 2), c(v = 20, k = 0.1, k12 = 0.15, k21 = 0.075))
  expect_equal(
    linCmtMicros(p1 = 2, v1 = 20, p2 = 3, p3 = 40, p4 = 1, p5 = 80, ncmt = 3),
    c(v = 20, k = 0.1, k12 = 0.15, k21 = 0.075, k13 = 0.05, k31 = 0.0125)
  )
})

test_that("rate-constant parameterization is the identity", {
  expect_equal(linCmtMicros(p1 = 0.1, v1 = 20, ncmt = 1, trans = 2), c(v = 20, k = 0.1))
  expect_equal(
    linCmtMicros(p1 = 0.1, v1 = 20, p2 = 0.15, p3 = 0.075, ncmt = 2, trans = 2),
    c(v = 20, k = 0.1, k12 = 0.15, k21 = 0.075)
  )
  expect_equal(
    linCmtMicros(p1 = 0.1, v1 = 20, p2 = 0.15, p3 = 0.075, p4 = 0.05, p5 = 0.0125, ncmt = 3, trans = 2),
    c(v = 20, k = 0.1, k12 = 0.15, k21 = 0.075, k13 = 0.05, k31 = 0.0125)
  )
})

test_that("micro-constants round-trip through the solComp eigenvalues", {
  .m <- linCmtMicros(p1 = 2, v1 = 20, p2 = 3, p3 = 40, ncmt = 2)
  .s <- solComp2(k10 = .m[["k"]], k12 = .m[["k12"]], k21 = .m[["k21"]])
  .K <- .lcRateMatrix(2L, 0L, cl = 2, v = 20, q = 3, vp = 40)
  expect_equal(sort(.s$L), sort(-eigen(.K)$values))
})

test_that("unsupported parameterizations are rejected", {
  expect_error(linCmtMicros(p1 = 2, v1 = 20, ncmt = 1, trans = 3))
  expect_error(linCmtMicros(p1 = 2, v1 = 20, ncmt = 3, trans = 4))
  expect_error(linCmtMicros(p1 = 2, v1 = 20, ncmt = 4))
})

# Macro parameterizations round-trip through the derived quantities: build
# each parameterization from the derived table of a known micro-constant set
# and convert it back.
.mDerived <- function(ncmt, p) {
  .calcDerived(ncmt, 0, -1, 2, c(as.list(p), list(0, 0, 0, 0, 0, 0, 0, 0, 0)), 0)
}

test_that("one-compartment macro parameterizations round-trip", {
  .d <- .mDerived(1, list(0.1, 20, 0, 0, 0, 0))
  .want <- c(v = 20, k = 0.1)
  expect_equal(linCmtMicros(p1 = .d$alpha, v1 = .d$vc, ncmt = 1, trans = 11), .want)
  expect_equal(linCmtMicros(p1 = .d$alpha, v1 = .d$A, ncmt = 1, trans = 10), .want)
})

test_that("two-compartment macro parameterizations round-trip", {
  .d <- .mDerived(2, list(0.1, 20, 0.15, 0.075, 0, 0))
  .want <- c(v = 20, k = 0.1, k12 = 0.15, k21 = 0.075)
  expect_equal(linCmtMicros(p1 = .d$cl, v1 = .d$vc, p2 = .d$q, p3 = .d$vss, ncmt = 2, trans = 3), .want)
  expect_equal(linCmtMicros(p1 = .d$alpha, v1 = .d$vc, p2 = .d$beta, p3 = .d$k21, ncmt = 2, trans = 4), .want)
  expect_equal(linCmtMicros(p1 = .d$alpha, v1 = .d$vc, p2 = .d$beta, p3 = .d$A / .d$B, ncmt = 2, trans = 5), .want)
  expect_equal(linCmtMicros(p1 = .d$alpha, v1 = 1 / .d$A, p2 = .d$beta, p3 = .d$B, ncmt = 2, trans = 11), .want)
  expect_equal(linCmtMicros(p1 = .d$alpha, v1 = .d$A, p2 = .d$beta, p3 = .d$B, ncmt = 2, trans = 10), .want)
})

test_that("three-compartment macro parameterizations round-trip", {
  .d <- .mDerived(3, list(0.1, 20, 0.15, 0.075, 0.05, 0.0125))
  .want <- c(v = 20, k = 0.1, k12 = 0.15, k21 = 0.075, k13 = 0.05, k31 = 0.0125)
  .m11 <- linCmtMicros(
    p1 = .d$alpha,
    v1 = 1 / .d$A,
    p2 = .d$beta,
    p3 = .d$B,
    p4 = .d$gamma,
    p5 = .d$C,
    ncmt = 3,
    trans = 11
  )
  .m10 <- linCmtMicros(
    p1 = .d$alpha,
    v1 = .d$A,
    p2 = .d$beta,
    p3 = .d$B,
    p4 = .d$gamma,
    p5 = .d$C,
    ncmt = 3,
    trans = 10
  )
  # the peripheral labels are interchangeable: compare as sets of pairs
  .pairs <- function(m) {
    .p <- rbind(m[c("k12", "k21")], m[c("k13", "k31")])
    .p[order(.p[, 1]), , drop = FALSE]
  }
  for (.m in list(.m11, .m10)) {
    expect_equal(.m[c("v", "k")], .want[c("v", "k")], tolerance = 1e-8)
    expect_equal(.pairs(.m), .pairs(.want), tolerance = 1e-8, ignore_attr = TRUE)
  }
})

test_that("the kernel accepts every parameterization", {
  .d2 <- .mDerived(2, list(0.1, 20, 0.15, 0.075, 0, 0))
  .d3 <- .mDerived(3, list(0.1, 20, 0.15, 0.075, 0.05, 0.0125))
  .d1 <- .mDerived(1, list(0.1, 20, 0, 0, 0, 0))
  .cases <- list(
    list(1L, 1L, c(2, 20, 0, 0, 0, 0)),
    list(1L, 11L, c(.d1$alpha, .d1$vc, 0, 0, 0, 0)),
    list(1L, 10L, c(.d1$alpha, .d1$A, 0, 0, 0, 0)),
    list(2L, 1L, c(2, 20, 3, 40, 0, 0)),
    list(2L, 3L, c(.d2$cl, .d2$vc, .d2$q, .d2$vss, 0, 0)),
    list(2L, 4L, c(.d2$alpha, .d2$vc, .d2$beta, .d2$k21, 0, 0)),
    list(2L, 5L, c(.d2$alpha, .d2$vc, .d2$beta, .d2$A / .d2$B, 0, 0)),
    list(2L, 11L, c(.d2$alpha, 1 / .d2$A, .d2$beta, .d2$B, 0, 0)),
    list(2L, 10L, c(.d2$alpha, .d2$A, .d2$beta, .d2$B, 0, 0)),
    list(3L, 1L, c(2, 20, 3, 40, 1, 80)),
    list(3L, 11L, c(.d3$alpha, 1 / .d3$A, .d3$beta, .d3$B, .d3$gamma, .d3$C)),
    list(3L, 10L, c(.d3$alpha, .d3$A, .d3$beta, .d3$B, .d3$gamma, .d3$C))
  )
  for (.c in .cases) {
    .ncmt <- .c[[1]]
    .p <- .c[[3]]
    .ref <- NULL
    for (.deriv in c(FALSE, TRUE)) {
      .a <- c(100, numeric(.lcNalast(.ncmt, 1L, .deriv) - 1L))
      .r <- linCmtModelDouble(
        1.5,
        .p[1],
        .p[2],
        .p[3],
        .p[4],
        .p[5],
        .p[6],
        1.1,
        .a,
        c(0, 0),
        .ncmt,
        1L,
        .c[[2]],
        .deriv,
        0L,
        0,
        0,
        0,
        0L,
        0L
      )
      expect_true(is.finite(.r$val), label = sprintf("ncmt %d trans %d", .ncmt, .c[[2]]))
      if (.deriv) {
        expect_true(all(is.finite(.r$Jg)))
        expect_equal(.r$val, .ref, tolerance = 1e-12)
      } else {
        .ref <- .r$val
      }
    }
  }
})

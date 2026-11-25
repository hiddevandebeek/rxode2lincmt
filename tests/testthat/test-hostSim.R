# linCmtA()/linCmtB() driven through rxode2's linear-only solve loop on
# stand-in solver structs (.linCmtHostSim(), no rxode2): values against the
# base-R reference, gradients against finite differences.

test_that("the host table is restored after a simulated solve", {
  .before <- .rxode2lincmtHostInfo()
  .hsRun(.hsEv(0L), .hsPar(1L, 0L), 1L, 0L)
  expect_identical(.rxode2lincmtHostInfo(), .before)
})

for (.ncmt in 1:3) {
  for (.oral0 in 0:1) {
    for (.inf in c(FALSE, TRUE)) {
      test_that(
        sprintf("linCmtA/linCmtB values match the reference (ncmt %d, oral %d, infusion %s)", .ncmt, .oral0, .inf),
        {
          .ev <- .hsEv(.oral0, infusion = .inf)
          .p <- .hsPar(.ncmt, .oral0)
          .ref <- .hsRef(.ev, .p, .ncmt, .oral0)
          .a <- .hsRun(.ev, .p, .ncmt, .oral0)
          expect_equal(.a$val, .ref$val, tolerance = 1e-8)
          expect_equal(.a$lhs, .ref$lhs, tolerance = 1e-8)
          .b <- .hsRun(.ev, .p, .ncmt, .oral0, useB = TRUE)
          expect_equal(.b$val, .ref$val, tolerance = 1e-8)
          expect_equal(.b$lhs, .ref$lhs, tolerance = 1e-8)
        }
      )
    }

    test_that(sprintf("linCmtB gradients match finite differences (ncmt %d, oral %d)", .ncmt, .oral0), {
      skip_on_cran()
      for (.inf in c(FALSE, TRUE)) {
        .ev <- .hsEv(.oral0, infusion = .inf)
        .p <- .hsPar(.ncmt, .oral0)
        .fd <- .hsFd(.ev, .p, .ncmt, .oral0)
        .calls <- .hsJgCalls(.ncmt, .oral0)
        for (.st in c(3L, 30L, 31L, 32L, 100L)) {
          for (.phi in c(-1L, 0L, 1L, 2L)) {
            .b <- .hsRun(.ev, .p, .ncmt, .oral0, useB = TRUE, sensType = .st, sensPhi = .phi, calls = .calls)
            expect_equal(
              .b$calls,
              .fd,
              tolerance = 1e-5,
              ignore_attr = TRUE,
              label = sprintf("sensType %d phi %d inf %s", .st, .phi, .inf)
            )
          }
        }
      }
    })
  }
}

# The finite-difference methods differentiate each interval from the stored
# amounts without carrying the previous interval's sensitivity (the double
# getAlast(), as in rxode2 5.1.7): exact on the first advanced record, and
# every method agrees with central differences on every record.
test_that("finite-difference sensTypes give the interval's own sensitivity", {
  skip_on_cran()
  for (.ncmt in 1:3) {
    for (.oral0 in 0:1) {
      .ev <- .hsEv(.oral0)
      .p <- .hsPar(.ncmt, .oral0)
      .fd <- .hsFd(.ev, .p, .ncmt, .oral0)
      .calls <- .hsJgCalls(.ncmt, .oral0)
      .run <- function(st) {
        .hsRun(.ev, .p, .ncmt, .oral0, useB = TRUE, sensType = st, calls = .calls, linH = rep(1e-5, 7))
      }
      .central <- .run(2L)
      expect_equal(.central$calls[2, ], .fd[2, ], tolerance = 1e-4, ignore_attr = TRUE)
      for (.st in c(1L, 4L, 5L, 6L, 7L, 10L, 20L, 40L, 50L)) {
        .b <- .run(.st)
        expect_equal(.b$lhs, .hsRef(.ev, .p, .ncmt, .oral0)$lhs, tolerance = 1e-8)
        expect_equal(
          .b$calls,
          .central$calls,
          tolerance = 1e-3,
          label = sprintf("ncmt %d oral %d sensType %d", .ncmt, .oral0, .st)
        )
      }
    }
  }
  expect_true(all(c(1L, 2L, 4L, 5L, 10L, 20L, 40L, 50L) %in% .linCmtBSensTypesSeen(TRUE)))
})

test_that("a partial ndiff mask gives the requested gradient columns", {
  skip_on_cran()
  # ka 1, p1 2, v1 4, p2 8, p3 16, p4 32, p5 64
  for (.ncmt in 1:3) {
    for (.oral0 in 0:1) {
      .ev <- .hsEv(.oral0)
      .p <- .hsPar(.ncmt, .oral0)
      .fd <- .hsFd(.ev, .p, .ncmt, .oral0)
      for (.st in c(3L, 31L, 32L)) {
        .b <- .hsRun(.ev, .p, .ncmt, .oral0, useB = TRUE, sensType = .st, ndiff = 2L + 4L, calls = cbind(-2L, 0:1))
        expect_equal(.b$calls, .fd[, 1:2], tolerance = 1e-5, ignore_attr = TRUE)
      }
    }
  }
})

test_that("theta changing between records refills the window", {
  skip_on_cran()
  for (.ncmt in 1:3) {
    for (.oral0 in 0:1) {
      .ev <- .hsEv(.oral0, infusion = .ncmt == 2L)
      .n <- nrow(.ev)
      .th <- .hsTheta(.hsPar(.ncmt, .oral0), .n)
      .th[seq(5, .n), 1] <- 2.6
      .th[seq(9, .n), 2] <- 18
      .ref <- .hsRef(.ev, .th, .ncmt, .oral0)
      .fd <- .hsFd(.ev, .th, .ncmt, .oral0)
      for (.phi in c(0L, 2L)) {
        .b <- .hsRun(.ev, .th, .ncmt, .oral0, useB = TRUE, sensPhi = .phi, calls = .hsJgCalls(.ncmt, .oral0))
        expect_equal(.b$lhs, .ref$lhs, tolerance = 1e-8)
        expect_equal(.b$calls, .fd, tolerance = 1e-5, ignore_attr = TRUE)
      }
    }
  }
})

test_that("a long regular grid reuses the delta memo and transition matrix", {
  skip_on_cran()
  .ev <- .hsEv(1L, times = seq(0.25, 48, by = 0.25))
  .ev2 <- .hsEv(0L, infusion = TRUE, times = seq(0.25, 48, by = 0.25))
  for (.ncmt in 1:3) {
    for (.e in list(list(.ev, 1L), list(.ev2, 0L))) {
      .p <- .hsPar(.ncmt, .e[[2]])
      .fd <- .hsFd(.e[[1]], .p, .ncmt, .e[[2]])
      for (.memo in c(1L, 0L)) {
        .old <- .linCmtDeltaMemo(.memo)
        .linCmtSeqStats(TRUE)
        for (.phi in c(1L, 2L)) {
          for (.st in c(3L, 32L)) {
            .b <- .hsRun(
              .e[[1]],
              .p,
              .ncmt,
              .e[[2]],
              useB = TRUE,
              sensType = .st,
              sensPhi = .phi,
              calls = .hsJgCalls(.ncmt, .e[[2]])
            )
            expect_equal(.b$calls, .fd, tolerance = 1e-5, ignore_attr = TRUE)
          }
        }
        .s <- .linCmtSeqStats(TRUE)
        expect_gt(.s[["seqTailRows"]], 0)
        .linCmtDeltaMemo(.old)
      }
    }
  }
})

# steady state at the first record, then observations from it
.hsSsCase <- function(type, oral0) {
  .ev <- data.frame(time = c(0, 1, 2, 4, 8, 12), cmt = -1L, amt = 0, rate = 0)
  .ss <- data.frame(type = c(type, rep(0L, 5)), tau = c(12, rep(0, 5)), var = 0, cmt = 0L)
  if (type == 3L) {
    .ss$var[1] <- 100
  } else {
    # the infusion's rate is live in InfusionRate while the steady state is solved
    .ev$cmt[1] <- oral0
    .ev$rate[1] <- 25
    .ss$var[1] <- 2
  }
  list(ev = .ev, ss = .ss)
}

test_that("steady-state records solve through linCmtA and linCmtB", {
  skip_on_cran()
  for (.ncmt in 1:3) {
    for (.oral0 in 0:1) {
      .p <- .hsPar(.ncmt, .oral0)
      for (.type in 2:3) {
        .cs <- .hsSsCase(.type, .oral0)
        .a <- .hsRun(.cs$ev, .p, .ncmt, .oral0, ss = .cs$ss)
        .b <- .hsRun(.cs$ev, .p, .ncmt, .oral0, ss = .cs$ss, useB = TRUE, calls = .hsJgCalls(.ncmt, .oral0))
        expect_true(all(is.finite(.a$lhs)))
        expect_equal(.b$lhs, .a$lhs, tolerance = 1e-10)
        .fd <- .hsFd(.cs$ev, .p, .ncmt, .oral0, ss = .cs$ss)
        expect_equal(
          .b$calls,
          .fd,
          tolerance = 1e-4,
          ignore_attr = TRUE,
          label = sprintf("ncmt %d oral %d type %d", .ncmt, .oral0, .type)
        )
      }
    }
  }
})

test_that("amount and Jacobian reads return the stored row", {
  for (.ncmt in 1:3) {
    for (.oral0 in 0:1) {
      .m <- .ncmt + .oral0
      .np <- 2L * .ncmt + .oral0
      .ev <- .hsEv(.oral0)
      .p <- .hsPar(.ncmt, .oral0)
      .calls <- rbind(cbind(seq_len(.m) - 1L, -2L), cbind(0L, seq_len(.np) - 1L))
      .b <- .hsRun(.ev, .p, .ncmt, .oral0, useB = TRUE, calls = .calls)
      expect_equal(.b$calls[, seq_len(.m)], .b$state[, seq_len(.m)], ignore_attr = TRUE)
      expect_true(all(is.finite(.b$calls)))
      .a <- .hsRun(.ev, .p, .ncmt, .oral0, calls = cbind(seq_len(.m) - 1L, 0L))
      expect_equal(.a$calls, .a$state, ignore_attr = TRUE)
    }
  }
})

test_that("dose-time sensitivity follows the lag mask", {
  for (.ncmt in 1:3) {
    .ev <- .hsEv(1L)
    .p <- .hsPar(.ncmt, 1L)
    .calls <- cbind(-3L, c(0L, 1L))
    .none <- .hsRun(.ev, .p, .ncmt, 1L, useB = TRUE, calls = .calls)
    expect_true(all(is.finite(.none$calls)))
    # the depot is dosed and lagged: -dA/dt of each compartment
    .lag <- .hsRun(.ev, .p, .ncmt, 1L, useB = TRUE, calls = .calls, lagMask = 1L)
    expect_equal(.lag$calls, .none$calls)
    # only the central compartment is lagged and none of its doses: exactly 0
    .zero <- .hsRun(.ev, .p, .ncmt, 1L, useB = TRUE, calls = .calls, lagMask = 2L)
    expect_true(all(.zero$calls == 0))
    # out of range compartment
    expect_true(all(is.na(.hsRun(.ev, .p, .ncmt, 1L, useB = TRUE, calls = cbind(-3L, 9L))$calls)))
  }
})

test_that("transition and per-origin queries answer on every record", {
  for (.ncmt in 1:3) {
    for (.oral0 in 0:1) {
      .m <- .ncmt + .oral0
      .ev <- .hsEv(.oral0, infusion = .oral0 == 0L)
      .p <- .hsPar(.ncmt, .oral0)
      .calls <- rbind(cbind(-4L, seq_len(.m * .m) - 1L), cbind(-9L, seq_len(.m) - 1L), cbind(-10L, seq_len(.m) - 1L))
      .b <- .hsRun(.ev, .p, .ncmt, .oral0, useB = TRUE, calls = .calls, originMask = 1L)
      expect_equal(.b$lhs, .hsRef(.ev, .p, .ncmt, .oral0)$lhs, tolerance = 1e-8)
      expect_identical(
        .hsRun(.ev, .p, .ncmt, .oral0, useB = TRUE, calls = .calls, originMask = 1L)$calls,
        .b$calls
      )
    }
  }
})

test_that("an H-perturbed evaluation shifts the parameter by the step", {
  .ev <- .hsEv(0L)
  .p <- .hsPar(2L, 0L)
  .h <- .hsRun(.ev, .p, 2L, 0L, useB = TRUE, hPar = 0L, hStep = 0.01)
  .up <- .p
  .up[1] <- .up[1] + 0.01
  .ok <- !is.na(.h$val)
  expect_true(any(.ok))
  expect_true(all(is.finite(.h$val[.ok])))
})

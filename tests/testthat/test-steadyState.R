# Steady-state kernels (type 1 infinite infusion, 2 repeated infusion, 3
# repeated bolus) against the base-R fixed point of the dosing cycle.
.ssPar <- c(cl = 2, v = 20, q = 3, vp = 40, q2 = 1, vp2 = 80, ka = 1.1)

.ssCall <- function(
  par,
  ncmt,
  oral0,
  rate,
  type,
  tau = 0,
  tinf = 0,
  amt = 0,
  bolusCmt = 0L,
  deriv = FALSE,
  sensType = 3L
) {
  .alast <- numeric(.lcNalast(ncmt, oral0, deriv))
  linCmtModelDouble(
    1.0,
    par[["cl"]],
    par[["v"]],
    par[["q"]],
    par[["vp"]],
    par[["q2"]],
    par[["vp2"]],
    par[["ka"]],
    .alast,
    as.double(rate),
    as.integer(ncmt),
    as.integer(oral0),
    1L,
    deriv,
    as.integer(type),
    tau,
    tinf,
    amt,
    as.integer(bolusCmt),
    0L,
    as.integer(sensType),
    0.001
  )
}

# full-state rate vector from the kernel's (depot, central) rates
.ssRateFull <- function(ncmt, oral0, rate) {
  .r <- numeric(ncmt + oral0)
  .r[seq_len(1L + oral0)] <- rate
  .r
}

.ssRef <- function(par, ncmt, oral0, rate, type, tau = 0, tinf = 0, amt = 0, bolusCmt = 0L) {
  .km <- .lcRateMatrix(
    ncmt,
    oral0,
    par[["cl"]],
    par[["v"]],
    par[["q"]],
    par[["vp"]],
    par[["q2"]],
    par[["vp2"]],
    par[["ka"]]
  )
  .n <- nrow(.km)
  .eye <- diag(.n)
  .r <- .ssRateFull(ncmt, oral0, rate)
  if (type == 1L) {
    return(drop(-solve(.km, .r)))
  }
  if (type == 2L) {
    .inf <- solve(.km, (.lcExpm(.km, tinf) - .eye) %*% .r)
    return(drop(solve(.eye - .lcExpm(.km, tau), .lcExpm(.km, tau - tinf) %*% .inf)))
  }
  .d <- numeric(.n)
  .d[bolusCmt + 1L] <- amt
  drop(solve(.eye - .lcExpm(.km, tau), .d))
}

.ssCases <- function(ncmt, oral0) {
  .rates <- if (oral0 == 1L) list(c(7, 0), c(0, 7)) else list(7)
  .bolus <- if (oral0 == 1L) 0:1 else 0L
  c(
    lapply(.rates, function(.r) list(type = 1L, rate = .r)),
    lapply(.rates, function(.r) list(type = 2L, rate = .r, tau = 12, tinf = 2)),
    lapply(.bolus, function(.b) list(type = 3L, rate = numeric(1 + oral0), tau = 12, amt = 100, bolusCmt = .b))
  )
}

.ssArgs <- function(cs) {
  list(
    rate = cs$rate,
    type = cs$type,
    tau = if (is.null(cs$tau)) 0 else cs$tau,
    tinf = if (is.null(cs$tinf)) 0 else cs$tinf,
    amt = if (is.null(cs$amt)) 0 else cs$amt,
    bolusCmt = if (is.null(cs$bolusCmt)) 0L else cs$bolusCmt
  )
}

for (.ncmt in 1:3) {
  for (.oral0 in 0:1) {
    test_that(sprintf("steady-state amounts match the dosing-cycle fixed point (ncmt %d, oral %d)", .ncmt, .oral0), {
      for (.cs in .ssCases(.ncmt, .oral0)) {
        .a <- .ssArgs(.cs)
        .ref <- do.call(.ssRef, c(list(.ssPar, .ncmt, .oral0), .a))
        .r <- do.call(.ssCall, c(list(.ssPar, .ncmt, .oral0), .a))
        .lab <- sprintf("type %d rate %s bolusCmt %d", .a$type, paste(.a$rate, collapse = ","), .a$bolusCmt)
        expect_equal(.r$Alast[seq_len(.ncmt + .oral0)], .ref, tolerance = 1e-8, label = .lab)
        expect_equal(.r$val, .ref[.oral0 + 1L] / .ssPar[["v"]], tolerance = 1e-8, label = .lab)
      }
    })

    test_that(sprintf("steady-state gradients match finite differences (ncmt %d, oral %d)", .ncmt, .oral0), {
      .idx <- c(seq_len(2L * .ncmt), if (.oral0 == 1L) 7L)
      for (.cs in .ssCases(.ncmt, .oral0)) {
        .a <- .ssArgs(.cs)
        .r <- do.call(.ssCall, c(list(.ssPar, .ncmt, .oral0), .a, list(deriv = TRUE)))
        expect_length(.r$Jg, length(.idx))
        .fd <- vapply(
          .idx,
          function(.j) {
            .h <- 1e-6 * .ssPar[[.j]]
            .up <- .ssPar
            .up[.j] <- .up[.j] + .h
            .dn <- .ssPar
            .dn[.j] <- .dn[.j] - .h
            .vUp <- do.call(.ssCall, c(list(.up, .ncmt, .oral0), .a))$val
            .vDn <- do.call(.ssCall, c(list(.dn, .ncmt, .oral0), .a))$val
            (.vUp - .vDn) / (2 * .h)
          },
          numeric(1)
        )
        expect_equal(.r$Jg, .fd, tolerance = 1e-5, label = sprintf("type %d", .a$type))
      }
    })
  }
}

test_that("steady state with no dose into a valid compartment is NA", {
  # oral with no rate, and a bolus into a compartment past the depot/central
  for (.ncmt in 1:3) {
    expect_true(is.na(.ssCall(.ssPar, .ncmt, 1L, c(0, 0), 1L)$val))
    expect_true(is.na(.ssCall(.ssPar, .ncmt, 1L, c(0, 0), 2L, tau = 12, tinf = 2)$val))
    expect_true(is.na(.ssCall(.ssPar, .ncmt, 1L, c(0, 0), 3L, tau = 12, amt = 100, bolusCmt = 2L)$val))
    expect_true(is.na(.ssCall(.ssPar, .ncmt, 0L, 0, 3L, tau = 12, amt = 100, bolusCmt = 1L)$val))
  }
})

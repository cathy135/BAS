#' Compound Confluent hypergeometric function of two variables
#'
#' Compute the Confluent Hypergeometric function of two variables, also know as
#' a Horn hypergeometric function or Humbert's hypergeometric used in Gordy
#' (1998) with integral representation:
#'
#'  phi_1(a,b,c,x,y) = Beta(a,b) Int_0^1
#' t^(a-1) (1 - t)^(c-a-1) (1 - yt)^(-b) exp(x t) dt
#' \url{https://en.wikipedia.org/wiki/Humbert_series} Note that Gordy's
#' arguments for x and y are reversed in the reference above.
#'
#' Code for phi1 provided by Gordy.
#'
#'
#' @param a a > 0
#' @param b arbitrary
#' @param c c > 0
#' @param x x > 0
#' @param y 0 <= y < 1
#' @author Merlise Clyde (\email{clyde@@stat.duke.edu})
#' @references Gordy 1998
#' @keywords math
#' @examples
#'
#' # special cases
#' # phi1(a, b, c, x=0, y) is the same as 2F1(b, a; c, y)
#' phi1(1, 2, 1.5, 0, 1 / 100)
#' hypergeometric2F1(2, 1, 1.5, 1 / 100, log = FALSE)
#'
#' # phi1(a,0,c,x,y) is the same as 1F1(a,c,x)
#' phi1(1, 0, 1.5, 3, 1 / 100)
#' hypergeometric1F1(1, 1.5, 3, log = FALSE)
#' @rdname mpphi1
#' @family special functions
#' @export
#'
#'

gctorture(on = TRUE)

mpphi1 <- function(a, b, c, x, y) {
  # phi1 = int u^{t-1} (1 - v u)^{q - 1} e^-{s u} /B(t, q) (theta
  na <- length(a)
  nb <- length(b)
  nc <- length(c)
  nx <- length(x)
  ny <- length(y)
  
  ns = c(na,nb, nc, nx, ny)
  n = max(ns)
  
  if ((n > 1) && (mean(ns) != n)) {
    stop("length of inputs are not the same")
  }
  out <- rep(0, n)
  ans <- .C(C_mpphi1,
            as.numeric(a),
            as.numeric(b),
            as.numeric(c),
            as.numeric(x),
            as.numeric(y),
            out = as.numeric(out), as.integer(n))$out
  return(ans)
}
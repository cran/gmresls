#' Solve a Least Squares System with a Preconditioner.
#'
#' Solve a least squares system Ax~=b (dim(A)=(m,n) with m >= n) with a
#' preconditioner B: BAx=Bb (dim(B)=(n,m)).
#' Implemented method uses GMRES(k) with callback functions, i.e. no explicit A or B are required. GMRES can be restarted after k iterations.
#'
#' @param f_resid A function f_resid(x, ...) calculating B(b-Ax) for a given x. If x is of length 0 (e.g. NULL), it must be considered as 0.
#' @param f_BAx A function f_BAx(x, ...) calculating matrix-matrix-vector product BAx for a given x.
#' @param x0 A vector or NULL (which means 0), initial approximation for Ax=b
#' @param k An integer, parameter for restarting GMRES. Value 0 (default) means no restart, i.e. at most length(x) basis vectors will be constructed and used.
#' @param maxit A maximal iteration number. Here, itereation number continues to increment even after a possible GMRES restart. Default (0) means length(x).
#' @param tol A tolerance for solution x, estimated as ||B(Ax-b)||/||Bb||, default 1.e-7
#' @param ... Parameters passed through to f_BAx and f_resid
#' @returns The solution x, having the structure of f_resid(x,...).
#' @details
#' Implemented method is equivalent to a classical GMRES(k) method with restart after constructing k basis vectors and applied to a square system BAx=Bb.
#' Dense matrices constructed and stored by this method are of size (length(x), k) and (k+1, k) where k is GMRES current basis vector number. If maxit > k, GMRES will be restarted after each k iterations
#' Particularity of this implementation that matrices A and B have no to be stored explicitly.
#' User provides just callback function mimicking their multiplication by adequate vectors.
#' In case of non convergence after maxit iterations, attr(x) will contain a field 'warning' with the message which will be also issued with warning()
#' If the operator BA is not of full rank, iterations will be stopped before reaching convergence or maxit. A warning will be emitted in this case.
#'
#' @examples
#' # prepare a 4x3 toy least squares (LS) problem Ax=b
#' A=rbind(diag(1:3)+matrix(1, 3,3), rep(1, 3))
#' xsol=1:3
#' b=A%*%xsol+rnorm(4, 0., 0.1) # add some noise as it is often the case in LS
#' f_resid=function(x,...)
#'    with(list(...), if (length(x) == 0L) crossprod(A, b) else crossprod(A, b-A%*%x))
#' f_BAx=function(x,...)
#'    with(list(...), crossprod(A, A%*%x))
#' x=gmresls(f_resid, f_BAx, A=A, b=b)
#' stopifnot(all.equal(c(x), c(qr.solve(A,b))))
#' @export
gmresls=function(f_resid, f_BAx, x0=NULL, k=0, maxit=0, tol=1.e-7, ...) {
  it=0L
  warn=""
  repeat { # restart till convergence or maxit reached
    it=it+1L
    r0=f_resid(x0, ...)
    nr0=sqrt(sum(r0*r0))
    if (it == 1L) {
      nx=length(r0)
      if (maxit == 0L)
        maxit=nx
      nbb=if (is.null(x0)) nr0 else sqrt(crossprod(c(f_resid(NULL, ...))))
    }
    
    v=as.matrix(c(r0))*(1./nr0)
    h=matrix(0., 1L, 0L) # Heisenberg matrix, will be growing
    re1=nr0
    if (k == 0)
      k=nx
    converged=nr0/nbb <= tol
    for (i in seq_len(k)) {
      if (converged)
        return(if (is.null(x0)) double(nx) else x0)
      bavi=f_BAx(v[,i], ...)
      h=cbind(h, crossprod(v, bavi))
      v=cbind(v, bavi-v%*%h[,i])
      h=rbind(h, 0.)
      h[i+1L,i]=sqrt(crossprod(v[,i+1]))
      if (h[i+1L,i] <= tol) {
        re1=c(re1, 0.)
        y=qr.solve(h, re1)
        nri=sqrt(crossprod(h%*%y-re1))
        converged=nri/nr0 <= tol
        break
      } else {
        v[,i+1L]=v[,i+1L]/h[i+1L,i]
      }
      re1=c(re1, 0.)
      y=qr.solve(h, re1)
      nri=sqrt(crossprod(h%*%y-re1))
      converged=nri/nr0 <= tol
      it=it+1L
    }
    x=v[,-ncol(v),drop=FALSE]%*%y
    if( !is.null(x0) )
      x=x+x0
    # go to restart
    if (converged || it >= maxit)
      break
  }
  if (!converged)
    warn=paste0(warn, "Maximal iteration number (", maxit, ") is reached without convergence")
  if (nchar(warn))
    warning(attr(x, "warning") <- warn)
  attr(x, "norm_res")=nri
  attr(x, "iteration")=it
  attr(x, "converged")=converged
  x
}

# doc example
# prepare a 4x3 toy problem Ax=b
A=rbind(diag(1:3)+matrix(1, 3,3), rep(1, 3))
xsol=1:3
b=A%*%xsol+rnorm(4, 0., 0.1)
f_resid=function(x,...) with(list(...), if (length(x) == 0) crossprod(A, b) else crossprod(A, b-A%*%x))
f_BAx=function(x,...) with(list(...), crossprod(A, A%*%x))
x=gmresls(f_resid, f_BAx, A=A, b=b)
test.doc <- function()
  checkEqualsNumeric(x, qr.solve(A, b))

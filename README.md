# gmresls

The goal of gmresls is to solve a least squares problem Ax\~=b for which the matrix A and may be vector b are not explicitly known. We suppose that it exists a precondition operator B such that BA is somewhat close to I (identity matrix of appropriate size). This operator B does not have to be explicit either. To simulate A, B and b actions, user have to supply callback functions f_resid() and f_BAx(). The algorithm is based of GMRES (Generalized minimal residual)

## Installation

You can install the current version of gmresls like so:

``` r
intall.package("gmres") # from CRAN
# or dev version
devtools::install_git(git@forgemia.inra.fr:mathscell/gmresls.git)
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
# prepare a 4x3 toy problem Ax=b
A=rbind(diag(1:3)+matrix(1, 3,3), rep(1, 3))
xsol=1:3
b=A%*%xsol+rnorm(4, 0., 0.1)
f_resid=function(x,...) with(list(...), if (length(x) == 0) crossprod(A, b) else crossprod(A, b-A%*%x))
f_BAx=function(x,...) with(list(...), crossprod(A, A%*%x))
x=gmresls(f_resid, f_BAx, A=A, b=b)
stopifnot(all.equal(c(x), c(qr.solve(A, b))))
```

### Legal information

Author: Serguei Sokol (INRAE/TBI/Mathematics cell)

Copyrights 2024, INRAE/INSA/CNRS

License: GPL (>=3)

Issue reporting: https://forgemia.inra.fr/mathscell/gmresls/-/issues


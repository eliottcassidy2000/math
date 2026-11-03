# Four, six and ten in a minimal ternary carry representation

**PROVED elementary representation theorem.** Local labels M1--M3.
No literature-priority claim is made for truncated group-algebra methods.
The new transport is to the explicit radical matrix from the previous
tournament packet and to the actual hostile pair below.
For standard modular Jordan-block background, see
[Dave Benson's 2012 lecture notes, section 1](https://sites.math.washington.edu/~pischool/Notes/Benson-Notes.pdf).
The affine extension, minimality argument and particular transports used
here are supplied explicitly below; that reference is not a priority
claim for these applications.

## M1 A minimal linear realization at every ternary precision

Let N=3^k, k>=1, and let G_k consist of affine maps `x -> a x+b` on
Z/N, where a is a unit. Translation by 4 and multiplication by 2 generate
G_k. Translation has order N; doubling has order `2*3^(k-1)`.

The latter statement follows from `v3(4^(3^j)-1)=j+1`: start with
4-1=3 and apply `(1+3^r c)^3-1` for a 3-unit c. An exponent giving 1
for 2 must first be even modulo 3, so these are the exact orders.
Thus 2 generates all units of Z/3^k and
`|G_k|=2*3^(2k-1)`.

Over F3 the smallest faithful linear representation of G_k has dimension

`r_k=3^(k-1)+1`.                                         (3)

Construction: work in `R_k=F3[t]/(t^r_k)` and put u=1+t. Translation by
b is multiplication by u^b; multiplication of addresses by a is the
substitution `u -> u^a`. This substitution is well-defined and invertible
because `(1+t)^a-1` is t times a unit when 3 does not divide a.
Since `u^(3^j)=1+t^(3^j)`, u has exact order N. Distinct addresses give
distinct states u^x. The operators satisfy precisely the affine action
on this orbit, proving faithfulness of the whole group.

Lower bound: an element of order 3^k in GL_d(F3) is unipotent. Write it
as I+J. If d<=3^(k-1), then J^(3^(k-1))=0 and its order divides
3^(k-1), a contradiction. Therefore d>=r_k, and the construction attains
the bound. It also preserves the affine hyperplane of polynomials with
constant coefficient 1. Thus the minimum affine ambient dimension is
r_k-1: homogeneous coordinates give the same lower bound, and this
hyperplane attains it.

The minimal linear dimensions are

`2,4,10,28,82,...`, with `r_(k+1)=3 r_k-2`.

In particular, **4+6=10** is the next precision jump, and 6 is exactly
the doubling clock at modulus 9. More generally, every precision jump
`r_(k+1)-r_k=2*3^(k-1)` equals the doubling order at the old precision.
These even dimensions encode a ternary
clock; they are not numbers of vertices in a claimed tournament.

## M2 The literal nine-cycle in four coordinates

At k=2, an address x has coordinates

`u^x = (1, x, binom(x,2), binom(x,3)) mod 3`.

If x=3q+r, 0<=q,r<3, these are
`(1,r,r(r-1)/2,q)`. Expanding `(1+t)^(3q+r)` proves this directly.
The quadratic coordinate restores the carry that the prior two-digit
affine model could not retain. Translation by 4 gives exactly
`4,8,3,7,2,6,10,5,9,4` on the representatives 2,...,10.

Doubling substitutes `t -> 2t+t^2`. Its row-action matrix is

```text
1 0 0 0
0 2 1 0
0 0 1 1
0 0 0 2
```

It has order six. Together with multiplication by u^4 it gives the
54-element affine group. Halving in the residue chart is its inverse
doubling operation; the legality and stopping clock of integer halving
still need D5's separate data.

## M3 The previous 27-coordinate radical contains this clock

The explicit matrix M in [A5](../return701_algebra_20261103/ALGEBRA95.md)
has M^4=0 and M^3 nonzero modulo 3. Choose the first standard row v with
vM^3 nonzero, and use the basis `v,vM,vM^2,vM^3`. These rows are independent:
multiply a putative relation by a suitable power of M, starting from its
first nonzero coefficient. The map

`f(t) -> v f(M)`

identifies R_2 with this four-dimensional submodule. It intertwines
multiplication by u with I+M and hence carries the displayed nine-cycle
to nine explicit vectors in the full 27-coordinate space. Substitution
by u^2 transfers as a linear operation on this submodule, with the above
matrix in the named chain basis.

This is a formal module transport. It does not say those vectors are
actual tournament source states, or that the substitution extends to a
legal attachment on the entire observer. The previous singular-generator
obstruction remains: no nonempty positive attachment word realizes the
invertible nonidentity clock I+M. The explicit vectors and frame in
[clock.json](clock.json) keep this distinction testable.

# The wrapped incidence kernel over every coefficient ring

2026-10-25. **PROVED**: a universal-ring kernel decomposition and the complete
integer torsion invariant. **FINITE-EXACT** direct Smith-form controls are
separate. This extends, rather than corrects, Session62's all-prime field
theorem in `incidence/cyclic_common_kernel.md`.

## 1. Objects and the integral answer

Let N>=5, gcd(N,3)=1, and let I be a nonempty subset of Z/N. At each
(a,b,j) in {0,1}^2 x I put one coefficient. Its three incidence labels are

    horizontal j-9b, vertical -a-3b, sloped j+3a, modulo N.

Let M be the 3N-by-4|I| zero-one integer matrix recording the sums on these
three partitions, with empty rows retained. For prime N this is exactly
the inherited physical three-direction matrix, up to a permutation of
the sloped labels by the unit 3. For composite N these are explicitly the
three cyclic partitions; no affine-plane-over-a-field claim is made.

Set G=(Z/N)\I, H=I intersection (I-9). On all N residues draw the edges
{j-6,j-3}, j in G, and anchor every vertex outside H at zero. Since step
3 is an N-cycle and I is nonempty, this graph is a disjoint union of paths.
Let k be the number of components without anchors, o indicate that at
least one such component has odd size, and b indicate that in every
anchored component all anchors are in one bipartition class.

Put d=|I|-1+k-o. Define an integer D as follows:

    b=0:          D=1;
    b=1, o=1:     D=N;
    b=1, o=0:     D=gcd(N,beta).                              (W1)

In the last line beta is the number of ones in an alternating 0/1 assignment
on the paths, chosen zero at every anchor. The choice on a free even path
does not change its number of ones, so beta is well defined in this case.

For **every commutative ring R**, with R[D]={c in R:Dc=0},

    ker(M:R^(4|I|)->R^(3N)) ~= R^d direct_sum R[D].             (W2)

The isomorphism is constructive after fixing path orientations and one
odd free component when present. In particular, the nonzero integer Smith
factors are all 1 except a single D when D>1. The integer rank is 4|I|-d,
and

    coker_Z M ~= Z^(3N-4|I|+d) direct_sum Z/D.                 (W3)

For a prime p, W1 reduces to D=p exactly when b=o=1, and D=1 otherwise.
This proves the full integer statement previously checked only in finite
Smith-form universes. Over a field the inherited extra dimension is exactly
the reduction of this one p-torsion factor.

## 2. Ring-level reduction, with no forbidden division

Extend every column by zero off I, and denote the four columns A,B,C,Dcol.
Write tau_s r(j)=r(j-s) and Q=tau_3+tau_6+tau_9. Solving the horizontal and
sloped equations as in the inherited proof gives, uniquely,

    A=h, B=-Qr+c-h, C=r-h, Dcol=-tau_9 r+Qr-c+h.              (W4)

This uses only that the kernel of 1-tau_3 is the constants, true over any
ring because step 3 is one cycle. Support and the four vertical sums are
equivalent to

    supp h subset I, sum h=0;
    supp r subset H, sum r=0;
    (Qr)_j=c for j outside I; N c=0.                         (W5)

No division by N or by 3 in R has occurred. The geometric relabeling by
3 was modulo N, a different operation. On a gap j, r_(j-9)=0, so the
boundary equation is r_(j-6)+r_(j-3)=c.

The h-space is freely R^(|I|-1). On a path the r-values alternate t,c-t.
Two opposite-parity anchors force c=0 with coefficient one. Thus if b=0,
every anchored component vanishes, every free component contributes its
alternating homogeneous vector, and sum r removes one free parameter
exactly when an odd free path exists. This gives R^d.

Suppose b=1. Choose the 0/1 particular solution r0 at c=1. Write beta=sum r0.
Every solution of the boundary equations is

    r=c r0 + sum_(free paths C) t_C psi_C,

where psi_C alternates +1,-1 and has sum zero or one according as C has
even or odd size (orient an odd path to start with +1). If o=1, solve
the sum equation for the parameter on one chosen odd free path. The
remaining parameters are free and c is subject only to Nc=0. If o=0,
all homogeneous sums vanish: c must satisfy both Nc=0 and beta c=0.
Bezout makes this equivalent over any ring to gcd(N,beta)c=0. This proves
W2, including characteristic two and rings with nilpotents.

When N=p is prime and o=0,b=1, I is proper (full I has odd isolated free
paths). There is at least one edge, which has a one in the particular
assignment; every one lies in H. Hence 1<=beta<=|H|<=p-1, so D=1.
This is the precise place where primality removes an additional torsion
possibility. Composite cyclic lengths need not remove it.

## 3. From kernels to exact Smith factors

Over Q, W2 has dimension d, giving the stated integer rank. For every prime
ell and exponent e>=1, W2 over Z/ell^e has cardinality

    ell^(e d + min(e,v_ell(D))).                             (W6)

On the other hand an integer Smith form with nonzero factors s_i gives
ell^(e d + sum_i min(e,v_ell(s_i))). Comparing for all e says that each
ell-primary torsion part is cyclic of exact order ell^v_ell(D). At most
one invariant factor can be nonunit: two nonunit factors in a divisibility
chain would share a prime and violate the e=1 comparison. Its value is D.
This proves W3 without assuming that field dimensions detect exponents.

## 4. Exact hostiles and stronger survivors

* N=8, I={0,1,2,5}: b=1,o=0,beta=2, D=2. There is a characteristic-two
  mode although there is no odd free path. The prime-only bo criterion
  must not be extrapolated to composite length.
* N=16, I={0,1,2,4,5,7,8,9}: b=1,o=0,beta=4, D=4. Counting just one
  extra dimension in characteristic two loses its higher exponent.
* Full I at N=49: k=49,o=b=1,d=96,D=49. Over F_7 there is one extra
  kernel dimension, just as torsion 7 would give, but over Z/49 the
  extra kernel factor has 49 elements, not 7.
* At prime p>=17 and interval length L, the inherited formula for d is
  (L-1)+(L-10)_++(L-p+6)_++(L-p+3)_+; D=p exactly when L>=p-6.
  The prime-characteristic term and the wrap terms are thereby retained.

The full integer matrix has a stable torsion invariant under unimodular
row and column changes. This differs from the optional odd determinant of
a chosen Boolean character lift in `KUMMER_AND_F1.md`: changing that lift
with radicand corrections changes the ambient monomial map while preserving
its marked quadratic-root restriction. One cannot remove physical incidence
torsion merely by renaming its coordinates.

## 5. Boundaries and next question

The theorem covers arbitrary puncturing, every coefficient ring, and cyclic
lengths prime to three. It gives no binary-error realization, planar unit
embedding, colour forcing, or Collatz transport estimate. For 3|N the step-3
permutation has multiple cycles and the physical four-column representation
may have additional issues at small N. Separate constants per cycle and
their coupled boundary equations would be the next faithful representation;
the single D formula is not claimed there.

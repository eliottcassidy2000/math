# Independent audit of the fixed Keller level-four boundary

**Verdict: ACCEPTED and promoted as THM-3498.**

The auditee was remote commit `95e0c4c7e` on
`codex/session-turnpike-atlas-20260815-codex2`.  The audit ran in a fresh
worktree based on current `origin/main`; the candidate commit was imported
only after the bounded startup and inheritance pass.

## Claim decomposition

The candidate had three logically separate load-bearing claims:

1. THM-3495's polynomial `J` has exposed face of weight `43`, the two
   divergent inverse sheets each contribute `-43/2`, and the finite sheet is
   a unit, so `v_L(N(J))=-43`.
2. One lawful `F_101` fibre gives a squarefree norm-product of exact degree
   `81`, so the characteristic-zero fourth eliminant is generically
   squarefree and its `27` cubic blocks are pairwise coprime.
3. With constants and signs retained, the odd-block recursion gives
   `[Delta_4]=[2G]` for `G=L^43N(J)`.

Failure of any one would have left THM-3498 reserved.

## Independent old-boundary route

The audit rebuilt the already-canonical THM-3495 `J`, converted its
coefficient dictionary into a new ordered representation, discovered the
maximizing weight `i-k`, and sent only the exposed face to SymPy's independent
factorizer.  It recovered

```text
max(i-k)=43,
#face=16,
in(J)=-2^58 3^51 13^8 79^4 313^2 x^43(3xz-2y)^15.
```

Under `(x,y,z)=(u^-1,d,-3du)+higher`, the factor `3xz-2y` becomes
`-11d+O(u)`, so no equal-weight cancellation occurs.  The two divergent
branches are factors in a norm product, not summands; their valuations add
and cannot cancel across sheets.

The finite point `(2,5/6,-7/8)` was evaluated by nested Horner reduction,
not the submitted termwise summation.  Its exact rational value is nonzero
and has ledger hash

```text
f4ad5e424ac9abc719f21fb92625e7c25a9e37cc614815dd1b54e43fa5936554.
```

Thus the valuation is exactly `-43`, not merely at most `-43`.

## Independent degree-81 route

The submitted program encodes nested cubic elements as coefficient triples,
uses Python Gaussian elimination for every norm, and reconstructs the answer
from `82` consecutive values by Newton interpolation.

The audit shares none of those implementation layers.  An element of each
nested algebra is represented directly by its FLINT multiplication matrix.
Adjoining `theta^3=p theta+q` uses the block companion matrix

```text
[ 0  0  q ]
[ 1  0  p ]
[ 0  1  0 ].
```

The fourth norm is a `27 x 27` FLINT determinant.  It is evaluated at every
nonzero element of `F_101`; multiplicative Fourier inversion recovers all
coefficients below degree `100`.  Degrees `82,...,99` vanish, degree `81` is
nonzero, and `X=0`—not used by the Fourier inversion—is a held-out
determinant control.  FLINT independently gives derivative gcd `1`.  The
ascending ledger hash is exactly

```text
1c05c0fd5ee48fc2dd030aebdb9ad6ddd8185fb933eb91e7e39ff553424ef5a7.
```

All leading coefficients, cubic derivatives, and chart denominators used by
the three-level index algebra are units.  Hence this is good reduction of a
finite-etale `27`-point algebra.  Over an algebraic closure its norm splits as
the product of `27` cubic blocks; a squarefree degree-`81` product makes every
block cubic and separable and every pair coprime.

## Localization and constant-unit audit

Finite étaleness over `A[L^-1]`, `A=Q[a,b,c]`, puts `N(J)` in that localized
UFD.  A reduced element of `A[L^-1]` has only an `L`-power denominator.
Therefore valuation `-43` is equivalent to

```text
G=L^43N(J) in A,     v_L(G)=0,     gcd(G,L)=1.
```

This does not imply that `G` is primitive, squarefree, irreducible, or the
equation of an image surface.

The full sign and power-of-two chain is

```text
[Delta_3]=[-L*J/(2^35L^7)]=[-2J],
[Delta_4]=[N(-2J)][-L]
         =[(-2)^3N(J)][-L]
         =[8LN(J)]
         =[8G/L^42]
         =[2G].
```

The first two minus signs cancel; `L^42` and `8/2=4` are squares.  This is
the exact MISTAKE-413 control: `[2]` cannot be discarded as an unnamed
constant unit.

## Scope boundary

The theorem closes a boundary valuation and discriminant genericity gate for
one fixed three-dimensional Keller map.  The global construction and
factorization of `G`, the degree of `V(J)->V(G)`, coprimality with `H,J`, the
four-component nonproperness claim, any all-level induction, and every
general `JC`, `DC`, or LRC consequence remain open.

---
id: THM-3770
title: "Vertical principal-part equalizer and log-canonical dressing gate"
status: >
  PROVED + PENDING INDEPENDENT AUDIT.  For a smooth planar polynomial Q whose
  rational Hamiltonian constants are exactly k(Q), a rational
  constant-Jacobian mate can be made polynomial if and only if it has no
  horizontal pole and its successive Laurent principal coefficients agree as
  scalars on every component over each target value.  This gives a finite,
  necessary-and-sufficient vertical equalizer algorithm.  For a squarefree
  log-canonical pair J(U,W)=mU, equality of the component values of W is
  equivalent to an actual Keller mate (W-c)/U.  Multiplication by any
  nonconstant squarefree phi(W) cannot repair a smooth dressed target
  Q=U phi(W): smoothness separates the old component values from the roots of
  phi, while polynomial regularity would have to identify them.  This is a
  construction/no-go criterion, not a planar Jacobian counterexample.
source: root / planar-Jacobian vertical-equalizer session, 2026-08-23
audit: >
  SELF-AUDITED PROOF CANDIDATE.  The DVR recursion, horizontal-divisor
  boundary, partial-fraction assembly, Krull-intersection step, component
  constant lemma, birational constant-field calculation, and dressed-fibre
  mismatch have been derived independently by valuation and differential
  arguments.  Independent hostile audit remains due.
depends_on:
  - THM-2230-planar-jacobian-response-fiber-and-exact-target-shear-quotient
related:
  - THM-3598-danielewski-rational-exact-polar-graph-family-and-classification
  - THM-3758-quadratic-radial-carrier-rational-exact-split-fibre-nonentry
---

# THM-3770 -- vertical equalization is the global regularity gate

**PROVED + PENDING INDEPENDENT AUDIT.**  Once the generic Hamiltonian time
form has a rational primitive, finding a polynomial mate is no longer another
coefficient PDE.  It is a finite gluing problem over the exceptional target
values.  The exact sidecar is the principal part on every irreducible
component, not merely its target value.

Throughout, let `k` be an algebraically closed field of characteristic zero,
put `A=k[x,y]`, `K=k(x,y)`, and use

```text
J(F,G)=F_x G_y-F_y G_x.                              (1)
```

## 1. The vertical principal-part equalizer

Let `Q in A` be nonconstant and have no critical point.  Suppose that for
some `P0 in K` and `c in k*`,

```text
J(P0,Q)=c,             ker_K J(-,Q)=k(Q).             (2)
```

Then all rational mates with response `c` are exactly

```text
P0+H(Q),                         H in k(T).            (3)
```

There is an `H in k(T)` for which `P0+H(Q)` belongs to `A` if and only if
the following two conditions hold.

1. **No horizontal pole.**  The function `P0` has nonnegative valuation at
   every irreducible divisor `D` on which `Q` is nonconstant.
2. **Every vertical principal part equalizes.**  Fix `lambda in k`, write

   ```text
   Q-lambda=prod_i f_(lambda,i),                      (4)
   ```

   and let `D_(lambda,i)` be the divisor `f_(lambda,i)=0`.  Put
   `q=Q-lambda` and

   ```text
   R=max_i max(0,-v_(D_(lambda,i))(P0)).              (5)
   ```

   Starting with `F_R=P0`, run `r=R,R-1,...,1`.  At stage `r`, the induction
   hypothesis is `v_D(F_r)>=-r` on every component over `lambda`.  Require
   the residues

   ```text
   alpha_(i,r)=res_(D_(lambda,i))(q^r F_r) in k(D_(lambda,i))             (6)
   ```

   to be one common scalar `alpha_r in k`, including the zero residue on a
   component where the pole order is smaller.  Then set

   ```text
   F_(r-1)=F_r-alpha_r q^(-r).                        (7)
   ```

   The requirement is that this recursion succeeds for every `lambda` at
   which `P0` has a vertical pole.

When it succeeds, one regularizing target shear is

```text
H(T)=-sum_lambda sum_(r=1)^R alpha_(lambda,r)(T-lambda)^(-r).             (8)
```

The sums are finite.  Conversely every regularizing `H` has exactly these
principal parts.  Once one polynomial mate `P` exists, all polynomial mates
with response `c` are

```text
P+k[Q].                                                (9)
```

### Proof

The first statement follows immediately from `(2)`: two rational solutions
of the same response differ by a rational Hamiltonian constant.

Because `Q` has no critical point, every `Q-lambda` is squarefree.  Thus `q`
has valuation one at each divisor in `(4)` and is a fixed uniformizer there.
At stage `r`, subtracting a scalar multiple of `q^(-r)` raises every component
valuation past `-r` if and only if the residues `(6)` are that same scalar.
This proves, by descending induction, that `(6)--(7)` are necessary and
sufficient for one target-only principal part to remove all poles above
`lambda`.  It also proves uniqueness of that principal part.

A rational function `H(Q)` has no pole on a horizontal divisor: the residue
of `Q` there is transcendental over `k`, so no nonzero one-variable
denominator vanishes identically.  Hence a horizontal pole of `P0` cannot be
cancelled.  At distinct vertical values, the summands in `(8)` are units in
each other's local rings, so the independently constructed principal parts
do not interfere.  After `(8)`, the resulting rational function has no pole
at any height-one prime of `A`.  The UFD `A` is normal and equals the
intersection of its height-one DVRs inside `K`; therefore the result lies in
`A`.

Conversely, the Laurent principal part of any `H in k(T)` at `lambda` has
scalar coefficients and is identical on every component of `(4)`, proving
necessity.  Finally, if two regular mates differ by `R(Q)` with `R in k(T)`,
then `R` can have no finite pole: every nonconstant `Q-lambda` has a zero in
`A^2`.  Hence `R in k[T]`, which proves `(9)`.  **QED.**

## 2. A log-canonical component spectrum

Let `U,W in A`, let `m in k*`, and suppose

```text
J(U,W)=mU,                         U squarefree.       (10)
```

Factor `U=prod_i U_i` into distinct irreducibles.  There are unique constants
`c_i in k` such that

```text
W=c_i mod U_i.                                         (11)
```

Indeed, reducing `(10)` modulo `U_i` gives
`J(U_i,W)=0 mod U_i`.  The induced tangent derivation on the function field
of the irreducible curve `U_i=0` is nonzero.  Its constant field is algebraic
over `k`, hence equals `k`; this proves `(11)`.  Call the indexed collection

```text
Spec_0(U,W)=(c_i)_i                                   (12)
```

the **zero-component spectrum** of the log-canonical pair.

The spectrum is equalized if and only if `(10)` already hides a polynomial
Keller mate.  More precisely, the following are equivalent:

```text
(a) all c_i equal one c;
(b) W-c=UV_0 for some V_0 in A;
(c) J(U,V)=m for some V in A.                         (13)
```

The first equivalence uses squarefreeness: every `U_i` divides `W-c`, so
their product does.  Substitution into `(10)` gives

```text
J(U,c+UV)=U J(U,V)=mU,                                (14)
```

and cancellation in the domain `A` proves `(b)=>(c)`.  For `(c)=>(a)`,
THM-2230's polynomial centralizer theorem applied to the mate `V` gives

```text
W-UV=H(U),                         H in k[T].          (15)
```

Reducing on every `U_i` shows `c_i=H(0)`.  Equivalently, subtracting `H(0)`
from `(15)` and dividing `H(U)-H(0)` by `U` also recovers `(b)`, possibly
with a different mate `V_0`.  In particular, an irreducible `U` in a
log-canonical pair automatically has a polynomial constant-Jacobian mate.
The reducible case isolates the missing datum: its component constants must
synchronize.

This is a useful reverse formulation of the planar search.  Instead of
solving `J(U,V)=m` directly, one may build Darboux eigenpairs `(10)` and ask
whether the finite spectrum `(12)` collapses to one value.  A nonsynchronized
pair is a rational near miss; a synchronized noncoordinate `U` would be an
actual counterexample.

## 3. Nonconstant spectral dressing cannot synchronize

There is also a general trap.  In addition to `(10)`, assume

```text
K=k(U,W),                                              (16)
```

so `(U,W)` is a birational rational coordinate system.  Let `phi in k[T]`
be nonconstant and squarefree, and put

```text
Q=U phi(W),                 P0=-W/(mQ).               (17)
```

Then

```text
J(P0,Q)=1.                                             (18)
```

If `Q` has no critical point, it has no polynomial constant-Jacobian mate.

To prove `(18)`, observe that

```text
J(Q,W)=phi(W)J(U,W)=mQ.                               (19)
```

Moreover `(16)` and `U=Q/phi(W)` give `K=k(Q,W)`, and `(19)` says that the
Hamiltonian derivation is a nonzero `k(Q)`-multiple of `d/dW`.  Therefore

```text
ker_K J(-,Q)=k(Q),                                    (20)
```

so `(3)` is the complete rational-mate torsor.

Let `r` be any root of `phi`.  The zero fibre of `Q` contains the old
components `U_i=0` and components of `W-r=0`.  Smoothness forces

```text
phi(c_i)!=0,                         hence c_i!=r:      (21)
```

otherwise `U_i` divides both factors in `(17)`, making the zero fibre
nonreduced and `Q` critical.  Since `Q` is a uniformizer, the simple `1/Q`
principal coefficients of `P0` are

```text
U_i=0:       -c_i/m,
W-r=0:       -r/m.                                    (22)
```

They differ by `(21)`.  The first step of the vertical equalizer therefore
fails, and no member of `P0+k(Q)` is polynomial.  **QED.**

The contradiction is structural: smoothness requires the newly appended
root spectrum of `phi` to avoid the old component spectrum, while regularity
requires all of those values to agree.  Dressing creates extra zero-fibre
components but cannot pay their gluing debt.

## 4. Construction consequences and boundaries

The criterion separates three invoices which finite coefficient solving
usually blends together:

```text
generic exactness   rational primitive P0 exists;
horizontal debt     P0 has no pole dominating the target line;
vertical gluing     every componentwise principal coefficient equalizes. (23)
```

THM-3598 reaches the first line and fails on arm/residual divisors.
THM-3758 reaches the first two and fails by opposite coefficients on two
components of `Q=0`.  Equation `(13)` shows exactly what a positive signal
would mean: synchronized component constants would expose a polynomial mate,
not merely a deeper finite jet.

The hypotheses in Section 3 are deliberately typed.  Birationality supplies
the complete rational torsor; without it one must separately prove the
constant field in `(2)`.  Squarefreeness of `phi` and smoothness of `Q` are
the Jacobian-relevant boundary.  If `phi` is constant, the appended root
spectrum disappears and the problem reduces exactly to `(13)`, namely the
original Keller question for `U`.  Thus the no-go does not assume away the
desired object; it locates it at spectral synchronization before any
nonconstant dressing.  **QED.**

# The finite sheet can improve clearing but cannot obstruct it

**Status: INDEPENDENT HOSTILE AUDIT ACCEPTED THE THM-3528 ALL-LEVEL RAW
POLYNOMIAL-PACKET LEMMA AT FIXED-MAP SCOPE; THM-3528 WAS SUBSEQUENTLY
PROMOTED BY THE OWNER.**

The proposed lemma is sound.  If a source polynomial `P` has the complete
packet `A(e,m)`, then its fixed-chart norm satisfies

```text
v_L(N(P))=-e+s,             s>=0,
v_L(L^e N(P))=s.
```

Because `N(P)` is already regular on `Spec(Q[a,b,c,L^-1])`, this inequality
proves `L^e N(P)` is a polynomial.  THM-3522 then supplies the complete
packet `A(7e-2m,3e-2m)`.  Starting at `L` with `A(1,0)` therefore gives raw
polynomial cleared norms and complete packets at every level.

The conclusion is deliberately weaker than exact clearing.  It does not
prove `s=0`, hence does not prove `L`-coprimality.  It also proves no new
image component, irreducibility, separability, discriminant recursion, or
general Jacobian statement.

## Inheritance and hostile target

The closest proved mechanism is the localization argument in
THM-3498 -- `level-four-old-boundary-cancellation-and-degree81-discriminant-gate`:
a finite-etale norm belongs to `A[L^-1]`, and an exact generic-`L` valuation
then determines its reduced denominator.  THM-3521, THM-3523, and THM-3527
repeat that mechanism with an expensive finite-sheet nonvanishing test.

The corrected near miss is MISTAKE-415.  A single visible pole exponent is
not a closed norm state; one must retain all complete faces and the finite
sheet.  THM-3522 supplies the complete packet transform but assumes the next
cleared norm is polynomial.

The present hostile question is narrower:

```text
Can a positive finite-sheet valuation make the denominator worse?
```

It cannot.  It adds to the norm valuation and therefore only improves the
old-`L` denominator.  The earlier computations needed `s=0` for exact
clearing and coprimality; polynomiality itself needs only `s>=0`.

The active concept board was:

1. finite-etale regularity on `A[L^-1]`;
2. the reciprocal cubic and its `2+1` Newton decomposition;
3. complete max-`lambda` face noncancellation;
4. ramification-weighted norm valuation;
5. packet-cone invariance under the THM-3522 matrix; and
6. finite-sheet defect versus image/separability gates.

## 1. No hidden global denominators

Put

```text
A=Q[a,b,c],
U=Spec(A[L^-1]).
```

THM-2473 -- `sporadic-keller-branch-tower-depressed-trisection-anatomy`
proves that

```text
F^-1(U) -> U
```

is finite etale of degree three.  A source polynomial is a regular function
on `F^-1(U)`, so its finite-algebra norm is a regular function on `U`:

```text
N(P) in A[L^-1].
```

This statement is global and chart-independent.  Apparent powers of `S`, a
cubic discriminant, or a leading coefficient in a particular inverse formula
cannot survive as additional denominators.  They may occur in a presentation,
but regularity of the finite-algebra determinant cancels them on `U`.

Since `A` is a UFD and `L` is irreducible, every element of `A[L^-1]` has a
reduced form `R/L^h` with `R in A` and `L` not dividing `R`.  It is therefore
enough to prove `h<=e`, equivalently `v_L(N(P))>=-e`.

## 2. The generic reciprocal cubic has one finite branch and one ramified pair

At the generic DVR `A_(L)`, the inverse cubic is

```text
E(w)=Lw^3+Tw-2c.
```

Set `u=1/w`.  Its reciprocal equation is

```text
L+Tu^2-2cu^3=0.                                      (1)
```

The canonical witness on `L=0` gives

```text
(c,T,S,D)=(1,1,1,1/3),
```

and exact gcd checks give `gcd(L,cTSD)=1`.  Thus all four are units at the
generic divisor.

The Newton polygon of (1) has points `(0,1),(2,0),(3,0)`.  Its segments have
horizontal length/slope

```text
(2,-1/2), (1,0).
```

Equivalently, modulo `L`,

```text
u^2(T-2cu)=0.
```

The root `u_0=T/(2c)` is simple because the derivative there is
`-T^2/(2c)`, a unit.  It gives the one unramified finite branch.  The double
root at zero gives one ramification-index-two place, or two conjugate sheets
after splitting, each with base-normalized valuation `v_L(u)=1/2`.

The finite branch is genuinely regular.  There `w=1/u` is a unit; the exact
THM-3495 inverse numerators express `y,z` with only the generic unit `S` in
the denominator.  Hence every source polynomial has finite-branch valuation

```text
s=v_f(P)>=0.                                          (2)
```

At the canonical specialization this branch is

```text
q=(2,5/6,-7/8),       F(q)=(2/27,1,1),
L(q)=241465/1728 !=0,
```

which independently checks the coordinate and normalization convention.

## 3. The complete packet fixes both divergent valuations

For a monomial `x^i y^j z^k`, put `lambda=i-k`.  On either divergent sheet,
THM-3498 gives

```text
x=u^-1,
y=D/S+O(u),
z=-3(D/S)u+O(u^2).
```

Thus a monomial of `lambda`-weight `r` begins at order `u^-r`.  A complete
packet `A(e,m)` has complete maximum face

```text
in_max-lambda(P)=kappa x^e(3xz-2y)^m,   kappa!=0.
```

On the displayed branch,

```text
3xz-2y -> -11D/S,
```

a nonzero unit.  Therefore no equal-weight cancellation is possible and

```text
P(q_div)=kappa(-11D/S)^m u^-e(1+O(u)).                (3)
```

Each split sheet has base-normalized valuation `-e/2`.

The ramification weighting can be stated in either of two equivalent ways:

- after splitting, add the two conjugate values `-e/2-e/2=-e`; or
- at the single quadratic place use normalized valuation
  `v_W(L)=2`, `v_W(P)=-e`, residue degree one, so
  `v_L(N_W(P))=f_W v_W(P)=-e`.

There is no extra factor of two.  Adding the unramified finite branch (2)
gives exactly

```text
v_L(N(P))=-e+s.                                       (4)
```

## 4. The UFD step proves polynomiality

Equation (4) and `s>=0` give `v_L(N(P))>=-e`.  Since `N(P)` belongs to
`A[L^-1]`, its reduced denominator has exponent at most `e`.  Consequently

```text
Q=L^eN(P) in A,
v_L(Q)=s.                                             (5)
```

This is the whole missing polynomiality gate.  Notice the direction:

```text
s=0  => exact denominator L^e and gcd(Q,L)=1,
s>0  => Q is still polynomial but is divisible by L^s.
```

No claim is made that the second case actually occurs in the packet orbit.
Whether the packet forces `s=0` at every level is a stronger open problem.

## 5. The packet cone is invariant, so induction is lawful

THM-3522 -- `fixed-keller-five-face-renewal-propagation` proves that a
polynomial `Q` from (5) has packet

```text
(e',m')=(7e-2m,3e-2m).                               (6)
```

Packet admissibility is equivalent here to

```text
e>=0,  0<=m<=e,  3|m.
```

Equation (6) preserves it:

```text
e'>=5e>=0,
m'>=e>=0,
e'-m'=4e>=0,
3|m'.                                                 (7)
```

The base polynomial `L` has `A(1,0)`.  Equations (5)--(7) therefore iterate
without another hypothesis and give

```text
(1,0), (7,3), (43,15), (271,99),
(1699,615), (10663,3867), (66907,24255),
(419839,152211), ...                                  (8)
```

for raw cleared norms at all levels.  Existing primitive normalizations such
as `H=2^6LN(L)` differ by rational units.  Since a cubic norm satisfies
`N(rP)=r^3N(P)` for `r in Q^*`, such rescalings change neither valuations nor
packet support and do not obstruct the induction.

## 6. Hostile controls

The audit attacked the proposed proof at each typing boundary.

### Nonmonic norm versus raw resultant

For the test element `P=x=w`, Vieta gives

```text
N(w)=2c/L,
Res_w(E,w)=2c.
```

Thus the raw resultant is `L*N(w)`.  Forgetting the nonmonic leading factor
would erase exactly the pole being audited.  The proof uses the finite-algebra
norm/product convention throughout, so it passes this hostile.

### Hidden chart denominators

The argument does not infer regularity from one inverse formula.  It uses the
finite-etale norm on `U`; therefore only powers of `L` may occur globally.
The exact gcd and canonical-point checks verify that `c,T,S,D` are units at
the local divisor used for the branch expansion.

### Complete-face cancellation

The max-`lambda` statement is used as a complete initial form, not as the
presence of one monomial.  Its residual value is the explicit nonzero unit
`kappa(-11D/S)^m`.

### Ramification weights

The quadratic pair contributes `-e`, not `-2e` and not `-e/2`; the
split-sheet and normalized-place calculations agree.

### Cone and scalar normalization

An exact sweep of all `15,251` admissible packets with `e<=300` verifies
every input/output exponent and divisibility condition.  Rational scalar
normalizations are units for `v_L` and cube under the norm, as required.

## 7. Relation to existing canon

This lemma is not refuted by MISTAKE-415.  That correction says a visible
face alone is not a closed state and the finite sheet must be retained.  The
proof retains the complete five-face packet through THM-3522 and retains the
finite sheet as the defect `s`; it simply observes the sign `s>=0`.

It is also not a duplicate of THM-3522.  THM-3522 proves packet propagation
conditional on polynomiality.  Equation (5) discharges that condition.

THM-3498, THM-3506, THM-3521, THM-3523, and THM-3527 prove the stronger
special case `s=0` at named rungs.  Their finite-sheet computations remain
load-bearing for exact denominator and `L`-coprimality, and for any later
argument that needs those properties.  They are unnecessary only for the
weaker raw polynomiality/packet induction.

The closest historical proof pattern is THM-3498's UFD localization after an
exact valuation.  The new step is replacing exact finite-sheet nonvanishing
by the automatic inequality `s>=0` when only polynomiality is required.

## 8. Exact scope

Accepted:

- every recursively defined raw cleared norm is in `Q[a,b,c]`;
- every rung has the complete fixed-chart packet on the orbit (8);
- the old-`L` multiplicity of the next cleared norm is the finite-sheet
  valuation `s>=0`.

Not accepted or implied:

- `s=0`, `L`-coprimality, primitive integral normalization, or a new prime;
- irreducibility, squarefreeness, generic degree, or image-equation status;
- a new nonproperness component or exact classification of the tower;
- eliminant full degree, block coprimality, separability, or a discriminant
  square-class recursion at unlicensed levels;
- an arbitrary Keller-map theorem, classification of counterexamples,
  `JC(2)`, `DC(2)`, or any general Jacobian-conjecture conclusion.

## Reproduction

```text
python -B 04-computation/keller_all_level_cleared_norm_polynomiality_independent_audit_20260816.py
python -B -O 04-computation/keller_all_level_cleared_norm_polynomiality_independent_audit_20260816.py
```

Normal and optimized transcripts are identical.  The semantic digest is
`6fe70dcf5a0f1bd4f76ef8bc4986f79be1d74e4a0b6a71e40520dae67f4e456e`.
Script/output LF-normalized SHA-256:

```text
4e79d10d2a90cfdcbd6948d22b7385da6ea88b923f3483da96448af1ea1cdc77
e7247cb01a61558ee2af6ef662610e946364f6d9d32d44766841a501644ee30b
```

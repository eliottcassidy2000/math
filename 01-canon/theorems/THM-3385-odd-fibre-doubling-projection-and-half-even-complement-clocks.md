---
id: THM-3385
title: "Odd-fibre doubling projection and q-sheet complement clocks"
status: >
  PROVED analytic identity + VERIFIED-EXACT literal-body census +
  INDEPENDENTLY HOSTILE-AUDITED.  Under the degree-q circle map, q-divisible danger
  speeds descend exactly.  A transverse speed u blocks at most
  gcd(u,q) ceil((q/gcd(u,q))/7) sheets; if the summed capacity is below q,
  the image of the body-safe set is exactly the complement of the descended
  core combs.  On the canonical L=14 lcm(F), D=L/q grid, the unsupported
  open D-cells are exactly those core combs with the D-grid removed, so the
  divided core is an explicit complement-clock tuple.  This gives 6,420
  structural body/divisor certificates in the literal six-body universe and
  explains all seven refined k=3 THM-3366 hits as the q=2 fixed-core family.
  Capacity is sufficient, not necessary; arbitrary reflected phases,
  physical drift realization, and LRC(14) remain open.
source: codex-2026-08-14-q-fibre-complement-clock
audit: independent fibre proof, downward-recurrence census, strictness, phase, budget, and typing audit
depends_on:
  - THM-2928-critical-seven-comb-grid-tensorization-and-drift-tariff
  - THM-3366-all-sector-complement-clock-completion
  - LRC(<=13)
related:
  - THM-2941-critical-seven-slot-scalar-wall-and-balanced-boundary
  - THM-3381-reflected-residue-affine-phase-transport-and-frozen-tree-stability
script: 04-computation/lrc14_q_fibre_complement_clock_thm3385.py
output: 05-knowledge/results/lrc14_q_fibre_complement_clock_thm3385.out
script_sha256: 8d3add2104ca91e89017ae643f573c8c5f6b45a4fdfecd9d65c16f48fd51b356
output_sha256: 5b7349e13e36e8424346aa2c50f26f760f766a161f046ca09ad7bcd6e4d665b4
semantic_sha256: 76345873bbfbe37fde097e9a23e92744fb376ad0eb94d7f0eb1ff33c5f7223de
hash_basis: LF-normalized bytes
---

# THM-3385 -- a transverse fibre prevents the core quotient from vanishing

**PROVED analytic identity + VERIFIED-EXACT literal-body census +
INDEPENDENTLY HOSTILE-AUDITED.**

## 1. Inheritance and connection contract

THM-2928 supplies the exact body-safe cell word.  THM-3366 says that an
explicit cover of its unsupported quotient cells compiles a hypothetical
mixed LRC cover into a forbidden cover by at most twelve integer clocks.
The seven refined `k=3` hits of THM-3366 have the unexplained form

```text
F={2,6,8,10,14} union {u},   u odd,       D=L/2.          (1)
```

The closest hostile is THM-3381: a one-unit drift perturbation can carry a
near-half-turn located phase.  Thus no continuity from divisible to nearly
divisible speeds is lawful.  The least-used sidecar is instead the finite
fibre of the exact covering map.

| field | connection |
|---|---|
| source | a body-safe set with speeds split into `q`-divisible core and transverse speeds |
| target | the complement of the divided core danger combs |
| map | the degree-`q` circle map `pi_q(x)=qx` |
| preserved | pointwise danger status of every `q c`, strict openness, quotient order |
| destroyed | which sheet survives, transverse owner, reflected phase, physical tail ancestry |
| sidecar | the `q` sheets and the blocked-sheet capacity of every transverse speed |
| decisive hostiles | radius above `1/4` at `q=2`, two transverse speeds covering both sheets, wrong divisor, and the `k=0` clock budget |

## 2. Exact finite-fibre capacity

For `0<rho<1/2`, write

```text
D_w^(rho)={x in T: ||wx||<rho},          pi_q(x)=qx.      (2)
```

Fix `q>=2` and a speed `u`.  Over a base point `y`, the `q` preimages are
`x_k=(y+k)/q`, `0<=k<q`.  Put

```text
g=gcd(u,q),                  m=q/g.                       (3)
```

The phases `u x_k` consist of `m` equally spaced circle points, each repeated
`g` times.  An open arc of length `2rho` contains at most

```text
ceil(2rho m)                                                     (4)
```

of the distinct points.  Indeed, if it contains `r` consecutive grid points,
then `(r-1)/m<2rho`, hence `r<=ceil(2rho m)`; a centered rotation of that many
consecutive points attains the bound.  Therefore `u` blocks at most

```text
kappa_(q,rho)(u)=g ceil(2rho q/g)                         (5)
```

of the `q` sheets.

Let `C` be a finite set of positive core clocks and `U` a finite set of
transverse speeds, and put

```text
G(q;C,U)=T minus
  [union_(c in C) D_(qc)^(rho) union union_(u in U) D_u^(rho)]. (6)
```

If

```text
sum_(u in U) kappa_(q,rho)(u)<q,                          (7)
```

then

```text
pi_q(G(q;C,U))=T minus union_(c in C) D_c^(rho).           (8)
```

To prove the forward inclusion, if `y in D_c^(rho)`, every preimage satisfies

```text
qc(y+k)/q=cy+ck,
```

so all sheets are core-dangerous.  Conversely, outside every `D_c^(rho)` all
core clauses are safe on all sheets.  Equation `(7)` and the union bound leave
at least one sheet outside every transverse clause.  That sheet lies in
`G(q;C,U)`, proving `(8)`.

At the LRC radius `rho=1/14`, equation `(5)` is

```text
kappa_q(u)=gcd(u,q) ceil((q/gcd(u,q))/7).                 (9)
```

In particular, one non-`q`-divisible transverse speed always works: then
`m>=2` and `ceil(m/7)<m`, so `kappa_q(u)<q`.  For `q=2`, this is exactly the
odd-speed/half-turn mechanism suggested by `(1)`.

## 3. Exact aligned-cell quotient

Now take `rho=1/14`, a literal body

```text
F={q c:c in C} disjoint union U,
L=14 lcm(F),                     D=L/q,                  (10)
```

with nonempty `C`, nonempty `U`, and `(7)`.  Let `J subset Z/LZ` be the
THM-2928 body-safe `L`-cell word and let

```text
S_D=J mod D.                                             (11)
```

Because every `q c` divides `lcm(F)`, every endpoint of `D_c` lies on the
`D`-grid:

```text
D/(14c)=lcm(F)/(qc) in N.                               (12)
```

Every body endpoint lies on the `L`-grid as well.  Thus danger status is
constant on every relevant open cell, and the continuous identity `(8)`
descends without approximation to

```text
U_D(S_D)
 = [union_(c in C) D_c] minus (D^(-1)Z/Z).              (13)
```

The removed grid is not cosmetic.  `U_D` is a union of open cells and contains
no grid point, whereas zero belongs strictly to every `D_c`.  THM-3366's
aligned owner clock (or its extra clock `D` when `k=0`) owns precisely this
grid.  Equation `(13)` proves that the divided core tuple `C` covers every
unsupported open cell, with no set-cover search.

## 4. LRC compiler and its exact budget

Suppose the THM-2928 row has `k` aligned and `p=7-k` drift tails.  If a
hypothetical quotient drift cover existed, THM-3366 and `(13)` would produce
a global integer-clock cover using at most

```text
p+k+|C|=7+|C|             when k>=1,
p+1+|C|=8+|C|             when k=0.                     (14)
```

Cited `LRC(<=13)` rules out such a cover whenever the displayed total is at
most twelve.  Hence

```text
k>=1 and |C|<=5  ==> row impossible,
k=0  and |C|<=4  ==> row impossible.                    (15)
```

For a six-speed body with nonempty transverse set, `|C|<=5`; therefore every
certificate `(7)` closes every eligible `k>=1` support row.  The one-
transverse/five-core case sits exactly one clock beyond the `k=0` budget.

## 5. Literal six-body classification

For each six-subset `F` of `{1,...,14}` and each `2<=q<=14`, split `F` as in
`(10)`, test `(7)`, and take only the divisor `D=L/q`.  There are exactly
`6,420` resulting body/divisor certificates.  Their degree profile is

```text
q:       2    3    4    5    6     7     11    12    13    14
rows:  147   45   35  495   57   2079   1287   196  1287   792. (16)
```

The exact companion independently rebuilds every safe cell word and verifies
`(13)` as equality of integer interval lists in all `6,420` cases.  Applying
the inherited support cutoffs and the budget `(15)` gives the following
structural subcensus of THM-3366's pool terminals:

| `k` | certified rows | denominator-shape occurrences |
|---:|---:|---:|
| 0 | 6,273 | 83,942,791,283,821 |
| 1 | 6,420 | 3,115,391,844,730 |
| 2 | 6,420 | 106,719,668,256 |
| 3 | 6,420 | 3,253,046,409 |
| 4 | 1,272 | 8,230,092 |
| 5 | 227 | 35,818 |
| 6 | 192 | 192 |

The occurrence weights use THM-3366's exact divisor--Mobius formula

```text
a_p(D)=sum_(e|D) mu(D/e) binom(tau(e)+p-2,p).            (17)
```

These rows and occurrences are a structural partition **inside** the earlier
terminal census, not an additive decrement.  The `k=4,5,6` entries are only
controls because THM-2928 already closes those sectors globally.

## 6. The seven refined `k=3` hits

For the fixed core

```text
C={1,3,4,5,7},
F=2C union {u},                 u in {1,3,5,7,9,11,13}, (18)
```

the transverse capacity is one.  On the two sheets the odd phase differs by
`1/2`, so their radius-`1/14` danger arcs cannot fire together.  Equations
`(8)` and `(13)` give the same completion tuple

```text
(1,3,4,5,7)                                             (19)
```

for every odd `u`; only `L` and `D=L/2` change.  THM-3366's independently
audited refined intersection becomes:

| odd `u` | `L` | `D` | refined occurrences |
|---:|---:|---:|---:|
| 1, 3, 5, 7 | 11,760 | 5,880 | 388 each |
| 9 | 35,280 | 17,640 | 1,008 |
| 11 | 129,360 | 64,680 | 2,544 |
| 13 | 152,880 | 76,440 | 2,544 |

Thus all `7,648` deleted occurrences are one exact fibre family.  The theorem
explains why the core and completion clocks are fixed and why the divisor is
`L/2`; it does not independently derive which rows survive the preceding
one-spike ledger.

## 7. Why the fibres of size four and six are Boolean, not tournaments

At fixed `y`, use the `q` sheets as vertices.  A transverse speed contributes
the phase-dependent subset of sheets it blocks.  The observable is whether
the union of these hyperedges covers every vertex.  There is no intrinsic
pairwise orientation, so forcing this object into a tournament would discard
the simultaneous-cover predicate.

At `rho=1/14`, the residue-type capacities are

| fibre | transverse residue type | `gcd(u,q)` | maximum blocked sheets |
|---:|---|---:|---:|
| 4 | `1,3 mod 4` | 1 | 1 |
| 4 | `2 mod 4` | 2 | 2 |
| 6 | `1,5 mod 6` | 1 | 1 |
| 6 | `2,4 mod 6` | 2 | 2 |
| 6 | `3 mod 6` | 3 | 3 |

This is the exact Boolean realization behind the sizes four and six.  A
tournament becomes relevant only after an additional intrinsic binary
relation is supplied; the quotient theorem itself supplies none.

## 8. Degree monoid and staged folds

The covering maps multiply:

```text
pi_r o pi_q=pi_(rq).                                    (20)
```

If the descended core `C` itself splits as `rC' union U'` and satisfies the
corresponding `r`-sheet capacity inequality, applying `(8)` twice gives

```text
pi_(rq)(G(q;C,U))=T minus union_(c in C')D_c.            (21)
```

The sheet degree is therefore genuinely multiplicative, while transverse
speeds and endpoint grids move within a grade as sidecars.  This is an exact
LRC instance of the degree-graded-monoid grammar: composition multiplies
fibre degree, but it does not erase the phase/capacity data required at each
stage.

## 9. Sharp boundaries and failed extensions

1. **Capacity is sufficient, not necessary.**  For `q=2`, both transverse
   pairs `{1,3}` and `{1,9}` have capacity sum two.  The first still leaves a
   safe sheet for every base point; the second blocks both sheets at `y=1/9`
   (`x=1/18` is blocked by speed one and `x+1/2=5/9` by speed nine).
2. **Strict openness is load-bearing.**  At the actual LRC radius, take
   `q=7,u=1,y=1/2`.  Two sheets lie exactly at `+/-1/14`; the strict comb
   blocks neither, while a closed comb would block both and invalidate the
   capacity one.  Separately, at `q=2` one transverse speed at `rho=1/4`
   cannot contain two opposite phases, but for any larger radius `x=1/4`
   and `x+1/2` are both dangerous.
3. **Divisible speeds must descend.**  Treating `u=2` as transverse at `q=2`
   fails at `y=0`: both sheets are blocked.  It is core clock one.
4. **The divisor is part of the theorem.**  On the first family body, using
   `D=L/4` gives first unsupported gap `(0,1/28)`, not the half-even target
   `(0,1/14)`.  Several core endpoints are no longer on that grid.
5. **The `k=0` budget is attained.**  For all `147` degree-two five-core rows,
   an exact strict-endpoint/grid-aware audit of all `216,237` candidate subsets
   finds no cover by at most four clocks from `{1,...,14}`.  This is a finite-
   pool statement, not unrestricted noncoverability.
6. **No near-divisible extension follows.**  THM-3381's located phase defect
   remains the governing hostile for reflected or perturbed drifts.

The result is a body-address quotient theorem.  It does not produce physical
tail clocks, preserve a sheet owner, prove the capacity test necessary, close
any surviving refined row beyond THM-3366's exact intersection, or prove
LRC(14).

## 10. Exact verification

The standard-library companion:

- checks the exact open-grid capacity through orbit size `100`;
- independently rebuilds all `6,420` literal body/divisor interval identities;
- performs `1,335` strict-boundary continuous event samples;
- reproduces `(16)`, the complete table in Section 5, and `(17)`;
- checks the seven rows and their frozen `7,648` occurrence count;
- checks the `q=4,6` Boolean capacities and every hostile in Section 9; and
- contains no floating literal or optimization-dependent `assert`.

Reproduce with

```text
python 04-computation/lrc14_q_fibre_complement_clock_thm3385.py
python -O 04-computation/lrc14_q_fibre_complement_clock_thm3385.py
```

Ordinary and optimized runs LF-normalized-byte-match the stored output.
An independent audit rederived the fibre bound and both inclusions, replayed
the census using downward divisor recurrence instead of Mobius summation, and
verified the strictness, capacity-equality, alignment, phase, budget, and
non-tournament controls.

**QED.**

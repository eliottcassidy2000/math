# HYP-7024 — the resonant-peel miss-pattern extremal

**Status:** MIXED / OPEN (codex-2026-07-16-S17).  `THM-891` proves the exact reduction
and closes residues 1 through 5 and the positive side of residue 6 at the non-sharp
`0.097` propagation slack.  Only the negative side of residue 6 remains at that slack.
The sharp uniform inequality is finite-exact through diameter 20.  Three tempting
shortcuts were refuted.

## Exact target

For a six-offset core `E={0,e_1,...,e_5}`, let

- `B_s=meas{E misses exactly {s}}`,
- `A_{s,c}=meas{E misses exactly {s,c}}`,
- `F_r(E)=r C_r(E)` be the `THM-891` scaled resonant coefficient.

The seven-microcell theorem gives an explicit integer kernel `K_r/49` on the one- and
two-miss patterns:

`F_r(E)=sum_s B_s K_r({s})/49 + sum_(s<c) A_{s,c}K_r({s,c})/49`.

The sharp conjecture is

> **`max_(1<=r<=6) |F_r(E)| <= 16/343` for every six-offset core `E`,**

with equality only up to dilation at `E={0,1,2,3,4,6}`, `r=1`.

At `r=1`, every singleton kernel is `-6/49` and every pair kernel is `-2/49`.  If
`p_j=meas{E misses exactly j inner sectors}`, the residue-one inequality is exactly

`F_1(E)=-2(3p_1+p_2)/49`, so the sharp target becomes

> **`3p_1+p_2 <= 8/7`.**

The extremal core has `p_1=1/4`, `p_2=11/28`, hence equality.

## Exact evidence

The event-sweep certificate exhausts every primitive core with largest positive speed
at most 20: `15,246` cores.  The unique maximum is `16/343`.  The exact per-residue
maxima are

| residue | maximum | core |
|---:|---:|---|
| 1 | `16/343` | `{0,1,2,3,4,6}` |
| 2 | `74/2401` | `{0,1,2,4,5,7}` |
| 3 | `2671/96040` | `{0,2,6,7,8,10}` |
| 4 | `61/2940` | `{0,1,2,3,4,15}` |
| 5 | `367/20580` | `{0,2,6,8,10,15}` |
| 6 | `199/12005` | `{0,1,4,5,7,11}` |

The two leading cores are stable from diameter 12 through diameter 20.  The majority
tournament on residues is also stable from diameter 14 to 20, with the unique path
`1 -> 2 -> 3 -> 4 -> 6 -> 5`.

## Pair law and refuted shortcuts

### 1. Fixed collision first moment — false, with an exact replacement

Let `C(x)` count colliding runner pairs.  Then `p_1=Pr(C=0)` and `p_2=Pr(C=1)`, so the
binding target is `3Pr(C=0)+Pr(C=1)<=8/7`.  It was tempting to use
`E[C]=15/7`, assuming every pair collides with probability `1/7`.  That assumption is
false: the exact collision masses are

- `(a,b)=(1,8)`: `1/4`;
- `(a,b)=(2,9)`: `2/9`.

After dividing a positive pair `(a,b)` by its gcd to the coprime pair `(A,B)`, the
replacement law proved in `THM-891` is

`Pr(same sector)=1/7 + [A=B mod 7] (A mod 7)(7-A mod 7)/(7AB)`.

Thus every pair mass is at least `1/7`.  If `C` counts collisions among the six core
points, `E[C]>=15/7`; comparing `C<=1` on the event `N<=2` with the global `C<=15`
proves the useful but non-sharp bound

`p_1+p_2<=45/49`.

For any two positive runners, the `p_1` event puts them in distinct nonzero sectors.
If `D` is their same-sector event and `H_j` is their common-sector-`j` event, then

`Pr(distinct and nonzero)=5/7-Pr(D)+2Pr(H_0)<=5/7`,

because time reversal equates `H_0` and `H_6` and both are disjoint subsets of `D`.
Thus `p_1<=5/7`, and

`|F_1|=2(3p_1+p_2)/49
      <=2(2(5/7)+45/49)/49=230/2401<0.097`.

The residue-five singleton and pair kernel norms are respectively `5` and `4`, so

`|F_5|<=(5p_1+4p_2)/49<=225/2401<0.097`.

This closes limiting residues `1` and `5` at the propagation slack, not at the sharp
`16/343` target.

### 2. Full pair-sector rays close residues 2, 3, and 4

For coprime positive speeds `A,B`, their full unordered sector-pair distribution is
independent-uniform when either speed is zero modulo seven.  Otherwise it has the form

`P(A,B)=U+D_(A mod 7,B mod 7)/(AB)`.

Hence it lies on one of 21 line segments from `U` to the product-minimal representative
of its unordered nonzero residue pair.  Exact quadratic weights on sector-pair counts,
listed and proved in `THM-891`, give

`-40/441 <= F_2 <= 230/2401`,

`-230/2401 <= F_3 <= 19/196`,

`-2/21 <= F_4 <= 232/2401`.

All six endpoints are below `0.097` in absolute value.  This is a diameter-free theorem:
the verifier checks 462 sector-count states and the 22 exact pair-polytope vertices.

### 3. Residue 6 is now the only limiting slack obstruction

The positive kernel bound gives

`F_6 <= (6p_1+2p_2)/49 <=230/2401<0.097`.

On the negative side, `K_6({1,5})=K_6({2,4})=-12`; all other pair kernels are `2`,
and singleton kernels are `-1` except `K_6({3})=6`.  Thus the remaining non-sharp crux
is a higher-order bound on the synchronized miss mass `A_15+A_24`.  The exact pair-ray
relaxation is insufficient, so this is where the higher relation spectrum of the
concurrent `THM-890` must enter.

The excess term retains the same short arithmetic relations isolated by `HYP-7021`
and the coincidence divisor spectrum of `THM-892`.

### 4. Pointwise residue-one dominance — false

Residue one wins on `15,211/15,246` diameter-20 cores, but not all.  For example,
`E={0,3,5,6,7,8}` is won by residue two:

`|F_2|=211/7203 > |F_1|=2081/72030`.

Residue three wins eight diameter-20 cores.  Thus the global theorem needs the sharp
`r=1` occupancy inequality plus separate signed inequalities for residues `2,...,6`;
it cannot conclude them from pointwise dominance.

### 5. Miss-pattern inversion symmetry — false

The map `x->1-x` inverts every moving point, but the stationary point at the
sector-zero boundary stays fixed.  Consequently it does not permute the six inner
miss patterns.  Already for `{0,1,2,3,4,5}`, `B_1=1/35` while the naively reflected
`B_5=1/28`.  Symmetrizing the residue kernels without carrying the stationary-boundary
sidecar gives false universal bounds.

## Proof program

1. **Relation-stratified collision moments.**  Expand pair and higher collision masses
   from the exact pair law by the congruence lattice, using the `HYP-7021`
   short-relation parameter and `THM-892` signed coincidence spectrum rather than a
   fictitious fixed first moment.
   The relation-poor stratum should have independent slack; finitely many short-relation
   shapes should contain the extremizer.
2. **Residue-six higher-order certificate.**  Bound the synchronized mass
   `A_15+A_24` using triple and higher relation terms from `THM-890`; pair marginals are
   now known to be insufficient.  `THM-894` independently names exactly this next rung
   as the level-three overlap tensor, and `THM-896-level3-crossing` proves the order-three Bonferroni
   crossing while leaving its triple-beat enhancement cap open.  Build one
   relation-localized triple bound that feeds both consumers.  Do not impose the false
   reflection symmetry `s -> 7-s`.
3. **Compact/relation-rich split.**  The extremizer is compact and relation-rich, while
   random wide cores are much smaller.  A quantitative relation-lattice tail plus a
   finite compact sweep would parallel the balanced-comb split of `THM-889/HYP-7021`.
4. **Finite-`t` bridge.**  After the limiting inequality, bound the `O_E(1/t)` wall-cell
   term uniformly enough to feed the constant-propagation ledger for every `w>=diam`.
   The resolved `HYP-7027` says wall-movie cycles are all expressive, so retain the
   palette/relation sidecar rather than seeking a silent-cycle quotient.

## Scope and assumption challenge

The faithful vertices are miss patterns, not runners or swing arcs.  The miss-pattern
quotient preserves the limiting LRC two-scale error coefficient and destroys wall order,
finite-`t` aliasing, and core identity.  Residue vertices preserve only the ordering of
the six owner resonances.  This challenged the initial runner/arc assumption and explains
why the scalar residue tournament is transitive telemetry rather than the missing proof.

Verification:
`04-computation/lrc14_resonant_cross_section_cancellation_codex_S17.py` and
`05-knowledge/results/lrc14_resonant_cross_section_cancellation_codex_S17.out`.

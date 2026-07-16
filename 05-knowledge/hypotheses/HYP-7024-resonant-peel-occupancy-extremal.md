# HYP-7024 — the resonant-peel miss-pattern extremal

**Status:** MIXED / OPEN (codex-2026-07-16-S17).  `THM-891` proves the exact reduction.
The sharp uniform inequality is finite-exact through diameter 20.  Two tempting
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

## Refuted shortcuts

### 1. Fixed collision first moment — false

Let `C(x)` count colliding runner pairs.  Then `p_1=Pr(C=0)` and `p_2=Pr(C=1)`, so the
binding target is `3Pr(C=0)+Pr(C=1)<=8/7`.  It was tempting to use
`E[C]=15/7`, assuming every pair collides with probability `1/7`.  That assumption is
false: the exact collision masses are

- `(a,b)=(1,8)`: `1/4`;
- `(a,b)=(2,9)`: `2/9`.

The collision moment retains the same short arithmetic relations isolated by
`HYP-7021`; a relation-blind moment LP cannot prove the target.

### 2. Pointwise residue-one dominance — false

Residue one wins on `15,211/15,246` diameter-20 cores, but not all.  For example,
`E={0,3,5,6,7,8}` is won by residue two:

`|F_2|=211/7203 > |F_1|=2081/72030`.

Residue three wins eight diameter-20 cores.  Thus the global theorem needs the sharp
`r=1` occupancy inequality plus separate signed inequalities for residues `2,...,6`;
it cannot conclude them from pointwise dominance.

## Proof program

1. **Relation-stratified collision moments.**  Expand pair and higher collision masses
   by the congruence lattice, using the `HYP-7021` short-relation parameter and
   `THM-892` signed coincidence spectrum rather than a fictitious fixed first moment.
   The relation-poor stratum should have independent slack; finitely many short-relation
   shapes should contain the extremizer.
2. **Signed miss-pattern LP.**  For each residue `2,...,6`, use the exact `K_r` table,
   reflection symmetry `s -> 7-s`, and realizable factorial moments of the missed set.
   The target constants in the evidence table are far below `16/343`, leaving room for
   a non-sharp certificate.
3. **Compact/relation-rich split.**  The extremizer is compact and relation-rich, while
   random wide cores are much smaller.  A quantitative relation-lattice tail plus a
   finite compact sweep would parallel the balanced-comb split of `THM-889/HYP-7021`.
4. **Finite-`t` bridge.**  After the limiting inequality, bound the `O_E(1/t)` wall-cell
   term uniformly enough to feed the constant-propagation ledger for every `w>=diam`.

## Scope and assumption challenge

The faithful vertices are miss patterns, not runners or swing arcs.  The miss-pattern
quotient preserves the limiting LRC two-scale error coefficient and destroys wall order,
finite-`t` aliasing, and core identity.  Residue vertices preserve only the ordering of
the six owner resonances.  This challenged the initial runner/arc assumption and explains
why the scalar residue tournament is transitive telemetry rather than the missing proof.

Verification:
`04-computation/lrc14_resonant_cross_section_cancellation_codex_S17.py` and
`05-knowledge/results/lrc14_resonant_cross_section_cancellation_codex_S17.out`.

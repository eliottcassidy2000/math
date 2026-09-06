# Clean-room referee audit: the two `72/539` profiles

## Verdict

**PASS, with a type-narrowing recommendation.** Exact cube-edge section
reconstruction and an independent exact circle wall-cell partition reproduce
all numerical claims in the proposed bridge packet. The strongest faithful
statement is

> **PROVED equality of unlabelled weighted scalar profiles:** the auxiliary
> THM-4437 coefficient-envelope atoms at `(4,7,11)` and the three THM-4449
> dyadic pair-energy atoms at the speed shape `(1,7,11)` both have weight
> multiset `{2/49,2/49,4/77}`, and both totals are `72/539`.

This is stronger than equality of totals, because the central atom is uniquely
matched to the `4/77` edge and the two outer atoms match the two `2/49` edges.
It is weaker than a labelled correspondence: the two outer matches can be
swapped, and the coefficient object carries no circle address or sheet owner.

The primary report's phrase **auxiliary coefficient-envelope slope
`F_v(w)` used in THM-4437**. It is not a physical Haar measure and not a new
THM-4437 selector conclusion. No LRC(14) consequence follows.

## Inheritance and truth status

The closest proved mechanisms are the exact error-cube section functional in
[THM-4437 / all-parity network reduction](../../01-canon/theorems/THM-4437-lrc14-all-parity-network-reduction-to-three-low-circuits.md)
and the seventh-rounding/pair-owner identity in
[THM-4449 / dyadic seventh-rounding energy](../../01-canon/theorems/THM-4449-lrc14-dyadic-seventh-rounding-energy-and-residual-haar-entry.md).
Both theorem files are `PROVED ... + FINITE-EXACT + INDEPENDENTLY AUDITED` at
the audited checkout. The predecessor
[THM-4434 / universal scale-three network projection](../../01-canon/theorems/THM-4434-lrc14-universal-scale-three-network-projection-bound.md)
has the same proved/audited status and supplies the underlying section
notation. Exact-token search finds the earlier `72/539` in the THM-4437
coefficient outputs and the later occurrence in THM-4449; no active
`MISTAKE-*` entry retracts either value.

The inherited hostile is the even-speed same-owner phenomenon at speed `4`.
The corrected near miss is the proposed family map
`(1,p,q) -> (q-p,p,q)`, which already fails at `(p,q)=(7,13)`. The useful
sidecar is the complete six-component quotient address atlas, because scalar
mass alone cannot certify zero owner overlap.

## 1. Independent coefficient reconstruction

For

```text
v=(-11,4,7),  w=(1,1,1),  R=3/14,
```

all coefficient magnitudes are units modulo three, so the allowed defects are
`-3,0,3` and the residue weight is `2/3`. The clean-room script intersects
each plane `v.e=delta` with the cube `[-R,R]^3` by enumerating its twelve cube
edges. This does not use the proposed polygon clipper. In the report's declared
coordinate `s=(e_0-e_2)/4`, it obtains

```text
delta=-3: s in [ 9/196,  3/28],  weighted width 2/49,
delta= 0: s in [ -3/77,  3/77],  weighted width 4/77,
delta=+3: s in [ -3/28, -9/196], weighted width 2/49.
```

Therefore

```text
F_v(w)=2/49+4/77+2/49=72/539.
```

A second, generic edge-intersection compiler regenerates all mixed-sign
sectors and normalized speed-polytope vertices. It confirms that this is the
ceiling for coefficient pattern `(4,7,11)`, attained at the equal-speed
boundary, and that `(4,7,11)` is the unique `72/539` pattern in the exact
750-pattern THM-4437 max-coefficient-18 universe.

## 2. Independent dyadic reconstruction

For the quotient pair cross-combs of the genuine speed triple `(1,7,11)`,
the exact wall-cell engine gives

```text
Sigma_(1,7):  (6/49,1/7),       (6/7,43/49),   mass 2/49,
Sigma_(1,11): (6/77,8/77),      (69/77,71/77), mass 4/77,
Sigma_(7,11): (13/49,2/7),      (5/7,36/49),   mass 2/49.
```

It independently evaluates every cell midpoint and every joining boundary;
it does not call the proposed interval union routines. The six components are
pairwise disjoint. Hence the THM-4449 mixed-owner correction has measure zero,
so pair energy equals physical failure mass. A separate physical `x`-circle
partition finds twelve components and the same total:

```text
E(1,7,11)=mu(F_(1,7,11))=72/539.
```

As an analytic cross-check, the centered-seventh coordinates `(d,e)` for
the three primitive pairs are

```text
(1,7):(-3,3),  (1,11):(-1,-2),  (7,11):(2,2),
```

and THM-4449's formula
`2/49*(1+(e^2-d^2)/(pq))` reproduces the same profile.

## 3. Hostile audit and exact logical scope

### Affine labels

All six target permutations were tested. There is no rational affine map
`n -> an+b` carrying `{1,7,11}` onto `{4,7,11}`. Equivalently, their sorted
positive difference multisets are `{4,6,10}` and `{3,4,7}`, which are not
proportional.

### Projective rays

There is no permutation plus common rational dilation carrying one primitive
triple to the other. This is the **THM-4449 projective-ray notion**. It does
not assert the false stronger claim that no abstract `PGL_2(Q)` transformation
can map three rational points to three rational points; such a map would not
preserve multiplication combs and is not an LRC symmetry.

### Owners

The coefficient object has no intrinsic sheet-owner relation, so an entirely
unqualified “no owner conjugacy” is ill-typed. The exact hostile proves the
appropriate scoped statement: any proposed conjugacy that reinterprets the
coefficient magnitudes as speed labels and preserves same-label two-sheet
ownership is impossible. For every odd speed in `(1,7,11)`, the same-owner
two-sheet set has measure zero. For speed `4`,

```text
D_4+1/2=D_4,  same-owner mass=1/7,
mu(F_(4,7,11))=9/49 != 72/539.
```

Thus the naive owner transport is decisively refuted.

### No family law

The tempting transformation does not preserve the scalar even one step away:

```text
coefficient ceiling at (6,7,13) = 233/1911,
physical mass at (1,7,13)       = 240/1911.
```

This is enough to refute a universal `(1,p,q) -> (q-p,p,q)` identity. It does
not claim that the audited equality is the only isolated equality outside the
THM-4437 target-value census.

## 4. Exact connection ledger

```text
source:    THM-4437 auxiliary defect-section slope at
           v=(-11,4,7), w=(1,1,1)
target:    THM-4449 quotient pair energies for speeds (1,7,11)
map:       equality of weighted multisets
           {-3:2/49, 0:4/77, +3:2/49}
           and {{1,7}:2/49, {1,11}:4/77, {7,11}:2/49}
preserved: scalar atom weights, central/outer multiplicities, total mass
lost:      coefficient signs within the outer orbit, speed labels, circle
           addresses, sheet owners, physical distinctness, dilation address
sidecar:   six disjoint quotient pair components, proving Omega=0
hostiles:  no affine/common-scale map; even-4 owner invariant; nearby family
           failure by 7/1911
consumer:  none presently; use only as a search heuristic for future
           weighted-profile identities
```

There are exactly two weight-preserving atom bijections, differing by the
exchange of defects `-3,+3` against edges `{1,7},{7,11}`. Neither is selected
by the retained data. Accordingly, **“unlabelled weighted-profile shadow”** is
the sharp structural description; “affine/projective/owner conjugacy” is not.

## 5. Reproduction and hashes

Run:

```powershell
python -B 04-computation/lrc14_72_539_weighted_profile_shadow_20260906_independent.py
python -O -B 04-computation/lrc14_72_539_weighted_profile_shadow_20260906_independent.py
```

Normal and optimized runs are byte-identical and end with
`PASS checks=1376867`.

The proposed packet's own `audit.py` was also rerun both normally and with
`-O`; those two transcripts are byte-identical, and their LF-normalization is
exactly the frozen proposed `audit.out`.

Clean-room raw-LF hashes:

```text
04-computation/lrc14_72_539_weighted_profile_shadow_20260906_independent.py
  0063ebf7fb63743323d849fc0fb24a8cbbe7991e451f30ceb3126a6e7e5d4485
05-knowledge/results/lrc14_72_539_weighted_profile_shadow_20260906_independent.out
  64a3d57de8c2dfe54915e3599e9198793048159c073107d8488d53b8e124ad6d
```

Primary-packet hashes independently reproduced:

```text
04-computation/lrc14_72_539_weighted_profile_shadow_20260906.py
  143cca9032ae9690ff33069d93a1011ad246027f1c9af2677c1078f00a0be91f
05-knowledge/results/lrc14_72_539_weighted_profile_shadow_20260906.out
  439843af0b3c9c7496c2419c50ed7576673524136024074c70760efd508f8acb
```

LRC(14) remains **OPEN**.

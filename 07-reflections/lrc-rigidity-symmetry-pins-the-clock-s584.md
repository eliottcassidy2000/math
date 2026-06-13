---
source: opus-2026-06-03-S584 (remote-control)
status: STRUCTURE — rigidity in LRC: the symmetry⟺witness-rigidity duality; the AP is lonely at j/n iff gcd(j,n)=1 (a rigid unit-orbit witness); local pinch rigidity at each clock point cascades globally via (Z/n)*; the non-units are the 2-adic seam
tags: [LRC, rigidity, symmetry, units, clock, pinch, cascade, fixed-point, duality, 2-adic, n14]
---

# Rigidity in LRC: symmetry pins the clock

**Prompt (user):** see where rigidity appears; local rigidity around a fixed point;
global rigidity that cascades, permeating between isomorphic copies in symmetric objects.

Rigidity appears as a **duality**: the *most combinatorially symmetric* config (the AP)
has the *most geometrically rigid* witness — and the symmetry group **cascades** the
local pinch-rigidity at one clock point across the whole witness orbit. The core fact is
exact and clean.

## 1. Where rigidity appears (the survey)

| rigidity | object | flexible ↔ rigid |
|---|---|---|
| **witness** (geometric) | the safe-set `S(V)` | positive-measure interval (loose) ↔ measure-0 point set (worry-set) |
| **pinch** (local) | the binding pair at `t*` | one binder (slidable) ↔ a straddle pinch (pins `t*`, `M`) |
| **combinatorial** | the runner tournament `Aut` | trivial `Aut` (asymmetric) ↔ vertex-transitive (regular, worry-set) |

The three line up by a duality (below): symmetric ⟺ rigid-witness.

## 2. The exact core: the AP is lonely at `j/n` iff `gcd(j,n)=1`

> **Lemma (rigid unit-orbit witness).** For the AP `{1,…,n-1}`, `t=j/n` is a witness
> (`‖vj/n‖ ≥ 1/n ∀v`) **iff `gcd(j,n)=1`.**
> *Proof.* `‖vj/n‖ < 1/n ⟺ n | vj`; some `v∈{1,…,n-1}` has `n|vj` ⟺ `n/gcd(n,j) ≤ n-1`
> ⟺ `gcd(n,j) > 1`. So no runner hits `0` iff `j` is a unit. ∎

So the AP's witness set is **exactly the `φ(n)` unit clock points** `(ℤ/n)^* · (1/n)` —
a finite, measure-zero, **rigid** orbit (n=14: `{1,3,5,9,11,13}/14`). The non-units
`{2,4,6,7,8,10,12}` *fail* — a runner lands on the observer there — and they are exactly
the residues sharing a factor with `n=2·7` (the **2-adic / composite seam**, S580).

## 3. Local rigidity around a fixed point

Each clock witness `t=u/n` is **pinned by a straddle pinch**: the binding pair
`(u^{-1}·1, u^{-1}·(n-1)) ≡ (a, n-a)` (sum `n`) has two runners moving *oppositely* as
`t` varies, so `t=u/n` is a strict local max of `min_i‖v_i t‖` — a rigid critical point,
not slidable. The pinch fixes both `t = m/n` and the value `M = 1/n` (S557: `r/s`,
`s=n`, `r=1`). *Local rigidity around the fixed point is the straddle pinch.*

## 4. Global rigidity that cascades (the symmetry permeates)

The AP is **unit-invariant**: `u·V ≡ V (mod n)` for every `u ∈ (ℤ/n)^*`. Consequence:

> **Scaling-invariance.** If `uV≡V`, the safe set `S(V)` is invariant under `t ↦ u^{-1}t`
> (`‖v·u^{-1}t‖` ranges over the same multiset as `‖vt‖` since `uV=V`). So the symmetry
> group acts on the witness set.

For the AP this action is **simply transitive on the unit clock orbit** (verified:
`{j : gcd(j,n)=1}` is permuted by every unit). So **the local pinch-rigidity at `t=1/n`
cascades to all `φ(n)` witnesses** — rigidity *permeates* through the orbit; the whole
witness set is one rigid `(ℤ/n)^*`-orbit. It cascades further **between isomorphic
copies** up the doubling tower (S579): `D:v↦2v` carries `AP_n`'s rigid clock into the
even layer of `AP_{2n}` (`t=1/n ↦ 1/2n`), so the rigid fixed point propagates across
levels.

## 5. The duality, verified

`clock-rigid %` (fraction of configs whose optimum sits on the clock) **rises
monotonically with the symmetry-group size `|G(V)|`**:
```
n=8:   |G|=1 → 7%,   |G|=2 → 13%,   |G|=4 → 24%
n=10:  |G|=1 → 8%,   |G|=2 → 22%,   |G|=4 → 36%
```
and the AP (`|G| = φ(n)`, the *full* unit group) is **always** `M=δ` with a rigid
clock witness. **Combinatorial symmetry ⟺ geometric witness-rigidity.** The worry-set is
the *maximally symmetric* tournament (the regular/rotational, S582) and therefore the
*maximally rigid* witness — the two are the same wall seen from opposite sides.

## 6. The rigidity proof-lens on LRC

> **Symmetry pins `M`.** A unit-invariant config's loneliness is *forced* onto the rigid
> unit-clock orbit, where `M = 1/n` exactly. It cannot dip below: the orbit value is
> pinned at `1/n`, and the only way to lose loneliness is for a witness `u/n` to fail —
> which happens iff `u` is a **non-unit**, i.e. `gcd(u,n)>1`, i.e. a runner is a multiple
> of `n/gcd` — exactly **C′** (a multiple of `n` breaks the clock). So:
> - **symmetric & all-clock-units-survive ⟹ `M=1/n`** (rigid, tight, lonely);
> - the only breakage is a **non-unit** clock point = the composite/2-adic seam (`n=14`:
>   the `2` and `7` directions) = the residual C′.

Rigidity therefore *unifies* the threads: the worry-set's symmetry makes its witness a
rigid unit-orbit; loneliness fails only at the non-unit (composite) directions, which is
C′ / the `2q` apex / the doubling degeneracy. **n=14's hardness is the non-unit
(`2,7`) breakage of an otherwise-rigid clock.**

## 7. Honest status

- **Exact/proved:** AP lonely at `j/n ⟺ gcd(j,n)=1`; the witness set is the unit clock
  orbit; scaling-invariance of `S(V)` under symmetry units; the straddle-pinch local
  rigidity.
- **Verified:** the symmetry⟺rigidity duality (clock-rigid% monotone in `|G|`, n=6..14);
  the AP's full-unit symmetry and rigid `M=δ`.
- **Lens (directional):** "symmetry pins `M`, breakage only at non-units = C′" is a
  reframing of C′ through rigidity, not an independent proof; the non-unit (composite)
  residual is the same open core.

**Artifacts:** `04-computation/lrc_rigidity_s584.py` (+`.out`),
`lrc_rigidity_scaling_s584.out`. Builds on S557 (pinch), S582 (regular = worry-set), S579
(doubling tower), THM-398 (C′), S580 (units/doubling). New: **HYP-2124**.

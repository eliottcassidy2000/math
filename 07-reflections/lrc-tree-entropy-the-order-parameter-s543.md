---
source: oracle-2026-06-01-S543
status: result (box-entropy of loneliness on the base-p tree; phase-transition order parameter; verified)
tags: [lonely-runner, entropy, box-dimension, p-adic-tree, phase-transition, order-parameter, odometer, thermodynamic-formalism]
---

# Entropy on the tree: the loneliness box-dimension is the order parameter, and the tight AP is the critical point

**Prompt (user):** attack entropy on the tree vigorously.

The recursive descent on the p-adic / base-`p` tree (a depth-`d` node = a base-`p`
interval `[j/p^d,(j+1)/p^d)`; the child map is `x -> px`, entropy `log p`) turns LRC
into a thermodynamic statement. The key structural fact up front: the **center =
shift** (S541, the `+1` odometer on `Z_p`) is **equicontinuous — zero entropy**. So
LRC is *blind* to the entropy of the shift; all the entropy lives **transverse** to
the center, in the `x -> p` descent. We attack *that*.

## Two box-entropies of the lonely set

Let `S = { t : min_i ||v_i t|| >= 1/n }` (the lonely/safe set). On the base-`p` tree,

```
N_touch(d) = #depth-d nodes meeting S        h_touch = lim (1/d) log_p N_touch
N_full(d)  = #depth-d nodes inside S          h_full  = lim (1/d) log_p N_full
```

**Verified (`lrc_tree_entropy_s543.py`, p=2,3, depth 7):**

| system | `|S|` | `h_full` | `h_touch` |
|---|---|---|---|
| generic (n=4,5,6) | `>0` | → **1** | → **1** |
| **tight AP** (n=4,5,6) | **0** | **−∞** (`N_full≡0`) | → **0** (`N_touch` BOUNDED) |

The generic lonely set has **box-dimension 1 in every base** (positive measure); the
tight AP's lonely set is a **finite set of wall points** — `N_touch` saturates
(`[2,2,2,...]` / `[2,4,4,4,...]`), so `h_touch = 0`, `h_full = -inf`. The dimension
is **base-`p` independent** (it is the archimedean transverse dimension; the tree is
just the counting scaffold, and the odometer contributes nothing).

So:

> **`h_touch` is the LRC order parameter.** `h_touch = 1`: lonely with positive
> measure (generic). `h_touch = 0`: lonely only on a finite wall-set (the tight
> extremal). `h_touch = -inf` (`N_touch -> 0`): NO lonely time — an LRC violation,
> never observed. **LRC ⟺ `h_touch >= 0`**, and the extremal AP sits exactly at the
> `h_touch: 1 -> 0` drop.

## The phase transition: dimension `D(theta)` and the critical threshold

Sweep the threshold `theta` (instead of fixing `1/n`). The lonely set
`S_theta = {min_i ||v_i t|| >= theta}` has box-dimension

```
D(theta) = 1     for  theta < theta_c        (positive measure)
         = 0     at   theta = theta_c        (the maximiser(s) only)
         = -inf  for  theta > theta_c        (empty),
```

a step function with critical point `theta_c = max_t min_i ||v_i t||` — exactly the
**loneliness radius / max collar** (S541). Verified:

```
 AP n=4,5,6 (tight):  theta_c = 1/n  EXACTLY     (|S|=0 at 1/n; dim drops to 0 there)
 generic:             theta_c > 1/n (often 1/2)  (|S|>0 at 1/n; dim 1)
```

So LRC is a **critical-threshold statement of an entropy phase transition**:

> **LRC(n) ⟺ `theta_c(v) >= 1/n` for every speed system** ⟺ the dimension-1
> (positive-entropy) phase of the lonely set reaches the threshold `1/n`. The
> extremal AP is the system whose critical point lands **exactly on** the threshold
> (`theta_c = 1/n`) — the conjecture is tight precisely at the phase boundary.

This reframes everything we have: the "fat collar" (S541) `= theta_c`; "wall-only /
tight" (S525) `=` the zero-entropy critical line `D=0`; "0 failures" (S520/S526) `=`
the dimension-1 phase always reaching `1/n`.

## Where the p-adic tree genuinely enters

The box-dimension is base-independent, but the **fine structure at the critical
threshold is `p`-adic for `p | n`**. At `theta = theta_c` the lonely set degenerates
to the wall points `t = a/n` and their descendants; these are exactly the S410
**moat endpoints** with `p`-adic depth `v_p(n)`, and the "debt export to child
vertices" is the self-similar generation of the critical set down the `p | n` tree
(S542). So: the **bulk** entropy (dim 1 vs 0) is archimedean and base-free; the
**critical set** (at `theta_c`, where LRC is decided) is the `p`-adic moat on the
`p | n` tree. The conjecture lives at the meeting of the two — the archimedean phase
boundary sitting on the `p`-adic moat.

## Vigorous summary / open (→ HYP)

- LRC `= ` "the lonely-set box-dimension is `1` at threshold `1/n`," i.e. the order
  parameter `h_touch(1/n) >= 0`, with the tight cases at the `1->0` critical drop.
- The center (odometer) has zero entropy: LRC's entropy is purely transverse, so a
  **thermodynamic formalism** (pressure `P(beta)` for the potential `-beta * log`
  collar over the `x->p` system) should have its phase transition at
  `beta_c <-> theta_c = 1/n` exactly at the AP. Computing the pressure/free-energy
  and its singularity at the AP is the next vigorous step.
- The critical-set multifractal spectrum at `theta_c` should be supported on the
  `p | n` moat (S410/S542) — a `p`-adic measure whose dimension encodes how tightly
  the AP fails to clear the threshold.

## Anchor
`04-computation/lrc_tree_entropy_s543.py` (+ `.out`): `h_touch` = 1 (generic) / 0
(tight) in p=2,3; `theta_c` phase transition (AP at `theta_c=1/n` exactly). Builds
on S541 (collar/center=odometer), S542 (p-adic tree/moat), S525 (wall-only),
S526/S520 (covering/0-failures).

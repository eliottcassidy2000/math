# The three axes of Erdős 592: tree-grid (tame/linear), Chang tower (wild/ordinal), and the bi-dyadic atom that bridges them

**Source:** kind-pasteur-2026-06-16-S3. Dispatch (terse): "think tree-grid dichotomy, Chang
tower-Ramsey, t-uniform bi-dyadic." These are three pieces of repo vocabulary from the
**Erdős 592** thread (the $1000 problem: which countable `α=ω^β` satisfy `α→(α,3)²`; smallest
open case `α=ω^(ω³)`). This reflection reads them as the **three coordinate axes** of the
finite calculus the lab has built (THM-453/460/464/465/469/471, HYP-2396), adds a framing that
organizes them, and re-verifies the keystone seam.

## The three axes

| axis | repo term | object | growth | status |
|---|---|---|---|---|
| **n** (CNF exponent, `ω^n`) | **tree-grid** | `R(n,2)` = least `t` with the `Q(n,t)` game UNSAT (triangle-free graph on the `t`-ary height-`n` tree leaves, hitting every binary subgrid) | **LINEAR**: `R(1,2)=3, R(2,2)=5`, conj. `R(n,2)=2n+1` (HYP-2396) | tame |
| **m** (tower rank, `ω^(ω^m)`) | **Chang tower** | `ω^(ω^m)→(·,3)²` | **WILD/tower** | `m=1` YES (Chang 1972), `m=2` YES (Schipperus), **`m=3` OPEN ($1000)**, `m≥4` NO |
| **b** (branching / grading) | **t-uniform bi-dyadic** | the binary (`b=2`) subgrid; `R_b(1)=R(3,b)` (THM-464) | the obstruction-**atom** | the `b=2` is *forced* (below) |

## The keystone: why "bi-dyadic" — the b=2 is forced by a unique-prime fact

The obstruction in the tree-grid game is the **binary** subgrid (`2^n` leaves), not a `b`-ary
one for `b>2`. Why 2? THM-469's core, re-verified here (`erdos592_three_axes_kps.py`): the
`v_p`-valuation level-sets of `ℕ` are **sum-free if and only if `p=2`** (`p=2`: every level is
sum-free since odd+odd=even raises the valuation; every odd `p`: level 0 already fails,
`1+1=2`). The game invariant is a *sum-free grading* (THM-469), and **`p=2` is the unique
prime supporting one** — so the grading, hence the subgrid, hence the whole obstruction, is
**binary**. "Bi-dyadic" = this `2` appearing twice: the **binary** subgrid (the geometry) *and*
the **2-adic** sum-free grading (the algebra) are the same `2`. The "t-uniform" qualifier is
the branching parameter `t` of `Q(n,t)` (the `b`-ary game `Q_b`), with the algebra living on
the `b=2` rung of the same dyadic ladder as the Sidon/`B_h` rungs (THM-446).

## The new framing: PRE-TOWER (n) vs TOWER (m), and the binary atom carried across

The two laddered axes look incomparable — one **linear** (`2n+1`), one a **tower** — but they
are the same problem on either side of the **first limit exponent**:

- **`ω^n`** (finite CNF exponent, `n<ω`) is the **pre-tower** regime. The CNF of the exponent
  `n` is a finite, limit-free number, so no genuine tower can form; the binary subgrid only
  needs **2 branches per level**, and the threshold is **linear** (`R(n,2)=2n+1`). Dodging a
  `2^n` grid in a `t`-ary height-`n` tree is a Kővári–Sós–Turán / Zarankiewicz-type density
  fight per level pair, burning ≈2 branches per level ⟹ `t≈2n` — the heuristic behind `2n+1`.
- **`ω^(ω^m)`** crosses the limit exponent `ω` into the **tower** regime. Now the CNF exponent
  `ω^m` is itself transfinite, and full-type sets are **stacked towers** (THM-460 B2). Chang's
  theorem (`m=1`) is the **base** of the tower; each additional CNF summand of `γ` (in
  `ω^(ω^γ)`) stacks one more level. The open `m=3` stacks **three** binary-subgrid-tower
  levels (THM-460 B3's shape grammar `BinM`: pairs at scale `e=0`, `M`-towers at limit scales).

> **The binary subgrid (`b=2`) is the universal obstruction-atom at every rank** — pairs in
> the pre-tower, `M`-towers-of-pairs across the limit — and it is *carried across the first
> limit exponent intact* because the `p=2` sum-free seam is what makes "binary" the right
> grading at all scales. So the same dyadic atom that gives the tame linear `2n+1` on the
> `n`-axis is what the wild `m`-tower is built by stacking.

This reframes the difficulty of `m=3`: it is **not** that a new obstruction appears at the
third summand, but that the *third stacking* of the same binary atom is the first the positive
machinery (Chang/Schipperus's Ramsey-over-forms uniformization, which works for `≤2` summands)
cannot uniformize — exactly where the lab's finite probes explode (THM-471: `m=2,M=1` already
UNSAT at `(2,2)`; the `m=3` enumerator's `j1`-march hits 3.5M leaves).

## A concrete consequence and a proof-strategy

If `R(n,2)=2n+1` for all `n` (HYP-2396), then `Q(n,t)` is UNSAT for `t≥2n+1`, so **strong
witnesses never exist** on the `n`-axis (a strong witness needs `Q(n,t)` SAT for *all* `t`).
The known negative relations (Specker `ω^n↛`, `n≥3`) must therefore use **non-strong**
witnesses — and the whole content of Erdős 592 lives in the gap between strong and non-strong
witnesses. So the `2n+1` conjecture is the precise statement "the tree-grid n-axis is tame; all
the difficulty is on the m-tower." **Proof strategy** (handoff): the linear lower bound is a
binary-subgrid construction (`Q(n,2n)` SAT, witnessed `n≤3`); the linear upper bound
(`Q(n,2n+1)` UNSAT) should be a Kővári–Sós–Turán count — every `t`-ary level pair forces a
`K_{2,2}`-density that triangle-freeness cannot sustain once `t>2n`. Settling it would close
HYP-2363 (strong witnesses) on the finite side.

## Status / honesty

This is **synthesis + framing**, not a solution to the $1000 problem (`m=3` stays open). The
keystone (`v_p` level-sets sum-free ⟺ `p=2`) is re-verified here and is THM-469's. The
three-axis organization and the **pre-tower/tower transition at the first limit exponent + the
binary atom carried across** are the new content (a navigation frame for a sprawling thread).
The `2n+1` heuristic and KST proof-strategy are a flagged route, not a proof. "Chang" =
C.C. Chang's 1972 theorem `ω^ω→(ω^ω,3)²` (the `m=1` base), per the repo's "first Chang number"
usage — *not* the set-theoretic Chang's conjecture. Cross-links: THM-453 (the tree-grid
ladder + tournament dictionary), THM-460 (the tower miniature / Chang shadow), THM-464/465
(the `b`-ary game + bi-dyadic witnesses), THM-469 (sum-free grading ⟺ `p=2`), THM-471 (the
shape enumerator), HYP-2396 (`R(n,2)=2n+1`), HYP-2556 (this synthesis), [[triangle_foundation]],
THM-446 (the dyadic `B_h` rungs). Source: Erdős Problem 592; Specker 1957; Chang 1972;
Schipperus 1999/2010; Džamonja–Koutsoukou-Argyraki–Paulson arXiv:2011.13218.

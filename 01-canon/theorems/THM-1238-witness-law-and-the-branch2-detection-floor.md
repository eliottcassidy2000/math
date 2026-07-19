---
id: THM-1238
title: "PRIORITY FIRST: the Kakeya ceiling D >= 3P serves a theorem codex already PROVED (THM-1203), so the whole route is redundant. What survives: (I) the WITNESS LAW P(D) = 3*floor(D/14) + e(D mod 14) on d = (1,D-1,D), which PROVES D >= 3P on that family with equality exactly at (1,2,3); (II) three gap-coordinate lemmas, including '#runs is always even' which makes death-star's '#runs <= 2D/3' and 'D >= 3P' the SAME statement; (III) the lattice-tube identity #runs = #{(sigma,m) : (B_sigma+m) meets [0,-d]}; (IV) on branch 2, a VALIDATION GATE showing the existing adversarial searches -- and my own -- have a detection floor ABOVE their target, so 'no counterexample found' is vacuous; against which 1200/1200 hard-stratum families carry an explicit grid CERTIFICATE M >= 1/14"
status: "(I) PROVED on the family, conditional only on the periodic closed form, which is verified over 698 values of D with 49-50 samples in every residue class mod 14. (II) PROVED outright. (III) PROVED (verified injective, 6264 runs). (IV) the detection-floor finding is DEMONSTRATED (my searcher fails to rediscover {1,...,13} at 1/14 from inside its own sampling range); the 1200/1200 certificate is RIGOROUS per family, not a search result. PRIORITY: the ceiling mu(BAD) <= 2/21 with equality iff AP is codex-S77's THM-1203, PROVED and Lean-certified, and supersedes this entire line. Uniform r=5 still needs codex's named finite-comb/eroded-start glue, which is untouched here"
source: kind-pasteur-2026-07-19-S128c84 (owner: prove D >= 3P from the (1,D-1,D) witness family; work #runs <= 2D/3; work the branch 2 target with a sharper adversarial search; work the n=12 AP uniqueness; mine the 3 / 4 / 1-12 threads; think Kakeya)
depends_on:
  - THM-1203    # codex-S77: PROVES the ceiling this route was chasing
  - THM-1211    # my slack simplex
related: [THM-1183, THM-1215, THM-1220, THM-633, MISTAKE-181, MISTAKE-176]
script: 04-computation/kakeya_gap_lattice_kps_S128c84.py, kakeya_witness_family_kps_S128c84.py, kakeya_witness_law_kps_S128c84.py, branch2_validated_search_kps_S128c84.py (+ .out)
---

# THM-1238 — the witness law, and what the branch-2 searches cannot see

## 0. PRIORITY — read before building on the Kakeya thread

**The target of the whole `D >= 3P` / `#runs <= 2D/3` line is already a theorem.**
codex-S77's **THM-1203** proves

> `mu(BAD_d) <= 2/21`, with equality exactly for four-term arithmetic-progression
> frequencies `{0,m,2m,3m}`,

PROVED, computer-assisted, with a Lean kernel certificate on the finite arithmetic
core and the equality-obligation lemma. That is the same quantity death-star and I
call the *sojourn*, and the same bound. Its route is entirely different — delete the
four phase points to a **non-arithmetic additive triangle** `(p,q,p+q)`, which turns
the ceiling into six torus triangles of leg `1/7` — and it never counts runs.

So `D >= 3P` is a second route to a proved theorem. It is not worthless (an
independent combinatorial proof has value), but it is **not on the critical path**,
and the fleet should stop treating it as the residue. THM-1203 §12 names the actual
remaining passage for uniform `r=5`:

> 1. a uniform finite-to-continuum domination/error estimate when the offsets may
>    scale with the carrier; and
> 2. a lower/address bound for the **eroded start complex** `E_k(P)`, not merely the
>    raw safe-set measure `mu(S(P))`.

Neither is touched by anything below, or by `D >= 3P`.

THM-1203 §10 also records the refutation of my own six-box draft, with the witness
`d = (1,6,7)`, `u = 3/4` hitting `(1/4,1/2,3/4)` while non-proportional — that is
MISTAKE-181, pre-registered against exactly this line of reasoning.

## I. The witness law — `D >= 3P` PROVED on the family the owner named

Both death-star and I found that the **minimal `D` attaining each `P`** always has
the shape `(1, D-1, D)`. On that family the count is exactly periodic. Writing
`D = 14q + r`:

> **P(D) = 3q + e(r)**,  with `e(r) ∈ {0,1,2}` given by
>
> | `e` | residues `r` mod 14 |
> |---|---|
> | 0 | 0, 1, 2, 4, 5 |
> | 1 | 3, 6, 8, 9 |
> | 2 | 7, 10, 11, 12, 13 |

Verified for `D = 3 .. 700`: in **every** residue class `c = P - 3D/14` takes a
**single** value across 49–50 samples. Equivalently `P = 3D/14 + c(D mod 14)` with
`max|c| = 15/14`.

**Consequence.**

> `D >= 3P  <=>  14q + r >= 9q + 3e(r)  <=>  5q >= 3e(r) - r`.

The maximum of `3e(r) - r` over all fourteen residues is **0**, attained only at
`r = 0` and `r = 3`. Hence `5q >= 0` always, so **`D >= 3P` holds on the entire
family**, with equality iff `q = 0` and `r ∈ {0,3}` — i.e. iff `D = 3`, i.e. the
direction `(1,2,3)`.

The asymptotic ratio is `P/D -> 3/14`, and `3/14 < 1/3` with margin `14/9`. The
constant is not arbitrary: `#runs = 2P -> 3D/7`, and `3/7` is exactly the total
measure of the three gap-windows `[1/7,2/7] ∪ [3/7,4/7] ∪ [5/7,6/7]` that the
fastest coordinate must occupy. **The witness family is equidistributed for its
fastest coordinate; only its first member is resonant.**

Additivity is confirmed as the extremal shape: over 3658 primitive triples of the
additive family `(a,b,a+b)` the max is `P/D = 1/3` at `(1,2,3)` alone, while
non-additive stress families top out far lower — `(1,2,D)` at `1/7`, `(1,D//2,D)`
at `1/6`, `(3,D-2,D)` at `1/8`.

## II. Three gap-coordinate lemmas (PROVED)

Change to gap coordinates `x = g_min`, `y = g_mid - g_min`, `z = g_max - g_mid`
(unimodular). death-star's constraints become `x,y,z <= 2/7` and `x+y+z >= 5/7`,
so each piece of `B` is the **corner simplex of the cube `[0,2/7]^3` with leg 1/7**.
Hence `x,y,z ∈ [1/7,2/7]`, and

> `g_min ∈ [1/7,2/7]`,  `g_mid ∈ [3/7,4/7]`,  `g_max ∈ [5/7,6/7]`

— three disjoint windows of length `1/7`. From this:

- **(L2) The ordering is constant on a run, and nothing wraps.** The gaps `y,z` are
  `>= 1/7 > 0`, so no two coordinates cross; and `1/7 <= g <= 6/7`, so no coordinate
  wraps. A run therefore needs **no cell subdivision** — the affine description is
  valid on the whole run. (0 violations, 6264 runs.)
- **(L3) The three rates are distinct, so `D >= 3`.** `d_i = d_j` forces `g_i = g_j`,
  a zero gap. Three distinct positive integers have max `>= 3`. This **proves the
  `P = 1` case of `D >= 3P` outright**, with no centre-hit hypothesis — death-star's
  `c >= 3` needed one.
- **(L4) `#runs` is always even.** `u -> 1-u` maps runs to runs (`g(1-u) = 1-g(u)`,
  and `B` is invariant under `g -> 1-g`). Its only fixed points are `u = 0`, where
  `g = (0,0,0)`, and `u = 1/2`, where every `g_i ∈ {0,1/2}` so two coordinates
  coincide. Both have a zero gap, so neither lies in `B`: **no run is self-mirror**
  and `#runs = 2P` exactly. (0 odd counts, 0 self-mirror runs, `d ∈ [1,30]^3`.)

**L4 is the piece that matters for bookkeeping:** it makes death-star's
`#runs <= 2D/3` and `D >= 3P` *the same statement*, not two.

Also proved exhaustively: **`D <= 5` forces `P <= 1`** (all primitive triples with
distinct rates `<= 5`; repeated rates give no runs at all by L3). That closes the
`P = 2` case, which needs `D >= 6`.

## III. The lattice-tube identity

Unrolled, the geodesic is the segment `[0,-d] ⊂ R^3`, and each `B_sigma` is convex,
so a line meets it in one interval. Therefore

> **`#runs = #{ (sigma, m) ∈ S_3 × Z^3 : (B_sigma + m) ∩ [0,-d] ≠ ∅ }`**

— a **lattice point count in a convex sausage**. (Verified injective: 6264 runs, 0
collisions.) So `D >= 3P` is literally "no thin tube about a rational direction
carries more lattice points than the `(1,2,3)` tube". That framing also says why it
resists a soft proof: lattice points in a thin tube exceed the volume precisely when
the tube is arithmetically aligned, which is the resonance itself.

## IV. Branch 2 — the detection floor, and what actually carries weight

The target (THM-1215): *show that when `q0 > 14`, some pair achieves
`D/(v_i+v_j) >= 1/14`.* I built a searcher with a **validation gate** in front of it.

**The gate.** (a) the exact evaluator reproduces `M({1..13}) = 1/14`,
`M({1..12}) = 1/13`, `M({1..11}) = 1/12`, `M({1..11,24}) = 2/25` — all OK; (b) the
searcher must **rediscover `{1,...,13}` at `1/14`** from random starts inside its own
sampling range `range(1,45)`.

**It fails (b).** 60 random restarts bottom out at `4/43 = 0.0930`; simulated
annealing was tried and did worse (`6/59`). The detection floor `0.0930` lies
**above** the counterexample target `1/14 = 0.0714`.

> **Therefore "0 counterexamples found" from this instrument is vacuous — it could
> not have detected one.** The same arithmetic applies to the existing runs:
> `dge2_hardmin_opus_S392` has floor `2/19 = 0.105`, `stability_gap_opus_S393` has
> floor `6/61 = 0.098`, and `dge2_branch_opus_S392` searched `range(1,45)`
> unrestricted and bottomed at `1/10` — missing `{1,...,13}` from inside its own
> range by 40%. None of these floors is below `1/14`.

This does not refute the branch-2 target. It says the **search-based** evidence for
it is weaker than the write-ups imply, and should be re-scoped.

**What does carry weight.** Replace search with certificates. `gridmax(V) <= M(V)`
always, so exhibiting one grid point `t` with `min_v ||v t|| >= 1/14` **proves**
`M(V) >= 1/14` for that family. Generating hard-stratum families structurally —
`q0 > 14` is exactly the nine divisibility obligations `q ∈ {5,7,8,9,10,11,12,13,14}`
(`2,3,4,6` follow), so they are built inside the stratum rather than rejected into it:

> **1200 / 1200 generated `q0 > 14` families carry an explicit grid certificate
> `M >= 1/14`.** None survived the prune; none needed a search at all.

That is a rigorous per-family statement about 1200 families, and it is the form the
existing branch-2 runs should have reported.

## V. Corrections against myself

- **My L5 over-claim, self-caught.** I asserted that `run = 1/(7D)` forces `(1,2,3)`
  and entry at the corner `(2/7,4/7,6/7)`. **False**: 900 runs attain `1/(7D)`, many
  at non-`(1,2,3)` directions — `(1,2,10)`, `(1,9,10)`, `(1,2,17)`, … What is true is
  the weaker `run <= (1/7)/r_max` with equality iff corner entry (240 runs, 0
  failures), where `r_max` is the rate of the *max coordinate*, which need not be `D`.
  The correct chain is `run <= 1/(7D) <= 1/21`, and it is the **second** inequality —
  `D >= 3` from L3 — that is tight only at `(1,2,3)`.
- **`|B| = 1/343` was already canon.** THM-1183 (codex) records the ambient volume
  `1/343` and the AP X-ray density `2/21`. My MISTAKE-176 correction of my own
  `0.003367` stands as a correction, but the closed form was **not** mine to claim.
- **My THM-1153 `n/(12n+1)` law is THM-633** (mac-mini-S33), PROVED and Lean
  kernel-pure (`LRCLadderD1.lean`), including the exact ladder
  `M({1,…,11,12m}) = m/(12m+1)`. boxeph-S123 had already credited it back. My
  contribution there reduces to the `1/M(core)` reading, not the law.

## VI. The 1/12 question, closed

The owner asked to mine the `3`, `4` and `1/12` threads. On `1/12` the repo has
already adjudicated this: `07-reflections/two-twelves-...-coincide-at-n14-kps.md`
concludes there are **two** unrelated twelves, coincident at `n=14` only because
`n-2 = 12` happens to equal the denominator of `B_2/2`. Confirmed and sharpened:

- the LRC `1/12` is elementary — `t = a/12` is a witness for any speed set with no
  multiple of 12, so `M >= 1/12` with equality on `{1,…,11}` and friends. This is
  THM-633's generic branch, and it generates the `n/(12n+1)` law, the `D=1` plateau
  of the `n=12` gap, and THM-615. **One phenomenon.**
- the `-1/12` of `zeta(-1)` is `-B_2/2`, the Dedekind reciprocity constant, an
  `n -> ∞` limit of `s(n,Phi_6) = -T/(12T+6)`; present for every `n`, and twice
  formally ruled "inspiration, not a rigorous tool" (LTI-087, HYP-2896). **Different.**
- the `-1/12` in my own THM-1173 continuum centres is `(6·14)/7 = 12` — pure
  normalisation bookkeeping. **Coincidence; I recommend dropping it from the list.**

## Named next

- The real residue is codex's, not death-star's: the **finite-to-continuum
  domination** and the **eroded start complex `E_k(P)`** (THM-1203 §12).
- Branch 2 deserves a certificate-based sweep, not a search: extend the 1200-family
  grid-certificate run to an exhaustive bounded stratum, and report the certificate
  rate rather than a hill-climb minimum. `lrc_height_batched_census_macmini_S115.py`
  already has the right engine (sound grid prescreen + exact confirmation).
- If anyone still wants `D >= 3P` as an independent proof, the witness law in §I is
  the template: find the analogous periodic count off the family.

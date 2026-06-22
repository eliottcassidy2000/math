# Node 1 progress: an elementary finite-V closure of the pure-dilation coordinated-growth core

*kind-mendel-2026-06-22-S3. A "real swing" at Open Node 1 (the finite-V_max Part A correction),
following the S1 audit (HYP-2847) and the S2 unification (decorrelation). Re-derives and sharpens
THM-523/THM-518/THM-525 for the specific "boundary core" the literature points to, with explicit
constants and a second argument that closes the regime the measure bound misses. All exact-verified
(`04-computation/lrc14_node1_decoupling_kindmendel.py`, `..._sheets_kindmendel.py`).*

## The target

The kps-S4 reflection identifies the hard residual as the **coordinated-growth core**
`{t,2t,…,12t,V}` and its CRT relatives — "where the cluster scale runs to infinity and the
margin sits exactly on the boundary." Its purest form (write `b` for `t`):
> `S = {b, 2b, …, 12b, V}`,  `V ≡ 0 (mod 14)` (the covering-forced parked runner), `gcd(b,V)=1`.

LRC(14) for `S` is `M(S) ≥ 1/14`, i.e. the lonely set `L(S)={τ:‖sτ‖≥1/14 ∀s∈S}` is nonempty.

**Constants (exact, this session).** Substituting `u=bτ mod 1`, the dilated part `{b,…,12b}` is
lonely iff `u ∈ G_12 := {u:‖ju‖≥1/14, j=1..12}`, a fixed scale-invariant set:
> `μ_12 := meas(G_12) = 6617/194040 ≈ 0.034101`, with `r_12 = 12` arcs; widest arc `=[1/14, ~1/13]`.

`{b,…,12b}` therefore has lonely measure `μ_12` (scale-invariant, verified) and exactly `12b` arcs.

## Proof 1 — the decoupling/comb floor (measure bound; needs `V ≫ b`)

`L(S) = meas(B ∖ D_V)`, `B={τ:bτ∈G_12}` (measure `μ_12`, `≤ b·r_12` arcs), `D_V={‖Vτ‖<1/14}`
(period `1/V`, density `1/7`). **Comb-on-arc lemma** (elementary): any arc of length `ℓ` meets
`D_V` in `≤ ℓ/7 + 1/(7V)` (`≤ ℓV+1` teeth, each `≤1/(7V)`). Summing over the `≤b·r_12` arcs:
> `L(S) ≥ (6/7)μ_12 − b·r_12/(7V)`,  **`> 0` ⟺ `V > κ·b`,  `κ := r_12/(6μ_12) = 388080/6617 ≈ 58.65`.**

This is THM-523/THM-518's decoupling floor specialized to the dilation, with the threshold made
explicit. **Verified exactly** over `b∈{1,2,3}` and many `V`: the inequality holds in every case,
and the lower bound is positive exactly when `V > κb` (e.g. `b=1`: positive for `V≥70`; the limit
floor is `(6/7)μ_12 ≈ 0.02923`).

## Proof 2 — sheet-counting (closes the comparable regime `V ≤ κb` the comb bound misses)

The measure bound is vacuous for comparable `V` (e.g. `b=3,V=14`: bound `= −0.34`), yet `L(S) > 0`
there too. The fix uses the **multiplicative structure**: write `τ = (n+u)/b`, `n∈{0,…,b−1}` (the
"sheet"), `u∈[0,1)`. Then `‖jbτ‖=‖ju‖` (depends only on `u`) and `‖Vτ‖=‖(V/b)(n+u)‖`.

Fix the widest arc `[u_lo,u_hi]⊂G_12` (width `w=6617/.. ≈0.00595`). On sheet `n`, as `u` sweeps the
arc, `Vτ` sweeps an interval `I_n` of length `ℓ=Vw/b`, with left endpoint
`p_n = frac((Vn+Vu_lo)/b)`. Since `gcd(V,b)=1`, `{frac(Vn/b):n} = {0,1/b,…,(b−1)/b}` are **`b`
equally-spaced points**. `I_n` fails to meet the safe set `[1/14,13/14]` only if `I_n ⊂ (−1/14,1/14)`
(the length-`1/7` danger arc), i.e. `p_n` lies in a window of length `1/7−ℓ`; that window holds
`≤ (1/7−ℓ)b + 1` of the `b` points. Hence
> **#good sheets `≥ b − (1/7−ℓ)b − 1 = (6/7+ℓ)b − 1 ≥ (6/7)b − 1`.**

Each good sheet `n` gives a `τ=(n+u)/b` with `u∈G_12` (all 12 dilated runners safe) **and**
`‖Vτ‖≥1/14` ⇒ `τ` is lonely. Since `V≡0 (mod 14)` is even, primitivity `gcd(b,V)=1` forces `b`
**odd**, so `b≥3` and `(6/7)b−1 ≥ 11/7 > 1`: **always ≥ 1 good sheet.** `L(S) > 0` for every
primitive pure-dilation core. **Verified:** good-sheet counts are `≥ 6b/7−1` (and `≥2`) for all
tested `b∈{3,5,13,…}`, `V` coprime.

## What this does and does not settle

**Settles (elementary, exact-verified):** the **pure-dilation** coordinated-growth core
`{b,…,12b,V}` is lonely for *all* primitive `(b,V)` — for large `V` with a quantitative measure
floor (Proof 1), and for comparable `V` by sheet-counting (Proof 2). The "hard core" of the
literature, in its purest form, is provably **not** a counterexample, and we see *why*: 12 speeds
sharing the base `b` create `b` equally-spaced sheets, of which one extra runner can spoil only
`~1/7`.

**Does not settle (the real Node-1 residual, now sharply localized):** the **CRT relatives** — 12
near-AP speeds that are *not* exact multiples of a single `b` (e.g. perturbations `{b+δ_j}` or
multi-modulus CRT clusters). There the sheet offsets are no longer exactly equally spaced, the
continuous slow-fast version reappears, and the finite-V error is the genuine `O(arcCount/V)` term.
This matches the unification (S2): Proof 2's exact equal-spacing is the *discrete, lossless* version
of the equidistribution Node 1 needs; perturbing the cluster reintroduces the discrepancy.

**Reusable lever:** the explicit `κ = r_12/(6μ_12)` and the decoupling floor give a clean,
formalizable (Lean-ready, elementary) certificate for any single-dominant-far config with `V>κ·b`.
The remaining work is (a) replacing exact equal-spacing by an **effective-equidistribution** bound
(Erdős–Turán) for perturbed clusters, and (b) the `V ≲ b` perturbed regime, where the bounded-ratio
analysis (THM-405) / proven LRC(≤13) should take over.

→ HYP-2847, HYP-2848, OPEN-Q-108, THM-518, THM-523, THM-525, THM-405; reflection
`lrc14-three-open-nodes-are-one-decorrelation-kindmendel-S2.md`.

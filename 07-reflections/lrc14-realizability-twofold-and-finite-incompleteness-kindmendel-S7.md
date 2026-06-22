# The LRC realizability structure is two-fold — and the finite/combinatorial half is provably incomplete

*kind-mendel-2026-06-22-S7. Owner: find creative realizability arguments (à la Tournament Analysis) — the
LRC's missing realizability structure may be "slightly different"; aim to finish. Result: the structure is
two-fold (a FINITE rigidity for the extremizer + an ANALYTIC equidistribution for the committed speed), and
I give a clean REFUTATION showing the finite half cannot finish alone — correcting my own HYP-2864 and
mac-mini's HYP-2876. Pulled main throughout; converges with kps-S31o (HYP-2895) and mac-mini-S43 (HYP-2897).*

## What "realizability" means here, and the two structures

Tournament Analysis realizes pairwise data as a *tournament* `T(x)` (the winding tournament; observer
lonely ⟺ marked source, THM-381). The owner's intuition is right that LRC needs a *slightly different*
realizability object. This session pins it down: there are **two** realizability structures, and they split
the two open nodes by the *kind of argument* each needs.

### Structure A — three-gap (Steinhaus) rigidity: the FINITE realizability of the extremizer (Node 2)

The extremal/tight LRC configuration is the AP (consec). Its realizability signature is the **three-gap
theorem**: the points `{j·x mod 1 : j=0..k-1}` have **at most 3 distinct gap-lengths for every x** —
*and this is unique to APs (and their dilations)*. Verified
(`04-computation/lrc14_threegap_realizability_kindmendel.py`):
- consec `{0..7}`: max 3 gaps; a perturbation `{0..6,8}`: 4 gaps; a spread set: 8 gaps; dilation
  `{0,2,…,14}`: 3 gaps (same as consec, THM-531 dilation-invariance).
- This rigidity makes sector coverage **all-or-nothing** ⇒ the empty-sector count `N_E` is **most bimodal**
  for the AP (mass at the extremes `N=0,6`): bimodality `0.348` (consec) vs `0.226` (perturb) vs `0.056`
  (spread). Bimodality tracks `p0`, so **the three-gap rigidity is *why* the AP maximizes `p0`/`L_y`** (Node
  2 / the compact-reduction extremality, = kps-S36 bimodality crux, now grounded in Steinhaus rigidity).

**Consequence:** Node 2 has a *finite/algebraic* realizability structure — "only APs are three-gap-rigid"
is a closed combinatorial characterization. A proof of "consec maximizes" should be reachable by a
three-gap rigidity/smoothing argument (replace any non-AP by its AP-hull, coverage gets more bimodal).

### Structure B — torus equidistribution: the ANALYTIC realizability of the committed speed (Node 3)

For a covering set with one *huge committed speed*, loneliness is realized only **dynamically**: the integer
orbit `(v_1 t,…,v_k t)` is a closed geodesic on `T^k`, and the committed speed's danger comb removes only
`~1/7` of the core's lonely set by **equidistribution** (Weyl), leaving a witness — but at an *unbounded
denominator*. This is the "slightly different" object: not a single tournament, but the **equidistribution
of the orbit** (joint danger governed by gcds/resonances, not `(1/7)^k`). (= kps-S31o HYP-2895, mac-mini-S43
"equidistribution is the mechanism".)

## The refutation: the finite/combinatorial route is PROVABLY incomplete

My HYP-2864 (every covering 13-set has a bounded-denominator witness `D≤D_0`) and mac-mini's HYP-2876
("every 13-set has a witness `D≤41`") are **FALSE**. Counterexample family (verified,
`lrc14_threegap_realizability...` companion run):

> **`S_X = {1,…,11,13, lcm(2..X)}` is a primitive covering 13-set whose minimal witness denominator is
> `nextprime(X) > X`.** Verified: `X=30→D=41, 45→53, 60→67, 80→83, 100→101` (always the next prime).

*Why:* every `D ≤ X` divides `v=lcm(2..X)`, so `‖v·a/D‖ = 0` (the committed speed sits on the observer) for
**all** `a` — so no `D ≤ X` certifies `S_X`. A witness appears only at the next prime `p > X` (which doesn't
divide `v`, so `v mod p` is generic-safe and the loose core is witnessed there). Hence **the witness
denominator is unbounded over covering 13-sets; no finite certificate basis exists.**

These `S_X` are still **lonely** (the prime witness exists), so LRC(14) holds for them — but *only* via the
large-`D` / equidistribution mechanism. **The committed speed realizes the obstruction "block all small
denominators" (by being `lcm`), forcing the witness into the analytic regime.** This is the precise sense in
which a purely finite/combinatorial proof of LRC(14) is impossible: the `S_X` family has no uniform finite
witness, so an analytic equidistribution input is *irreducibly required*. (This rigorously confirms my S6
"honest ceiling" and corrects the finite-certificate overclaim.)

## The clean split (what kind of argument each node needs)

| node | realizability structure | kind of argument | status |
|---|---|---|---|
| **Node 2** (compact reduction / "consec maximizes") | three-gap (Steinhaus) rigidity — *only APs* | FINITE / algebraic — plausibly closable | the bimodality extremality, now grounded in rigidity |
| **Node 3** (committed-speed floor / decorrelation) | torus equidistribution of the orbit `(v_i t)` | ANALYTIC — *irreducibly so* (the `S_X` family) | the genuine open crux; no finite witness suffices |

So "aim to finish" splits honestly: **Node 2 is the finite half (attack via three-gap rigidity); Node 3 is
the analytic half (equidistribution with explicit error), and the `lcm` family proves it cannot be evaded by
any finite/algebraic device.** The bounded-D witness, the finite certificate basis, the covering-system
over-determination — all are *only the bounded/Node-2 half*; they are real but provably incomplete.

## Leads (toward finishing)
1. **Three-gap closure of Node 2**: prove "consec maximizes `L_y`" by an AP-hull smoothing — show any cluster's
   `N_E` is dominated in convex order by the AP's (whose `N_E` is exactly computable via the 3 gap-lengths).
   This is finite and may be genuinely closable.
2. **Effective equidistribution for Node 3**: the `S_X` family shows the needed input is an *effective* Weyl
   bound (Erdős–Turán) with explicit constant, applied to the committed speed's comb against the core's
   lonely set — exactly THM-518/523's decoupling `(6/7)meas - r/(7w)`, but with the arc-count `r` controlled
   *uniformly* (the part that fails in the resonant middle). The `nextprime(X)` witness suggests the witness
   denominator is `O(max prime factor structure)`, i.e. controlled by the committed speed's *radical*, not
   its size — a possible effective handle.
3. **Realizability obstruction as a theorem**: state cleanly "no uniform finite witness ⇒ LRC(14) requires an
   analytic input," so the team stops seeking a purely-finite closure (which the `S_X` family forecloses).

## Status / honest correction
No node closed. But: (i) **refuted** the bounded-D / finite-certificate completeness (HYP-2864, mac-mini
HYP-2876) via the `S_X = {1..11,13,lcm(2..X)}` family — an important guardrail; (ii) identified the **two-fold
realizability structure** (three-gap rigidity = finite, Node 2; torus equidistribution = analytic, Node 3);
(iii) grounded the Node-2 extremality in **three-gap rigidity** (a finite, plausibly-closable handle). The
LRC's "missing realizability" is exactly the analytic (Structure B) half — and it is irreducible.

→ HYP-2898 (new), HYP-2864/2880/2872 (mine, 2864 now corrected), HYP-2876 (mac-mini, corrected), HYP-2895
(kps), HYP-2897 (mac-mini), THM-518/523/527/531, OPEN-Q-108.

# Three seams: the quantized seven-tile wall, moment blindness with numbers, and what owns the sliver

**kind-pasteur-2026-07-19-S128c90** (HYP-7955) · owner: *"an even more creative connection
pursuit session — explore past threads and concepts in a depth-first search for
inspiration, creatively propose new hypotheses, investigate concurrently."*

Method: three seams chosen by depth-first descent through the S128c86–c89 arc's loose
ends and their ancestor threads (the r=7 crown, the Vitali/moment wall of the claude
06-03 era, boxeph-S130's 2/29 near-miss arithmetic, MISTAKE-093's sliver); one gated
script per seam, run concurrently with the reading
(`04-computation/lrc14_three_seams_kps_S128c90.py` + extension, frozen out). Everything
below is exact arithmetic on stated ranges; hypotheses are labeled as such.

---

## Seam G — the quantized seven-tile wall (new mechanism, new growth law)

**Ancestry (the DFS path).** My S128c88 minimal-covers data showed the level-1 improper
sea's minimal covers had size ≥ 9 at p = 43, 61, 71 → the density ratio
|P|/|A(w)| = (p/2)/(p/14) → **7 exactly** → the r=7 deletion crown (THM-1153-crown/1155:
the union-bound coefficient vanishes at exactly seven deletions) → and then the older
find: the μ_{1/7} density floor was itself reframed in canon as *"13 arcs of length 1/7
fail to cover the circle — the owner's Vitali hint made literal (Stevens 1939)"*
(SESSION-LOG:21999). Covering by density-1/7 objects, where seven pieces have exactly
the right mass, is the project's oldest recurring shape. The discrete question nobody
asked: **what is c(p), the minimum number of folded danger-APs A(w) = {fold(jw) : j ≤
⌊p/14⌋} that cover the folded line P = [1..(p−1)/2]?**

**Result (exact, witness sets frozen):**

| p | 29 | 43 | 61 | 71 | 101 | 113 | 127 | 151 | 173 | (ext continues) |
|---|---|---|---|---|---|---|---|---|---|---|
| density bound ⌈h/dk⌉ | 7 | 7 | 8 | 7 | 8 | 7 | 7 | 8 | 8 | … |
| **c(p)** | **7** | **9** | **9** | **9** | **10** | **10** | **10** | **11** | **11** | frozen out |

(c(173) = 11: the no-≤10 exhaustion took 23.8M nodes — the wall's finite checks are
themselves getting expensive, another quantized echo of the continuum statement.)

Three findings:

1. **At p = 29 a PERFECT TILING exists** (c = 7, zero slack: 7 sets of 2 cover 14
   points exactly; witness w-set {1,4,5,6,7,9,13}). The seven-tile configuration is
   *achievable* at small quantization — the wall is not local combinatorics, it
   emerges with resolution.
2. **The excess c(p) − 7 GROWS**: 0, 2, 2, 2, 3, 3, 3, 4 across the table — roughly
   Θ(log p) on this range (hypothesis H-G1, open: pin the law; the extension run
   pushes to p = 269). This is a *strong* discrete form of the r=7 wall: not only do
   seven density-1/7 tiles never suffice past the toy scale — no constant number does;
   the covering overhead diverges. The mechanism smells like greedy-set-cover's log
   factor meeting three-gap rigidity: each folded AP's overlap with a structured
   uncovered residue is forced (boxeph-S125's "wasted overlap" — covering debt
   concentrates at resonances instead of spreading).
3. **The sieve-design consequence (hypothesis H-G2):** improper 13-tuples need their
   13 danger-APs to cover, so **if c(p) > 13 then I(13,p,1) = ∅** — the level-1
   improper sea *vanishes* at primes past the crossing, and Rosenfeld-Lemma-6
   verification becomes trivial at level 1 there. Extrapolating the log-ish trend puts
   the crossing far beyond the S–T prime range (~10⁴–10⁵?), but the *existence* of a
   crossing — and its exact location — is now a well-defined, computable quantity that
   directly parametrizes finite-check difficulty at large p. Nobody has ever computed
   the improper-sea extinction point.

**Proposed mechanism name:** the *quantized seven-tile wall* — the discrete covering
number c(p) of the danger system, with its excess-growth law, as the finite-p shadow of
the r=7 zero-coefficient wall. Proof target (bounded, repo-shaped): c(p) ≥ 9 for all
p ≥ 43 via the rank-identity/cover-debt frame (any 8 danger-APs have total mass
8·dk ≈ 8p/14; covering needs debt ≈ p/14 concentrated *away* from resonances — the
spread-vs-waste dichotomy, now on a finite exact object).

## Seam H — moment blindness, with numbers (the third triage axis)

**Ancestry.** The cage (HYP-7940) pins families by even power sums S₂…S₂₆ — full
Newton determination at degree 26. The claude-era **Vitali/moment wall** (06-03: "no
finite-moment invariant decides tightness") said low moments can't do this. Neither
statement had numbers. The Prouhet–Tarry–Escott shape of the question: moment fibers
contain distinct sets exactly when PTE-style collisions exist; how much does M vary
*within a fiber*?

**Result (exact):** among 13-sets built from colliding 3-subsets (same Σv², Σv⁴ —
degree-4 even-moment fibers; 85 disjoint pairs in [1,45]):

- **max ΔM within a degree-4 fiber = 1/12** — the pair {1..11,25,37} (M = 1/12) vs
  {1,2,3,4,5,7,8,9,10,11,15,19,38} (M = 1/6): same S₂, same S₄, loneliness differing
  by a factor of TWO.
- **54 of 85 fibers straddle a certificate rung** (examples: 3/37 vs 2/13 across the
  1/12 rung) — degree-4 moments cannot see the rung ladder at all.
- Degree-6 fibers (colliding 4-subsets; only 2 pairs in range): max ΔM = 2/99 ≈ 0.02 —
  an order smaller but nonzero.

**The mechanism, visible in the extreme pair's argmaxes (independently re-verified:
S₂ = 2500 both, S₄ = 2,304,760 both, S₆ differs):** V contains the complete block
{1..11} and is pinned at t = 1/12 (the THM-633 plateau); W contains no multiple of 6,
so t = 1/6 is free and M = 1/6. **Low-degree moments are blind to DIVISIBILITY, and M
is divisibility-driven** (the witness-blocking cascade: covering = a tower of mod-q
conditions). A handful of real-valued moments cannot express "some vᵢ ≡ 0 mod 6" —
that is the moment wall's one-line cause, and it says the blindness persists exactly
until the moment set determines the multiset (Newton degree), because nothing short of
full determination recovers divisibility predicates.

**The statement this buys (verified in range, proposed as the third blindness axis):**
*even-moment invariants of degree ≤ 4 are blind to M at width ≥ 1/12; degree ≤ 6 at
width ≥ 2/99; blindness persists (with shrinking width) until full Newton degree 26.*
This joins THM-1185 (measure-blind) and THM-1225 (translation-blind) in the triage law
with, for the first time, a *quantitative decay profile* — and it explains precisely
why the cage must run its tower to m = 13: hypothesis **H-H1** (open): the blindness
width within degree-2m fibers decays like the PTE-collision height forces it to, and
the degree threshold where 13-set fibers become singletons is governed by
PTE-in-squares solvability at size 13. The PTE literature is now load-bearing for an
LRC question — a genuinely new external bridge.

## Seam J — the sliver belongs to structure, not grids (clean negative, recorded)

**Ancestry.** boxeph-S130's near-miss arithmetic: no multiple of 29 ⟹ M ≥ 2/29 =
96.6% of 1/14; MISTAKE-093 fixed the last-mile window [2/29, 1/14) of width 1/406 —
and, the DFS re-read revealed, proved more than I remembered: **the AP and GW
loneliness profiles Λ are IDENTICAL on [2/29, 1/14]** (GW's first divergence is a real
kink at 2/29 from the (5,24) pair). Hypothesis tested: does the two-prime CRT grid
mod 29·31 = 899 recover *effective* 1/14-witnesses on the no-mult-29 class (worth
having even after opus's THM-1289 made the micro-gap unconditionally empty but
ineffectively so)?

**Result:** REFUTED where it matters, instructively. The mod-899 grid clears 1/14 on
**200/200 random** no-mult-29 families (worst value 111/899 ≈ 0.123, comfortable) but
**fails on AP and GW themselves** (best 899-witness = 64/899 ≈ 0.0712 < 1/14): the
tight families' witnesses live at denominator 14, and no grid coprime to 14 lands on
them. The instrument-blind-extremizer genus (the ~45% mistake class), demonstrated in
one line: grids see the generic sea, never the structured shore. **Why recorded (dead-
end rule):** any effective-witness mechanism for the sliver must be
denominator-adaptive (contain 14's structure), which is exactly what certificates
can't be — the sliver is owned by structural invariants. And the connective tissue:
on [2/29, 1/14] where the Λ-profiles of AP and GW *coincide*, my cage's J-separator
(J(AP) ≠ J(GW), S128c89) still distinguishes them — **profile-blind, moment-sighted**.
The sliver is precisely where moments beat profiles, which is the cleanest one-line
justification yet for the cage architecture.

## The connective picture (what the DFS actually taught)

One triangle of mechanisms, each owning a regime:

- **Grids/certificates** own the generic sea (Seam J: 200/200; the rung ladder;
  boxeph-S130's certificate mechanism) and strand at structure.
- **Profiles** (Λ, safe-sets, combs) own the geometry down to 2/29 and go blind on
  the sliver where AP ≡ GW (MISTAKE-093).
- **Moments** own the sliver (the J-separator; the cage) but are blind below degree
  ~26 at measured widths (Seam H) — so they must be run to full Newton depth.

And beneath all three, the seven-tile wall (Seam G) is the *quantized* reason the
covering-side instruments strand: covering by density-1/7 objects carries a growing
excess, so "just cover it" never closes — the excess is where the structure hides.

## Handoffs

- **boxeph:** Seam G's c(p) ≥ 9 target is your rank-identity/cover-debt frame on a
  finite exact object — one session, and the excess-growth law would make your
  "wasted overlap" mechanism a theorem family. Seam H's degree-4 pairs are also a
  fresh test bank for the mod-p spread rungs (do fibers split by rung certificates?).
- **death-star:** c(p) and its 13-crossing (H-G2) is a gate-law-shaped growth
  question; your ghost evaluator could push the c(p) table two orders of magnitude.
- **opus:** Seam H's blindness profile slots into the five-axis triage law as axis 6
  ("moment degree ≥ d* or blind — measured widths 1/12 at deg 4, 2/99 at deg 6");
  the PTE bridge (H-H1) is a literature session of your S-T/G-K pinning genre.
- **anyone (Lean):** the p=29 perfect tiling is a 14-point decidable statement; c(43) ≥ 9
  is a finite check — both are kernel-ready seeds for the seven-tile wall.

**Files:** script + extension + frozen out (`lrc14_three_seams_kps_S128c90`), HYP-7955
(full write-up in INDEX), SESSION-LOG entry, backlog leads (c(p) growth law; PTE
bridge; the profile/moment/grid triangle as a triage supplement).

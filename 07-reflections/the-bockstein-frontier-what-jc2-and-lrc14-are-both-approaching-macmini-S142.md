# The Bockstein frontier: what JC₂ and LRC(14) are both approaching

**Instance:** mac-mini-2026-07-27-S142 (owner: push the JC₂ and LRC14 frontier; examine how our
thinking changed over the repo's course; boldly suggest new approaches; feel for what is most
fundamental). **Type:** synthesis + one exact computation (the χ ledger, frozen out) + three
typed approach proposals. Stands on opus-20260726 (spectral vs geometric rank), boxeph-S205
(nullcone rigidity), kps-S131b (no-cancellation grammar), CURRENT-FRONTIER.md (the rank-11/12
atlas), and the era maps; does not repeat them.

## 1. How the thinking actually changed (both problems, compressed to the moves that mattered)

**LRC(14), four object-corrections deep:** statistics of families (means, moments — died at the
Vitali wall) → the max-type object M and the covering reduction (THM-523) → exact certificate
coordinates (the (D,s) rung lattice, duty calculus, the plateau/strata closures) → and now, in
the 07-22..26 window, a **relation-code geometry**: a hypothetical counterexample is pinned to
sparse relation codes of rank ≥ 11 with a two-anchor star atlas, CRT clock compatibility, and
the crux living on the **joint 7⊗13 spectrum mod 91** (kps's mixed unit character). Each
correction moved the invariant one level up: statistic → certificate → lattice → **cocycle**.
The current frontier documents say it plainly: the live LRC objects are word-currents, holonomy
groupoids, endpoint cocycles — H¹ language — with a **first-Bockstein sidecar** (THM-2337/2356)
on the primitive-91 edge.

**JC₂, the same ladder in one week:** algebraic geometry of Keller maps → the engine-dimension
lemma and the golden corner (the Newton-polygon corner as a continued-fraction/Fibonacci
extremal — Diophantine, not algebraic) → the planar strata closures (THM-2102/2110/2113/2118/
2127/2129) → and the JC side now runs on **flux (the Keller one-form) on the trigonal spectral
cover with pole-descent** — again H¹ language, with geometric (S₃/dyadic) coefficients.

**The convergence is not decorative.** opus-20260726 verified the shape on both sides: the two
programs are each asking *"does a signed 1-cocycle on a branched cover fail to be a
coboundary?"* — differing only in coefficient type (spectral: the mod-91 character group;
geometric: S₃ monodromy and the dyadic tail). The rank-1 seed of the common nullcone frame is
proved functional-agnostically (THM-1840); GMC(2) — the third sibling — was actually PROVED and
kernel-formalized this week (THM-2022, GMC2Main unconditional): the first atom of the family to
fall to the repo's own machinery, and it fell on the *spectral* side.

## 2. The most fundamental thing we are approaching (bold, typed as a thesis)

> **The Bockstein thesis.** Every frontier the repo now touches is an instance of one question:
> *is the first Bockstein β of a canonical mod-N cocycle class nonzero?* Cancellation-
> completeness = coboundary-ness; a "one signed bit" is literally a mod-2 lift obstruction.
> - **Rédei's theorem is the solved prototype**: H ≡ 1 mod 2 is a nonvanishing mod-2 class,
>   PROVED by an involution with odd fixed locus — the repo's founding move is a Bockstein
>   computation avant la lettre. LEM-020's 2-adic tower is the Bockstein ladder; the S123
>   mirror palindrome is the involution that computes β at depth 1 on Wall-B data.
> - **JC₃ fell** because its class had geometric coefficients and dimension 3 gave the class
>   room to be nontrivial — the collision realized it (and the anatomy was still
>   parity-lawful: odd fibers over the symmetric locus; escape carried 2-adically, det = −2,
>   ℤ[1/2] inverse data, Mersenne fibers at the bad prime).
> - **LRC(14)'s crux** is the nonvanishing of the 7⊗13 character line — and canon already
>   carries its Bockstein sidecar (THM-2337/2356). **JC₂/DC₁/GMC(3)** are the atoms where the
>   corresponding class has no room (dimension/floor) and the conjecture is that it vanishes.
> The solved/open boundary in this repo tracks exactly: *solved where β is computable by an
> involution or where a collision realized the class; open where neither mechanism yet has a
> handle.* What we are building, across every problem, is a **computability theory for β** —
> when can complete cancellation be certified, when can it be refuted by construction.

## 3. Three bold approaches (typed, with referee paths)

**(1) Bockstein-force the 91-line.** Aim the involution technology (Rédei/LEM-020/the S123
palindrome — which localizes parity obligations onto 2-adic fixed loci) at the rank-11
residual's CRT clock discharge: compute the THM-2059 histogram dot products **mod 2** under the
mirror pairing, and separately compute the THM-2337/2356 Bockstein sidecar explicitly on the
primitive-91 edge. Success = the mixed-character nonvanishing converts from analysis to parity
(a forced-odd count). This is the repo's founding mechanism aimed at its current deepest edge;
the instruments exist and are mine.

**(2) The rank-12 box first.** The geometric-rank law (opus-20260726: collisions live at the
maximal rank the ambient allows — JC at 3 = dim, GMC at 4) predicts that if LRC(14) fails at
all, it fails in the **rank-12 maximal-minor finite box**, not the rank-11 stars. The box is
finite and named-undecided in CURRENT-FRONTIER. Discharge it exhaustively and the geometric
channel is CLOSED, leaving only the spectral channel — converting the conjecture's status into
a structural asymmetry ("no room for a collision; only Diophantine rigidity remains").

**(3) The χ-selection law — with this session's honest negative, AND a referee caveat added
in-session (S143): the fiber-count STRATIFICATION (3/2/1 on exactly those strata) is NOT yet
referee-solid — a fast fiber-cubic counter proved UNSOUND at degenerate strata (it undercounts
(0,1,0): the 0/0 branch of the y-parametrization resolves to a genuine point; S131's Gröbner
count of 3 stands), so the LEDGER BALANCE below is a hand-derivation pending a corrected exact
count map; the disc-cube, χ(K) = 1, and χ(Σ) = 1 computations are solid.** Computed tonight
(frozen out): the Jelonek quartic K is quadratic in c with **disc_c(K) = (b² − 12a)³ — a
perfect cube**; K is a double cover of the plane branched on a parabola, degenerating over
a = 0, and **χ(K) = 1 exactly**. The deep stratum Σ (1-point fibers) is the line {a=b=0} plus
the curve {a=0, bc=1}, **χ(Σ) = 1**. The stratified Euler ledger of the Keller map balances
perfectly: 1 = 3·(1−χ(K)) + 2·(χ(K)−χ(Σ)) + 1·χ(Σ) = 3·0 + 2·0 + 1, and every load-bearing
object — the engine cubic (THM-1340), the Jelonek surface, the deep stratum — has **χ = 1**.
BUT: the same numerology balances formally in dimension 2 (d = 3, χ(A) = 1, χ(Σ) = 1), so
**χ-accounting alone is NOT the JC₂ obstruction**. Approach (3) therefore refines to: extract
the actual mechanism of boxeph-S144(B)'s descent block (why THIS collision cannot fit in ℂ²)
and generalize it to all geometric ranks — the true JC₂ target is the general descent-blocker,
not the Euler ledger. (Negative results that kill one's own fresh idea in the same session are
the culture working; filed as such.)

## 4. Scope honesty

§1 is sourced narrative; §2 is a THESIS (organizing, falsifiable in its predictions — e.g. it
predicts the THM-2337/2356 sidecar computation is decisive-shaped, and that the rank-12 box is
empty); §3(3)'s computation is exact and frozen; nothing here modifies canon. The week-sweep
agent was killed by a disk-full incident mid-run; this synthesis stands on the primary
documents read directly (CURRENT-FRONTIER, the opus-20260726 reflection, the INDEX digest) —
a fuller week narrative remains available to any session that re-runs the sweep.

**Files:** `04-computation/jc_chi_ledger_macmini_S142.py` (+ frozen .out), HYP-9040.
**Cross-links:** opus-20260726 (spectral/geometric), boxeph-S205, kps-S131b, THM-2022 (GMC(2)
proved), THM-1840, THM-2051..2059 (the atlas), THM-2337/2356 (the Bockstein sidecar),
THM-1310/1315/1340, LEM-020, HYP-7965 (palindrome), THM-523 (the first object correction).

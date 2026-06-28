# The Joukowski/De Moivre bridge: the LRC covering bound's Lee-Yang circle IS the tournament's real-rooted class; the ideal is the 7th cyclotomic, the dip is deviation from 7-fold symmetry

*mac-mini-2026-06-27-S73. The owner asked to explore deeply for the niche connection that helps close
LRC(14) and to test bold ideas as I go. The frontier (S66-S72) reduced the covering bound to **bounding the
dip** (the off-circle / odd / phi^4 content of the miss-count PGF), with cap = binomial exact (THM-577). The
niche connection: the two halves of the WHOLE project are the two classical "zeros on a curve" theorems, and
the **Joukowski/De Moivre map `w = z + R^2/z`** is the bridge between them. Tested and (mostly) verified.
Continues [[the-lee-yang-circle-web-radius-coverage-off-circle-dip-and-the-compression-hierarchy]],
[[real-rootedness-of-independence-polynomial]]; feeds codex HYP-3153/HYP-3132/HYP-3147, HYP-3127 (Asano).*

## The two halves are the two zeros-on-a-curve theorems
- **Tournament side:** `I(Omega(T), x)` (independence polynomial of the odd-cycle conflict graph; `H = I(Omega,2)`)
  is **REAL-ROOTED** — zeros on the negative real axis — because `Omega` is claw-free, by Chudnovsky-Seymour
  (the Heilmann-Lieb class). PROVED for `n <= 8` (claw-free), conjectured all `n` ([[real-rootedness-of-independence-polynomial]]).
- **LRC side:** the miss-count PGF `G_N(z) = sum_t q_t z^t` (degree 6 over the **7 sectors mod 7** — the apex
  prime) is **CIRCLE-ROOTED** (Lee-Yang): zeros near `|z| = R`, coverage `p0 = q0 = q6 R^6` (S72).

Real axis and circle are the two stable-polynomial loci. The map that takes one to the other is **Joukowski**:
`w = z + R^2/z` sends `|z| = R` (`z = R e^{i th}`) to `w = 2R cos th` (**real**, in `[-2R, 2R]`). This is
*exactly* the "k=8 De Moivre biquadratic resolvent" of S70 — and now we see WHY it deserves the name.

## VERIFIED EXACT: the apex-7 ideal IS the 7th cyclotomic, Joukowski image = the de Moivre angles
The perfectly-symmetric (lambda = 0, on-circle) coverage PGF is the uniform `1 + z + z^2 + ... + z^6`, whose
zeros are the 6 nontrivial **7th roots of unity** (on the unit circle), and whose Joukowski image is
```
   {2 cos(2 pi j / 7) : j = 1,2,3} = {-1.8019, -0.4450, 1.2470}   (EXACT)
```
These are the **de Moivre angles**. So the S70 "De Moivre resolvent" is literally the resolvent whose ideal
roots are roots of unity. **cap = binomial = the 7-fold cyclotomic ideal; the dip = the deviation from
perfect 7-fold symmetry.** (`lrc_joukowski_resolvent_macmini_S73.py`.)

## VERIFIED ROBUST: consec is the cyclotomic-closest extremizer
For consec runners `{0,...,k-1}` the Joukowski resolvent root `j=3` stays **pinned at `~ -1.80` = 2cos(6pi/7)`**
while the other two drift with `k` — the drift is the dip. And the extremality is rock-solid:
- Over **400 random sets at each of n = 7, 8, 9**, ZERO beat consec's coverage `p0`. consec is the **max-coverage,
  min-off-circle** set at every n.
- **Scale-invariance** falls straight out: `2*{0..7} = {0,2,4,...,14}` has *identical* `R, p0, Im` to `{0..7}`
  — the `x2` dilation face (`14 = 2*7`, the H2 multiplicative correction kps flagged in HYP-3084).

The off-circle deviation `Im(w)` and the coverage `p0` are **negatively associated** (off-circle => less
coverage = the dip >= 0), but only loosely as a scalar correlation on random sets (`-0.03..-0.42`); the clean,
robust statement is the **extremality** (consec = max p0, min Im, 0/400 exceptions), not a tight correlation.

## The reframing this buys (toward the proof)
The open core of the covering bound is **dip >= 0 / consec maximizes coverage**. Via the bridge:
> **dip >= 0  <=>  the Joukowski resolvent stays as real-rooted as possible (zeros nearest the circle)  <=>
> a Lee-Yang/stability statement for `G_N`.**

And real-rootedness/stability is **exactly the tournament side's PROVED property** (Chudnovsky-Seymour for
claw-free `Omega`). The Joukowski map carries the LRC circle onto the real axis where that theorem lives. So
the covering-bound dip is not a new analytic beast — it is the **real-rootedness defect** of the De Moivre
resolvent, the same object class the project already controls on the tournament side.

**The apex lands in PROVEN territory.** Chudnovsky-Seymour real-rootedness of `I(Omega,x)` is unconditional
exactly while `Omega` is claw-free, i.e. while 3 vertex-disjoint odd cycles can't fit — `n <= 8` (three
3-cycles need 9 vertices). The apex is the prime **7 < 9**: the winding tournament at the LRC optimum is the
circulant on 7 vertices (S57), whose `Omega` is claw-free, so its `I(Omega,x)` is **provably** real-rooted —
no conjecture. The Joukowski bridge therefore lands the apex-7 resolvent on the real axis precisely where the
real-rootedness theorem is a theorem, not a conjecture.

This dovetails with codex's concurrent Asano program: codex finds Asano contraction **blocked on the R/tail
block** because it has *interior Lee-Yang zeros* (HYP-3128/3132) — that is precisely the **off-circle**
(dip > 0) here, and the apex/tip block being Lee-Yang-safe is the **on-circle** (cyclotomic) part. My bridge
gives the geometry: the obstruction is the gap between `G_N` and the 7th-cyclotomic ideal, measured by `Im(w)`.

## Ferromagnetic merge (S73b): tournament = AFM, LRC = FM, dip = frustration — and the decoupling the niche sets expose
This bridge is the physics bridge between the project's two magnet pictures, which already existed separately:
- **[[the-antiferromagnetic-tournament]]** (opus-S15): the tournament IS a frustrated spin system; `I(Omega,2)`
  is the hard-core (repulsive) partition function; regular/Paley = the Néel/AFM ground state = the H-maximizer;
  transitive = the fully-polarized ferromagnet. So the tournament side is **ANTI-ferromagnetic / repulsive**,
  which is exactly why its zeros are **real** (Heilmann-Lieb, the repulsive-gas locus).
- **[[the-two-maps-lee-yang-extremality-tournament-real-axis-lrc-circle-kps]]** (kps-S31ah): the two partition
  functions are the two Lee-Yang loci — tournament `I(Omega,x)` real-axis, LRC `Q(z)` circle (the original
  ferromagnetic-Ising Lee-Yang theorem). My Joukowski `w=z+R^2/z` is **the explicit map between kps's two
  loci** — they are not two parallel facts but one object seen on the circle and on its real shadow.
So: **TOURNAMENT = anti-ferromagnet (real axis); LRC coverage = ferromagnet (Lee-Yang circle); the dip / off-
circle deviation = frustration; the AP is the unfrustrated FM configuration.** (φ⁴ HYP-3099: the dip = the
quartic coupling λ that the off-circle motion turns on.)

**BUT the niche-set test (`lrc_ferromagnetic_niche_sets_macmini_S73.py`) tempers this — coverage and circle-
tightness DECOUPLE.** Fixed at 13 runners (= LRC(14)), ranking by circle spread `max|r|/min|r|`:
```
 primes      spread 1.247  (TIGHTEST circle / least frustrated)   p0=0.479
 GW          1.341                                                 p0=0.676
 AP {1..13}  1.361  (NOT tightest)                                 p0=0.705  (MAX coverage)
 squares 2.26 ; random 2.77 ; Fibonacci 3.25 ; Sidon 3.34 ; dyadic 8.9 (MOST frustrated, |Im|=14.8)
```
The **AP maximizes coverage `p0`** (the LRC quantity) but is **not** the tightest Lee-Yang circle — the
**primes** are. So "AP = tightest circle" (kps MAP 2, my opening) is really "**AP = max coverage**"; circle-
tightness is a *separate* functional that the primes win. The reusable caution: **the Lee-Yang circle captures
the GEOMETRY, not the LRC coverage coefficient** — `p0 = q6·R^6` decouples because `q6` (the all-missed mass)
varies too (primes have larger `R` but smaller `q6`, netting lower `p0`). The off-circle↔coverage link is a
loose heuristic, not a law; do not over-invest in circle-tightness as the proof handle. Two more honest reads:
- **Commensuration, not cooling:** consec's circle is tightest at **k=8 = 7 runners on 7 sectors** (spread
  1.14), then *loosens* to 1.36 by k=13 — a commensurate-filling effect, not "more runners = colder".
- **dyadic/geometric is the spin-glass:** one far element (the multi-far of codex's program) blows the spread
  to 8.9 and `|Im|` to 14.8 — the maximally frustrated, maximally off-circle configuration. The far element
  "pushes a zero out" (codex HYP-3131) = a single dominant frustrated bond.

## Honest status
- **EXACT / verified:** Joukowski sends circle->real; uniform PGF -> the de Moivre angles `2cos(2 pi j/7)`
  (the cyclotomic ideal = the cap); consec is max-coverage & min-off-circle (0/400 random exceptions, n=7,8,9);
  `x2` scale-invariance; the `j=3` root pinned at `2cos(6pi/7)`.
- **Loose:** the `Im <-> p0` scalar correlation (a crude proxy; the extremality is the real content).
- **STRUCTURAL / BOLD, NOT a proof:** the identification "dip >= 0 = resolvent real-rootedness defect = the
  tournament-side stability class (Chudnovsky-Seymour), bridged by Joukowski." It REFRAMES the open dip as a
  Lee-Yang/real-rootedness problem and pins its ideal (the 7th cyclotomic); it does not close it. LRC(14)
  remains open. The next move is to ask whether `G_N`'s single-runner factors admit a Grace-Walsh-Szego /
  Asano argument that keeps the Joukowski image real-rooted — i.e. prove the resolvent is in the stability
  class, not merely near it.

The reusable lesson: **coverage lives on a circle, tournament H lives on the real axis, and Joukowski
`w = z + R^2/z` is the same map that makes the apex-7 resolvent "de Moivre" — its ideal is the 7th roots of
unity, and the dip is how far the runners push the coverage off perfect 7-fold symmetry.**

Related: HYP-3154 (this), HYP-3152 (Lee-Yang circle), HYP-3132 (k=8 biquadratic = the resolvent), HYP-3147
(Worpitzky = odd/off-circle), HYP-3153 (codex packet), HYP-3127/3128 (Asano), THM-577 (cap = binomial),
[[real-rootedness-of-independence-polynomial]] (Chudnovsky-Seymour), OPEN-Q-108.

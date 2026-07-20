# The novelty audit, the Mersenne-depth law, and the first working of the six targets

**Instance:** opus-2026-07-20-S422 (HYP-8190). Owner: work the six targets; audit
whether Target 1 (and the rest) is novel and meaningful GIVEN the JC disproof is one
day old and external; decide what deserves rigorous writeups + Lean; sweep for
attackable threads including the ones that appear only once or twice.
Scripts + frozen outs: `targets_multithread_opus_S422.py`,
`jacobian_inverse_2adic_mersenne_law_opus_S422.py` (extends death-star S59m to K=16).

## 1. THE NOVELTY AUDIT (the owner's direct question, answered honestly)

**Target 1 (Zhao transport).** Mathematically CROLLARY-GRADE: every expert knows the
universal VC/IC forms fell with JC; the equivalences are literature. What is genuinely
valuable: (a) FIRST + EXPLICIT + MACHINE-VERIFIED witness — explicit artifacts after an
existence-level shock get cited for decades (the community will produce one within
weeks-to-months; SPEED is the edge); (b) our unique adds that no fast-follower will
have: the witness inherits the torus equivariance, the 2-adic ladder (now: the
Mersenne-depth law, section 2), and the sporadicity frame — we can CHARACTERIZE the
witness, not just exhibit it. VERDICT: do it now, as a short sharp paper + Lean-verified
identities; claim service + structure, not depth. CAUTION LEARNED THIS SESSION (both
demonstrated symbolically, frozen in the .out): the naive substitution gadget does NOT
preserve Keller-ness (det picks up 89 off-section monomials), and the naive symmetric
potential P = <y,H(x)> is NOT Hessian-nilpotent even for toy H — the BCW cubic
reduction and the dBvdE symmetrization both have non-obvious twists that MUST be pinned
from the sources (van den Essen book Ch. 6; de Bondt-van den Essen 2004-05; Zhao 2007)
before any transport run. S421's "pin the dimension" understated this; corrected in the
ledger. The machinery (Hessians, nilpotency checks, Delta^m towers) is now built and
validated as code.

**The rest of the JC thread:**
- THM-1300 (explicit Dixmier phi on A3): corollary-grade math, first-artifact value
  high. Lean the 18 Weyl identities (cheap, finite, high credibility). Merge with T1
  into ONE paper: "Explicit consequences of the Jacobian counterexample".
- THM-1345 (equivariant JC2 + Poisson bridge): **GENUINE standalone theorem,
  independent of the disproof — the highest paper-grade item we hold.** Full rigorous
  writeup; Lean the hyperbolic leading-coefficient lemma (medium project, worth it).
- THM-1330 (Keller monoid ideal) + S414b sporadicity + THM-1310/1335 anatomy: novel
  framing; the classification program is the medium-term paper; literature will move
  here fast — stake the monoid-ideal statement in the writeup above.
- LRC THM-1288 (refutation of published S-T Conjecture 7.1): INDEPENDENT of JC news,
  already complete, low effort to write — short arXiv note NOW. Well-regarded per the
  owner's criterion: it corrects a published claim with an explicit construction.
- Shear/Faulhaber/OEIS layer: OEIS submissions appropriate (kp's four + klein's batch
  + S420 candidates); not paper-grade alone.

**Lean priority queue:** (L1) F verification (det = -2 + collisions) — hours, the
foundation; (L2) Dixmier identities; (L3) the transported witness identities once the
constructions are pinned; (L4) THM-1345 hyperbolic lemma. NOT the full equivalence
chains (months, low marginal value). klein's LRC sorry-ledger infrastructure carries
over.

## 2. NEW RESULT: the Mersenne-depth law of the 2-adic inverse ladder (Target 2)

Extending death-star's formal-inverse computation from K=12 to K=16 (script above):
min v_2 per degree d = 1..16:
-1,-1,-3,-3,-4,-4,-7,-7,-8,-8,-10,-10,-11,-12,-15,-15.

**LAW (data-supported at all four available Mersenne degrees): the ladder reaches FULL
dyadic depth v_2 = -d exactly at d = 1, 3, 7, 15 = 2^k - 1** — predicted from the
K=12 data before the run, confirmed by the -3 jump at d=15 (-12 -> -15). Off-Mersenne
deficits epsilon(d) = d + v_2: 1,1,1,2 at d = 5,9,11,13 — exact epsilon-law OPEN
(neither binary-weight nor zeros-count fits all four; more terms needed, K=20+ next).
Also new: the even-degree pairing (v(2j) = v(2j-1)) BREAKS at d = 13,14 and re-locks
at 15,16. Reading: the inverse's dyadic descent is clocked by the torus doubling
lambda -> lambda^2; full depth at Mersenne degrees rhymes with the S420 Mersenne
shear-row sums — the counterexample's "2" is being spent one bit per degree except
where a Mersenne resonance spends it all. Feeds Target 2's hunt for a p-adic invariant
separating A_1 from A_3.

## 3. Target 6 delivered: first unit-distance facts about our metagraphs

Embedder validated on the Moser spindle (residual 4e-13). Then:
- **G_4 (arc-flip metagraph, V=4, E=5) IS a unit-distance graph** (embedded, min pair
  distance 1.000 — all pairs at >= unit).
- **G_5 (V=12, E=30) is NOT: it contains K_4** (no four points pairwise unit in the
  plane). The tournament carriers fail unit-distance at the FIRST nontrivial size, by
  clique obstruction. Next moves: max unit-distance subgraphs of G_5; whether K_4-free
  spanning subgraphs (spine/ribs only?) embed; E_n duals.
The analogy cards ("unit-distance cyclotomic norm", lowest retention in the technique
index) are now attached to actual theorems about our objects.

## 4. The comprehensive thread census (the "don't miss the niche ones" sweep)

File-count census of famous-problem mentions (git grep, whole repo):
- **MAJOR HIDDEN THREADS (need dedicated mining sessions):** Goldbach 151 files (!),
  Sidon 114, Kakeya 104, Collatz 93, twin primes 57, perfect numbers 55, RH 32,
  Hadwiger-Nelson 30. These are mostly codex-era lens/carrier work — the unit-distance
  precedent says each may hide a card->math conversion. ONE mining session each to
  triage which have substance.
- **GENUINELY NICHE (1-3 files — the owner's request):**
  - **Casas-Alvero conjecture (2 files):** a degree-n polynomial sharing a root with
    each derivative is (x-a)^n. OPEN, p-adic partial results in literature — DIRECTLY
    adjacent to our new 2-adic ladder + derivative-tower machinery. Attackable and
    fashionable-adjacent post-JC (polynomial rigidity). PROMOTE to target list.
  - **Sendov conjecture (2 files):** zeros vs critical points in discs — same
    derivative-geometry family.
  - **cap set (3), abc (3: one S116 survey script), union-closed (18), sunflower (18),
    Lehmer (12), Erdos-Straus (11).**
- **NEW BRIDGE SPOTTED:** Lehmer's Mahler-measure problem vs the S420 shear dial —
  Pascal-type grids dial PISOT numbers (x^s = x^(s-1)+1), exponential grids dial
  radicals: "which Mahler measures arise as shear-growth constants of relational
  grids?" is a well-posed new question connecting our shear calculus to a famous open
  problem. Logged as a thread seed.

## 5. Session verdicts in one line each

T1: constructions must be pinned from sources — negatives demonstrated, machinery
built. T2: Mersenne-depth law found and confirmed at d=15. T3: fleet's lane (kp/
mac-mini), monoid-ideal statement to be staked in the joint paper. T4: untouched this
session (finite computation spec stands). T5: fleet-active, untouched. T6: G_4 yes,
G_5 no (K_4) — delivered. Audit: one joint explicit-consequences paper + THM-1345
standalone + THM-1288 arXiv note + OEIS batch; Lean queue L1-L4.

Cross-links: HYP-8185 ledger (updated) · THM-1300/1310/1330/1335/1345 · death-star
S59m (2-adic base run) · S420 Mersenne shear sums · S421 targets · MISTAKE-198 genus
(the epsilon-law is left OPEN rather than jet-fit — lesson applied).

# The Jacobian counterexample through the repo's eyes — doubling ladders, observers, and the staircase with no bottom

**death-star-2026-07-19-S59m** (HYP-8075, THM-1300; owner prompt: "investigate deeply
connections to our large body of work and the dixmier conjecture. tournaments and their
recursion is reminiscent of the n vs 2n relationship between jacobian and dixmier").
Everything cited as *verified* below was computed exactly this session or by
kind-pasteur S128c97; the analogies are labeled as analogies. Grading is honest:
§1–§3 are mathematics, §4–§6 are structural resonances with testable offshoots,
§7 is the meta-read.

## 1. The two conjectures were one ladder, and the ladder just snapped asymmetrically

The classical lattice of implications is a **two-rail ladder with rungs in both
directions**: DC_n ⟹ JC_n (same level, Weyl → polynomial; the constructive transfer
is THM-1300 §1), and JC_{2n} ⟹ DC_n (Tsuchimoto, Belov-Kanel–Kontsevich — *descend
a doubling*: prove on the 2n-dimensional symplectic side, harvest on the
n-dimensional quantized side). The owner's counterexample severs BOTH rails at every
level n ≥ 3 — but the failure propagates **only upward and across, never down**:
identity-padding pushes ¬JC from 3 to every n ≥ 3, the same-level rung carries it to
¬DC₃, padding to ¬DC_{n≥3}. What survives is exactly the bottom: JC₂, DC₂, DC₁.
The doubling rung JC₄ ⟹ DC₂ is now *vacuous* (its hypothesis is false), so the
bottom's only live rungs are DC₂ ⟹ JC₂ and DC₁'s isolation.

The repo has lived inside exactly this grammar all year: **refutation kills a rung
and everything ladder-reachable from it; the base cases survive and become the
frontier** (LRC settled ≤ 13, open at 14; per-rung gate laws; the K-ladder's residual
corner). The JC/DC event is the same phase transition executed on a 90-year-old
tower: overnight, "prove it for all n" became "the open problem lives at n ≤ 2" —
a *floor* problem, which is the shape of problem this repo knows how to love.

## 2. n vs 2n is phase-space doubling — and the repo's Mode B is its inverse move

Why 2n? Because A_n is the quantization of the symplectic ℂ^{2n}: n positions X_i,
n momenta D_j, the Weyl relation as the symplectic form's shadow. The BKK bridge
runs through characteristic p, where A_n ⊗ F_p acquires a huge center — a polynomial
ring in the p-th powers X_i^p, D_j^p — which is **a fresh copy of the 2n-dimensional
space**; the endomorphism's restriction to the center is a Keller map in dimension
2n, and JC_{2n} (if true) forces it to be an automorphism, dragging DC_n along.

The repo's recursion has both moves natively:

- **Mode A (n → n−1, hypotenuse removal)** is the same-level rung: one observer
  removed, fast time scale — the analogue of DC_n ⟹ JC_n (stay at level n,
  change category: Weyl → polynomial; tournament → subtournament).
- **Mode B (n → n−2, both legs)** is the (un)doubling rung: the slow time scale,
  the Cayley-Dickson descent — the analogue of JC_{2n} ⟹ DC_n. Mode B removes a
  *conjugate pair* (source column + sink row = the two legs); BKK removes a
  conjugate pair at scale (positions + momenta). In both, the double-step is the
  one that changes what KIND of object you hold (CD: algebra halves and loses a
  property; BKK: symplectic-classical becomes quantized).
- The GF(2) **Cut ⊕ Cycle decomposition** (base-path arcs = cut space, tiles =
  cycle space) is the repo's own canonical "positions ⊕ momenta" split of one
  edge space into two complementary halves whose pairing carries all the
  structure. The tiling model fixes the cut half and lets the cycle half vary —
  exactly "fix configuration, quantize momentum."

Labeled honestly: these are structural isomorphisms of *ladder grammar*, not of
content. But note they pointed correctly: the owner's hint ("recursion ~ n vs 2n")
predicted that the productive thing to do with the counterexample was to RUN THE
LADDER CONSTRUCTIVELY — which yielded the explicit A₃ self-embedding (THM-1300 §1),
the session's hardest artifact.

## 3. The doubling is INSIDE the counterexample too (verified)

The map itself is ℂ*-equivariant (weights (1,−1,−2) → (−2,−1,1), THM-1300 §2), and
its non-injectivity is powered by **λ ↦ λ² on one torus orbit** — the fiber is
1 (fixed branch) + 2 (doubled orbit). The squaring map on ℂ* is the multiplicative
n-vs-2n; the counterexample to a conjecture *about* a doubling ladder is itself an
equivariant package around a doubling map. Three repo echoes, each checkable:

- **Rédei-shaped fibers**: 3 = 1 + 2 = fixed + swapped pair is precisely the odd
  count "identity + involution-orbits" that powers Rédei's theorem (inshat's
  1 + 2·#3-cycles) and the repo's SC-spine vs paired-classes census
  (V_merged = (A000568 + SC)/2 — the same "count = fixed + halved-pairs" Burnside
  arithmetic). kind-pasteur saw the ℤ/2; the full torus upgrades it: the parity
  mechanism is the residual finite part of a continuous symmetry.
- **The dyadic staircase**: the formal inverse exists at every finite 2-adic layer
  and dies at each one (v₂ ladder −1,−1,−3,−3,…,−10,−10, unbounded, paired by the
  λ↔−λ parity; THM-1300 §3). Compare the repo's newest lesson (S59l): "the D-gate
  structure lives entirely at the finest 2-adic layer." Here the OBSTRUCTION to
  inverting lives at the 2-adic layers all at once: det = −2 means the map
  degenerates exactly at p = 2, and the inverse's staircase has no bottom. A
  residual with a *name and a gradient* — the repo's typed-residual discipline
  applied to a century-old conjecture's corpse.
- **Sparsity = grading**: the inverse's absurd sparsity (7 terms where 91 could
  live) is the weighted grading acting — the same phenomenon as the repo's
  band-limitedness results (Walsh degree ≤ 2⌊(n−1)/2⌋): a hidden symmetry showing
  up as a support constraint before you know the symmetry is there.

## 4. The observer principle, read into the Weyl tower (analogy, with a testable)

The repo's CD tower ties algebra dimension 2^j to tournament size n = 2^j + 1 —
"the +1 is the observer." The Weyl tower ties phase dimension 2k to the algebra
A_k. Composing the two dictionaries: **A_k should shadow tournaments on n = 2k+1
vertices** — A₁ ↔ n=3 (the 3-cycle, the atomic odd cycle, where Rédei parity is
born), A₂ ↔ n=5, A₃ ↔ n=7 — the repo's H7/seven-wall home ground. The refutation
lives at A₃; the surviving open bottom DC₁ lives at n=3, the 3-cycle. If the
dictionary has any teeth, the repo's 3-cycle/atomic-parity technology is aimed at
the SURVIVING conjecture, not the dead ones. Testable offshoot (backlog): the
fiber-parity statistics of Keller maps mod p (kind-pasteur's stream) vs the repo's
odd-cycle parity census — is "1 + 2·(pairs)" the universal fiber shape of
equivariant Keller maps, i.e. is every counterexample's fiber an odd
"observer + doubled pairs" count? (The verified one is.)

## 5. The tree series whose tail never dies (analogy, with a computation attached)

The Bass–Connell–Wright inversion formula writes the formal inverse of a Keller
map as a sum over rooted trees — the same Catalan-shaped combinatorics as the
repo's cluster integrals (THM-438: Paley cluster integrals ARE Catalan numbers)
and path-sum cancellations (the OCF's signed collapse). JC-as-was said: for unit
Jacobian, the tree series *telescopes to a polynomial* — a perfect cancellation
claim, the exact genre of the repo's "too-clean cancellation" reflexes. The
counterexample says the cancellation is a lie at n ≥ 3: the tail of the tree
series survives at every degree (measured: sparse, dyadically descending, forever).
The repo's residual discipline gives the right epitaph: *the conjecture was the
claim that a certain typed residual vanishes; the residual is real, sparse,
graded, and 2-adically unbounded — now measure it, don't mourn it.*

## 6. What the repo can actually DO here (filed, priced)

- (a) **Fiber-parity census of equivariant Keller maps** (with kind-pasteur's
  mod-p stream): is the odd "1 + 2k" fiber universal under ℂ*-equivariance?
  Sub-hour once their Groebner degree lands.
- (b) **The weighted-torus classification**: which weight vectors (a,b,c) admit
  equivariant non-injective Keller maps? The verified one uses (1,−1,−2) with
  component weights the REVERSAL (−2,−1,1) — an antidiagonal flip, the repo's
  hypotenuse move. Enumerate small-weight candidates by the same
  validate-then-emit generator discipline as the LRC pipeline.
- (c) **Lean**: the verification stack is kernel-shaped already — det JF = −2 is
  a polynomial identity (linear_combination-able in principle, large), the three
  point-evaluations are rational arithmetic (decide/norm_num), and R1/R2 are 18
  polynomial identities. A `JacobianCounterexample.lean` with the collision +
  det certificate would make the repo host a machine-checked refutation of a
  90-year conjecture. Hours-scale, honest stub filed, not owed.
- (d) **The exact-poly engine** (pmul/padd/pdiff/adjugate/flatness/Laurent/
  truncated-inverse, zero dependencies, built this session under a pip-less
  sandbox) joins the engineering shelf next to mod_rank — the same
  verify-exactly-or-refuse ethos, now proven on a live mathematical emergency.

## 7. The meta-read

A conjecture that stood since 1939 fell to a *checkable object* — three
polynomials any exact engine verifies in under a second. Both in-repo
verifications agreed before any literature could (searches find nothing yet:
the object outran the indexes). That is the repo's entire epistemology playing
out in the large: claims are cheap, certificates are the currency, and the
right response to big news is neither belief nor doubt but a dependency-free
script and a typed statement of what exactly is now known. The Dixmier transfer
ran the same way: not "the literature says DC falls too" but an explicit φ, 18
exact identities, and a one-line module argument — the difference between
knowing THAT the ladder transmits failure and HOLDING the failed object at the
other end.

## Cross-links

THM-1300 (the artifact) · kind-pasteur S128c97 / HYP-8070 (verification,
σ-equivariance, Groebner + mod-p streams) · HYP-8075 (this session) ·
everything-is-the-triangle (Mode A/B) · merged-metagraph-invariants (Burnside
1 + 2k arithmetic) · S59l dyadic-collapse datum (the 2-adic layer lesson) ·
THM-438 (Catalan cluster integrals) · the LRC certificate pipeline (the
verify-then-emit discipline reused verbatim).

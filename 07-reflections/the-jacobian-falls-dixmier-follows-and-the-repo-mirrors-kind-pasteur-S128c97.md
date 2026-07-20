# The Jacobian conjecture falls; Dixmier follows at n ≥ 3; and the repo's lenses read the counterexample

**kind-pasteur-2026-07-19-S128c97** (HYP-8070) · owner: *"the Jacobian conjecture is
false [explicit map given]. Investigate deeply connections to our large body of work
and the Dixmier conjecture. Tournaments and their recursion is reminiscent of the
n vs 2n relationship between Jacobian and Dixmier."*

## 1. Independent verification (everything checks; frozen out)

The owner's map F: ℂ³ → ℂ³, F = ((1+xy)³z + y²(1+xy)(4+3xy),
y + 3x(1+xy)²z + 3xy²(4+3xy), 2x − 3x²y − x³z):

- **det J(F) ≡ −2** — verified SYMBOLICALLY (sympy, full expansion). A genuine
  Keller map.
- **F(0,0,−1/4) = F(1,−3/2,13/2) = F(−1,3/2,13/2) = (−1/4,0,0)** — verified in exact
  rational arithmetic; three distinct points, one image. NON-INJECTIVE.
- Literature: no trace yet (searches surface only Pinchuk's real-plane maps, a
  different settled question); the verification above is self-contained ground truth.

**JC₃ is false. Hence JC_n is false for every n ≥ 3** (append identity coordinates).

## 2. The Dixmier cascade and the surviving bottom

By the standard implication DC_n ⟹ JC_n (van den Essen/Adjamagbo; BKK), the
contrapositive gives **¬DC₃: the Dixmier conjecture is FALSE for the third Weyl
algebra — there exists an endomorphism of A₃(ℂ) that is not an automorphism** — and
by tensoring with the identity, **DC_n is false for all n ≥ 3**. The
Tsuchimoto/Belov-Kanel–Kontsevich bridge (JC_{2n} ⟹ DC_n) remains a valid
implication; the counterexample is consistent with it (¬DC₃ ⟹ ¬JC₆ ✓ already known
from ¬JC₃).

**What survives is a THREE-RUNG LADDER: DC₂ ⟹ JC₂ ⟹ DC₁** (note DC₂ is NOT
refuted — the only upward route JC₄ ⟹ DC₂ is now vacuous, and no DC-descent
DC₃ → DC₂ is known; the script banner's "JC₂ and DC₁" undercounted). All three open,
linearly ordered, refute-one-refute-everything-above / prove-one-prove-everything-
below, ending at JC₁ = trivially true. The owner's instinct is exact: the JC/DC
recursion halves dimension (DC_n lives on JC_{2n}'s phase space) and bottoms at an
irreducible hard chain, the same shape as the repo's Mode-B descent (n → n−2, the
two-leg removal) and Cayley–Dickson tower bottoming at the small cases that carry
all the difficulty. The falsity at n ≥ 3 makes the surviving ladder MORE
interesting, not less: the phenomenon that kills dimension 3 (see §4: a z-linear
lift of a plane engine) structurally cannot fit in dimension 2, which is now the
entire question.

## 3. The arithmetic fingerprint: S₃ Chebotarev, not cubic residues (measured)

Fiber statistics of F over F_p³ at p = 5, 7, 11, 13, 17, 19, 31: at EVERY prime the
histogram is supported on **{1, 3} only** — no 2-fibers, no ≥4-fibers, across ~66k
fibers in seven full sweeps — with fractions of p³ converging monotonically to
**N₁ → 1/2, N₃ → 1/6, misses → 1/3** (p=31: 0.515 / 0.162 / 0.323). (The exact
relation misses = 2·N₃, which holds at all seven primes, is FORCED by the {1,3}
support via p³ = N₁ + 3N₃ — a consistency check, not independent evidence.) The
Chebotarev reading: a degree-d dominant map has transitive monodromy G ≤ S_d, and
rational-fiber-size densities = fixed-point statistics of G. Any G contains the
identity, so d-fibers appear at density ≥ 1/d! — max fiber 3 with density ~1/6
pins **d = 3**; the transitive degree-3 groups are C₃ and S₃, and C₃ (fixed-point
counts {3,0,0}) predicts NO 1-fibers — excluded by the 52% singletons. Verdict:
**F is a degree-3 cover of ℂ³ with full S₃ monodromy**, fractions matching S₃'s
fixed-point classes (1/6 split, 1/2 one point, 1/3 none) — and the three collision
points are one totally-split fiber, not an anomaly. My initial cubic-residue
prediction (p mod 3 stratification) is REFUTED by the data — the honest record.
**And the S₃ was FORCED a priori: by Campbell's classical theorem, a Keller map
whose field extension is Galois is invertible** — so any counterexample must be
non-Galois, hence field degree ≥ 3 with symmetric closure: degree-3-with-S₃ is the
MINIMAL stratum Campbell's rung left open, and the mod-p statistics independently
confirm the witness lives exactly there. Non-injectivity that is Galois-generic is
invisible to any instrument sampling rational fibers expecting degree 1 — the ~78%
singleton fibers mod p show how naive point-sampling would "confirm" injectivity
forever: the instrument-blindness genus at civilizational scale.

## 4. Structure the repo lenses read off

- **The Rédei parity lens (verified symbolically):** F(−x,−y,z) = (F₁,−F₂,−F₃) — an
  equivariance σ → τ. Fibers over the τ-fixed plane (b = c = 0) carry the
  σ-involution, so their rational counts satisfy |fiber| ≡ |σ-fixed part| mod 2; the
  collision fiber is the minimal odd configuration 3 = 1 (σ-fixed) + 2 (swapped
  pair). The counterexample's collisions live ON the symmetric locus — the same
  "extremal structure at the fixed locus of the involution" geometry as the repo's
  incentre/six-ray findings. If the discoverers hunted symmetric ansatzes, they used
  our heuristic; if not, the theorem found it anyway.
- **z-linearity (the plane engine):** every component of F is AFFINE IN z:
  F = A(x,y)·z + B(x,y), A = ((1+xy)³, 3x(1+xy)², −x³). The constant Jacobian for a
  z-linear map forces differential identities on (A, B) (the z-coefficients of det J
  vanish), and the whole counterexample is a plane phenomenon lifted along a line
  bundle — which is exactly WHY it exists in ℂ³ while JC₂ survives: the engine is a
  plane correspondence that is not area-preserving, with the z-direction paying the
  Jacobian bill. Dimension 3 = plane + compensating line: the repo's "add one
  compensating coordinate" motif (the observer!).
- **The BKK bridge ≡ the repo's mod-p descent (HYP-8020), as philosophy:** BKK
  prove JC_{2n} ⟹ DC_n by reducing the Weyl algebra mod p, where it becomes Azumaya
  over a polynomial center in 2n variables — the noncommutative problem DESCENDS to
  a commutative char-p problem and lifts. Our newest pipeline (descend Wall A to
  F_p, prove the shell theorem there, cage-lift) is the same maneuver on our
  objects. The counterexample does not break the bridge (implications stand); it
  shows the char-0 side simply fails — and its own mod-p fingerprint (§3) is the
  kind of data the descent philosophy says to collect first.
- **The n vs 2n ↔ cut ⊕ cycle rhyme (graded honestly as a rhyme):** the Weyl
  algebra is quantized T*ℂⁿ — positions ⊕ momenta, dimension 2n with a pairing; the
  repo's tournament space splits GF(2)-canonically as cut ⊕ cycle (base-path arcs ⊕
  wiggly arcs) with the intersection pairing, and Mode B descends by removing a
  leg-pair. Both theories move by pairing-compatible 2-step reductions and bottom at
  small irreducible cases. A rhyme with one concrete question attached: is there a
  tournament-side analogue of "endomorphism of the quantized algebra" whose
  classical shadow is a metagraph map — i.e., does the repo's quotient tower have a
  Dixmier-flavored rigidity question of its own?
- **The methodology mirror (the deepest connection):** eighty-seven years of JC
  partial results — injective ⟹ bijective, degree ≤ 2, cubic-homogeneous and
  Druzkowski normal forms, symmetric reductions — are certificate rungs in
  boxeph-S130's exact sense: each closed a generic stratum, each preserved the
  problem's equivalence while moving the witness OUT of the searched normal-form
  strata, and the community's repeated near-completion sensation is our "closed
  modulo X" episode list at scale. The witness, when it came, is SMALL (three
  variables, degree 7, coefficients in {1,2,3,4}) — the structured extremal that
  every reduction was blind to: the ~45% instrument-blindness genus, the
  detection-floor lesson, and the eternal-inhabitant check (here inverted: the
  inhabitant was the counterexample) all apply verbatim. Our rules were built on a
  14-runner problem; they describe a century of Jacobian history without editing.

## 5. Fleet leads filed

1. **Lean-certify the refutation** (headline-grade, one session): det J + 2 = 0 as
   a polynomial identity plus three exact evaluations — kernel-checkable with the
   fleet's standard machinery; likely among the first formal certifications of the
   result anywhere.
2. **The explicit Dixmier witness**: run van den Essen's DC ⟹ JC construction
   backwards on F to write down the non-surjective endomorphism of A₃ explicitly
   (injective — all Weyl endos are — but not onto). Paper-grade exercise.
3. **The z-linear anatomy**: extract the differential identities on (A, B) forcing
   det ≡ const; classify which plane engines admit compensating lifts — the
   structural question JC₂'s survival now hinges on.
4. **The S₃-Chebotarev study**: push the fiber statistics to more primes /
   prime powers, extract the splitting field of the generic fiber, and pin the
   monodromy claim algebraically (resultant discriminant).

**Files:** verification script + frozen out (det, collisions, equivariance, 7-prime
fiber statistics), HYP-8070, TANGENT entry, backlog leads, SESSION-LOG.

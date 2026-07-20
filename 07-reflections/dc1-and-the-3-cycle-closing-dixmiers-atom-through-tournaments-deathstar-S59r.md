# DC₁ and the 3-cycle — reading Dixmier's surviving atom through the repo's founding object

**death-star-2026-07-20-S59r** (HYP-8160; owner: "see how the DC₁ frontier can be
creatively closed by checking analogies with tournaments"). **Honesty banner,
first and loudest: DC₁ is OPEN and this reflection does NOT close it.** DC₁
(Dixmier's conjecture for the first Weyl algebra A₁) is, via BKK, tied to JC₂ —
the graveyard of false proofs. What follows is a graded analogy hunt with
*verified teeth*, not a resolution. Grades: **bridge** (mechanism-level),
**rhyme** (shape-level), **stretch** (evocative). Credits up front: opus-S418
(the cyclic 3-tournament inside the counterexample's Galois fiber; Rédei-sign =
discriminant character), mac-mini-S137 (the golden/worst-approximable degree
corner of JC₂), boxeph-S144 (the dim-2 no-go).

## 0. What DC₁ is, and why the repo has any right to look at it

A₁ = ⟨p, q | [p,q] = 1⟩, the Weyl algebra / quantized symplectic plane. DC₁:
every endomorphism (every pair P, Q with [P,Q]=1) is an automorphism. A₁ is
simple, so endomorphisms are injective; DC₁ = "injective ⟹ surjective." It is
the classical shadow's quantization: JC₂ ({P,Q}=1 ⟹ coordinate) is DC₁ at
ℏ→0, and **JC₂ ⟹ DC₁** (BKK, n=1). DC₁ is the surviving Dixmier atom — the
whole doubling tower's floor (S59m). And the repo's **observer principle**
(A_k ↔ tournaments on n = 2k+1 vertices) puts DC₁ ↔ **A₁ ↔ n = 3**: the
3-vertex tournament, the first odd cycle, the repo's founding object. So the
frontier of a 90-year quantization problem lands, under the repo's own
dictionary, on the 3-cycle. That is the license.

## 1. The bridge with teeth: A₁'s weight triple IS the oriented 3-cycle

Represent A₁ on ℂ[x]: p = d/dx, q = x·, N = qp. Verified exactly (operator
matrices on {1,…,x⁶}):

  [p, q] = 1,   [N, q] = +q,   [N, p] = −p.

Read as a graph on three vertices {N, q, p}: N is the **Euler/grading**
operator (the "observer"), q the weight **+1** raiser, p the weight **−1**
lowerer, and the central relation [p,q]=1 the edge between them. The signs
(+ on N→q, − on N→p) orient a triangle — **the n = 3 tournament**. The classical
Poisson shadow is identical ({p,q}=1, {H,q}=q, {H,p}=−p, H=pq), so the triple is
exact at both ℏ=0 and ℏ=1. This is standard sl₂/Heisenberg representation
theory wearing a tournament coat — I file it as a **dictionary entry, not a
theorem** (MISTAKE-197 discipline: textbook facts stay textbook). The value is
the identification, next.

## 2. The identification that unifies three threads: observer = ℏ = the +1

The observer principle's "+1" (A_k has 2k generators; the tournament has 2k+1
vertices) is, concretely, the vertex N — the grading operator the symplectic
**pair** (p,q) gains to become the **triangle** {N,p,q}. And the central
relation is [p,q] = **1** — the ℏ. So:

  **observer vertex  =  ℏ  =  the conserved +1 of S59p/S59q.**

The +1 I tracked through the JC counterexample (u = 1+xy, the Yagzhev X, the
formal-group denominator, "conserved across every normal-form change; absent in
2D, which is why 2D is safe") is the SAME +1 that quantizes the Poisson plane
to the Weyl algebra and promotes the symplectic pair to the 3-cycle. This is
the content of the old wild guess W8 ("observer = ℏ"), now with an exact bracket
computation under it. **Bridge** (the arithmetic is real; the claim that this
*explains* anything about DC₁ is still to be earned).

## 3. The creative strategy: DC₁ as a two-lens problem

The tournament corpus's sharpest transferable lesson (mac-mini-S129):
**parity-protected closures are robust; bare-uniqueness beliefs (the JC/DC
genus) are fragile.** Rédei's theorem is the archetype — #Hamiltonian paths is
ALWAYS odd (1 + 2·#3-cycles), a quantity that structurally cannot vanish. So
the tournament-native attack on DC₁ is: **find the Rédei-parity of an
endomorphism, and pair it with a bound.**

The candidate parity is the **fiber-degree**. Via BKK, an A₁-endomorphism
reduces mod p to a Keller map of the char-p center (a symplectic plane); its
degree = generic fiber size. Verified teeth (T3): the two n=3 tournaments have
1 and 3 = 1+2 Hamiltonian paths, and the JC counterexample's fiber over the
axis is exactly **1 (fixed) + 2 (orbit) = 3**, the Rédei-shape (THM-1305, opus
S413/S418). So known counterexamples carry an **odd** fiber-degree — a genuine
Rédei shadow. But odd is *necessary, not sufficient*: 3 is odd and > 1, so
parity alone permits the higher-dimensional counterexamples. The second lens is
where DC₁'s dimension helps:

  **Lens A (tournament / Rédei):** fiber-degree is odd (parity-protected).
  **Lens B (Poisson filtration / THM-1345):** in dim 2 the leading-form
   obstruction {P_A,Q_B}=0 descends through the weight filtration and bounds
   the degree — the mechanism that fails at higher dim (where the +1 has room
   to build a coupled unit cascade, S59q).
  **Together:** odd ∧ bounded ⟹ degree 1 ⟹ automorphism ⟹ DC₁.

This is a **strategy, graded speculative**, not a proof — Lens B's "bound" is
exactly the open AMS-hard content of JC₂ (THM-1345 §5). But the framing is
honest and it says something structural: **DC₁ may hold where DC₃, DC₄, … fail
because dimension 2 supplies only one unit (no crossing), while the Rédei
parity is dimension-blind.** The failure is dimensional; the parity is not; DC₁
lives exactly where the dimensional lens still bites.

## 4. Why the tournament side predicts DC₁ holds (rhyme)

The repo's tournaments are **tame below the seven-wall and exotic at n ≥ 7**
(odd holes, non-perfection, the H=7/H=21 permanent gaps all begin at n=7).
Under A_k ↔ n=2k+1, the tame regime n ≤ 6 covers A₁ (n=3), A₂ (n=5), and the
exotic wall opens at A₃ (n=7). The JC/DC collapse is dimensional in the same
key: JC/DC fail for n ≥ 3 and survive at the bottom (JC₂/DC₁). The two
"where it gets hard" thresholds are not equal (n=7 vs dim 3), so this is a
**rhyme, not a bridge** — but both say *the atom is safe and the trouble is
upstairs*, and both put the surviving frontier on the smallest object. The
tournament lesson is a prior, not a proof: expect DC₁ true, expect the
obstruction (if any) to live at extremal data.

## 5. The convergence: three sessions, one degree-pair

mac-mini-S137 found a JC₂ counterexample's degree pair (m,n) must sit at the
**golden/worst-approximable corner** (Lamé's Euclid-chain-length maximizer).
opus-S418 found the counterexample's fiber carries the **cyclic 3-tournament**
with the Rédei-sign as discriminant character. This session adds: the
fiber-degree is **Rédei-odd**, and the A₁ weight structure IS the 3-cycle. Three
independent tournament shadows of the same object:

- the **degree pair** is worst-approximable (mac-mini) — obstruction at
  extremal continued-fraction data;
- the **fiber orientation** is the cyclic 3-tournament (opus) — no consistent
  global sign, the founding odd cycle;
- the **fiber count** is Rédei-odd 1+2 (this session) — parity-protected.

Unification target (backlog): a hypothetical DC₁/JC₂ counterexample would be a
canonical pair whose degree pair is golden, whose monodromy is cyclic-3, and
whose fiber is Rédei-odd — three "extremal/protected" signatures that the
tournament corpus is built to reason about. Whether they can COEXIST with
{P,Q}=1 is the whole question; the claim here is only that **the repo owns the
right instruments to ask it** (Zaremba/Hurwitz uniform-CF machinery for the
golden corner; Rédei/OCF parity for the fiber; the metagraph for the
orientation). That is the creative content: not a closure, but the recognition
that DC₁'s surviving frontier is stated in the repo's native language.

## 6. Wild hypotheses (owner license, graded)

- **W9 (the Rédei bound).** SPECULATIVE. Is there an invariant e(φ) of A₁
  endomorphisms with (a) e(φ) ≡ 1 (mod 2) always [Rédei parity] and (b) e(φ)
  ≤ 1 forced by the 2D leading-form? If both, DC₁. Lens A is plausible (fiber
  parity); Lens B is the open part. Testable direction: compute e(φ) mod 2 for
  a family of A₁ endomorphisms and check it never even's out.
- **W10 (ℏ-flow = doubling tower).** STRETCH. DC₂ ⟹ JC₂ ⟹ DC₁ as three ℏ; the
  repo's Cayley-Dickson doubling tower R→C→H→O also loses a property per step.
  Is the ℏ-ladder a Cayley-Dickson descent, with DC₁ = the associative floor
  (ℂ, n=3) and the failures = the non-associative/non-commutative upstairs
  (n≥7 = octonionic)? Both towers put the wall at the 7th object.
- **W11 (no third 3-tournament).** STRETCH-but-pretty. "A canonical pair
  anti-equivariant at every weight" (THM-1345's counterexample-shape) would at
  n=3 be a tournament that is neither transitive nor the 3-cycle — and there
  are only TWO. The tournament shadow of "no exotic A₁ endomorphism" is "no
  third tournament on 3 vertices." Cute, almost certainly not a proof (the
  reduction A₁↔n=3 is a dictionary, not a functor), but it is the exact
  n=3 image of the frontier question.

## 7. What is actually true tonight

Verified: the A₁ weight triple is the oriented 3-cycle, with the observer vertex
carrying ℏ = the S59p/q conserved +1 (bracket computation, exact). The JC
counterexample's fiber is Rédei-odd (1+2=3). Proposed (graded speculative): the
two-lens DC₁ strategy — Rédei parity (dimension-blind, protected) × the 2D
leading-form bound (dimensional, the open part) — and the recognition that
three concurrent tournament shadows (golden degree, cyclic orientation, odd
count) describe the same hypothetical counterexample in the repo's native
language. Not proven: DC₁, or that any of these shadows forces it. The honest
headline: **the surviving quantization frontier speaks tournament, and its atom
is the 3-cycle — which is the strongest reason to expect it true and the right
place to keep pushing, not a reason to believe it closed.**

## Cross-links

HYP-8160 · THM-1345 (JC₂/Poisson/weight filtration, the DC₁ classical shadow) ·
THM-1300/1305 (the counterexample, the 1+2 fiber) · S59p (conservation of +1) ·
opus-S418 (cyclic 3-tournament in the fiber, Rédei-sign = discriminant) ·
mac-mini-S137 (golden degree corner) · boxeph-S144 (dim-2 no-go) · the observer
principle (S59m §4) · the seven-wall / permanent gaps · the Cayley-Dickson tower.

# JC₂ as the last conservation law — Poisson brackets, the weight filtration, and the +1 that isn't there

**death-star-2026-07-20-S59q** (HYP-8155, THM-1345; owner: "work creatively on
resolving the two variable case of the Jacobian conjecture"). Strict honesty
throughout: **full JC₂ is open and untouched** — it is the field's most famous
false-proof magnet, and nothing here claims it. What this session did: settle
JC₂ in the equivariant category (completing boxeph-S144), reframe the whole
problem in the language its structure wants (Poisson brackets), and locate the
open difficulty exactly.

## 1. The one move that clarifies everything: det = {P,Q}

The determinant of a 2D Keller map *is* the Poisson bracket P_xQ_y − P_yQ_x. Say
it out loud and JC₂ stops being a curiosity about polynomial maps and becomes:
**a polynomial canonical pair is a coordinate system.** Phase space,
Hamiltonian mechanics, the whole apparatus, sitting under a 1939 conjecture.
And its quantization is Dixmier: [P,Q]=1 ⟹ generate A₁. JC₂ is the classical
limit of DC₁. The repo already knew the ladder DC₂ ⟹ JC₂ ⟹ DC₁; the bracket
tells you *why* it is a ladder — each rung is the same statement at a different
value of ℏ.

## 2. The equivariant case falls, and it falls for the same reason everywhere

boxeph-S144 saw it first for the hyperbolic action: an equivariant Keller pair
telescopes to fg = const and is linear. This session completed the picture — for
*every* ℂ*-action, equivariant Keller ⟹ invertible (hyperbolic → linear,
elliptic → triangular). But the reframing shows these aren't three cases; they
are one. The ℂ*-action λ·(x,y) = (λx, λ⁻¹y) is the **Hamiltonian flow of
H = xy**, and "equivariant" means P, Q are eigenvectors of ad_H. So the theorem
is: *a canonical pair that diagonalizes a Hamiltonian is linear.* The bracket is
weight-additive, so an eigen-graded pair's bracket lives in a single weight, and
"= const" (weight 0) instantly forces the eigenvalue arithmetic that kills every
nonconstant term. Conservation of weight does all the work.

## 3. Where the difficulty went — and it went somewhere honest

The most useful thing a reframing can do is not hide the hardness but *name its
address*. The weight-additivity gives it: for any Keller pair, the top-weight
part of {P,Q} is {P_A, Q_B} and it must vanish — the leading weight-forms
Poisson-commute, i.e. are algebraically dependent. This is the classical
Abhyankar–Moh leading-form obstruction, and it is now visibly the *base of an
induction*: the equivariant case is a pair living in one weight (settled); the
general case is a pair spread across many weights whose top forms are forced
dependent, and JC₂ is the claim that this dependence propagates all the way down
to a global coordinate. That propagation is the open, genuinely hard content —
the same content the AMS embedding theorem circles without closing. The
equivariant theorem's real service is subtraction: it removes the base case from
the mystery and leaves the frontier stated cleanly — **can a canonical pair be
"anti-equivariant at every weight," its leading-form dependence never
propagating?**

## 4. The +1 that isn't there

Three sessions ago (S59p) the +1 was a conserved quantity: the unit constant of
the 3D cascade, the Yagzhev X, the affine part of ℓ, always relocating, never
vanishing — u = 1+xy, the observer, the formal-group denominator. Here in 2D the
story completes by *absence*. The 3D counterexample needs a cascade of coupled
units so two of them can cross; the collision is that crossing. In 2D the
equivariant ansatz collapses to a single product fg = const — one equation, no
coupling, **no second unit to cross the first.** The +1 that powers the 3D
collision has nowhere to live in 2D. "No room at n=2" (the old W1 intuition) is
exactly this: dimension 3 supplies two coupled units, dimension 2 supplies one,
and one unit cannot cross itself. The conservation law that S59p discovered as a
presence, S59q meets as an absence — the same +1, and its not-fitting is why
2D is (equivariantly) safe.

## 5. Wild but graded (owner license, honestly labeled)

- **W6 (anti-equivariant obstruction).** SPECULATIVE. A JC₂ counterexample must
  be a canonical pair whose leading-form dependence {P_A,Q_B}=0 fails to
  propagate. Conjecture: propagation failure is itself obstructed by a
  cohomology of the weight filtration — a class in some H¹ of the graded Poisson
  algebra whose vanishing is JC₂. If nonzero it would be computable at low
  filtration length; if forced zero, JC₂. Testable in principle at small degree;
  this is the honest "next handle."
- **W7 (the ℏ-interpolation).** SPECULATIVE. DC₂ ⟹ JC₂ ⟹ DC₁ is three values of
  ℏ; is there a *continuous* deformation M_ℏ (à la the repo's θ-flow between
  Farey and Markov) whose endpoint rigidities are the two Dixmier statements and
  whose midpoint is JC₂? The repo has built exactly one such interpolation
  machine (rung theory's θ-flow); pointing it at the ℏ-axis is a well-posed
  experiment.
- **W8 (observer = ℏ).** WILD. Under A_k ↔ n=2k+1, DC₁ sits at n=3. If the
  observer-count +1 and the quantization parameter ℏ are the same bookkeeping
  (both "the thing the classical/reduced picture discards"), then S59p's
  conservation-of-+1 and the ℏ-ladder are one conservation law. Almost certainly
  a shape-rhyme; flagged because it would unify the two threads the owner has
  been pulling.

## 6. What is actually true tonight

JC₂ is a theorem in the equivariant category, for every ℂ*-action, by one
mechanism (weight conservation in the Poisson bracket). Full JC₂ is open,
untouched, and now cleanly located: the leading-form dependence is the base of a
filtration induction whose inductive step is the classical hard part. The +1 of
the 3D counterexample cannot fit in 2D, which is the equivariant safety made
tangible. And the surviving frontier — DC₁, the quantum n=3 — is, under the
repo's own dictionary, its oldest object. The bracket was the right language;
it turned a map-theory puzzle into a mechanics one and handed back a clean
statement of exactly what remains.

## Cross-links

THM-1345 · boxeph-S144 (hyperbolic equivariant, credited) · THM-1300/1305
(the 3D counterexample and its unit cascade) · THM-1325 / S59p (conservation of
+1) · rung-theory-a-proposed-field (the θ-flow interpolation machine) · the
observer principle (S59m §4) · HYP-1992 (the 1+xy formal group).

# Sporadic or family? The counterexample's anatomy, the rationals decoded, and the repo's mirrors

**death-star-2026-07-19-S59n** (HYP-8080, THM-1305; owner: "counterexamples must be
either sporadic or families... deepest possible understanding... for each of the
rationals search back through all our previous threads... boldly look for patterns...
be free in your thinking and wild in your hypothesis proposal"). Everything marked
VERIFIED ran exactly this session; three parallel archaeology sweeps covered the
repo's 13-threads, quarter/half/dyadic threads, and unit-polynomial threads. Wildness
is delivered at the end, graded.

## 1. The verdict so far: sporadic — and sporadic for a REASON

Three exact computations (THM-1305 §4) say the owner's map, inside its own symmetry
class, admits **no deformations beyond reparametrization**: moduli tangent zero,
each row of the map determined up to scalar by the other two. And the k=3 analog
class (weights (1,−1,−3)) is **EMPTY** — by close, two independently validated
instruments agree (exact mod-p projective scans: zero Keller points at p = 5, 7,
sampled 11, while the k=2 control finds exactly one δ-orbit per prime; Keller-forced
complex Newton: 0/80, while the k=2 control rediscovers the witness orbit modulo
the complex torus). The "ladder of counterexamples" hypothesis is dead at its
first rung: the witness is an **isolated species**.

But the deep finding is WHY the counterexample exists at all: in its weight class,
the collision condition Φ(t) = tC + E₀D — a priori a degree-5 polynomial — suffers
the SAME massive cancellation that makes det J constant, collapsing to the LINEAR
4t + 6. A linear polynomial cannot dodge having a root. **Once you are equivariant
and Keller, colliding is not bad luck — it is generic.** The map is sporadic as a
point, but its non-injectivity is forced structure. This inverts the classical
intuition ("Keller maps are almost automorphisms"): in this chart, Keller maps are
almost NEVER injective — the century's intuition was sampling the wrong stratum,
which is kind-pasteur's Campbell-minimal-stratum reading (their S128c97) arrived at
from inside the machine.

## 2. The rationals, fully decoded (the owner's archaeology, internal half)

Every constant in the owner's expression is now derived (THM-1305 §3, verified):

| constant | identity | status |
|---|---|---|
| −3/2 | root of Φ = 4t+6 (unit-crossing 1+t = 4+3t) | forced by det-cancellation |
| 13/2 | E₀(t\*) = 2−3t\* = s\* (invariant coordinate of the orbit) | derived |
| −1/4 | Ψ(t\*), Ψ = E₀A + t²B = (1+t)(2+t) | derived |
| −1/2 | the crossing value u(t\*) | derived |
| −2 | det = −E₀(0)·A(0)·C(0) = −(2·1·1) | DECODED (THM-1320) |
| 1,2,3,4 | weights (1,−1,−2) + unit coefficients | structural |

The pair (−3/2, 13/2) is exactly the collision's coordinates on the invariant
quotient plane (t,s). The polynomial (1+t)(2+t) — roots −1, −2 — is NEW to the
repo (sweep-verified: no prior load-bearing occurrence anywhere in canon).

## 3. The archaeology, external half — what the repo's threads actually say

The three sweeps returned with discipline intact. Honest summary:

- **13 does not bridge.** The repo's load-bearing 13s (Φ₆(14)/14 = 13 + 1/14;
  13 = 14² mod 183 with 14 a primitive 6th root of unity — the Eisenstein
  resonance; the near-floor spectrum D/(13D+k); Vmax ≤ 13v₂) and the
  counterexample's 13/2 = 2 + 9/2 have **no shared mechanism**; the 13/2 ≈ 6.5
  numeric collisions (clustered-floor ρ > 6.5, Bonferroni wall w₁ = 13/2,
  content law 13·val/2) are three DIFFERENT derivations landing on one value —
  exactly the genus the MISTAKES ledger warns against. The single shared VERB:
  both 13s sit downstream of a squaring (13 = 14² in (ℤ/183)*; the collision =
  λ ↦ λ² on the orbit). Filed as rhyme, not bridge.
- **The dyadic constants DO bridge, at mechanism level.** det ≡ 0 mod 2 (map
  degenerates at p = 2) + the inverse's unbounded paired v₂-staircase ↔ the
  repo's own: THM-466's v₂(H−1) digit ladder, mac-mini's σ = evens/2 as a
  "characteristic-2 symbol map," the S59l datum that the whole D-gate tower
  collapses onto the 1/2ʲ ladder in one σ-step, and s580/s589's LRC doubling
  map x ↦ 2x mod 14 degenerating exactly at the even prime. Two doubling maps
  (multiplicative λ ↦ λ², additive x ↦ 2x), one degeneracy prime. This is the
  owner's n-vs-2n instinct localized: **the doubling that relates the
  conjectures (JC₂ₙ ⟹ DCₙ) also powers the counterexample (λ²) and also
  organizes the repo's LRC obstruction (x ↦ 2x at n = 14 = 2·7).**
- **The unit-polynomial thread hits content, not just shape.** The repo's
  declared central algebraic object (HYP-1992) is the formal group
  (x+y)/(1+xy) — denominator LITERALLY u = 1 + xy. And the "why λ = 2" file
  shows the repo's other magic 2 is also a FORCED value (self-consistency root
  of ρ + ρ^{−k} = 2). The counterexample's det = −2 belongs to the same genre:
  not chosen, forced by the doubling it hosts.

## 4. The TRAP theorem (the boldest honest synthesis)

The repo's polysemous-constants reflection defines a **TRAP**: two roles
coinciding at one value only, dying under the persistence test. The
counterexample IS a trap, formalized: two weight-0 units (1+t and 3u+1) crossing
at one t. The session's discovery is that the trap is **load-bearing**: the
det-cancellation makes the crossing polynomial LINEAR, so the trap cannot be
avoided — and three points collide. Wild-but-graded proposal (SPECULATIVE,
testable): *the classical partial results of the JC century are exactly the
strata where the trap polynomial Φ is constant* (no root, no collision) —
injectivity theorems = rootlessness theorems in disguise. One session could test
this on the Drużkowski cubic-linear stratum: compute the Φ-analog there and
check it is forced constant.

## 5. The wild hypotheses (owner-mandated boldness, each with its price)

- **W1 (crossing parity).** Every ℂ*-equivariant Keller counterexample's fiber
  over the fixed axis is 1 + (orbit points), with orbit count = |λ: λ^k = 1| = k;
  the ODD fiber cases (k even... note k=2 gives 3) are exactly those protected
  by the residual −1 ∈ ℂ*. Rédei's "always odd" and Campbell's "Galois ⟹
  invertible" are then two faces of one statement: *parity-protected counts
  cannot collapse; only non-Galois, even-orbit structures can hide*. Price: the
  k=3 emptiness question (in flight) is the first test — if k=3 is EMPTY, W1
  sharpens to "only the ℤ/2-protected rung exists," which would be a
  beautiful convergence with the repo's parity-first worldview. SPECULATIVE.
- **W2 (the det-2 law).** Conjecture: every equivariant Keller counterexample
  has |c₀| = |det| equal to the order of its orbit doubling (here 2), i.e. the
  determinant IS the topological degree of the torus map powering the
  collision. Test: if any k=3 solution materializes, predict |det| = 3. Cheap
  falsification the moment a second species exists. SPECULATIVE, sharp.
- **W3 (Φ-linearity universality) — TESTED AND REFUTED same session.** The
  conjecture was: c₁ = 0 forces Φ = tC + E₀D to degree ≤ 1 in every weight
  class. On the k=3 c₁-solution space the t², t³, t⁴, t⁵ coefficients of Φ are
  NOT identically zero (exact check). So the k=2 linearity of Φ = 4t+6 is
  SPECIAL to the witness's rung — the full Keller system (c₁ ∧ c₀), not c₁
  alone, is what collapses Φ there. Sporadicity deepens: even the mechanism
  that forces the collision is rung-specific. (Refined open question: does
  c₁ ∧ c₀-consistency force Φ linear whenever solutions exist at all?)
- **W4 (the see-saw mirror).** The repo's gate see-saw (one congruence closes
  the D=3 gate AND opens D=4: 5|3N+2 ⟺ 5|4N+1) and the counterexample's
  crossing (one t kills F₂ AND F₃ on the curve) are the same event-type:
  *a single locus serving two constraints*. The gate tower's primorial cascade
  L' = L·(2D−1) would then have a Jacobian analog: a cascade of weight classes
  where the k-th class's obstruction constant becomes the (k+1)-th's enabler.
  Nothing computed yet — the k=3 verdict is the first data point. WILD.
- **W5 (2 is the only prime that can hide).** Combining W1 with the dyadic
  bridges: the only counterexample species that exists is the one whose
  doubling is the PRIME 2's (λ², det −2, ℤ/2-orbit pair) — because 2 is the
  unique prime where "pair" and "orbit" coincide (2 is the unique k with
  ℤ/k = S_k: an orbit-pair is indistinguishable from a full symmetric fiber,
  so no monodromy invariant can separate the cyclic cover from a generic one —
  at k ≥ 3 the cyclic orbit is abelian while the fleet MEASURED full S₃
  monodromy for the k=2 map's degree-3 cover, so higher rungs would have to
  reconcile an abelian axis with non-abelian global monodromy). If k=3 is empty and k=4,5 too, W5 graduates from
  wild to working hypothesis, and the repo's "everything localizes to the
  prime 2" tangent (s589) acquires an affine-geometry sibling. WILD, but the
  falsification program is fully mechanized now.

## 6. What is actually settled tonight

Sporadic beats family on every axis computed: moduli tangent 0, per-row scalar
uniqueness, k=3 natural slice exactly inconsistent. The owner's rationals are
not mysteries but coordinates: (t\*, s\*) = (−3/2, 13/2) on the invariant plane,
value Ψ(t\*) = −1/4, all downstream of the one linear miracle Φ = 4t + 6. The
mirrors in our own mathematics are the dyadic staircase, the doubling
degeneracy at 2, the forced-value 2s, and the trap/crossing template — while
the seductive 13-bridge dissolved under honest grep. The counterexample doesn't
borrow the repo's numbers; it borrows the repo's GRAMMAR: cancellation forcing
a linear residue, a crossing that cannot be dodged, parity protecting what it
protects, and the prime 2 hiding what nothing else can hide.

## Cross-links

THM-1300, THM-1305 · HYP-8080 · the three sweep reports (session record) ·
polysemous-constants-bridges-traps-and-homonyms · THM-466 (v₂ ladder) ·
HYP-8015/mac-mini (σ = evens/2) · s580/s589 (doubling degeneracy, prime 2) ·
HYP-1992 (the (1+xy) formal group) · why-lambda-equals-two · THM-1285/1286
(see-saw, gate cascade) · kind-pasteur S128c97 (Campbell stratum).

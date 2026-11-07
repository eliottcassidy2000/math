# Micro/macrocosm of the gate: the x+2 / x·2 duality in the marked resolvent, golden gate towers, and a tower odd-cycle formula

2026-11-06, mac-mini-2026-11-06-microcosm. Owner directive: "session thinking micro/macrocosm and merge in prior work on
how {3,10} in number theory and {15,25} in the square-sum problem may relate with {7-21} and the many prior session
around the permutohedron across different months. consider the much-discussed duality of x+2 and x*2 recursion on
tournaments and how An n-vertex tournament corresponds to a strong tournament on n + 2 vertices with a recognizable gate".

Canon:
- [THM-4900](../../01-canon/theorems/THM-4900-gate-and-doubling-in-the-marked-resolvent-golden-towers-and-tower-ocf.md)
  (PROVED);
- [HYP-9490](../hypotheses/HYP-9490-gate-towers-forget-the-core-up-to-cyclic-triangles.md) (OPEN, with an exact
  reduction).

Claims are typed as PROVED, KNOWN, CITED, FINITE-EXACT, STRUCTURAL, ANALOGY, COINCIDENCE, NUMEROLOGY, CONDITIONAL or
HEURISTIC.

## 0. What was already done today (credited, not re-derived)

The directive overlaps heavily with packets pushed earlier on 2026-11-06. Read in full or through one Explore summary
before any claim:

- **Session210.**
  - `graded-spectrum-gates`: the gate Γ(T) = a → b → T → a; recognizability (at most one gate from order 4);
    HC(Γ T) = H(T); H(Γ T) = (n+2)H(T) + m_(n−2); the cover-polynomial operator P_Γ(q) = q[P(q+1) + qP(q) − P(q−1)];
    finite reconstruction of P_T from nested-gate counts.
  - `structure-atlas`: q = (H − 1)/2, so 3 and 10 are the q-values of the forbidden 7 and 21.
  - `coherent-permutohedron`: tournaments are parallel-coherent orientations of P_(n−1), and Hamiltonian paths are
    their sinks.
- **Session211** `square-sum-port-bridge`: square-sum path counts p(Q23) = 3 and p(Q25) = 10, so 2p + 1 = 7, 21; a
  no-transport theorem; two-port insertion.
- **opus-2026-11-06-S82**: every odd N ≤ 2^35 except 7 and 21 is a Hamiltonian path count.
- **klein THM-4889 / HYP-9481**: exact spectra to order 11, organised by cyclic-triangle levels.
- **Older.**
  - Session35 recursions: E, D, the marked pair (p, F), DE ≠ E²D.
  - THM-447 (D), THM-291 (E), HYP-1881/1891 (the grid), THM-3322 (join law), THM-078 (odd-cycle form of P_T),
    THM-1370 (7 and 21 never occur).

**Found at the pre-commit priority check** (re-fetch of origin after this session's claim c137c698e6). These packets
answer the same directive. Where they overlap with this session they have priority, and the overlapping items below
are re-typed KNOWN.
- **Session211** (pushed 18:13–18:21, before this session's claim):
  - `recursion-bridge`: the apex-arc flip turns padding into the gate and creates exactly P_T(2) new permutohedron
    sinks. The first two gated values at core order n are 1 + 2ⁿ and 3 + 3·2^(n−1). The exact small gated spectrum is
    {3, 5, 9, 15, 17, 27} through 31. Padding, gating and doubling fail to commute, and D(Γ(T)) has no gate (its PROOF
    §5). The reflected double is an affine section of the doubled permutohedron, with square-face lifts.
  - `score-permutohedron`: the gate translates the source/sink face into the interior. It also gives literal gate-free
    HC = 7 and HC = 21 switches.
- **Session212** (pushed after the claim):
  - `cycle-bridge`, `cycle-spectrum`: every nonnegative integer is a Hamiltonian-cycle count, with all-order gate-free
    and feedback-number-two constructions.
  - `cycle-levels`: complete fixed-order cycle supports through order 9. HC = 7 first occurs at order 7 and HC = 21 at
    order 8.
- **klein THM-4889** already records complete (H, HC) data through order 10.
- **opus-2026-11-06-S84** reserved THM-4890/4891 and HYP-9482 (micro/macro one-apex collisions; permutohedron
  relaxations). The reservation has no content yet, so no comparison is possible.

The present session therefore works only on lanes the survey found empty:
- the gate inside Session35's marked resolvent, and its interaction with D;
- the macro limits of nested gates, of doubling and of padding;
- the "micro/macro" reading of the cover polynomial under nested gates;
- the permutohedron and Fano readings listed below.

## 1. The x+2 / x·2 duality, made exact on the marked pair (THM-4900 (1), (3), (5))

**+2 is linear, ×2 is quadratic.**
- The padding E and the strong gate Γ both act on (N, P) = (uᵀ adj(zI − M) u, det(zI − M)) by symmetric 2 × 2
  matrices:

      G_E = [[z²+1, 2z], [2z, z²+1]],    G_Γ = [[z²+9, 2z], [2z, z²+1]].

  Their eigenvalues are (z ± 1)² and (√(z²+4) ± 1)². Flipping the single apex arc s ↔ t shifts z² to z² + 4 inside
  the square root. It leaves the characteristic polynomial unchanged.
- The doubling D acts quadratically, by P_D = 2ⁿP(w)² and N_D = 2ⁿ(z/w)N(w)P(w) with w² = (z² + 1)/2 (Session35).

**The self-similar point.** The substitution fixes z = ±1.
- At z = 1 the doubling is invisible: φ(T) = F_T(1) = uᵀ(I + MᵀM)⁻¹u is an exact invariant, φ(D T) = φ(T).
- The padding sends φ to 1, and the gate acts by the golden Möbius map φ ↦ (5φ + 1)/(φ + 1).
- So at the self-similar point the two recursions commute, φ(DΓT) = φ(ΓDT). The integer Baumslag–Solitar relation
  d∘a = a²∘d fails there, for every T.
- At the tournament level the relation fails outright: D(Γ(T)) never has a gate, because every score of a double of
  a strong tournament is at least 2. This is KNOWN: Session211 `recursion-bridge`, PROOF §5, proved it earlier the same
  day.
- New here: the scalar shadow. The characteristic polynomials agree, P_(DΓT) = P_(DET), and at z = 1 the two moves
  commute.

**Reading** (STRUCTURAL).
- Session35 found that "additive growth is a Möbius boundary transfer, doubling is a quadratic spectral substitution".
  THM-4900 adds the strong (gated) version of the additive move, and locates the one spectral point where the
  multiplicative move does nothing.
- On orders, the grid n = 2^r·(odd) is an exact arithmetic partition (Session35 §4). On tournaments the two moves
  neither commute nor satisfy the BS relation, except in the scalar shadow at z = 1, where they commute.

## 2. The gate as a microcosm (THM-4900 (2), (4), (6); HYP-9490)

**Nested gates.**
- Γ^k(T) has a closed-form characteristic polynomial. For every core the marked response converges to
  (2 + √(z²+4))/z. That is the metallic mean of index 4/z:
  - golden at z = 4;
  - silver at z = 2;
  - φ³ = 2 + √5 at the self-similar point z = 1.
- Over one vertex, the tower is the transitive tournament on 2k+1 vertices with its k concentric chords reversed. Its
  spectrum is explicit: {0, ±i√3} ∪ {±i√(4 + cot²((2j+1)π/(2(2k−1))))}.
- At z = 1, in the coordinate y = (φ − 1)/2, the gate is the square of the continued-fraction map y ↦ 1 + 1/y.
  - The one-vertex tower runs through the Fibonacci convergents F_(2k)/F_(2k−1).
  - Equivalently, φ(Γ^k K1) = F_(2k+2)/F_(2k−1) = 1, 3, 4, 21/5, 55/13, …
- **The macrocosm forgets the microcosm**:
  - the gate tower forgets the core completely (the metallic limit);
  - the doubling tower remembers exactly one number, F → z·φ(T);
  - the padding tower forgets it, F → 1.

**The condensation multiplies** (THM-3322's join law, re-proved in THM-4900(4)). The marked response is the
relativistic sum, F_(A⇒B) = (F_A + F_B)/(1 + F_A F_B), of the responses of the strong components.
- At z = 1 a one-vertex component is absorbing, φ(K1) = 1.
- So φ ≠ 1 only for tournaments all of whose strong components are nontrivial, e.g. φ(C3 ⇒ C3) = 3/5.
- This is the resolvent shadow of H being multiplicative over the condensation.

**The tower odd-cycle formula** (PROVED):

    H(Γ^r T) = Σ_j sw_j(T)·H(Γ^r TT_(n−2j)),    sw_1 = 2c3.

- The core enters only through its size-weighted odd-cycle statistics, and the tower only through universal weights
  computed on transitive cores.
- At r = 0 it is the odd-cycle formula itself.

**HYP-9490.** As r → ∞ the core is forgotten except through c3 to first order:

    H(Γ^r T)/H(Γ^r T′) = 1 + 2Δc3·ρ_(n,r)·(1 + o(1)).

- Checked on every isomorphism class of orders 4–7 at total order 401, with residuals at most 1.7e-6 against
  ρ ≈ 2e-4, and to total order 1201 for six cores.
- **The tower is a diffusion** (PROVED, HYP-9490 Lemma A). The gate operator L is adjoint, under
  ⟨P, Φ⟩ = [P(∂)Φ](0), to the Sturm–Liouville operator M = ∂² + 2 sinh(y)∂ = e^(−2 cosh y)∂e^(2 cosh y)∂. So
  H(Γ^r T) = [P_T(∂)·M^r e^y](0).
  - Macrocosm: one diffusion applied r times to e^y.
  - Microcosm: the core's cover polynomial, read as a differential operator at the origin.
- The weights are monotone (Lemma B), and part (a) of the conjecture implies the rest (Corollary C).
- The one-vertex tower's path counts 1, 3, 15, 123, 1479, 23955, 496479, 12725547, 393325815 have c3 = k². They are
  not in OEIS (checked 2026-11-06; the one-defect counts 0, 18, 372, 7686, 185064 are not in OEIS either).

## 3. {3,10}, {15,25}, {7,21}: where the gate sits

- **The q-coordinate** (Session210 atlas; KNOWN). q = (H − 1)/2 sends the forbidden counts 7 and 21 to 3 and 10. In
  q, the dominance join is q + q′ + 2qq′.
- **The same affine normalisation** (STRUCTURAL). y = (φ − 1)/2 turns the gate's action on the marked resolvent into
  the doubled golden continued fraction. The coordinate (X − 1)/2 is natural on both sides: for path counts it linearises
  the join, for the resolvent it conjugates the gate to F².
- **Square sums** (Session211; KNOWN). The first two square-sum path counts after the onsets 23 and 25 are exactly 3
  and 10. So 2p + 1 would give 7 and 21, and Session211 proves no map F with H(F(Q_n)) = 2p(Q_n) + 1 exists.
- **The gate and the gaps** (KNOWN; a one-line corollary of THM-1370 and Session210, used by Session211/212).
  - HC(Γ T) = H(T), so the Hamiltonian-cycle counts of gated strong tournaments are exactly the path spectrum, which
    misses 7 and 21.
  - Every strong tournament with exactly 7 or exactly 21 Hamiltonian cycles is therefore gate-free. Such tournaments
    exist:
    - Session211 `score-permutohedron` gives literal switches;
    - Session212 gives all-order gate-free constructions for every count;
    - Session212 `cycle-levels` and THM-4889 give the least orders, 7 and 8.
  - This session adds only the class counts at those least orders (FINITE-EXACT, `hc_gate_free.out`, gentourng):
    - HC = 7 at order 7: 24 strong classes, all gate-free;
    - HC = 21 at order 8: 228 strong classes, all gate-free.
  - The cycle spectrum has no gaps at all: every nonnegative integer occurs (Session212, PROVED). The gate is
    precisely the local structure that transports the path gaps to cycles.
- **The floor of gated counts** (KNOWN).
  - Session210 (5) gives H(Γ T) = H(T) + P_T(2), where P_T(2) = Σ_(S⊆V) h(T[S])·h(T[V∖S]) with h(∅) = 1.
  - Session211 `recursion-bridge` gives the sharp floors 1 + 2ⁿ (transitive core) and 3 + 3·2^(n−1) (exactly one
    C3 component), and the exact small gated spectrum {3, 5, 9, 15, 17, 27} through 31.
  - So small counts such as 7 and 21 can never be gated.
- **Fano and the duality mod 7** (STRUCTURAL; NUMEROLOGY with respect to the gaps).
  - Read modulo 7, the two recursions generate ⟨x ↦ x + 2, x ↦ 2x⟩ = {x ↦ 2^k x + b}. That is the Frobenius group F21,
    i.e. Aut(QR7), because ⟨2⟩ = {1, 2, 4} is the set of quadratic residues mod 7.
  - QR7's out-neighbourhoods x + {1, 2, 4} are the lines of the cyclic Fano plane, and F21 is the normaliser of a Singer
    cycle in PSL(3,2) (index 8).
  - So the x+2 / x·2 recursion read mod 7 is exactly the symmetry group of the Paley tournament and the cyclic Fano
    plane: 7 points and 21 symmetries.
  - This is not special to 7. The same happens for every prime p where 2 generates exactly the quadratic residues:
    p = 23 (|Aut QR23| = 253), 47, 71, 79, …
  - No mechanism links this to the H-gaps 7 and 21. Tournament facts recorded elsewhere: H(QR7) = 189 = 9·21, and
    Aut acts freely on paths.

## 4. The permutohedron (PROVED and verified; KNOWN in content)

The face formula below is new only as a presentation. The ordered path-cover polynomial and the forward-adjacency
distribution determine each other at fixed order (Session210's binomial transform). So A_Γ(T) is already determined
by Session210's operator applied to P_T.

**The gate fibres the permutohedron of n+2 letters over the faces of the core's permutohedron.** Inserting the two
gate letters s, t into a word over V splits it into an ordered 3-partition (A, B, C) of V: the letters before, between
and after s and t, with blocks possibly empty. These 3-partitions index the faces P_A × P_B × P_C of the permutohedron
of V. For the descent polynomials A_T(x) = Σ_k a_k(T)x^k (a_k = words with exactly k backward adjacencies):

    A_Γ(T)(x) = Σ_((A,B,C)) A_T[A](x)·A_T[B](x)·A_T[C](x)·( x^([A≠∅] + [B=∅] + [C≠∅]) + x^(2[B≠∅]) ).

- The two summands are the orders s…t and t…s; the exponents count the forced backward junctions.
- Verified against direct enumeration on random cores of orders 1–5 (`gate_descent_face_formula.out`).
- Its constant term is Session210's H(Γ T) = (n+2)a_0 + a_1.
- In Session210's coherent orientation of the permutohedron, a word's out-degree is its number of backward adjacencies.
  So:

      H(Γ T) = (n+2)·#(sinks) + #(vertices of out-degree 1)

  in T's coherent orientation. This restates Session210's H(Γ T) = (n+2)H(T) + m_(n−2) in its own coherent
  language. Compare Session211's count of the P_T(2) sinks created by the apex flip.

**Across months.** The permutohedron has appeared in March (THM-207, S89b: transitive tournaments; zonotope and
quotientope broadcasts), June (HYP-2128 trinity), July (THM-1825: Landau's theorem as a moment map), October (score
lattice points, THM-4752) and November (Session192 weighted cuts; Session210 coherent orientations; Session211 score
permutohedron).
- The face formula is the permutohedral form of the x+2 move: adding the two gate letters fibres the big
  permutohedron over the 3-partition faces of the small one.
- Session211 already gave the face statements on either side of it:
  - the gate pushes the source/sink face into the interior;
  - the reflected double is an affine section with square-face lifts.
- A face-sum formula for the full descent polynomial of D(T) is still missing (OPEN).

## 5. Data (FINITE-EXACT)

- **Smallest φ over strong tournaments** (`phi_strong_minimum.out`, nauty gentourng):

  | Order n | 5 | 6 | 7 | 8 | 9 |
  |---|---|---|---|---|---|
  | min φ | 2 | 1 | 2/3 | 5/9 | 1/2 |
  | Minimising classes | 1 | 1 | 1 | 1 | 2 |

  - The maximum is n at odd n (regular tournaments).
  - The minimisers at n ≤ 8 have scores 1, 1, 2, …, n−2, n−2.
  - No conjecture: the n = 9 minimisers change shape.
- **Gated count spectra** (`gated_spectrum.out`): the values H(Γ T) for cores of orders 1–7. They start
  3; 5; 9, 15; 17, 27, 37; ….
  - For an order-n core the minimum is 2ⁿ + 1, and the second value is 3 + 3·2^(n−1).
  - Both agree with Session211's sharp floors and with its exact small spectrum through 31 (KNOWN; the table only
    extends the range).

## 6. Ledger

| Claim | Type |
|---|---|
| Gate linear on (N, P), eigen-coordinates, P_Γ = P_E | PROVED |
| Closed-form towers, metallic limit, explicit one-vertex spectrum | PROVED |
| z = 1 doubling invariance of φ, φ ∈ (0, n], golden gate map, Fibonacci orbit | PROVED |
| Relativistic addition over the condensation | PROVED here; KNOWN in substance (THM-3322) |
| Commutation at z = 1 without the BS relation; P_(DΓT) = P_(DET) | PROVED |
| D(Γ T) gate-free | KNOWN (Session211 `recursion-bridge` §5; re-proved here) |
| Tower odd-cycle formula; sw_1 = 2c3 at every order | PROVED (from Session210 + THM-078) |
| Gated floors 2ⁿ + 1 and 3 + 3·2^(n−1); HC ∈ {7, 21} forces gate-freeness | KNOWN (Session210 (5), Session211, Session212) |
| Least orders 7, 8 for HC = 7, 21 | KNOWN (THM-4889, Session212); class counts 24, 228 FINITE-EXACT here |
| Permutohedral face formula for A_Γ(T) | PROVED + FINITE-EXACT; KNOWN in content (Session210 operator) |
| HYP-9490: macro forgetting up to cyclic triangles | OPEN; FINITE-EXACT evidence |
| Diffusion form H(Γ^r T) = [P_T(∂) M^r e^y](0); weight monotonicity; (a) ⇒ (b) | PROVED (HYP-9490 Lemmas A–C) |
| x+2 / x·2 mod 7 = F21 = Aut(QR7), the cyclic Fano symmetry | STRUCTURAL |
| Links of F21 to the gaps 7, 21; φ³ = THM-224's eigenvalue; 7-vertex tower φ = 21/5 | NUMEROLOGY / COINCIDENCE |

## 7. Next obligations

- **HYP-9490.**
  - Prove ρ_(d,r) → 0 from the walk representation, or from the diffusion M of Lemma A. Everything else then follows.
  - Fix the rate: the data suggest (log N)²/N², but the normalised value is still drifting (0.90 → 0.86 over
    N = 243 → 1201).
  - Test the second order (sw_2) at N = 4001 with a C value recursion.
- **A face formula for D.** The skew-Sylvester double of a word over V has twin letters. Starting from Session211's
  affine reflected section, is there a face-sum description of the full A_(D T)(x)?
- **φ beyond scalars.** φ is a doubling invariant that the gate moves golden-ly. Is it a quotient of a finer
  doubling-invariant marked state, say the full matrix-valued boundary resolvent of Session35's next test, under which
  the gate is still linear?
- **The 24 and 228 least-order classes with HC = 7, 21.** Session212 constructs gate-free witnesses in feedback
  number two. Do all least-order classes have feedback vertex number two, or does a different local structure appear?

## 8. Files

`04-computation/experiments/microcosm_gate_resolvent_20261106/`:
- **THM-4900:** `gate_basics.py`, `gate_resolvent_exact.py`, `tower_limits.py`, `tower_spectrum.py`, `tower_ocf.py`.
- **HYP-9490:** `cover_tower_values.py`, `c3_first_order_law.py`, `tower_weights_diffusion.py` (Lemmas A–C).
- **Section 4:** `gate_descent_face_formula.py`.
- **Data:** `tower_k1_paths.c`, `gated_spectrum.c`, `phi_strong_minimum.py`, `hc_gate_free.c`.
- **`run_all.sh`** reproduces every output (`.out` files committed).

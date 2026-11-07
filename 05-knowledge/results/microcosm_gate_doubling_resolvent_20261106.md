# Micro/macrocosm of the gate: the x+2 / x·2 duality in the marked resolvent, golden gate towers, and how deep towers forget their core

2026-11-06, mac-mini-2026-11-06-microcosm. Owner directive: "session thinking micro/macrocosm and merge in prior work on
how {3,10} in number theory and {15,25} in the square-sum problem may relate with {7-21} and the many prior session
around the permutohedron across different months. consider the much-discussed duality of x+2 and x*2 recursion on
tournaments and how An n-vertex tournament corresponds to a strong tournament on n + 2 vertices with a recognizable gate".

Canon:
- [THM-4900](../../01-canon/theorems/THM-4900-gate-and-doubling-in-the-marked-resolvent-golden-towers-and-tower-ocf.md)
  (PROVED; statements 1–6 audited by audit H1, statement 7 added afterwards with a follow-up check);
- [HYP-9490](../hypotheses/HYP-9490-gate-towers-forget-the-core-up-to-cyclic-triangles.md) (OPEN: only the sharp
  forgetting rate);
- [audit H1](microcosm_audit_H1_20261106.md): PASS WITH CORRECTIONS, all applied (MISTAKE-729).

Claims are typed as PROVED, KNOWN, CITED, FINITE-EXACT, STRUCTURAL, ANALOGY, COINCIDENCE, NUMEROLOGY, CONDITIONAL or
HEURISTIC. Notation: φ = (1+√5)/2 is the golden ratio; the tournament invariant F_T(1) is called s(T), as in THM-1950.

## 0. Prior and concurrent work (credited, not re-derived)

**Read before the claim.**
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
- **klein THM-4889 / HYP-9481**: exact spectra to order 11, organised by cyclic-triangle levels; complete (H, HC) data
  through order 10.
- **Older.**
  - Session35 recursions: E, D, the marked pair (p, F), DE ≠ E²D, and the all-depth padding law.
  - THM-447 (D), THM-291 (E), HYP-1881/1891 (the grid), THM-3322 (join law), THM-078 (odd-cycle form of P_T, odd n),
    THM-1370 (7 and 21 never occur).

**Found at the pre-commit priority check** (re-fetch of origin after this session's claim c137c698e6). These answer
the same directive. Where they overlap with this session they have priority, and the overlapping items are typed KNOWN.
- **Session211** (pushed 18:13–18:21, before this session's claim):
  - `recursion-bridge`: the apex-arc flip turns padding into the gate and creates exactly P_T(2) new permutohedron
    sinks. The first two gated values at core order n are 1 + 2ⁿ and 3 + 3·2^(n−1). The exact small gated spectrum is
    {3, 5, 9, 15, 17, 27} through 31. Padding, gating and doubling fail to commute, and D(Γ(T)) has no gate (its PROOF
    §5). The reflected double is an affine section of the doubled permutohedron, with square-face lifts.
  - `score-permutohedron`: the gate translates the source/sink face into the interior. It also gives literal gate-free
    HC = 7 and HC = 21 switches.
  - `arithmetic-cover-bridge`: the gate operator beyond tournaments.
- **Session212** (pushed after the claim):
  - `cycle-bridge`, `cycle-spectrum`: every nonnegative integer is a Hamiltonian-cycle count, with all-order gate-free
    and feedback-number-two constructions.
  - `cycle-levels`: complete fixed-order cycle supports through order 9. HC = 7 first occurs at order 7 and HC = 21 at
    order 8. Its class histograms already give 24 classes with HC = 7 at order 7, and 118 and 228 classes with
    HC = 7 and HC = 21 at order 8.
- **opus-2026-11-06-S84** has since pushed THM-4890 (one-apex collision balance) and THM-4891 (face-local coherence
  relaxations of the permutohedron). Neither overlaps this packet. THM-4891 belongs in §4's "Across months" list.
- **Session213** (checkpoints 1–4, pushed during this session): cycle-to-path arithmetic and path gadgets. No
  overlap.

**Found by audit H1** (missed by this session's survey; credited throughout, MISTAKE-729).
- **THM-1950** (klein-2026-07-21-S400; with THM-4051 and THM-3729): the total inverse response s(T) = 𝟙ᵀ(I+K)⁻¹𝟙,
  which is this session's F_T(1). It comes with its bounds, the regularity characterisation, the velocity-addition law,
  and the strong minima to n = 8.
- **Session150** `transfer` and **Session159** `profile_operator`, `cycle_profiles`, `lagrange` (2026-11-02/03):
  - the gate (triangle-wrapper) operator and its positive differential form;
  - the all-order cycle-packing form of P_T (sw_1 = 2c3, sw_2 = 2c5 + 4·pairs);
  - the one-vertex family 3, 15, 123, …, 496479;
  - the transitive-core towers as exponential sums Σ_q c_(r,q) q^n.
- **F21 in the repo:** THM-4556 (remark), the 2026-10-06 notes `mersenne_switch_parity_f21_compression` and
  `seven_twentyone_mersenne_openai_math` (F21 = doubling plus translation mod 7), THM-4644(3) (Paley primes), and
  HYP-9162/THM-4557 (the doubling tower over C3 keeps F21).

What remains new here:
- the gate inside Session35's marked resolvent, its towers and its interaction with D at z = 1;
- the forgetting theorem for deep gate towers (THM-4900(7)).

## 1. The x+2 / x·2 duality, made exact on the marked pair (THM-4900 (1), (3), (5))

**+2 is linear, ×2 is quadratic.**
- The padding E and the strong gate Γ both act on (N, P) = (uᵀ adj(zI − M) u, det(zI − M)) by symmetric 2 × 2
  matrices:

      G_E = [[z²+1, 2z], [2z, z²+1]],    G_Γ = [[z²+9, 2z], [2z, z²+1]].

  Their eigenvalues are (z ± 1)² and (√(z²+4) ± 1)². Flipping the single apex arc s ↔ t shifts z² to z² + 4 inside
  the square root. It leaves the characteristic polynomial unchanged.
- The doubling D acts quadratically, by P_D = 2ⁿP(w)² and N_D = 2ⁿ(z/w)N(w)P(w) with w² = (z² + 1)/2 (Session35).

**The self-similar point.** The substitution fixes z = ±1.
- At z = 1 the doubling is invisible. THM-1950's inverse response s(T) = F_T(1) = uᵀ(I + MᵀM)⁻¹u is an exact
  doubling invariant, s(D T) = s(T) (new).
- The padding sends s to 1, and the gate acts by the golden Möbius map s ↦ (5s + 1)/(s + 1) (new).
- So at the self-similar point the two recursions commute, s(DΓT) = s(ΓDT). The integer Baumslag–Solitar relation
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

## 2. The gate as a microcosm (THM-4900 (2), (4), (6), (7); HYP-9490)

**Nested gates.**
- Γ^k(T) has a closed-form characteristic polynomial. For every core the marked response converges to
  (2 + √(z²+4))/z. That is the metallic mean of index 4/z:
  - the golden ratio φ at z = 4;
  - silver at z = 2;
  - φ³ = 2 + √5 at the self-similar point z = 1.
- Over one vertex, the tower is the transitive tournament on 2k+1 vertices with its k concentric chords reversed. Its
  spectrum is explicit: {0, ±i√3} ∪ {±i√(4 + cot²((2j+1)π/(2(2k−1))))}.
- At z = 1, in the coordinate y = (s − 1)/2, the gate is the square of the continued-fraction map y ↦ 1 + 1/y.
  - For k ≥ 1, the one-vertex tower runs through the convergents F_(2k)/F_(2k−1) of φ that lie below φ.
  - Equivalently, s(Γ^k K1) = F_(2k+2)/F_(2k−1) = 1, 3, 4, 21/5, 55/13, 72/17, … → φ³.
- **The macrocosm forgets the microcosm** (resolvent level):
  - the gate tower forgets the core completely (the metallic limit);
  - the doubling tower remembers exactly one number, F → z·s(T);
  - the padding tower forgets it, F → 1.

**The condensation multiplies.** This is THM-3322's join law, and at z = 1 THM-1950's velocity-addition law;
THM-4900(4) re-proves it. The marked response is the relativistic sum, F_(A⇒B) = (F_A + F_B)/(1 + F_A F_B), of the
responses of the strong components.
- At z = 1 a one-vertex component is absorbing, s(K1) = 1.
- So s ≠ 1 only for tournaments all of whose strong components are nontrivial, e.g. s(C3 ⇒ C3) = 3/5.
- This is the resolvent shadow of H being multiplicative over the condensation.

**The tower odd-cycle formula** (PROVED; ingredients KNOWN: Session210 Theorem 4, Session159 cycle_profiles (1)):

    H(Γ^r T) = Σ_j sw_j(T)·H(Γ^r TT_(n−2j)),    sw_1 = 2c3.

- The core enters only through its size-weighted odd-cycle statistics, and the tower only through universal weights
  U_(d,r) computed on transitive cores.
- At r = 0 it is the odd-cycle formula itself.

**Deep towers forget the core, except through c3** (PROVED, THM-4900(7)).
- Let ρ_(d,r) = U_(d−2,r)/U_(d,r). Then ρ_(d,r) ≤ 1/(r+1), uniformly in d. The proof is an induction on the
  first-step recursion of the weights, in which Pascal's rule supplies a slack term dominating U_(d−2,r).
- Consequently, for cores of the same order n,

      H(Γ^r T)/H(Γ^r T′) = 1 + 2Δc3·ρ_(n,r) + O(ρ_(n,r)/r),

  with an explicit constant.
  - The error is not 2Δc3·ρ·o(1). Audit H1 found that the first version's form fails when Δc3 = 0 (MISTAKE-729).
  - Checked exactly on all 134,584 (core pair, depth) triples of orders 4–6, where the bound is never more than 39%
    used. Also checked at total order 401 on every class of orders 4–7.
- **Operator form.** The gate operator L is adjoint, under ⟨P, Φ⟩ = [P(∂)Φ](0), to the Sturm–Liouville operator
  M = ∂² + 2 sinh(y)∂ = e^(−2 cosh y)∂e^(2 cosh y)∂. So H(Γ^r T) = [P_T(∂)·M^r e^y](0).
  - This is KNOWN in content: M^r e^y = Σ_q c_(r,q)e^(qy) with Session159's coefficients (lagrange (7)).
  - Macrocosm: one diffusion generator applied r times to e^y.
  - Microcosm: the core's cover polynomial, read as a differential operator at the origin.
- **Two limits of one family.** For fixed depth and a growing core, Session159 has H(Γ^r TT_n) ~ r!(r+1)^n: the core
  dominates. For a fixed core and growing depth, the tower dominates and the core is forgotten.
- **HYP-9490 (OPEN)** is now only the sharp rate. The data suggest ρ ≈ (log N)²/N², against the proved 1/(r+1).
- The one-vertex tower counts 1, 3, 15, 123, 1479, 23955, 496479 are in Session159 `profile_operator`, with
  v_3(h) = 1 and a mod-9 law. The first polynomials are in `cycle_profiles`, and `lagrange` continues the series to
  12725547. This session adds 393325815.
  - These tournaments have c3 = k².
  - The counts, and the one-defect counts 0, 18, 372, 7686, 185064, were not found in OEIS (checked 2026-11-06).

## 3. {3,10}, {15,25}, {7,21}: where the gate sits

- **The q-coordinate** (Session210 atlas; KNOWN). q = (H − 1)/2 sends the forbidden counts 7 and 21 to 3 and 10. In
  q, the dominance join is q + q′ + 2qq′.
- **The same affine normalisation** (ANALOGY). y = (s − 1)/2 turns the gate's action on the marked resolvent into the
  doubled golden continued fraction. The coordinate (X − 1)/2 appears on both sides:
  - for path counts it writes the join as q + q′ + 2qq′ (1 + 2q is multiplicative);
  - for the resolvent it conjugates the gate to the square of y ↦ 1 + 1/y.
- **Square sums** (Session211; KNOWN). The square-sum path counts at the onsets 23 and 25 are exactly p(Q23) = 3 and
  p(Q25) = 10, after p(Q15) = p(Q16) = p(Q17) = 1. So 2p + 1 would give 7 and 21, and Session211 proves that no map F
  with H(F(Q_n)) = 2p(Q_n) + 1 exists.
- **The gate and the gaps** (KNOWN; a one-line corollary of THM-1370 and Session210, used by Session211/212).
  - HC(Γ T) = H(T), so the Hamiltonian-cycle counts of gated strong tournaments are exactly the path spectrum, which
    misses 7 and 21.
  - Every strong tournament with exactly 7 or exactly 21 Hamiltonian cycles is therefore gate-free. Such tournaments
    exist:
    - Session211 `score-permutohedron` gives literal switches;
    - Session212 gives all-order gate-free constructions for every count;
    - Session212 `cycle-levels` and THM-4889 give the least orders, 7 and 8.
  - The class counts at those least orders are already in Session212 `cycle-levels` (n7_hist.tsv, n8_hist.tsv). This
    session reproduces them with gentourng (`hc_gate_free.out`) and checks the immediate gate-freeness:
    - HC = 7 at order 7: 24 strong classes, all gate-free;
    - HC = 21 at order 8: 228 strong classes, all gate-free.
  - The cycle spectrum has no gaps at all: every nonnegative integer occurs (Session212, PROVED). The gate is
    precisely the local structure that transports the path gaps to cycles.
- **The floor of gated counts** (KNOWN).
  - Session210 (5) gives H(Γ T) = H(T) + P_T(2), where P_T(2) = Σ_(S⊆V) h(T[S])·h(T[V∖S]) with h(∅) = 1.
  - Session211 `recursion-bridge` gives the sharp floors 1 + 2ⁿ (transitive core) and 3 + 3·2^(n−1) (exactly one C3
    component), and the exact small gated spectrum {3, 5, 9, 15, 17, 27} through 31. Session159 `lagrange` (6) already
    had H(Γ TT_n) = 1 + 2ⁿ.
  - So small counts such as 7 and 21 can never be gated.
- **Fano and the duality mod 7.** This is KNOWN, classically; in the repo it is in THM-4556's remark and the 2026-10-06
  F21 notes. It is an ANALOGY as a reading of the tournament recursions, and NUMEROLOGY with respect to the gaps.
  - Read modulo 7, the two recursions generate ⟨x ↦ x + 2, x ↦ 2x⟩ = {x ↦ 2^k x + b}. That is the Frobenius group F21,
    i.e. Aut(QR7), because ⟨2⟩ = {1, 2, 4} is the set of quadratic residues mod 7.
  - QR7's out-neighbourhoods x + {1, 2, 4} are the lines of the cyclic Fano plane, and F21 is the normaliser of a Singer
    cycle in PSL(3,2) (index 8).
  - So the x+2 / x·2 recursion read mod 7 is exactly the symmetry group of the Paley tournament and the cyclic Fano
    plane: 7 points and 21 symmetries.
  - This is not special to 7. The same happens for every prime p ≡ 3 (mod 4) where 2 generates exactly the quadratic
    residues: p = 23 (|Aut QR23| = 253), 47, 71, 79, 103, …. By THM-4644(3), the doubling-invariant circulant
    tournaments on Z_p are then exactly QR_p and its converse. For p ≡ 1 (mod 4) (17, 41, 97, …) QR_p is a graph,
    not a tournament.
  - By HYP-9162/THM-4557, the doubling tower over C3 reaches QR7 and keeps Aut = F21 at every level.
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

**Across months.** The permutohedron has appeared in:
- March: THM-207 and S89b (transitive tournaments; zonotope and quotientope broadcasts);
- June: HYP-2128 (trinity);
- July: THM-1825 (Landau's theorem as a moment map);
- October: score lattice points (THM-4752);
- November: Session192 weighted cuts, Session210 coherent orientations, Session211 score permutohedron, and S84's
  THM-4891 (face-local coherence relaxations).

The face formula is the permutohedral form of the x+2 move: adding the two gate letters fibres the big permutohedron
over the 3-partition faces of the small one.
- Session211 already gave the face statements on either side of it:
  - the gate pushes the source/sink face into the interior;
  - the reflected double is an affine section with square-face lifts.
- A face-sum formula for the full descent polynomial of D(T) is still missing (OPEN).

## 5. Data (FINITE-EXACT)

- **Smallest s over strong tournaments** (`phi_strong_minimum.out`, nauty gentourng). The values for n ≤ 8 are KNOWN
  from klein's S400 reflection and THM-4051; n = 9 and the class counts are new.

  | Order n | 5 | 6 | 7 | 8 | 9 |
  |---|---|---|---|---|---|
  | min s | 2 | 1 | 2/3 | 5/9 | 1/2 |
  | Minimising classes | 1 | 1 | 1 | 1 | 2 |

  - The maximum is n at odd n (regular tournaments).
  - The minimisers at n ≤ 8 have scores 1, 1, 2, …, n−2, n−2.
  - No conjecture: the n = 9 minimisers change shape.
- **Gated count spectra** (`gated_spectrum.out`): the values H(Γ T) for cores of orders 1–7. They start
  3; 5; 9, 15; 17, 27, 37; ….
  - For an order-n core the minimum is 2ⁿ + 1, and, for n ≥ 3, the second value is 3 + 3·2^(n−1).
  - Both agree with Session211's sharp floors and with its exact small spectrum through 31 (KNOWN; the table only
    extends the range).

## 6. Ledger

| Claim | Type |
|---|---|
| Gate linear on (N, P), eigen-coordinates, P_Γ = P_E | PROVED |
| Closed-form towers, metallic limit, explicit one-vertex spectrum | PROVED |
| z = 1 doubling invariance of s, golden gate map, Fibonacci orbit | PROVED |
| s = THM-1950's inverse response; s ∈ (0, n], = n iff regular; velocity addition at z = 1; strong minima n ≤ 8 | KNOWN (THM-1950, klein S400, THM-4051); re-proved |
| Relativistic addition over the condensation | PROVED here; KNOWN in substance (THM-3322 at every z; THM-1950 at z = 1) |
| Commutation at z = 1 without the BS relation; P_(DΓT) = P_(DET) | PROVED |
| D(Γ T) gate-free | KNOWN (Session211 `recursion-bridge` §5; re-proved here) |
| Tower odd-cycle formula | PROVED (Session210 Theorem 4 + Session159 (1)) |
| sw_1 = 2c3 (and the whole size-weighted form) at every order | KNOWN (Session159 cycle_profiles (1)); re-proved |
| Deep towers forget the core except through c3: ρ ≤ 1/(r+1), first-order law with O(ρ/r) error, universality in N | PROVED (THM-4900(7)) |
| Operator form H(Γ^r T) = [P_T(∂) M^r e^y](0) | PROVED; KNOWN in content (Session159 lagrange (7)) |
| HYP-9490: the forgetting rate is Θ((log N)²/N²) | OPEN; FINITE-EXACT evidence, normalised value still drifting |
| Gated floors 2ⁿ + 1 and 3 + 3·2^(n−1); HC ∈ {7, 21} forces gate-freeness | KNOWN (Session159, Session210 (5), Session211, Session212) |
| Least orders 7, 8 for HC = 7, 21, with class counts 24, 228 | KNOWN (THM-4889; Session212 cycle-levels histograms); reproduced here |
| Permutohedral face formula for A_Γ(T) | PROVED + FINITE-EXACT; KNOWN in content (Session210 operator) |
| x+2 / x·2 mod 7 = F21 = Aut(QR7), the cyclic Fano symmetry | KNOWN (classical; THM-4556 remark, 2026-10-06 F21 notes, THM-4644(3)); ANALOGY as a reading |
| Links of F21 to the gaps 7, 21; φ³ = THM-224's eigenvalue; 7-vertex tower s = 21/5 | NUMEROLOGY / COINCIDENCE |

## 7. Next obligations

- **HYP-9490.**
  - Improve ρ ≤ 1/(r+1) to o(1/N): track where U_(d−2,r)'s mass sits in the induction.
  - Fix the rate: the data suggest (log N)²/N², but the normalised value is still drifting (0.90 → 0.86 over
    N = 243 → 1201).
  - Test the second order (sw_2) at N = 4001 with a C value recursion.
- **A face formula for D.** The skew-Sylvester double of a word over V has twin letters. Starting from Session211's
  affine reflected section, is there a face-sum description of the full A_(D T)(x)?
- **s beyond scalars.** s is a doubling invariant that the gate moves golden-ly. Is it a quotient of a finer
  doubling-invariant marked state under which the gate is still linear? A candidate is the full matrix-valued
  boundary resolvent of Session35's next test.
- **The 24 and 228 least-order classes with HC = 7, 21.** Session212 constructs gate-free witnesses in feedback
  number two. Do all least-order classes have feedback vertex number two, or does a different local structure appear?

## 8. Files

`04-computation/experiments/microcosm_gate_resolvent_20261106/`:
- **THM-4900 (1)–(6):** `gate_basics.py`, `gate_resolvent_exact.py`, `tower_limits.py`, `tower_spectrum.py`,
  `tower_ocf.py`.
- **THM-4900 (7):** `tower_weights_diffusion.py` (operator form, monotonicity, universality sandwich) and
  `tower_forgetting_bound.py` (ρ ≤ 1/(r+1), the explicit first-order error, the audit's Δc3 = 0 example).
- **HYP-9490:** `cover_tower_values.py`, `c3_first_order_law.py`.
- **Section 4:** `gate_descent_face_formula.py`.
- **Data:** `tower_k1_paths.c`, `gated_spectrum.c`, `phi_strong_minimum.py`, `hc_gate_free.c`.
- **Audit H1:** `audit_H1/` (`sh run_audit.sh`).
- **`run_all.sh`** reproduces every output (`.out` files committed).

# Fresh re-marking two-port quotients, strong-atom line responses, and divisor arrangements

2026-10-26, mac-mini-2026-10-25-pyth. Canon produced: [HYP-9416](../hypotheses/HYP-9416-fresh-re-marking-two-port-quotient-joint-classes-are-minimal-and-line-responses-classify-strong-atoms.md) (OPEN; FINITE-EXACT for n ≤ 7). Not yet independently audited.

**Owner prompt.**
- **User's own text.** Cuboids, spheroids and simplices in each dimension, and their vertices, points and intersections; the triangle/circle instance s(30) = 42.
- **Pasted block, read as part of the request:**
  - classify strong-atom joint-port responses;
  - find the smallest attachment data closed under repeated complement and re-marking, with the 81/27 pair as hostile control, tracking the Frobenius phase;
  - weave the Vitali-atom concepts;
  - merge the 1-2-4 / {2,5} / {7,21} / power–doubling / aliquot / abc threads, with microcosm/macrocosm patterns.

**Concurrent sessions.** The same request was evidently given to several sessions today. Their packets are the primary record, and this note is a complement:
- **[vitali-response](strong_atom_continuation_vitali_20261026/README.md):**
  - the aligned-pair carrier C1, closed under complement, re-marking and assembly;
  - unbounded alignment memory C4;
  - the operation-kernel quotient for named, guarded re-marking, C5;
  - the complete strong atlas through order 7;
  - the Vitali history ledger.
- **[aliquot124](aliquot_power_124_20261026/README.md):** the simplex/spheroid/cuboid test and the unique perfect prime-square fibre.
- **[Session69](session69_prime_ports_power_aliquot_20261026/README.md):** the closed multi-root interface, power rows and abc maps.
- **opus-2026-10-26-S60:** THM-4805 / HYP-9415, the arithmetic merge.

Scripts are in `04-computation/experiments/joint_port_bisimulation_20261026/` and `04-computation/experiments/cuboid_simplex_aliquot_20261026/`.

## 1. The fresh re-marking quotient (HYP-9416)

**The open question.** vitali-response C5 answers the minimal-continuation question for **named** guarded re-marking: its guard signatures separate every presentation through order 7. It notes that "if the observation keeps only distributions under fresh uniform re-marking, a smaller stochastic quotient may exist". This section computes that quotient exhaustively for n ≤ 7.

**Setup.**
- **States:** two-port classes `[C; s, t]`.
- **Dynamics:** V1 (fixed ports) takes the successor multiset over Hamilton paths `s → t`, followed by complement. V2 resets the ports to the ends of each path.
- **Quotient:** the coarsest partition refining an observation that is stable under successor multisets.
- **Validation:** the gluing code reproduces the published 289/253 hostile (`validate_289_253.out`).

**Findings (table in HYP-9416).**
- **Two-sided attachments admit no compression on strong states.** Nine two-sided gluings `k*h*k'`, with contexts of order at most 2, already separate all 92, 946 and 14162 strong joint classes at n = 5, 6, 7.
- **Closure turns one-sided observers into joint ones.**
  - A one-sided gluing `h*k` or `k*h` depends only on the exit-rooted or entrance-rooted class of h. So no one-sided observation separates the 7, 20, 70 one-sided collision pairs.
  - Closing under V1, i.e. repeated fresh re-marking and complement, separates every one of them.
  - This is a dynamic version of 289/253: the joint orbit is invisible to every one-sided observer, but it propagates through the dynamics into one-sided observations.
- **The bare class observable compresses.**
  - Under V1, tracking only the unrooted class at every step needs 42, 607, 11033 cells for the 92, 946, 14162 strong states.
  - So the "smaller stochastic quotient" exists for this observable, and only for observables that ignore attachment.
  - Under V2 the class partition is already closed.

**Answer to "the smallest attachment data closed under repeated complement and re-marking".**
- **Strong atoms, n ≤ 7:**
  - with named operations, it is the full presentation (vitali-response C5);
  - with fresh re-marking and any two-sided attachment observer, it is the joint two-port class, with no further compression;
  - with only one-sided attachment observers, closure forces the joint class anyway.
- **Beyond one run:** vitali-response C4 shows that no constant number of extra bits suffices for multi-run alignments.

## 2. Strong-atom line responses (HYP-9416 (3))

`F_A(h)` is the number of presentations of the strong atom A whose rooted complement is h.

| n | atoms | distinct targets | strong targets | atoms with blue presentations |
|---|---|---|---|---|
| 3 | 1 | 1 | 0 | 1 |
| 4 | 1 | 5 | 2 | 1 |
| 5 | 6 | 44 | 32 | 6 |
| 6 | 35 | 547 | 484 | 9 |
| 7 | 353 | 9985 | 9516 | 79 |

**Classification.** In every order the histograms `F_A` are pairwise distinct. They stay distinct after forgetting the ports, that is, keeping only the multiset of target free classes. So the line response determines the strong atom through order 7. The joint data of vitali-response's `strong-atlas.json` is the same object, recomputed here independently. Its totals agree with their table: 1, 5, 50, 903, 30773 strong presentations, and largest joint target fibres 1, 1, 2, 3, 8 (`crosscheck_atlas.out`).

## 3. The 81/27 control and the Frobenius phase

The pair is `T₊ = C3[A₊,A₊,A₊]` and `T₋ = C3[A₊,A₊,A₋]`, with `A₊ = C3→1` (`control_81_27.out`).

**What agrees:**
- 67 labelled cyclic triples in both;
- `H = 462915` in both.

**What differs:**
- `|Aut|` is 81 and 27, so `tc` is 5715 and 17145.
- In both tournaments every presentation has its own partner class; the partner map is injective, confirming HYP-9414 (1) at m = 3.
- Only 74 partner classes are shared.

**The Frobenius phase.** For the period-3 decoration, the entrance block relative to the odd child splits T₋'s 462915 paths equally, 154305 per phase. Each phase family contains 5715 partner classes. Injectivity makes the three families disjoint.

**Why this pair is the right control.** Any child datum that keeps only cyclic triples, or child types up to converse, sees the two decorations as the same. It therefore cannot predict a single partner. The datum must keep:
- the child type;
- its stabilizer in the quotient, i.e. the Frobenius period;
- the port's phase.

This agrees with vitali-response F1 and with Session68's gcd/lcm phase count.

## 4. The Vitali weave, in one statement

vitali-response's HISTORY.md already threads the Vitali atom from March to October:
- **THM-169 / THM-170:** the ±2 lambda switch, and the witness masks 82468 / 84707 with H = 141 vs 143;
- **the recovery principle:** "a target descends through an observation exactly when it is constant on each observation fibre";
- **the endpoint, half-edge, joint-root, reflection, stabilizer and phase hostiles:** 15/11, 11/13, 289/253, 9/11, 81/27 and gcd(d,e).

The one formal object behind all of these is the **operation-stable quotient** of an observation: the coarsest refinement through which every continuation descends. That is C5 for named operations and HYP-9416 for fresh ones. Each hostile is a pair inside one fibre of a too-coarse observation.

The arithmetic threads below are the same pattern, one level down:
- the divisor shape forgets the prime weights;
- the aliquot value forgets the number;
- the power-row address forgets which abc map was used (1+8=9 arises from both 72 and 9, per aliquot124).

## 5. Cuboids, spheroids, simplices: vertices, points, intersections

aliquot124's CONNECTION_BOARD already makes the basic reading precise:
- the divisor box is a cuboid;
- in the squarefree case it is a Boolean cube with its top corner removed;
- what remains is a simplex boundary, a sphere;
- for 30 this is a circle, weighted to 42.

The new piece here is the **arrangement count dictionary** (`divisor_shapes.out`). It realises "vertices, points, intersections" in every dimension, for a squarefree n with k primes:

| geometric object (dimension k−1) | count | arithmetic meaning | with weight d | with weight μ(d) |
|---|---|---|---|---|
| regions of the k facet hyperplanes of a (k−1)-simplex in R^(k−1) | 2^k − 1 | proper divisors (the missing all-negative region is n; the bounded region is 1) | s(n) | −μ(n) = (−1)^(k+1) |
| intersection points (simplex vertices) | k | maximal proper divisors n/p | – | – |
| regions of k great spheres on S^(k−1) (orthants) | 2^k | all divisors (cube vertices, orthoplex facets) | σ(n) | 0 |
| vertices of that arrangement (±e_i) | 2k | the 2k facets {p ∤ d}, {p \| d} of the divisor cube | – | – |
| regions in RP^(k−1) | 2^(k−1) | complementary pairs {d, n/d} | σ(n) (as d + n/d) | – |

The region counts are realised exactly by an explicit construction for k ≤ 8. In general the barycentric sign vector `(−,…,−)` is impossible because the coordinates sum to 1, and every other sign vector is realised.

**What the geometry loses.**
- `s(pqr) = 1 + e₁ + e₂`. A single linear form in the elementary symmetric functions does not fix the primes: `s(170) = s(182) = 154`.
- The value alone does not fix the shape either: `s⁻¹(42) = {30, 41²}`, a triangle and a two-step chain.

## 6. The arithmetic threads, against prior work

- **1-2-4, {2, 5}, {7, 21}: nothing new beyond prior verdicts.**
  - THM-1370: H ≠ 7, 21 at every order.
  - HYP-9362 (d): PG(2,2) and PG(2,4) are labelled numerology. The next value, 73 = |PG(2,8)|, is attained (THM-4744; HYP-2516's falsification of the base-2^k pattern).
  - aliquot124 and Session69 reject a map between untouchables and the H-gaps: 7 = s(8) and 21 = s(18) are touchable, while H = 5 occurs.
  - {2, 5} as the only untouchable primes is conditional on the distinct-prime Goldbach statement, which is CITED as verified to 4·10¹⁸.
  - 7 = 2 + 5 is noted as numerology only.
- **Your table** is s⁻¹ (with the correction s(9) = 4). Complete fibres are in Session66 `small_fibers.json` and `aliquot_preimages.out`.
- **Twin primes (confirming Session69).** For v > 3, v has both a prime-square and a 2q preimage iff (v−3, v−1) are twin primes. The only self-hit is v = 6 (aliquot124 Theorem A).
- **A second microcosm/macrocosm (new here).**
  - 31 = σ(2⁴) = σ(5²) is the only value 1 + p + … + p^(a−1), with a ≥ 3, shared by two prime bases below 10¹² (p < 2000). This is Goormaghtigh's 31, so s⁻¹(31) contains 32 and 125.
  - It is the macrocosm of 7 = 1+2+4 = s(8), which has a unique preimage.
- **The two sixes and abc (new remark).**
  - The power chain of 3 and the doubling chain of 1 come within 1 of each other only at the Gersonides pairs (2,3), (3,4), (8,9).
  - The resulting triples 1+2=3, 1+3=4 and 1+8=9 all have radical 6. The last is the first abc hit, with quality 1.2263.
  - The remaining contact 1+1=2 (radical 2) is the perfect-number normalisation σ(n)/n = 2.
  - So the "two sixes" are rad = 6 of the 2–3 contacts and the first perfect number. Session69's binary bijection (6 ↦ 1+2=3, 72 ↦ 1+8=9) is the doubling-row form of the first and last of these.
- **Power/doubling.**
  - Each odd non-power b₀ gives a grid {2ⁱ b₀ʲ}.
  - 6 is the first n whose square lies in a later doubling row: 36 = 9·4.
  - The general condition is n = 2ᵃu with u^(k−1) > 2ᵃ.
  - All of this is PROVED in aliquot124 chains §2 and Session69-power; `power_doubling.out` reproduces it.

## 7. Open problems

1. All-order proofs of HYP-9416 (1)–(3), and a structural reason why contexts of order 2 suffice.
2. The fresh-re-marking quotient for multi-run attachments (k > 1), and whether some compression survives there that is invisible to C4's guards.
3. Whether the class-observable compression ratio tends to 1.

## Reproduction

- `bisim3.py n` for n = 3..7: about 70 s at n = 7.
- `bisim_detail.py n`, `bisim_detail2.py n`.
- `atom_responses.py n`.
- `control_81_27.py`: about 50 s.
- In `cuboid_simplex_aliquot_20261026/`: `divisor_shapes.py`, `aliquot_preimages.py 300`, `power_doubling.py`, `perfect_repunit_abc.py`.

Requirements: nauty (`gentourng`, `labelg`), sympy and numpy.

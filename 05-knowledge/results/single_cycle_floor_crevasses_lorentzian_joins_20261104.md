# The single-cycle floor and the crevasses: the bottom band of the strong spectrum lies in the uniquely Hamiltonian stratum, 7 and 21 are holes of that stratum, and one-arc joins (session29's Lorentzian law) cover [23, 84705] with a 3-adic bias

**Session:** opus-2026-11-04-S73. Worktree `C:/tmp/wt-s73`, cut from main after S72's checkpoint.

**Owner prompt.**
- Start a session on the single-cycle floor and the crevasses.
- Don't aim to fill the spectrum naively, one ladder at a time; be more creative with groupings and constructions.
- Use the attached papers as inspiration, combined with the repository's prior work:
  - Bannwart, *The real Betti realization of very effective Hermitian K-theory*, arXiv:2610.20688;
  - Grogan, *Prime pairs at the prime number theorem threshold*, 2026;
  - Ding-Kim, *Grassmann-Plücker functions for orthogonal matroids*, arXiv:2601.19830.

**Status.**
- **KNOWN.** The single-cycle floor itself:
  - Busch-extremal strong tournaments are Moon's interval chains (Moon-Yang 2007; THM-4744 (4), THM-4745);
  - interval chains have the upset score list (1, 1, 2, ..., n-2, n-2) (THM-4744 (2)), and upset tournaments each have a
    unique Hamiltonian cycle (Brualdi-Li, JCT B 35 (1983); Poet-Shader, EJC 5 (1998) #R24);
  - the same sentence appears in Savchenko (arXiv:2311.13377, Introduction), right after he cites Burzio-Demaria
    (1990) for the score-list characterization;
  - the repository already had the case of at most two blocks (THM-4655 (1e)).

  A short independent cut proof is given here.
- **KNOWN in the repository.** The one-arc join law H(J) = x_h y_h + x_R y_R - x_e y_e, with its Lorentzian form:
  - session29 `ADDENDUM.md`, sections 2 and 5, and `LORENTZ_MOMENT_COROLLARIES.md`;
  - the C1 case of `polarized_path_calculus_20261103`.

  Added here: the statement in klein's root-response notation for arbitrary rooted pairs, and the cycle formula
  hc(J) = s_a(A) e_b(B), an elementary generalization of THM-4862 (b).
- **PROVED (elementary).** The real-point decomposition H = H(core) x (product)^2 for self-converse tournaments, with its
  parity crevasses. The component pairing is THM-4743 (1b), and H = 3 at even order also follows from THM-791.
- **FINITE-EXACT.**
  - The cycle ceiling and the two-cycle floor through order 10.
  - The single-cycle band [f(n), f_2(n)) is the bottom part of the uniquely Hamiltonian stratum; that stratum extends
    higher, up to 467 at order 10.
  - The crevasse localization, through order 10.
  - Coverage of [23, 84705] by the one-arc-join family, with its residue structure.
  - The Pfaffian null result.
- **PRE-REGISTERED** (HYP-9457):
  - P2, P4 and P5 were genuine tests and are confirmed;
  - P6 was a genuine test and is refuted (the 3-adic bias is real but weaker than predicted);
  - P1 and P3 were sanity checks. P1 is implied by THM-4744 (4) and the cut proof; P3's one-arc-cut half is classical
    (Shelor 1977).
- **CLOSED, KNOWN.** HYP-9457 (b): every uniquely Hamiltonian tournament has a one-arc cut. Shelor's 1977 complete
  insertion construction gives a vertex of in-degree 1, which is itself a one-arc cut (session103;
  `operation_closed_synthesis_20261029` section II; klein THM-4862 (b), (f)).

Independently audited, PASS WITH CORRECTIONS, all applied (MISTAKE-705): section 10.

## 0. What was known

**THM-1370 (the impossibility half).** 7 and 21 never occur. The proof uses:
- multiplicativity over strong components;
- the insertion lemma;
- Moon's strong-deletion theorem;
- exhaustive strong spectra through order 6, with f(7) = 25.

It calls the gaps crevasses of the strong floor.

**THM-4744 / THM-4745 (the floor).**
- Moon's interval chains IC(p) have H = prod (2^(p_i - 1) + 1).
- They are the extremal strong tournaments (Moon-Yang), with E(m) = 1, 1, 2, 2, 1, 5, 3, 1, 9, 4, 1 floor classes at m = 3..13.
- Every strong tournament with H < 625 has order at most 12.

**THM-4846 (the uniquely Hamiltonian class).** It cites:
- Douglas 1970, for the characterization of tournaments with exactly one Hamiltonian cycle;
- Garey 1972, for their count F(2n - 6).

**Shelor 1977** (cited in the repository since session103). Every uniquely Hamiltonian tournament of order n >= 5 arises
from one of order n - 1 by inserting a vertex that beats every old vertex except one parent.

**Session29 (2026-10-22).** Its `ADDENDUM.md` sections 2 and 5 already prove the one-arc join law of section 4, with the
same case split, and its Lorentzian form.

**klein (same day).**
- THM-4862: one-exit law; readiness digit = Hamiltonian cycles of the entry attachment.
- THM-4865: the sole-exit tower is tribonacci, and its root response B_t(T) = sum over A of h(T[A]) s_t(T[V - A]).
- HYP-9453 and HYP-9454: the near-join and C2 sweeps.
- Their sweeps from bases of order <= 10 reach every odd number in [69, 8,079,323] (THM-4865 (e)). The repository's banks realise every odd number up to 99,999,999 except 7 and 21 (HYP-9453).

**THM-4094.** Completeness reduces to the prime lanes.

**THM-4863 (S72, Theorem H).** At odd order H = G + 2hc, where G is determined by the deck. So the strong layer is the hc >= 1 layer, by Camion.

**This note.** It groups tournaments in four ways:
- by their number of Hamiltonian cycles, where the bottom of every strong layer (through order 10) lies in the uniquely
  Hamiltonian stratum;
- by one-arc joins, using session29's Lorentzian law with the cycle formula added, in a Goldbach-style representation
  count (an ANALOGY);
- by converse symmetry, the "real points";
- by Pfaffians.

Program [`scf_crevasses_20261104.c`](../../04-computation/experiments/scf_crevasses_20261104.c). Outputs:
- [census](../../04-computation/experiments/scf_crevasses_20261104_census.out), orders 3-10, every class;
- [chains](../../04-computation/experiments/scf_crevasses_20261104_chains.out);
- [joins](../../04-computation/experiments/scf_crevasses_20261104_joins.out).

## 1. The single-cycle floor (CITED; independent proof)

**Claim.** Every interval chain IC(p) has exactly one Hamiltonian cycle.

**Proof.**
1. The vertices are x_0 < ... < x_N. The only backward arcs are the k block arcs x_(l_(i+1)) -> x_(l_i), and the blocks
   [l_i, l_(i+1)) partition the cuts j | j + 1. So each cut is spanned by exactly one backward arc.
2. A Hamiltonian cycle crosses every cut at least once in each direction. It therefore crosses each cut backward exactly
   once and forward exactly once, and it uses all k back arcs.
3. Its forward arcs then tile [0, N). Their tails are forced: every vertex except the back-arc tails x_(l_2), ..., x_N.
   A tiling is determined by its left endpoints, so the arc set is determined.
4. Camion gives existence.

QED.

**Single-cycle floor.** With Moon-Yang, every strong tournament with H = f(n) is uniquely Hamiltonian. This is classical:
interval chains are upset tournaments (Poet-Shader 1998; Brualdi-Li 1983). The proof above works verbatim for all upset
chains, with parts >= 1; the audit checked all 2,046 labelled compositions through order 12.

**Checks.**
- Every floor class at orders 3-10 has exactly one Hamiltonian cycle. The floor-class counts are 1, 1, 2, 2, 1, 5, 3,
  1, matching THM-4744.
- All 1,596 interval chains of order <= 17 have exactly one Hamiltonian cycle and the product value. Their number at
  order n is F(n - 2), so the total is F(17) - 1. They form a Fibonacci subfamily of Garey's F(2n - 6).

**Deck form at odd order.** By THM-4863 Theorem H, G = H - 2hc is determined by the deck. The statement G_min = f(n) - 2
(= 1, 7, 23, 73 at n = 3, 5, 7, 9) is equivalent to the cycle ceiling of section 2.

## 2. The bottom band and the cycle ceiling (FINITE-EXACT through order 10; HYP-9457 P1-P4 pre-registered)

**Definitions.**
- f_k(n) = min{H : hc = k}. We call f_2 the "two-cycle floor". In S50's note, "second floor" means the second-lowest
  strong count, a different quantity.
- The single-cycle band at order n is [f(n), min over k >= 2 of f_k(n)). By definition it is the bottom part of the
  uniquely Hamiltonian stratum; the stratum extends higher.
- U(n) is the uniquely Hamiltonian stratum (hc = 1).

| n | f(n) | band end f_2(n) | U(n): classes / distinct H / range | band values missing at this order |
|---:|---:|---:|---|---|
| 5 | 9 | 13 | 3 / 2 / [9, 11] | none |
| 6 | 15 | 23 | 8 / 4 / [15, 23] | **21** |
| 7 | 25 | 41 | 21 / 11 / [25, 51] | none |
| 8 | 45 | 69 | 55 / 18 / [45, 107] | 47, 55, 67 |
| 9 | 75 | 115 | 144 / 48 / [75, 227] | 77, 79, 83 |
| 10 | 125 | 205 | 377 / 91 / [125, 467] | 127, 129, ... (13 values) |

**Cycle ceiling.** H(T) >= f(n) + 2(hc(T) - 1) for every strong class at orders 3-10. Equality holds only at the floor
classes. At odd order this is the deck form G >= f(n) - 2.

**Two-cycle floor.** The k-cycle floors lie far above the ceiling. The two-cycle floor f_2(n) = 13, 23, 41, 69, 115,
205 (n = 5..10) satisfies f_2(n) < f(n + 1) in every case.

**Pre-registered results.**
- P1 (sanity check): f(10) = 125 is attained by one class, which is uniquely Hamiltonian. This was implied by THM-4744 (4)
  and the cut proof.
- P2: no violation of the ceiling at order 10; equality only at the floor. CONFIRMED.
- P3 (sanity check): 377 = F(14) uniquely Hamiltonian classes, all with a one-arc cut. Classical: Garey's count, and
  Shelor's in-degree-1 vertex.
- P4: f_2(10) = 205, inside [175, 215] and below 225. CONFIRMED.

## 3. The crevasses are holes of the uniquely Hamiltonian stratum (FINITE-EXACT + THM-1370)

The bottom band at every order 5-10 lies in U(n), and through order 10 every hole of the strong spectrum in
[f(n), f(n + 1)) lies in the band. So, through order 10, the crevasse region is a question about U(n):

    U(4) = {5},    U(5) = {9, 11},    U(6) = {15, 17, 19, 23},    U(7) = {25, ..., 39, 45, 47, 51}, ...

- **7** lies between the order-4 stratum {5} and the order-5 floor 9. It is prime, and no smaller value provides it as a
  product.
- **21** is the unique hole of U(6) = {15, 17, 19, _, 23}, below f(7) = 25, and 21 = 3 x 7 is not a product of
  attained values.
- The band holes at orders 8-10 (47, 55, 67; 77, 79, 83; 127, ...) are all filled at other orders or as products
  (55 = 5 x 11). Only 7 and 21 survive.

**Typing.** This is a localization, not a new impossibility proof; the proof is THM-1370's. The finite content, through
order 10, is that the region where crevasses can arise lies in the bottom band of each order, and that this band lies in
the uniquely Hamiltonian stratum. `exclusion_7_21_20261103` re-proves the 7/21 exclusion, with a core-counting
refinement through H <= 21.

## 4. One-arc joins are Lorentzian inner products (KNOWN in the repository: session29; the cycle formula is added here)

**Setup.** For tournaments A, B, vertices a in A and b in B, let J = J(A, a; B, b) be A => B with the arc a -> b
reversed. Write p(X) for the number of Hamiltonian paths of an induced subtournament, with p(empty) = 1. Define:

    x(A, a) = (h(A), R_a(A), e_a(A)),    R_a(A) = sum over X subset V(A) - a of p(X) s_a(A - X)

(klein's root response B_a, THM-4865), and y(B, b) = x(B^op, b) = (h(B), R'_b(B), s_b(B)).

**Theorem (one-arc join law).**

    H(J)  = x_h y_h + x_R y_R - x_e y_e,
    hc(J) = s_a(A) e_b(B)            (|A| + |B| >= 3).

The form is a Lorentzian inner product of signature (2, 1).

**Proof.**
1. Only one arc, b -> a, leads from B to A, so a Hamiltonian path crosses back at most once.
2. Paths that never cross back are an A-path followed by a B-path, minus those that used a -> b:
   h_A h_B - e_a s_b.
3. Paths that cross back once have the form A_1 B_1 (b -> a) A_2 B_2, where:
   - A_1 is a path on a subset of A - a;
   - B_1 is a path ending at b;
   - A_2 = A - A_1 is a path starting at a;
   - B_2 = B - B_1.

   Every junction is an available A -> B arc. The count factorizes as R_a(A) R'_b(B).
4. A Hamiltonian cycle must use b -> a exactly once. What remains is an A-path from a followed by a B-path to b, giving
   s_a e_b. (The junction can only fail when |A| = |B| = 1.)

QED.

**Prior art.** The H formula, its proof and its Lorentzian form are session29's (`ADDENDUM.md` sections 2 and 5;
`LORENTZ_MOMENT_COROLLARIES.md`), and polarized_path_calculus C1 is a special case. The cycle formula hc = s_a e_b
generalizes THM-4862 (b).

**Checks.** The program checks H on three runs of 3,000 random joins that share one seed. The audit checked H and hc
exhaustively on all 138,905 rooted class pairs with |A| + |B| <= 9.

**Special cases.**
- B = one vertex gives klein's T_1 = h - e + B.
- B = a transitive chain gives the near-join h + 2^k B (HYP-9453).
- TT_k and TT_m, rooted at source and sink, give Moon's NT_(k+m) with 1 + 2^(k+m-2).
- Interval chains are iterated one-arc joins: cut at the first shared vertex.

**Consequences for uniquely Hamiltonian tournaments.**
- A join with s_a(A) = e_b(B) = 1 is uniquely Hamiltonian.
- Conversely, every uniquely Hamiltonian tournament has a one-arc cut, and the cut can always be taken to be a single
  vertex. Shelor's 1977 complete insertion construction gives a vertex of in-degree 1, which is classical (HYP-9457 (b)
  is CLOSED, KNOWN). The pieces with s_a = 1 are classified in `polarized_path_calculus_20261103`, with class series
  x^3/(1 - 3x + x^2).

## 5. Goldbach-style covering (Grogan-inspired ANALOGY; FINITE-EXACT; HYP-9457 P5 confirmed, P6 refuted)

The join values are Q(x, y) for x, y in the set X of realizable rooted triples. The y-set {x(B^op, b)} equals X because
B -> B^op is a bijection on tournaments. Covering by this family becomes a representation problem, by analogy with
N = p + q.

| both sides of order <= | distinct triples | odd N with no representation |
|---:|---:|---|
| 5 | 42 | 7, 21, 47, 67, 83, ... (reach about 50) |
| 6 | 254 | only 7, 21 below 2,395 |
| 7 | 2,732 | **only 7, 21 below 84,707** |

**The comet.** At orders <= 7, the minimum representation count r(N) over odd N in each dyadic window is:

| window | [128, 256) | [512, 1024) | [2048, 4096) | [4096, 8192) | [8192, 16384) | [16384, 32768) |
|---|---:|---:|---:|---:|---:|---:|
| min r(N) | 34 | 34 | 94 | 160 | 198 | 152 |

After that the range edge sets in: with a finite piece set, r(N) must vanish beyond the range. Within the range the
counts grow and the exceptional set is {7, 21}. The resemblance to Hardy-Littlewood is an ANALOGY.

**The singular series.** r(N) depends on N mod 3, not on N mod 5 or mod 7. The mean ratio of residue 0 to the other
residues is:

| sides of order <= | window | ratio |
|---:|---|---:|
| 5 | [375, 750) | 2.0 |
| 6 | [5000, 10000) | 2.1 |
| 7 | [25000, 50000) | 1.37 |

The mod 5 deviations are <= 0.04 and the mod 7 deviations <= 0.01.

P6 had predicted a ratio in [1.6, 3.0] at orders <= 7. It is REFUTED: the 3-adic factor is real but weakens as the
pieces grow.

**Mechanism (the audit's analysis, HEURISTIC).**
- Triples equidistributed mod p give a generic local factor 1 + p/(p^3 - 1): 1.115, 1.040 and 1.020 at p = 3, 5, 7.
  The residue-0 excesses observed mod 5 and mod 7 (about 1.05 and 1.01) are of this generic size.
- The extra 3-adic bias comes from triples ≡ (0, 0, 0) mod 3: 238 of 2,732 (8.7 %), against an expected 3.7 %.
- It does not come from H being divisible by 3: only 32.6 % of the triples have h ≡ 0 mod 3.

This is the opposite of a crevasse. Multiples of 3 are over-represented, and both crevasses 7 and 21 sit far below the
scale where r(N) is large.

**Relation to known coverage.** The section is a statement about this family, not about the spectrum. The spectrum is
already known to cover every odd number up to 99,999,999 except 7 and 21 (HYP-9453). klein's sweeps from bases of
order <= 10 reach [69, 8,079,323]. The family here reaches [23, 84705] from pieces of order <= 7. Its counts are those
of a product set: about |X|^2 representations spread over about max(h)^2 targets.

## 6. Real points (Bannwart-inspired; PROVED + FINITE-EXACT)

**The decomposition.** Real Betti realization takes fixed points of complex conjugation. Here the involution is
tournament converse, and the real points are the self-converse tournaments.
- An anti-automorphism reverses the condensation. Components pair as C_i with C_(k+1-i) ≅ C_i^op, and a middle
  component, if present, is self-converse.
- By multiplicativity and H(C^op) = H(C):

      H(T) = H(M) x (prod over pairs of H(C_i))^2,    |T| ≡ |M| (mod 2).

So the real spectrum is {m s^2 : m a strong self-converse value or 1, s in the spectrum}: a self-dual core times
hyperbolic squares.

The component pairing is THM-4743 (1b). The resemblance to L-theory, where Bannwart's real realization of ko lands, is an
ANALOGY.

**Parity crevasses (PROVED, for self-converse T only).**
- H = 5 never occurs at odd order: s = 1 would force a strong odd-order core with H = 5, but strong odd orders give 1,
  3 or >= 9. For tournaments in general this is false: NT_4 plus a dominated vertex has order 5 and H = 5.
- H = 3 never occurs at even order: it is neither a square nor an even-order strong value. This also follows from
  THM-791 (H ≡ 1 mod 4 under a fixed-point-free involutive anti-automorphism; THM-4818 (A)).

**Observed real gaps.** The real spectrum has far more gaps than the full one. Below 200 there are 38 odd gaps at order
9 and 67 at order 10 (census). They sort by parity of the core.

## 7. Pfaffians (Ding-Kim-inspired; FINITE-EXACT null result, plus an analogy)

A tournament is a ±1 skew-symmetric matrix S, and |Pf(S)| is an isomorphism invariant at even order. At orders 4-10,
|Pf| mod 4 agrees with H mod 4 in 4/4, 30/56, 3,608/6,880 and 4,868,460/9,733,056 classes. That is chance level, and the
same holds against hc parity. There is no Pfaffian law for H at this resolution.

The structural analogy that survives is Cayley's identities, which express principal and almost-principal minors of a
skew-symmetric matrix through its Pfaffians; Ding-Kim build on them. The parallel is Theorem H (THM-4863): at odd order
the deck (principal data) determines H up to the one spanning term 2hc (at even order it determines H). That term is invisible to the deck, and the one-arc join law of section 4 makes it multiplicative,
hc = s_a e_b. This is an INTERPRETATION, not a theorem.

## 8. Typing

| item | type |
|---|---|
| single-cycle floor (Busch-extremal ⇒ uniquely Hamiltonian) | KNOWN (Moon-Yang; upset tournaments: Brualdi-Li 1983, Poet-Shader 1998; THM-4655 (1e)); short proof here |
| interval chains: one Hamiltonian cycle, all n | KNOWN (upset tournaments); independent cut proof; checked on 1,596 chains through order 17 |
| cycle ceiling; two-cycle floor f_2 = 13, 23, 41, 69, 115, 205; f_2(n) < f(n+1) | FINITE-EXACT through order 10 (P2, P4 pre-registered) |
| bottom band lies in the uniquely Hamiltonian stratum; every uniquely Hamiltonian tournament has a one-arc cut (a single vertex) | band: FINITE-EXACT through order 10; cut: KNOWN (Shelor 1977), HYP-9457 (b) CLOSED |
| crevasses as holes of U(n) (7 between U(4), U(5); 21 the hole of U(6)) | FINITE-EXACT localization through order 10; impossibility is THM-1370 |
| one-arc join law H = Lorentzian form | KNOWN (session29 ADDENDUM sections 2 and 5) |
| cycle formula hc = s_a e_b | PROVED (generalizes THM-4862 (b)); audit-checked exhaustively on 138,905 pairs |
| family coverage [23, 84705] at orders <= 7 | FINITE-EXACT (P5 confirmed); not new spectrum information |
| 3-adic bias | FINITE-EXACT; P6 refuted as stated (1.37); mechanism HEURISTIC; Hardy-Littlewood ANALOGY |
| real-point decomposition and parity crevasses (self-converse T) | PROVED (elementary; pairing THM-4743 (1b); H = 3 also THM-791); L-theory ANALOGY |
| Pfaffian | FINITE-EXACT null |
| Cayley analogy | INTERPRETATION |

## 9. Pre-registration record

- HYP-9457 reserved in `9db3423759`; predictions P1-P6 committed in `4d5099bc95`, before `census 10` and `joins 7 7`.
- Census 3-9, chains 17 and joins 5 5 / 6 6 were exploratory.
- The real-point and Pfaffian items were exploratory.

## 10. Audit

Independent adversarial audit: [`AUDIT.md`](../../04-computation/experiments/scf_crevasses_20261104_audit_S73/AUDIT.md).
Verdict: **PASS WITH CORRECTIONS.**

**Reproduction.** The auditor wrote their own code and reproduced:
- the full census through order 10;
- the join law and the cycle formula, exhaustively on 138,905 rooted class pairs;
- the coverage, the window minima and the mod-3 ratios;
- the real-point decomposition on all 8,784 self-converse classes of order 10.

**Majors, all applied.**
- **M1.** The one-arc join law and its Lorentzian form were session29's; now KNOWN.
- **M2.** HYP-9457 (b) is classical: Shelor 1977, as recorded in session103 and operation_closed_synthesis. It is now
  closed, and the cut is a single vertex.
- **M3.** "The bottom band equals the uniquely Hamiltonian stratum" was false as worded; it is corrected throughout.

**Minors, all applied.**
- the chain result is classical (upset tournaments);
- the chain count is 1,596;
- P1 and P3 are now sanity checks;
- the scope and credit of (E);
- the mod-3 mechanism replaced;
- the coverage context;
- "through order 10";
- the "two-cycle floor" rename;
- the Cayley direction and the odd-order qualifier;
- the converse-bijection wording;
- the credit to polarized_path_calculus and exclusion_7_21.

**Logged as MISTAKE-705.**

## 11. Reproduction

```
gcc -O2 -fopenmp -o scf 04-computation/experiments/scf_crevasses_20261104.c
./scf census 10        # 4.5 min on 12 threads (order 10: 9,733,056 classes)
./scf chains 17        # 7 s
./scf joins 7 7 200000 # under 1 s
```

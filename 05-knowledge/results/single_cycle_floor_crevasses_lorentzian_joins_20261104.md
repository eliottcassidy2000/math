# The single-cycle floor and the crevasses: the bottom band of the strong spectrum is the uniquely Hamiltonian stratum, its members are one-arc (Lorentzian) joins, 7 and 21 are holes of that stratum, and Lorentzian joins cover [23, 84705] with a 3-adic singular series

**Session:** opus-2026-11-04-S73. Worktree `C:/tmp/wt-s73`, cut from main after S72's checkpoint.

**Owner prompt.**
- Start a session on the single-cycle floor and the crevasses.
- Don't aim to fill the spectrum naively, one ladder at a time; be more creative with groupings and constructions.
- Use the attached papers as inspiration, combined with the repository's prior work:
  - Bannwart, *The real Betti realization of very effective Hermitian K-theory*, arXiv:2610.20688;
  - Grogan, *Prime pairs at the prime number theorem threshold*, 2026;
  - Ding-Kim, *Grassmann-Plücker functions for orthogonal matroids*, arXiv:2601.19830.

**Status.**
- **CITED.** The single-cycle floor itself:
  - Busch-extremal strong tournaments are Moon's interval chains (Moon-Yang 2007; THM-4744 (4), THM-4745);
  - interval chains have the minimum n - 2 cyclic triples (THM-4744 (2));
  - strong tournaments with n - 2 cyclic triples have exactly one Hamiltonian circuit (Burzio-Demaria, J. Graph Theory 14 (1990) 663-672, as stated in Savchenko, arXiv:2311.13377).

  A short independent cut proof is given here.
- **PROVED (elementary).**
  - The one-arc join law: H(J) = x_h y_h + x_R y_R - x_e y_e and hc(J) = s_a(A) e_b(B). It contains klein's tower step and near-join (THM-4865) as special cases.
  - The real-point decomposition H = H(core) x (product)^2 for self-converse tournaments, with its parity crevasses.
- **FINITE-EXACT.**
  - The cycle ceiling and the second floor through order 10.
  - The bottom band equals the uniquely Hamiltonian stratum, and all 1 + 1 + 3 + 8 + 21 + 55 + 144 + 377 of its classes are one-arc joins.
  - The crevasse localization.
  - Lorentzian join coverage of [23, 84705], with the representation comet and its residue structure.
  - The Pfaffian null result.
- **PRE-REGISTERED** (HYP-9457): P1-P5 confirmed; P6 refuted (the 3-adic bias is real but weaker than predicted).

Independently audited: section 10.

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

**klein (same day).**
- THM-4862: one-exit law; readiness digit = Hamiltonian cycles of the entry attachment.
- THM-4865: the sole-exit tower is tribonacci, and its root response B_t(T) = sum over A of h(T[A]) s_t(T[V - A]).
- HYP-9453 and HYP-9454: the near-join and C2 sweeps.
- Their three sweeps together give every odd number in [69, 702189].

**THM-4094.** Completeness reduces to the prime lanes.

**THM-4863 (S72, Theorem H).** At odd order H = G + 2hc, where G is determined by the deck. So the strong layer is the hc >= 1 layer, by Camion.

**This note.** It does not fill the spectrum ladder by ladder. It groups tournaments in four ways:
- by their number of Hamiltonian cycles, where the bottom of every strong layer is the uniquely Hamiltonian stratum;
- by a two-sided join that is a Lorentzian inner product of rooted triples, which turns covering into a Goldbach-type problem;
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

**Single-cycle floor.** With Moon-Yang, every strong tournament with H = f(n) is uniquely Hamiltonian. As typed in the
status line, this is classical by combination.

**Checks.**
- Every floor class at orders 3-10 has exactly one Hamiltonian cycle. The floor-class counts are 1, 1, 2, 2, 1, 5, 3,
  1, matching THM-4744.
- All 1,601 interval chains of order <= 17 have exactly one Hamiltonian cycle and the product value. Their number at
  order n is F(n - 2), a Fibonacci subfamily of Garey's F(2n - 6).

**Deck form at odd order.** By THM-4863 Theorem H, G = H - 2hc is determined by the deck. The strong floor is
therefore G_min + 2 with G_min = f(n) - 2 = 1, 7, 23, 73 at n = 3, 5, 7, 9; this is equivalent to the cycle ceiling
of section 2.

## 2. The bottom band and the cycle ceiling (FINITE-EXACT through order 10; HYP-9457 P1-P4 pre-registered)

**Definitions.**
- f_k(n) = min{H : hc = k}.
- The single-cycle band at order n is [f(n), min over k >= 2 of f_k(n)).
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

**Second floor.** The k-cycle floors lie far above the ceiling. The second floor f_2(n) = 13, 23, 41, 69, 115, 205
(n = 5..10) satisfies f_2(n) < f(n + 1) in every case.

**Pre-registered results.**
- P1: f(10) = 125 is attained by one class, which is uniquely Hamiltonian. CONFIRMED.
- P2: no violation of the ceiling at order 10; equality only at the floor. CONFIRMED.
- P3: 377 = F(14) uniquely Hamiltonian classes, all one-arc joins. CONFIRMED.
- P4: f_2(10) = 205, inside [175, 215] and below 225. CONFIRMED.

## 3. The crevasses are holes of the uniquely Hamiltonian stratum (FINITE-EXACT + THM-1370)

The bottom band at every order 5-10 is entirely uniquely Hamiltonian, so the crevasse region is a question about U(n):

    U(4) = {5},    U(5) = {9, 11},    U(6) = {15, 17, 19, 23},    U(7) = {25, ..., 39, 45, 47, 51}, ...

- **7** lies between the order-4 stratum {5} and the order-5 floor 9. It is prime, and no smaller value provides it as a
  product.
- **21** is the unique hole of U(6) = {15, 17, 19, _, 23}, below f(7) = 25, and 21 = 3 x 7 is not a product of
  attained values.
- The band holes at orders 8-10 (47, 55, 67; 77, 79, 83; 127, ...) are all filled at other orders or as products
  (55 = 5 x 11). Only 7 and 21 survive.

**Typing.** This is a localization, not a new impossibility proof; the proof is THM-1370's. The new content is that the
region where crevasses can arise is exactly the bottom band of each order, and that this band is the uniquely
Hamiltonian stratum.

## 4. The creative grouping: one-arc joins are Lorentzian inner products (PROVED; checked)

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

QED. The formula is checked by direct counting on 3 x 3,000 random joins.

**Special cases.**
- B = one vertex gives klein's T_1 = h - e + B.
- B = a transitive chain gives the near-join h + 2^k B (HYP-9453).
- TT_k and TT_m, rooted at source and sink, give Moon's NT_(k+m) with 1 + 2^(k+m-2).
- Interval chains are iterated one-arc joins: cut at the first shared vertex.

**Consequences for the bottom band.**
- A join with s_a(A) = e_b(B) = 1 is uniquely Hamiltonian.
- Conversely, every uniquely Hamiltonian class at orders 3-10 (all 610 of them) has a cut crossed backward by exactly
  one arc, so it is such a join. This is pre-registered P3 at order 10.

So, through order 10, the bottom band is the image of the Lorentzian form on pairs of rooted tournaments with
s_a = e_b = 1. This is a two-parameter numeration of the crevasse region, not a one-parameter ladder. Whether every
uniquely Hamiltonian tournament has such a cut is HYP-9457 (b); Douglas's characterization is the natural place to prove
it.

## 5. Goldbach-style covering (Grogan-inspired; FINITE-EXACT; HYP-9457 P5 confirmed, P6 refuted)

The join values are Q(x, y) for x, y in the set X of realizable rooted triples. X is closed under converse, so both
sides draw from the same set. In this two-sided form, completeness is a representation problem, like N = p + q.

| both sides of order <= | distinct triples | odd N with no representation |
|---:|---:|---|
| 5 | 42 | 7, 21, 47, 67, 83, ... (reach about 50) |
| 6 | 254 | only 7, 21 below 2,395 |
| 7 | 2,732 | **only 7, 21 below 84,707** |

**The comet.** At orders <= 7, the minimum representation count r(N) over odd N in each dyadic window is:

| window | [128, 256) | [512, 1024) | [2048, 4096) | [4096, 8192) | [8192, 16384) | [16384, 32768) |
|---|---:|---:|---:|---:|---:|---:|
| min r(N) | 34 | 34 | 94 | 160 | 198 | 152 |

After that the range edge sets in. The representation counts grow, as Hardy-Littlewood predicts for N = p + q, and the
exceptional set is {7, 21}.

**The singular series.** r(N) depends on N mod 3, not on N mod 5 or mod 7. The mean ratio of residue 0 to the other
residues is:

| sides of order <= | window | ratio |
|---:|---|---:|
| 5 | [375, 750) | 2.0 |
| 6 | [5000, 10000) | 2.1 |
| 7 | [25000, 50000) | 1.37 |

The mod 5 deviations are <= 0.04 and the mod 7 deviations <= 0.01.

P6 had predicted a ratio in [1.6, 3.0] at orders <= 7. It is REFUTED: the local factor at 3 is real but weakens as the
pieces grow. The likely mechanism is that small tournaments have H divisible by 3 disproportionately (C_3 factors).

This is the opposite of a crevasse. Multiples of 3 are over-represented, and both crevasses 7 and 21 sit far below the
scale where r(N) is large.

**Relation to klein.** klein's three sweeps (orders <= 9 bases) reach [69, 702189]. The single two-sided family here
reaches [23, 84705] from pieces of order <= 7. Its counts are those of a product set: about |X|^2 representations
spread over about max(h)^2 targets.

## 6. Real points (Bannwart-inspired; PROVED + FINITE-EXACT)

**The decomposition.** Real Betti realization takes fixed points of complex conjugation. Here the involution is
tournament converse, and the real points are the self-converse tournaments.
- An anti-automorphism reverses the condensation. Components pair as C_i with C_(k+1-i) ≅ C_i^op, and a middle
  component, if present, is self-converse.
- By multiplicativity and H(C^op) = H(C):

      H(T) = H(M) x (prod over pairs of H(C_i))^2,    |T| ≡ |M| (mod 2).

So the real spectrum is {m s^2 : m a strong self-converse value or 1, s in the spectrum}: a self-dual core times
hyperbolic squares. This is the shape of L-theory, which is where Bannwart's real realization of ko lands.

**Parity crevasses (PROVED).**
- H = 5 never occurs at odd order: s = 1 would force a strong odd-order core with H = 5, but strong odd orders give 1,
  3 or >= 9.
- H = 3 never occurs at even order: it is neither a square nor an even-order strong value.

**Observed real gaps.** The real spectrum has far more gaps than the full one. Below 200 there are 38 odd gaps at order
9 and 67 at order 10 (census). They sort by parity of the core.

## 7. Pfaffians (Ding-Kim-inspired; FINITE-EXACT null result, plus an analogy)

A tournament is a ±1 skew-symmetric matrix S, and |Pf(S)| is an isomorphism invariant at even order. At orders 4-10,
|Pf| mod 4 agrees with H mod 4 in 4/4, 30/56, 3,608/6,880 and 4,868,460/9,733,056 classes. That is chance level, and the
same holds against hc parity. There is no Pfaffian law for H at this resolution.

The structural analogy that survives is Cayley's identities, which Ding-Kim use to recover Pfaffians from principal and
almost-principal minors. The parallel is Theorem H (THM-4863): the deck (principal data) determines H up to the one
spanning term 2hc. That term is invisible to the deck, and the one-arc join law of section 4 makes it multiplicative,
hc = s_a e_b. This is an INTERPRETATION, not a theorem.

## 8. Typing

| item | type |
|---|---|
| single-cycle floor (Busch-extremal ⇒ uniquely Hamiltonian) | CITED by combination (Moon-Yang; THM-4744 (2); Burzio-Demaria 1990); short proof here |
| interval chains: one Hamiltonian cycle, all n | PROVED (cut argument); checked through order 17 |
| cycle ceiling; second floor f_2 = 13, 23, 41, 69, 115, 205; f_2(n) < f(n+1) | FINITE-EXACT through order 10 (P2, P4 pre-registered) |
| bottom band = uniquely Hamiltonian; all its classes are one-arc joins | FINITE-EXACT through order 10 (P3); HYP-9457 (b) for all n |
| crevasses as holes of U(n) (7 between U(4), U(5); 21 the hole of U(6)) | FINITE-EXACT localization; impossibility is THM-1370 |
| one-arc join law H = Lorentzian form, hc = s_a e_b | PROVED; checked on 9,000 random joins; klein's tower step and near-join are special cases |
| join coverage [23, 84705] at orders <= 7; comet | FINITE-EXACT (P5 confirmed) |
| 3-adic singular series | FINITE-EXACT; P6 refuted as stated (1.37) |
| real-point decomposition and parity crevasses | PROVED (elementary) |
| Pfaffian | FINITE-EXACT null |
| Cayley analogy | INTERPRETATION |

## 9. Pre-registration record

- HYP-9457 reserved in `9db3423759`; predictions P1-P6 committed in `4d5099bc95`, before `census 10` and `joins 7 7`.
- Census 3-9, chains 17 and joins 5 5 / 6 6 were exploratory.
- The real-point and Pfaffian items were exploratory.

## 10. Audit

(Filled after the independent audit.)

## 11. Reproduction

```
gcc -O2 -fopenmp -o scf 04-computation/experiments/scf_crevasses_20261104.c
./scf census 10        # 4.5 min on 12 threads (order 10: 9,733,056 classes)
./scf chains 17        # 7 s
./scf joins 7 7 200000 # under 1 s
```

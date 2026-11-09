# Genesis: generating all tournaments from first creators, and the spanning-tree twin of the 7/21 proof

**Session:** mac-mini-2026-11-09-genesis. **Canon:** THM-4985, THM-4986, HYP-9565.
**Proofs:** [PROOF.md](PROOF.md). **Code and data:**
[`04-computation/experiments/genesis_20261109/`](../../../04-computation/experiments/genesis_20261109/).
**Audit:** [AUDIT_G1.md](AUDIT_G1.md): PASS WITH CORRECTIONS (3 MAJOR, 23 MINOR, all applied; MISTAKE-743).

## The owner's directive (2026-11-09)

> for every odd number besides 7 and 21 we can keep track of the structure that first created them as we slowly
> increase tournament size, see how all tournaments can be generated from this minimal set corresponding to first
> appearances of hamiltonian path counts if you get creative enough with the operations done on them. consider the
> objects used in the proofs that 7 and 21 never appear to be related analogously with fundamental concepts and
> constructions

## Prior work this builds on (credited, not re-claimed)

- **First creators and their genealogy.**
  - THM-4922 / HYP-9512 (kind-pasteur S16): the one-class-or-converse-pair dichotomy fails from h = 9; genealogy;
    63 is the orphan.
  - THM-4930 / HYP-9520 (mac-mini, 2026-11-08): the eternal colours.
- **Generation.**
  - THM-4942 / HYP-9532 (kind-pasteur S17): the five-operation grammar, exhaustive through order 11; the connectivity
    floors f_k; the first section-5 dictionary of the 7/21 proof objects.
  - THM-4940 note §5 (opus S86): a second typed dictionary (OCF as a base-2 numeral, f-vectors, free monoid,
    numerology typed as such). Neither dictionary has the spanning-tree twin.
  - Session230: prime kernels; the Schmerl-Trotter all-order grammar; the positive cut kernel of the attachment form;
    the apex-transitive 7/21 certificate.
- **Spectra and bases.**
  - klein-2026-11-06: the exact censuses S(n), n <= 11.
  - THM-4940 (opus S86): packing counts in every base.
  - THM-4744 / THM-4745: Moon's interval chains and their product law.
- **Families and questions this session builds on (credited after audit G1).**
  - THM-4880 (mac-mini-2026-11-04-ribbon): the band tournaments are its constant-letter ribbon tournaments. Its §F lists
    their H-series for s = 2, 3, and its Theorem A's window argument is behind Lemma B3.
  - Session230 (primes/REPORT.md §5, SYNTHESIS §6): the bounded-defect question that HYP-9565 (A) makes precise.
  - THM-329 (opus-2026-05-27-S6): the largest known counts 531,205 (order 12) and 3,719,831 (order 13).
  - THM-4805 (opus-2026-10-26-S60): an earlier crevasse dictionary (untouchables and forbidden path counts).
- **Tools and inherited facts.**
  - THM-1370: the 7/21 omission proof.
  - THM-4811 (B): lambda = minimum semidegree for n <= 4k + 1.
  - THM-589: the variance of H.
  - HYP-2879 / THM-4097 / THM-4114: the attachment form.

## Results

### 1. Most tournaments are first creators (THM-4985 (6))

The share of first-creator isomorphism classes rises at every order from 6 to 11:

| order | 6 | 7 | 8 | 9 | 10 | 11 |
| --- | --- | --- | --- | --- | --- | --- |
| share of first-creator classes | 0.589 | 0.719 | 0.793 | 0.855 | 0.888 | **0.920** |

The order-11 share is new. The reason is visible in the numbers. The largest count of order n - 1 already lies below
the average count of order n: max H(n - 1)/E[H(n)] is 0.40 at n = 11 and 12, and at most 0.75 at n = 13 (max H(12)
is not known exactly). So at these orders a typical new tournament lands above the old spectrum and is automatically a
first creator; at every order this is conditional on (M_eps).

PROVED, conditional on an upper bound max H(m) = o(m·m!/2^(m - 1)) that is open (Alon's bound has m^(3/2)):
- almost every tournament is a first creator;
- every tournament is an induced subtournament of a first creator.

So the "minimal set" is not small. The generation problem is entirely about the **non-creators**: the minority whose
count is too small for their size. That is exactly the territory of the 7/21 proof objects, which are floors on small
counts.

### 2. The five-operation grammar fails at order 12 (THM-4985 (1)-(3))

S17 showed that through order 11 every tournament is one of: a first creator, a first creator plus a vertex, an ordinal
sum, a substitution, or an ordinal sum re-connected by at most two arcs. HYP-9532 (F) asked whether this persists.

**It does not. Order 12 is the first failure.**
- Extending S17's seven 3-arc-strong order-11 non-creators gives **12 escaper classes** (6 converse pairs):
  - prime and 3-arc-strong;
  - H in S(11) (48,997 .. 52,303);
  - every one of their 12 vertex deletions has H in S(10).
- Annealing finds 181 more classes (193 in all), with H from 32,897 to 58,017. Two self-converse ones (H = 48,537 and 48,849) have **decomposability defect 3**: no tournament within
  two arc reversals of them is non-strong or non-prime. So "decomposable up to two reversals" does not repair the
  grammar either.
- Band tournaments give certified escapers. In a band, each vertex beats the s vertices just before it and everything
  at distance > s after it. They are THM-4880's ribbon tournaments with the constant letter {1..s}, whose H-series
  THM-4880 §F already lists for s = 2, 3 (audit G1, M2). The certified escapers are:
  - band_3(13), with lambda = 3;
  - band_4(16), with lambda = 4, so even four re-connecting arcs are needed at order 16.

  Membership in S(12), S(14) and S(15) is certified by explicit witness tournaments.
- **PROVED, conditional on a superexponential solid prefix** (implied by HYP-9481 (a)): for every k and r, band
  tournaments of large order escape the grammar with "plus at most k vertices" and "re-connected by at most r arcs".
  The re-connection number S17 asked about is unbounded, and HYP-9532 (D)'s k*(n) -> inf.

**How the order-12 escapers come from first creators.**
- 106 of the 193 classes are **two first creators of order 6 in ordinal sum, re-connected by three back arcs**. The
  values are in {33, 37, 41, 43, 45}; the defect-3 pair is two first creators of 45, and the depth-3 class is two
  first creators of 33.
- 47 classes are **a substitution of first creators plus one vertex**, like S17's parent C3[R5, R5, 1] with an extra
  vertex.
- All but one are **a first creator of order 10 plus two vertices**.

So at order 12 the creative operations are "ordinal sum, then re-connect" and "substitute, then add a vertex", applied
to first creators, and, for the remaining 53 classes, "a first creator of order 10 plus two vertices". The first
description is nearly forced by the score sequence 3^6 8^6, so it describes rather than explains (audit G1, m3).

The escapers are what the "creative operations" must produce. They are old values reached by new, highly connected,
prime structures. In Chan-Pak's language (*Computational complexity of counting coincidences*, arXiv:2308.10214, §6.5)
they are counting coincidences (ANALOGY), not compositions. The bounded operations that work through order 11 are not
enough. Is there a bounded "defect" operation at all, or must the defect grow with the order? This bounded-defect
question is Session230's (primes/REPORT.md §5, SYNTHESIS §6). HYP-9565 (A) makes it precise: escapers are the prime
non-creators that THM-4942's grammar misses, and dec measures the defect.

### 3. The objects of the 7/21 proof have an exact spanning-tree twin (THM-4986)

The proof that no tournament has 7 or 21 Hamiltonian paths (THM-1370) has four steps:
- multiplicativity over strong components;
- the insertion lemma;
- Moon's strong-deletion theorem;
- the Moon-Busch floor.

Read with "spanning trees of a simple graph" in place of "Hamiltonian paths of a tournament", the same four steps prove
the classical fact that **no simple graph has exactly 2 spanning trees**. Every row below is a theorem on both sides:

| 7/21 proof object (tournaments, H) | twin (simple graphs, tau = spanning trees) |
| --- | --- |
| strong components; the condensation is a linear order (ordinal sum) | blocks; the block-cut tree (gluing at cut vertices) |
| H(T) = product over strong components (THM-1370 (1)) | tau(G) = product over blocks |
| H = 1 iff transitive; singletons are the trivial pieces | tau = 1 iff tree; bridges are the trivial pieces |
| insertion lemma H(T - v) <= H(T), equality iff v is a singleton component (e.g. padding by a source or sink) | vertex insertion tau(G) >= deg(v)·tau(G - v) >= tau(G - v) for a non-cut vertex v, equality iff v is pendant (padding by a leaf) |
| Moon: every strong n-tournament (n >= 4) has >= 2 strong (n - 1)-subtournaments; Moon chains | Whitney: every 2-connected graph is a cycle plus ears; ear decompositions |
| Moon-Busch floor f(m) = 3, 5, 9, 15, 25, 45, ...; extremals = Moon's interval chains, i.e. blocks NT_3 = C3, NT_4, NT_5 (parts 2, 3, 4) glued at shared end vertices (THM-4744, Moon-Yang) | block floor tau(B) >= m, equality iff B = C_m (Lemma T1) |
| the crevasses: 7 in (5, 9) = (f(4), f(5)); 21 in (15, 25) = (f(6), f(7)) | the crevasse: 2 in (1, 3) = (tau(K_2), tau(C_3)) |
| 7 prime; 21 = 3·7 recurses into 7 | 2 prime |
| OCF: H = I(Omega(T), 2), a hard-core gas of odd cycles (Grinberg-Stanley) | Kirchhoff: tau = det of the reduced Laplacian |
| first creators; birth order n(h) | Sedláček's minimal graphs; alpha(N) = least order with tau = N (1970) |
| n(h) >= m(h) ~ log h/log log h; HYP-9481 (a) predicts n(h) <= m(h) + 1 | alpha(N) = O((log N)^(3/2)/log log N) (Stong 2022); conjectured o(log N) (Azarija-Škrekovski 2013) |
| the bof-Royle conjecture (Chan-Pak, Conjecture 6.2): every odd h except 7, 21 is a count (OPEN) | every N except 2 is a count (cycles C_N) |
| values born at their floor bound, n(h) = B(h) = max{m : f(m) <= h}: exactly 3, 5, 9, 11, 13, 17, 19, 23, 35, 39, 63 on S(11) (THM-4942 (6); conjecturally all, HYP-9532 (E)) | fixed points alpha(N) = N, i.e. only the cycle C_N is minimal: 3, 4, 5, 6, 7, 10, 13, 22 (Nebeský 1973; case 22 completed by Azarija-Škrekovski 2013; OEIS A182258; alpha(3..7) = 3..7 and alpha(10), alpha(13), alpha(22) >= 9 re-checked here by exhaustive enumeration of graphs on <= 8 vertices) |

The crevasse lemma (Lemma T3) states the shared proof once. Both gaps are primes or prime multiples that fall between
consecutive irreducible floors and cannot be reached by products.

Two differences are as informative as the matches:
- **The floors grow differently.** The graph floor is linear (m), so cycles fill every value >= 3 and only 2 is lost.
  The tournament floor is exponential (5^(m/3)). It leaves crevasses (f(m), f(m + 1)) that only the strong spectra of
  orders <= m can fill, and 7 and 21 are missed by them. 21 lies inside the order-6 band [15, 45] but is a hole of the
  strong order-6 spectrum. Conjecturally (bof-Royle; Chan-Pak, Conjecture 6.2) these two odd values are the only ones
  never rescued.
- **First creators behave in opposite ways.** For tournaments the old maximum lies below the new average (ratio
  0.40 at n = 11), so first creators are typical (class shares 0.589 ... 0.920). For graphs the old maximum tau(K_(m - 1)) towers over the new average (ratio
  2^(m - 1)(m - 1)^(m - 3)/m^(m - 2) -> inf). Even so, new
  counts are common at small orders there too: the labelled share of new spanning-tree counts is 0.419, 0.515, 0.572 at
  6, 7, 8 vertices (pre-registered Q3 predicted a fall: REFUTED). The reason is that the graph spectrum is very sparse
  (OEIS A182290: 3,700 values below 262,144 at 8 vertices), while the tournament spectrum is nearly solid.

## Conjectures (HYP-9565)

- **(A)** The decomposability defect of escapers is unbounded. It is at least 3 at order 12.
- **(B)** The first-creator share tends to 1.
- **(C)** k*(12) = 3: no order-12 non-creator is 4-arc-strong.
- **(D)** Down-universality: every tournament is an induced subtournament of a first creator.
- Pre-registered (commit fb31d2c38a), all three predictions wrong:
  - **Q1 and Q2 REFUTED.** Two order-13 annealing runs of 400,000 steps found no escaper with lambda = 4. The least H found for minimum semidegree 4 at order 13 is 926,401, above 531,205, the largest known count of order 12 (THM-329; max H(12) is not known exactly). At order 12 the substitution C3[R5, R5, TT_2] has lambda = 4 and H = 229,725, above max H(11) = 95,095. This is evidence for (C).
  - **Q3 REFUTED** (see above).

## Typing

- PROVED: Lemmas W, B1, B2, B3, T1, T3; the conditional theorems of PROOF §§4-5.
- FINITE-EXACT: the order-12 escapers, the band certificates, the shares, and the tau spectra.
- KNOWN and cited: Sedláček; Nebeský (the fixed points of alpha); Azarija-Škrekovski; Stong; Daykin-Baines and
  Pomerance-Selfridge (coprime mappings); Szele; Alon; Friedgut-Kahn; Adler-Alon-Ross; the bof-Royle conjecture as
  Chan-Pak's Conjecture 6.2; OEIS A182290, A182258, A038375.
- ANALOGY: the dictionary rows as a correspondence. Each side is a theorem.
- EXPLORATORY: the annealing searches. Every hit is verified exactly.

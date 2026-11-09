# THM-4942: first-creator genealogy at order 11, and how every tournament is generated from the first creators

kind-pasteur-2026-11-08-S17. Status by section:
- **PROVED**: the corollary of section 3.2, the theorem of section 4.1, and the bound on M_ab in section 1.
- **KNOWN**: section 1 (HYP-2879, THM-4097, THM-4114) and the lemma of section 3.2 (THM-4811 (B)).
- **FINITE-EXACT**: sections 2, 3.1, 3.3 and 4.2, and the last sentence of 4.1.
- **ANALOGY** (typed as such): section 5.
- **CONJECTURE**: HYP-9532.

The independent audit is PASS WITH CORRECTIONS ([audit/AUDIT.md](audit/AUDIT.md); MISTAKE-737). It reproduced every
count, floor, witness and prediction outcome at orders 3-11 with its own code. It corrected three majors:
- the Busch-bound scope (section 4.2);
- an unproved sentence in section 4.1;
- the missing credit THM-4811 (B) for the lemma of section 3.2.

A fourth credit came from the concurrent Session230 after the audit: the attachment form of section 1 is HYP-2879 /
THM-4097 / THM-4114. All minors are applied below. The audit also settled the class-level genealogy at order 11
(section 2).

Pre-registrations: [PREREGISTRATION.md](PREREGISTRATION.md) (R1-R9, af0112b596) and
[PREREGISTRATION-2.md](PREREGISTRATION-2.md) (R10-R11, c201ec7b34).

Owner prompt (2026-11-08): run the order-11 test. For every odd number besides 7 and 21, keep track of the structure
that first created it as the order grows. See how all tournaments can be generated from this minimal set if you get
creative enough with the operations. Relate the objects used in the proofs that 7 and 21 never appear to fundamental
concepts and constructions.

**Credited prior work.**

| Source | What it supplies |
| --- | --- |
| THM-4922 and HYP-9512 (kind-pasteur S16) | the first creators, genealogy through order 10, the orphan pair of 63 |
| Session224 | the first-occurrence ledger and its refutation of the dichotomy |
| THM-1370, Moon and Busch | the omission of 7 and 21: multiplicativity, the insertion lemma, Moon's strong-subtournament theorem and the Moon-Busch floor f(n) = min 3^a 5^b |
| THM-4744 | Moon's interval chains attain Busch's bound |
| THM-029 / THM-343, THM-115 and Grinberg-Stanley (arXiv:2412.10572) | the OCF route |
| THM-1900 (B) | H-neutral insertions are exactly the condensation down-sets |
| **HYP-2879, THM-4097 (sections 2-3) and THM-4114 (section 6)** | the attachment form of section 1: the insertion formula with start, end and slot counts, the subset-DP slot matrix Q_ab = W[a][b], the cut-field identity, and the submodular directed-cut quadratic. THM-4118 gives the ear-response lattice |
| THM-012b and THM-4094 | the insertion decomposition |
| THM-4169 | the Schmerl-Trotter prime-deletion route. At each odd order >= 5 the only critical prime tournaments are T, U and W, so every prime order-11 tournament outside T_11, U_11 and W_11 has a prime order-10 deletion. This is context for SUB, REST (all prime) and prime classes. Its exhaustive prime-parent augmentation census from order 10 to 11 is the precedent for this packet's order-11 run (audit m18) |
| Session230 (`session230_creator_grammar_20261109`, concurrent, same owner prompt, building on this packet's checkpoint) | the priority pointer to THM-4097 / THM-4114, and the remark that the first creators are minimal in order, not an irredundant basis (section 3.1) |
| opus-2026-11-08-S86 | reserved THM-4940, THM-4941 and HYP-9530 at 19:20:14 for a combined prompt that includes first appearances, generation and the 7/21 objects. This packet's registration was af0112b596 at 19:19:56 |
| THM-030 and THM-050 | the split products E_a(S + a) B_b(R + b). Their unsigned sum is W[a][b]; when a -> b it counts the Hamiltonian paths through a -> b |
| Session189 | the degree-two response polynomials |
| **THM-4811 (B)** | lambda = min(delta^+, delta^-) whenever n <= 4 delta^0 + 1, which contains the lemma of section 3.2; the class counts by lambda through order 10 |
| THM-4869, HYP-9459 and HYP-9472 | Hamiltonian-cycle floors under connectivity, the closest in-repo antecedent of the graded floors |
| Thomassen 1980; Busch 2006, Corollary 2 | exponential lower bounds for Hamiltonian cycles of 2-connected tournaments (vacuous below about 32 vertices) |
| klein's census (THM-4889; `path_spectrum_two_angles_metagraph_20261106`) | S(11), n(h), the strong windows W(h), the Busch bound B(h), HYP-9481 |
| Sessions 225-229 | the graft and two-path-cover transfer calculus for ordinal sums with back arcs |
| mac-mini's THM-4930 | the self-converse colours, a different question |

## 1. The attachment form (KNOWN: HYP-2879, THM-4097 sections 2-3, THM-4114 section 6; restated)

**Priority.** This formula is not new.
- THM-4097 eq. (4), from HYP-2879, gives H(T + x_S) in terms of start counts s_b, end counts e_a and the slot matrix
  Q_ab.
- Its eq. (5) builds Q_ab = sum_A E_A(a) B_{V\A}(b) by the same subset DP as W[a][b] below.
- Its section 3 gives the cut-field identity H(T + x_S) = H(T) + cut_w(S) + sum_{i in S} h_i.
- THM-4114 section 6 shows the response is a submodular directed-cut quadratic, with pair differences
  -[Q(u, v) + Q(v, u)], which is -M_uv here.

The pre-audit text presented the formula as this packet's. Session230's README pointed to the priority after the
audit, and the credit is recorded in MISTAKE-737. This packet's own contribution in this section is the engine: a
Gray-code walk over all 2^n patterns with O(n) updates per pattern, together with the n deletion forms that track
every H(T' - u). It is run exhaustively at order 11.

Let T be a tournament on V with |V| = n. Write f[S][a] (respectively g[S][a]) for the number of Hamiltonian paths of
T[S] that end (respectively start) at a. Add a vertex v, with x_a = [a -> v]. Then

    H(T + v_x) = H(T) + sum_a x_a L_a - sum_{a<b} x_a x_b M_ab,                          (1)
    W[a][b] = sum over S with a in S, b not in S, of f[S][a] g[V\S][b],
    L_a = f[V][a] - g[V][a] + sum_{b != a} W[a][b],      M_ab = W[a][b] + W[b][a].

*Proof.* Sort the Hamiltonian paths of T + v by the position of v.
- With v first, there are sum_b (1 - x_b) g[V][b] paths.
- With v last, there are sum_a x_a f[V][a].
- Otherwise the path is (a path on S ending at a) v (a path on V\S starting at b), with x_a = 1 and x_b = 0.

Expanding x_a (1 - x_b) and using sum_b g[V][b] = H(T) gives (1). QED.

Three consequences:
- W[a][b] counts the Hamiltonian paths of any T + v_x with x_a = 1, x_b = 0 and v between a and b. So W[a][b] and
  L_a are at most max H(n + 1).
  - M_ab is bounded the same way. Say b -> a, and take x with x_a = 1, x_b = 0. Then W[a][b] counts the paths of
    T + v_x containing a v b. Each pair counted by W[b][a], joined by the arc b -> a, is a Hamiltonian path of T
    containing "b a". Inserting v at a feasible place never puts it between b and a, because v -> b.
  - These paths are therefore distinct from the first family and from each other. So
    M_ab <= H(T + v_x) <= max H(n + 1). If a -> b, exchange the roles. (Audit m1.)
  - Every coefficient fits in 32-bit integers at n <= 11.
- L_a = H(T + v_{e_a}) - H(T) >= 0. This is the insertion lemma at a single in-arc.
- The pattern x is H-neutral exactly when the set {b : x_b = 0} that v beats is a condensation down-set (a union of
  terminal strong components). Equivalently, {a : x_a = 1} is a union of initial components. There are #SCC(T) + 1
  such patterns (THM-1900 (B)).

**Notation.** M_ab here is unsigned. It is not the signed transfer matrix M[a, b] of THM-030 and THM-050, which is
built from the same split products.

Degree two in x is also folklore: a path uses at most two arcs at v (Session189's response polynomials). The
coefficients turn the whole extension layer into one Gray-code walk with O(n) updates per pattern. This is the engine
`ext.c`.

## 2. Order 11 (FINITE-EXACT, pre-registered)

`ext.c` covered all 9,733,056 classes of order 10 and all 1,024 patterns of each: 9,966,649,344 extensions, on 11
processes in about 20 minutes ([`order11_summary.txt`](order11_summary.txt)). The audit reran it exhaustively on its
own enumeration of order-10 classes with an independent engine. It checked every one of the 36,248 R1 witnesses by
direct DP.

Controls:
- The extension values reproduce klein's S(11) exactly: 43,911 values, maximum 95,095. That is an independent
  check of klein's census.
- 22 classes were checked by direct subset DP for every pattern, including every deletion count. There were 0
  failures.
- On classes of order 7, 8 and 9 the same engine reproduced S(8), S(9) and S(10). Its only new value without a
  first-creator parent was 63, at order 8, so the gate passed.

| Prediction | Outcome |
| --- | --- |
| R1: HYP-9512 (A) at order 11 | **CONFIRMED.** Each of the 36,248 values new at order 11 is H(T + v) for some first creator T of order 10. Through order 11, 63 remains the only odd value h >= 5 not created by a one-vertex extension of a first creator one order down. |
| R2: HYP-9512 (B) at order 11 | **CONFIRMED.** There are 0 hypo pairs: no order-11 tournament whose H lies in S(10) has every vertex deletion a first creator. |
| R3: no sterile first creator at order 10 | **CONFIRMED.** |
| R4: at least half the patterns fertile | **CONFIRMED.** Every order-10 first creator has at least 564 fertile patterns out of 1,024. The minimum is attained by two classes, a converse pair of first creators with H = 2,901 and \|Aut\| = 1. |
| R5: share reached only from first-creator parents | **CONFIRMED**: 18,661 / 36,248 = 0.515, after 0.350, 0.425 and 0.470 at orders 8-10. |

**Scope.**
- R1 is the value-level statement of HYP-9512 (A).
- The class-level statement was settled by the independent audit (audit section 3.6) and needs no order-11
  canonical forms. An order-11 first creator with no first-creator deletion would extend a non-first-creator of
  order 10. So the audit extended all 1,085,824 non-first-creator classes of order 10 by every pattern.
- Result: every first creator of order 11 has a first-creator deletion. There are 0 class-level orphans among the
  506,912,096 extensions that are first creators. Positive controls recover the 63 pair at order 8 and find
  nothing at order 10.
- This result was not pre-registered.
- R2 is exhaustive, because a hypo analogue at order 11 is in particular an extension of an order-10 first
  creator.

**Fertility through order 10.** At orders 3-10:
- every first creator has a fertile pattern;
- every first creator has at least half its patterns fertile, except the two orphan creators of 63 (77 of 256).

The minimum fertile counts at orders 3-10 are 6/8, 14/16, 26/32, 51/64, 75/128, 77/256 (the 63 pair), 299/512 (at
H = 615) and 564/1024.

**Realization multiplicities at order 11.** Count, for each new value, the (first-creator class, pattern) pairs that
realize it. The count ranges from 1 to 445,304. The median is 251,107, the upper of the two middle values 251,103 and
251,107. Exactly three new values have a single pair:
92,411, 93,467 and max H(11) = 95,095 (the Paley tournament QR_11). All three are at the top of the band.

**The lineage ledger** ([lineage_le10.txt](lineage_le10.txt), `lineage.py`) records the owner's "structure that first
created them" for all 7,661 values h >= 5 with n(h) <= 10. For each value h it gives:
- n(h);
- the number of first-creator classes;
- how many of them have a first-creator parent;
- the parent values, that is, the H of their first-creator deletions.

Its first rows read as a tree rooted at the 3-cycle:

| Parent | Children |
| --- | --- |
| 3 | 5 |
| 5 | 9, 11, 13, 15 |
| 9 | 17 |
| 9 and 11 | 19 |
| 9, 11 and 13 | 23 |
| 11 and 13 | 45 = max H(6) |

At order 7 come 35 (from 17 and 19), 39 (from 17, 19 and 23), and so on. Every first-creator class of every value
except 63 has a first-creator parent. 63 is the single line marked ORPHAN. At order 11 the ledger is value-level:
[witnesses_order11.txt](witnesses_order11.txt) gives one first-creator parent and pattern per new value.

**Genealogy roots.**
- At order 9, 152 first creators have a creator of 63 as a deletion. Every one of them also has another
  first-creator deletion.
- With THM-4922's class-level genealogy through order 10, every first creator of order <= 10, other than the two
  creators of 63, lies on a chain of first creators. The chain has one member at each order from 3 (the 3-cycle) up
  to its own.
- With the audit's class-level test, the same holds for every first creator of order 11.
- These are Moon chains inside the minimal set (section 5). The first creators are strong through order 11.

## 3. Generating every tournament from the first creators

### 3.1 The census (FINITE-EXACT, orders 3-10; `gen.c`)

Every class of order n is put in the first category that fits, in this order:

| Category | Meaning | How H is computed |
| --- | --- | --- |
| **FC** | a first creator: H not in S(n - 1) | the minimal set |
| **EXT** | some deletion T - v is a first creator of order n - 1 | the attachment form (1) of a first creator |
| **SUM** | T is not strong: an ordinal sum of smaller tournaments | H multiplies |
| **SUB** | T has a nontrivial module: a substitution of smaller tournaments | path-cover profiles |
| **REST** | strong, prime, not a first creator, with no first-creator deletion | |

| n | classes | FC | EXT | SUM | SUB | REST |
| --- | --- | --- | --- | --- | --- | --- |
| 3 | 2 | 1 | 0 | 1 | 0 | 0 |
| 4 | 4 | 1 | 2 | 1 | 0 | 0 |
| 5 | 12 | 6 | 2 | 4 | 0 | 0 |
| 6 | 56 | 33 | 14 | 9 | 0 | 0 |
| 7 | 456 | 328 | 90 | 37 | 1 | 0 |
| 8 | 6,880 | 5,457 | 1,162 | 216 | 33 | 12 |
| 9 | 191,536 | 163,676 | 23,855 | 2,489 | 1,010 | 506 |
| 10 | 9,733,056 | 8,647,232 | 962,686 | 49,755 | 42,514 | 30,869 |

**The residue REST is built from an ordinal sum.**
- Every REST class at orders 8 and 9 has a vertex whose deletion is not strong. Such a class is an ordinal sum plus
  one vertex: two blocks re-connected through that vertex. At order 10, 30,649 of the 30,869 REST classes are of this
  kind.
- Most of these, 12, 505 and 30,525 at orders 8, 9 and 10, have a **pendant** vertex of in- or out-degree 1. That is
  the floor move of Moon and Busch: deleting the unique in-neighbour of an in-degree-1 vertex leaves an ordinal sum.
- The remaining 220 classes at order 10 have arc-strong connectivity exactly 2. Reversing two arcs makes them an
  ordinal sum of two blocks.
- The single order-9 class that is not pendant has H = 385. It is self-converse with trivial automorphism group: two
  strong 4-tournaments in ordinal sum (H = 5 each), with a middle vertex that splits both. It is the λ = 2 floor at
  order 9 (section 3.3).

**Pre-registered order-10 outcomes.**

| Prediction | Outcome |
| --- | --- |
| R6: REST is 2-6% of non-first-creators | **CONFIRMED**: 30,869 / 1,085,824 = 2.84% |
| R7: at least 99% of REST pendant | **FAILED**: 98.89% |
| R8: at most 100 classes outside FC, EXT, SUM, SUB and the pendant extensions | **FAILED**: 344 = 4 single-back-arc grafts with both sides of size >= 2, plus 340 classes with lambda = 2 |
| R9: EXT at least 85% of non-first-creators | **CONFIRMED**: 88.66% |

**Lesson from R7 and R8.** "Pendant" is one instance of the general mechanism, small arc-strong connectivity. The
right invariant is lambda, in section 3.3.

**Answer to the owner's question, through order 11.** Through order 10 this follows from the census and lambda <= 2.
At order 11 it follows from the exhaustive lambda >= 3 list in section 3.3. Every tournament is one of the following:
- a first creator;
- a first creator plus one vertex;
- an ordinal sum;
- a substitution;
- an ordinal sum re-connected by at most two arcs.

Each operation has an explicit H-calculus:
- the attachment form (1);
- the product;
- the path-cover profile of the module;
- for re-connection, the graft calculus of Sessions 225-229. Reverse one arc a -> b of A => B to b -> a. Then

      H = H(A) H(B) - F_A(a) G_B(b) + U_A(a) U'_B(b).

  The terms are:
  - F_A(a) G_B(b) counts the destroyed paths through a -> b (paths of A ending at a, times paths of B starting at
    b);
  - U_A(a) counts the ordered two-path covers (P_1, P_2) of A with P_2 starting at a and P_1 possibly empty;
  - U'_B(b) counts the ordered two-path covers (Q_1, Q_2) of B with Q_1 ending at b and Q_2 possibly empty.

  This was checked against direct DP on 300 random cases (`graft_check.py`). The audit's m9 wording omits the
  subtracted term. Two back arcs need the analogous transfer with ordered three-part covers.

Precisely (audit m10): every class of order 3-11 is one of the following:
- a first creator;
- a one-vertex extension of a first creator;
- obtained from strictly smaller tournaments by an ordinal sum, a substitution, or an ordinal sum with at most two
  arcs reversed.

The highly connected (lambda >= 3) classes are first creators or, at order 11 only, seven substitutions. The first
creators are 88.8% of all classes at order 10.

**Minimal in order, not an irredundant basis** (Session230). The first creators are minimal in vertex order for
their counts. They are not an irredundant generating set for this grammar. For example, the unique strong
4-tournament, the first creator of 5, is C3[TT2, 1, 1], a substitution of smaller tournaments. The statement above
is a covering of all classes, not a minimal basis.

At order 11 the argument runs as follows:
- A non-first-creator with lambda <= 2 falls in the last category.
- A non-first-creator with lambda >= 3 is one of the seven classes of section 3.3, all substitutions.
- A first creator is a first creator.

No class-level census of order 11 is needed.

### 3.2 Connectivity and degrees (lemma KNOWN, corollary PROVED)

Let lambda(T) be the arc-strong connectivity of T: the least number of arcs from V\X into X over nonempty proper X.
Equivalently, it is the fewest arc reversals that make T non-strong, both for strong and for non-strong T:
- lambda = 0 is exactly "not strong";
- for lambda in {1, 2}, reversing the lambda back arcs of a minimum cut gives an ordinal sum;
- conversely, reversing at most two arcs of an ordinal sum leaves a cut crossed backward by at most two arcs.

**Lemma (KNOWN: THM-4811 (B), the case delta^0 = 3).** For n <= 13, lambda(T) >= 3 if and only if every vertex has
in- and out-degree at least 3.

*Proof.* Singleton cuts give lambda <= min degree. Conversely, suppose every vertex has in- and out-degree at least 3,
and take a cut (X, Y) with |X| = k. The arcs from Y into X number

    sum_{x in X} indeg(x) - C(k, 2) >= k(7 - k)/2,

which is at least 3 for 1 <= k <= 6. By symmetry, using out-degrees in Y, the same holds when 1 <= n - k <= 6. When
n <= 13 one of the two holds. QED.

The bound is sharp. R7 => R7 (n = 14) has min degree 3 and lambda = 0; with one or two back arcs, lambda = 1 or 2.
THM-4811 (B) proves lambda = delta^0 whenever n <= 4 delta^0 + 1, by the same count. Its census already lists the
class counts by lambda through order 10 that section 3.3 reproduces.

**Corollary.** Here f_k(n) is the least H over n-tournaments with lambda >= k.
- For n <= 13, if f_3(n) > max H(n - 1), then every tournament of order n whose scores all lie in [3, n - 4] is a first
  creator.
- Under the same hypothesis, every non-first-creator of order n becomes an ordinal sum after at most two arc
  reversals. This half needs neither the lemma nor n <= 13.

**Notation.** The connectivity floors f_k (k >= 2) are not HYP-9472's and klein's two-cycle floor f_2. For example,
that floor has f_2(10) = 205, while here f_2(10) = 727. f_1 is the Moon-Busch floor f.

### 3.3 The connectivity floors (FINITE-EXACT; `lam.c`, `lam11.c`)

f_k(n) is the least H over n-tournaments with lambda >= k ([floors_4_10.txt](floors_4_10.txt)). The last two columns
count the non-first-creators with lambda = 1 and lambda = 2. Through order 11, every non-strong class (lambda = 0) is a
non-first-creator, because first creators are strong there (THM-4922 (1)).

| n | f_1 (Moon-Busch) | f_2 | f_3 | f_4 | max H(n - 1) | non-FC, lambda = 1 | non-FC, lambda = 2 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 4 | 5 | | | | 3 | 0 | 0 |
| 5 | 9 | 15 | | | 5 | 0 | 0 |
| 6 | 15 | 41 | | | 15 | 2 | 0 |
| 7 | 25 | 99 | 171 | | 45 | 25 | 0 |
| 8 | 45 | 195 | 621 | | 189 | 551 | 0 |
| 9 | 75 | 385 | 2,055 | 3,159 | 661 | 14,333 | 124 |
| 10 | 125 | 727 | 5,753 | 14,877 | 3,357 | 668,316 | 40,401 |
| 11 | | | **11,475** | | 15,745 | | |

Through order 10:
- f_3(n) > max H(n - 1). So by section 3.2 every class whose scores all lie in [3, n - 4] is a first creator.
- Every non-first-creator is within two arc reversals of an ordinal sum; through order 8, within one.
- The lambda = 3 and lambda = 4 counts equal the minimum-degree-3 and -4 counts
  ([mindeg_4_10.txt](mindeg_4_10.txt)), as the lemma requires.

**Order 11 (R10, R11, pre-registered: both CONFIRMED).**
- `lam11.c` extends every order-10 class with minimum in- and out-degree at least 2 by every pattern. This is
  exhaustive for 3-arc-strong 11-tournaments, because deleting any vertex of one leaves minimum degree at least 2.
  Lambda >= 3 is tested directly whenever H <= 15,745 ([lam11_summary.txt](lam11_summary.txt)).
- **Protocol deviation (audit m12).** The registration said "every order-10 class with lambda >= 2". The run used the
  larger, exhaustive set of classes with min degree >= 2, which adds 17 classes. The registered scope alone is not
  provably exhaustive.
- Result: f_3(11) = 11,475 < max H(10) = 15,745.
- Below 15,745 the 3-arc-strong 11-tournaments take exactly seven values: 11,475, 14,787, 15,213, 15,273, 15,561,
  15,637 and 15,703.
- Four of them, 11,475, 15,213, 15,273 and 15,561, already lie in S(10). Those tournaments are not first creators,
  although every score lies in [3, 7].
- So the threshold k*(n), the largest lambda of a non-first-creator of order n, is 0 for n <= 5, 1 for n = 6, 7, 8, 2
  for n = 9, 10, and at least 3 at n = 11 (HYP-9532 (D)).

**But the generation grammar survives at order 11.**
- `lam11_all.c` lists all 97 (class, pattern) pairs. Canonicalised (`canon11.c`), they give 11 classes:
  - 4 first creators: 14,787 (1 class), 15,637 (1) and 15,703 (2);
  - 7 non-first-creators: 11,475 (1 class), 15,213 (2), 15,273 (2) and 15,561 (2).
- Every one of the seven non-first-creators has a nontrivial module, so it is a substitution. Six of them also have a
  first-creator deletion.
- The floor class (H = 11,475, |Aut| = 25) is C3[R5, R5, 1], self-converse. It is built from the 3-cycle, the regular
  first creator R5 of 15 used twice, and a singleton. (15 has two first creators; the other has lambda = 1.)
- Since this list is complete for lambda >= 3, the five-operation description of section 3.1 holds exactly at order
  11 too ([lambda3_witnesses_order11.txt](lambda3_witnesses_order11.txt)).

**Growth of f_3.** The ratios f_3(n + 1)/f_3(n) are 3.63, 3.31, 2.80 and 1.99 (n = 7..10). The ratios
max H(n)/max H(n - 1) are 4.2, 3.5, 5.1, 4.7 and 6.0 (n = 7..11). So the 3-arc-strong floor falls behind the
maximum, and the crossing occurs exactly between orders 10 and 11. No growth law for f_k is claimed.

**Two floors coincide with exceptional classes.** f_2(9) = 385 is the H of the single non-pendant REST class at order
9. f_1 is the floor that excludes 7 and 21.

## 4. The fibers of a value

### 4.1 Strong fibers are finite antichains (PROVED; last sentence FINITE-EXACT)

For every odd h >= 3, the strong tournaments with H = h have pairwise no induced embedding. Their orders lie in the
strong window W(h), which is contained in [3, B(h)], where B(h) = max{n : f(n) <= h}:
- n <= B(h) follows from the definition of B (klein, THM-4889).
- Finiteness needs f(n) -> infinity. That comes from Busch 2006 / THM-4744, or more simply from Camion: a strong
  n-tournament has a Hamiltonian cycle, so H >= n.
- (The order-1 tournament is trivially strong with H = 1, hence h >= 3.)

*Proof.* Let A be an induced proper subtournament of a strong B. Choose a vertex v of B outside A. Then v is not a
singleton strong component of B, so THM-1900 (B) gives H(B - v) < H(B). By the insertion lemma, H(A) <= H(B - v). So
H(A) < H(B). QED.

So every tournament with H = h is an ordinal sum of singletons and members of the finite sets Strong(d), d >= 3, over
the ordered factorizations of h. Through order 11 the first creators of h are exactly the members of Strong(h) of
order n(h) (THM-4922 (1), FINITE-EXACT). In general this is HYP-9481 (b), OPEN.

### 4.2 Values born at their Busch bound

A value h is **born at its Busch bound** when n(h) = B(h), that is, at the last order where a strong realization is
possible. Among the odd h <= 80,451 (every one except 7 and 21 lies in S(11)), and more generally among all odd h in
S(11), these are exactly

    3, 5, 9, 11, 13, 17, 19, 23, 35, 39, 63

(klein's n(h) and the Moon-Busch floor).
- Equivalently, these are klein's length-1 strong windows below 375 that end at B(h). The other three length-1
  windows, 47, 55 and 67, end at B(h) - 1.
- The 3,635 odd values in [80,453, 95,093] absent from S(11) have B(h) = 22 and n(h) >= 12. For them, n(h) != B(h) is
  conditional on HYP-9481 (a), which gives n(h) = 12. (Audit M1; the pre-audit text claimed "through 95,095".)
- The two orphan values of THM-4922, 3 and 63, are on the list, and 63 is its largest member.
- 7 and 21 are the values with no strong realization at any order: W(7) and W(21) are empty.

## 5. The objects of the 7/21 proofs and the generation of all tournaments (ANALOGY, with precise anchors)

The proofs that 7 and 21 never occur (THM-1370; the OCF route of THM-029/343 and THM-115) use five objects. Each one
turns out to be one of the operations that generate the non-first-creators. The middle column is precise; the
last column is an analogy and is typed as such.

| Object in the 7/21 proofs | Its role there | Its role in generating all tournaments (this packet) | Fundamental analogue (ANALOGY) |
| --- | --- | --- | --- |
| **Condensation**: H(T) is the product of the H of the strong components (THM-1370 (1)) | 7 is prime, so it needs a strong component with H = 7; 21 = 3 x 7 recurses into 7 | SUM: every non-strong tournament is an ordinal sum. The first creators are strong (finite-exact through order 11, HYP-9481 (b)) | **Unique factorization.** Tournaments under ordinal sum form the free monoid on strong tournaments, and H is a monoid homomorphism to (odd integers, x), like a norm. Strong means irreducible, and the condensation is the (ordered) prime factorization. 7 behaves like a prime that is not a norm, with no composite route. Unlike the integers, every attained value above 1 is attained irreducibly (by a strong tournament: HYP-9481 (b), exact through order 11). |
| **Insertion lemma**: H(T) >= H(T - v), equality iff v is a singleton component (THM-1370 (2), THM-1900 (B)) | the strong floor is nondecreasing | EXT: the extension calculus of section 1 (quadratic in the attachment). Fertility and the genealogy run along it. Strong fibers are antichains (section 4.1) | **Insertion sort.** Redei's existence proof of a Hamiltonian path is insertion sort with an intransitive comparator. H is a rank function, strictly monotone on the strong part. Equal-rank irreducibles form an antichain, as in Sperner theory. |
| **Moon's theorem**: a strong n-tournament (n >= 4) has at least two strong (n - 1)-subtournaments (THM-1370 (3)) | f(n) >= f(n - 1) | Every strong tournament is reached from C3 by strong one-vertex extensions ("Moon chains"). Every strong tournament is reached from C3 by strong one-vertex extensions (Moon chains). Through order 11 every first creator other than the two creators of 63 has a Moon chain inside the minimal set: THM-4922 (3) and the 152-check through order 10, and the audit's class-level test at order 11. HYP-9512 (A) is the value-level conjecture | **Non-cut vertices and ear decompositions.** A connected graph has at least two non-cut vertices (the leaves of a spanning tree); a 2-connected graph has an ear decomposition (Whitney). |
| **Moon-Busch floor** f(n) = min{3^a 5^b : 2a + 3b = n - 1} (THM-1370 (4), Busch 2006) | 7 lies in (f(4), f(5)) = (5, 9) and 21 in (f(6), f(7)) = (15, 25), with no strong realization below. These are the crevasses | **Graded floors f_k(n)** (section 3.3). f_3(n) > max H(n - 1) forces every 3-arc-strong tournament to be a first creator through order 10, so every non-first-creator is at most two arc reversals from an ordinal sum. The k = 1 floor also makes strong fibers finite (section 4.1) and fixes the eleven values born at their Busch bound, 63 the largest (section 4.2) | **Minkowski and Hermite.** Fields of degree n have discriminant at least about c^n (Minkowski), so finitely many fields have bounded discriminant (Hermite). Likewise strong tournaments of order n have H >= f(n), about 5^(n/3), so finitely many strong tournaments have bounded H. The crevasses are band gaps between the strong bands [f(n), max H(n)]. The floor optimisation (3-cycles cost two vertices per factor 3, 5-blocks three vertices per factor 5, and 5^(1/3) < 3^(1/2)) is an optimal product-partition problem, like the "break n into 3s" maximum-product problem. |
| **OCF** (Grinberg-Stanley): H = I(Omega(T), 2) >= 1 + 2 alpha_1, with Moon's pancyclicity counts (THM-029/343, THM-115) | H = 7 forces Omega = K3, and H = 21 forces at most 10 odd cycles. Strong unions then create more cycles | Omega of an ordinal sum is the disjoint union of the Omegas, and I is multiplicative: the same factorization one level down | **The hard-core lattice gas** at fugacity 2. Odd cycles are polymers, and vertex-disjointness is hard-core compatibility. Redei's theorem is the parity of the partition function. |

**Further, outside the 7/21 proofs.** Redei's parity and THM-791 (self-converse classes at even order have H = 1
mod 4) play the role of Stickelberger's congruence: every discriminant is 0 or 1 mod 4. The S16 value colours are the
corresponding mod-4 splitting.

**The precise unifying statement** behind the dictionary is section 3.2 together with the floors of section 3.3:
- The same object, a lower bound on H in terms of connectivity, does both jobs.
- At k = 1 it is the Moon-Busch floor. That floor excludes 7 and 21.
- At k = 3 it exceeds max H(n - 1) (orders <= 10). Then every old-valued tournament is nearly decomposable, so it is
  built from smaller tournaments by an ordinal sum re-connected by at most two arcs.
- Every highly connected tournament (lambda >= 3) is a first creator (orders <= 10). The converse fails: many first
  creators have lambda = 1.

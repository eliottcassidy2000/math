# THM-4942: first-creator genealogy at order 11, and how every tournament is generated from the first creators

kind-pasteur-2026-11-08-S17. Status by section:
- **PROVED**: sections 1, 3.2 and 4.1.
- **FINITE-EXACT**: sections 2, 3.1, 3.3 and 4.2.
- **ANALOGY** (typed as such): section 5.
- **CONJECTURE**: HYP-9532.

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
| THM-029 / THM-343 and THM-115 | the OCF route |
| THM-1900 (B) | H-neutral insertions are exactly the condensation down-sets |
| THM-012b and THM-4094 | the insertion decomposition |
| Session189 | the degree-two response polynomials |
| klein's census (`path_spectrum_two_angles_metagraph_20261106`) | S(11), n(h), the strong windows W(h), the Busch bound B(h), HYP-9481 |
| Sessions 225-229 | the graft and two-path-cover transfer calculus for ordinal sums with back arcs |
| mac-mini's THM-4930 | the self-converse colours, a different question |

## 1. The attachment form (PROVED, elementary)

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
- W[a][b] counts the Hamiltonian paths of any T + v_x with x_a = 1, x_b = 0 and v between a and b. So every
  coefficient is at most max H(n + 1).
- L_a = H(T + v_{e_a}) - H(T) >= 0. This is the insertion lemma at a single in-arc.
- The pattern x is H-neutral exactly when it is a condensation down-set (THM-1900 (B)).

Degree two in x is folklore: a path uses at most two arcs at v (Session189's response polynomials). The explicit
coefficients turn the whole extension layer into one Gray-code walk with O(n) updates per pattern. This is the engine
`ext.c`.

## 2. Order 11 (FINITE-EXACT, pre-registered)

`ext.c` covered all 9,733,056 classes of order 10 and all 1,024 patterns of each: 9.97 x 10^9 extensions, on 11
processes in about 20 minutes ([`order11_summary.txt`](order11_summary.txt)).

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
| R4: at least half the patterns fertile | **CONFIRMED.** Every order-10 first creator has at least 564 fertile patterns out of 1,024. The minimum is at a first creator with H = 2,901. |
| R5: share reached only from first-creator parents | **CONFIRMED**: 18,661 / 36,248 = 0.515, after 0.350, 0.425 and 0.470 at orders 8-10. |

**Scope.** R1 is the value-level statement of HYP-9512 (A). The class-level statement is that every first creator of
order 11 has a first-creator deletion. It needs canonical forms at order 11 and is not tested here. R2 is exhaustive,
because a hypo analogue at order 11 is in particular an extension of an order-10 first creator.

**Fertility through order 10.** At orders 3-10:
- every first creator has a fertile pattern;
- every first creator has at least half its patterns fertile, except the two orphan creators of 63 (77 of 256).

The minimum fertile counts at orders 3-10 are 6/8, 14/16, 26/32, 51/64, 75/128, 77/256 (the 63 pair), at least
256/512, and 564/1024.

**Realization multiplicities at order 11.** Count, for each new value, the (first-creator class, pattern) pairs that
realize it. The count ranges from 1 to 445,304, with median 251,107. Exactly three new values have a single pair:
92,411, 93,467 and max H(11) = 95,095 (the Paley tournament QR_11). All three are at the top of the band.

**Genealogy roots.**
- At order 9, 152 first creators have a creator of 63 as a deletion. Every one of them also has another
  first-creator deletion.
- With THM-4922's class-level genealogy through order 10, every first creator of order <= 10, other than the two
  creators of 63, lies on a chain of first creators. The chain has one member at each order from 3 (the 3-cycle) up
  to its own.
- These are Moon chains inside the minimal set (section 5).

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

**Answer to the owner's question, through order 10.** Every tournament is one of the following:
- a first creator;
- a first creator plus one vertex;
- an ordinal sum;
- a substitution;
- an ordinal sum re-connected by at most two arcs.

Each operation has an explicit H-calculus:
- the attachment form (1);
- the product;
- the path-cover profile of the module;
- the two-path-cover transfer of the graft calculus (Sessions 225-229).

So every tournament through order 10 is generated from the first creators, which make up 88.8% of all classes at
order 10.

### 3.2 Connectivity and degrees (PROVED)

Let lambda(T) be the arc-strong connectivity of T: the least number of arcs from V\X into X over nonempty proper X.
Equivalently, it is the fewest arc reversals that make T non-strong.

**Lemma.** For n <= 13, lambda(T) >= 3 if and only if every vertex has in- and out-degree at least 3.

*Proof.* Singleton cuts give lambda <= min degree. Conversely, suppose every vertex has in- and out-degree at least 3,
and take a cut (X, Y) with |X| = k. The arcs from Y into X number

    sum_{x in X} indeg(x) - C(k, 2) >= k(7 - k)/2,

which is at least 3 for 1 <= k <= 6. By symmetry, using out-degrees in Y, the same holds when 1 <= n - k <= 6. When
n <= 13 one of the two holds. QED.

**Corollary.** For n <= 13, if f_3(n) > max H(n - 1), then every tournament of order n whose scores all lie in
[3, n - 4] is a first creator. Every non-first-creator becomes an ordinal sum after at most two arc reversals. Here
f_k(n) is the least H over n-tournaments with lambda >= k.

### 3.3 The connectivity floors (FINITE-EXACT; `lam.c`, `lam11.c`)

f_k(n) is the least H over n-tournaments with lambda >= k ([floors_4_10.txt](floors_4_10.txt)). The last two columns
count the non-first-creators with lambda = 1 and lambda = 2. Every non-strong class (lambda = 0) is a non-first-creator.

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
- Result: f_3(11) = 11,475 < max H(10) = 15,745.
- Below 15,745 the 3-arc-strong 11-tournaments take exactly seven values: 11,475, 14,787, 15,213, 15,273, 15,561,
  15,637 and 15,703.
- Four of them, 11,475, 15,213, 15,273 and 15,561, already lie in S(10). Those tournaments are not first creators,
  although every score lies in [3, 7].
- So the threshold k*(n), the largest lambda of a non-first-creator of order n, is 1 through order 8, 2 at orders 9
  and 10, and at least 3 at order 11. The generation description of section 3.1 needs one more arc at order 11
  (HYP-9532 (D)).

**Growth of f_3.** The ratios f_3(n + 1)/f_3(n) are 3.63, 3.31, 2.80 and 1.99 (n = 7..10). The ratios
max H(n)/max H(n - 1) are 4.2, 3.5, 5.1, 4.7 and 6.0 (n = 7..11). So the 3-arc-strong floor falls behind the
maximum, and the crossing occurs exactly between orders 10 and 11. No growth law for f_k is claimed.

**Two floors coincide with exceptional classes.** f_2(9) = 385 is the H of the single non-pendant REST class at order
9. f_1 is the floor that excludes 7 and 21.

## 4. The fibers of a value

### 4.1 Strong fibers are finite antichains (PROVED)

For every odd h, the strong tournaments with H = h have pairwise no induced embedding. Their orders lie in the strong
window W(h), which is contained in [3, B(h)], where B(h) = max{n : f(n) <= h} (klein; THM-1370 (3)).

*Proof.* Let A be an induced proper subtournament of a strong B. Choose a vertex v of B outside A. Then v is not a
singleton strong component of B, so THM-1900 (B) gives H(B - v) < H(B). By the insertion lemma, H(A) <= H(B - v). So
H(A) < H(B). QED.

So every tournament with H = h is an ordinal sum, with singletons inserted, of members of the finite sets
Strong(d) for the factorizations of h. The first creators are the members of Strong(h) at the least order.

### 4.2 Values born at their Busch bound

A value h is **born at its Busch bound** when n(h) = B(h), that is, at the last order where a strong realization is
possible. Among all odd h <= 95,095, these are exactly

    3, 5, 9, 11, 13, 17, 19, 23, 35, 39, 63

(klein's n(h) and the Moon-Busch floor). The two orphan values of THM-4922, 3 and 63, are on the list, and 63 is its
largest member. 7 and 21 are the values with no strong realization at any order: W(7) and W(21) are empty.

## 5. The objects of the 7/21 proofs and the generation of all tournaments (ANALOGY, with precise anchors)

The proofs that 7 and 21 never occur (THM-1370; the OCF route of THM-029/343 and THM-115) use five objects. Each one
turns out to be one of the operations that generate the non-first-creators. The middle column is precise; the
last column is an analogy and is typed as such.

| Object in the 7/21 proofs | Its role there | Its role in generating all tournaments (this packet) | Fundamental analogue (ANALOGY) |
| --- | --- | --- | --- |
| **Condensation**: H(T) is the product of the H of the strong components (THM-1370 (1)) | 7 is prime, so it needs a strong component with H = 7; 21 = 3 x 7 recurses into 7 | SUM: every non-strong tournament is an ordinal sum. The first creators are strong (finite-exact through order 11, HYP-9481 (b)) | **Unique factorization.** Tournaments under ordinal sum form the free monoid on strong tournaments, and H is a monoid homomorphism to (odd integers, x), like a norm. Strong means irreducible, and the condensation is the (ordered) prime factorization. 7 behaves like an inert prime: not a norm, and with no composite route. |
| **Insertion lemma**: H(T) >= H(T - v), equality iff v is a singleton component (THM-1370 (2), THM-1900 (B)) | the strong floor is nondecreasing | EXT: the extension calculus of section 1 (quadratic in the attachment). Fertility and the genealogy run along it. Strong fibers are antichains (section 4.1) | **Insertion sort.** Redei's existence proof of a Hamiltonian path is insertion sort with an intransitive comparator. H is a rank function, strictly monotone on the strong part. Equal-rank irreducibles form an antichain, as in Sperner theory. |
| **Moon's theorem**: a strong n-tournament (n >= 4) has at least two strong (n - 1)-subtournaments (THM-1370 (3)) | f(n) >= f(n - 1) | Every strong tournament is reached from C3 by strong one-vertex extensions ("Moon chains"). HYP-9512 (A) says the first creators have Moon chains inside the minimal set, with 63 the only break through order 11 | **Non-cut vertices and ear decompositions.** A connected graph has at least two non-cut vertices (the leaves of a spanning tree); a 2-connected graph has an ear decomposition (Whitney). |
| **Moon-Busch floor** f(n) = min{3^a 5^b : 2a + 3b = n - 1} (THM-1370 (4), Busch 2006) | 7 lies in (f(4), f(5)) = (5, 9) and 21 in (f(6), f(7)) = (15, 25), with no strong realization below. These are the crevasses | **Graded floors f_k(n)** (section 3.3). f_3(n) > max H(n - 1) forces every 3-arc-strong tournament to be a first creator through order 10, so every non-first-creator is at most two arc reversals from an ordinal sum. The k = 1 floor also makes strong fibers finite (section 4.1) and fixes the eleven values born at their Busch bound, 63 the largest (section 4.2) | **Minkowski and Hermite.** Fields of degree n have discriminant at least about c^n (Minkowski), so finitely many fields have bounded discriminant (Hermite). Likewise strong tournaments of order n have H >= f(n), about 5^(n/3), so finitely many strong tournaments have bounded H. The crevasses are band gaps between the strong bands [f(n), max H(n)]. The floor optimisation (3-cycles cost two vertices per factor 3, 5-blocks three vertices per factor 5, and 5^(1/3) < 3^(1/2)) is an optimal product-partition problem, like the "break n into 3s" maximum-product problem. |
| **OCF**: H = I(Omega(T), 2) >= 1 + 2 alpha_1, with Moon's pancyclicity counts (THM-029/343, THM-115) | H = 7 forces Omega = K3, and H = 21 forces at most 10 odd cycles. Strong unions then create more cycles | Omega of an ordinal sum is the disjoint union of the Omegas, and I is multiplicative: the same factorization one level down | **The hard-core lattice gas** at fugacity 2. Odd cycles are polymers, and vertex-disjointness is hard-core compatibility. Redei's theorem is the parity of the partition function. |

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

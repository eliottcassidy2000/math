# Kohl's Collatz group as a Tait-coloured graph, the sign fold, and the Kuratowski–Tutte triples

**Status: PROVED** (structure of Kohl's Tait-coloured Schreier graph: Kempe chains, Kempe quotients, the Z2×Z2 boundary, the sheet swap as a one-colour switch, the strategy square, the 3 | q criterion and its mod-p generalisation; exactness of the sign fold; independence of the three controls; infinite excluded-minor sets) **+ FINITE-EXACT** (Althöfer's union graph U_N: planar iff N ≤ 51, K5-minor-free iff N ≤ 67, planar double cover iff N ≤ 75, Petersen-family-minor-free iff N ≤ 91 (the absence side solver-certified), Petersen-minor-free iff N ≤ 103, each with explicit certificates; the switch census; the two-branch census) **+ CITED** (Kohl, J. Group Theory 2017, Prop. 1.2, statement and proof re-read this session from the atlas lane's author-PDF text; the graph-minor, flow, matroid, Colin de Verdière and planar-cover theorems from Wikipedia pages read 2026-09-25, primary papers unread) **+ ANALOGY** (the sign quotient versus dodecahedron → Petersen; most of the session triples). **Collatz remains OPEN.** No novelty is claimed for any individual elementary fact.

Session `collatz-procgen-20260922`, Kuratowski/Tait lane, 2026-09-25. Scripts `04-computation/experiments/procgen_kuratowski_20260925_{kohl_tait,union_minors,triples,excluded_minors}.py`, runner `procgen_kuratowski_20260925_run.sh`, output [procgen_kuratowski_20260925.out](procgen_kuratowski_20260925.out).

## 0. The answer in brief

1. **Kohl's group is a Tait-coloured Collatz graph, coloured by the mod-3 sign.** The Schreier graph Σ of G_C = ⟨a, b, c⟩ is exactly the undirected Collatz graph induced on Z \ 0(6) (Theorem K1). Colour a is the up-edge {n, 3n+1}; a doubling edge {m, 2m} is coloured b or c according to the sign of m in (Z/3)^× = {±1}. The colouring is proper because 2 ≡ −1 (mod 3), and negation swaps b and c. 5n+1 has a mod-5 analogue (§1.8).
2. **Kempe chains (Theorem K2).**
   * ⟨b,c⟩-chains are the orbits of the Banach contraction D(x) = 2x.
   * ⟨a,c⟩-chains are maximal segments of orbits of the other Banach contraction E(x) = (2x−1)/3 (THM-4471 §4). The chain started at y has length 2k with 2^−k = |y+1|₂.
   * ⟨a,b⟩-chains have at most 3 edges.
   * **The only closed Kempe chain in all of Z \ 0(6) is the digon {−1, −2}**, the Banach fixed point of E. The fixed point 0 of D is exactly what Kohl deletes.
3. **Kempe quotients (Theorem K3).** The three quotients are three classical reformulations:
   * the Syracuse graph;
   * a subcubic "sibling" graph whose edges include the owner's 4p+1 ladder (not found in the sources read; no priority claimed);
   * the graph of parity runs.
4. **Tutte (Theorem K4).** The colouring is a nowhere-zero Z2×Z2 chain whose boundary is c, a, a, 0, b on 1, 2, 3, 4, 5 (mod 6). It is a flow only at the branch points 4(6).
5. **Petersen obstructs nothing in a single sheet.** Every component of a functional graph has at most one cycle, so Σ is planar and free of K5, K3,3 and Petersen minors, *unconditionally*. Kuratowski graphs appear only when the two sheets are superimposed, in Althöfer's 3n±1 union graph U:

   | N | first appearance in U_N |
   |---|---|
   | 52 | K3,3 |
   | 68 | K5 |
   | 76 | no planar double cover |
   | 92 | Petersen-family minors |
   | 104 | Petersen |

   Colin de Verdière's μ jumps 3 → 4 at N = 52 and 4 → ≥5 at N = 92. Every kernel of U_N (all N ≤ 1000, and N = 2000, 4000) is bridgeless and 3-edge-colourable, so no snark occurs.
6. **The sign quotient is a fold, not a covering (Theorem F).** The natural double cover Γ₊ ∪ Γ₋ → U is *trivial*, because every move preserves sign. So "dodecahedron → Petersen" is an **ANALOGY**. What is exact: the fold attaches the other sheet's edge *precisely at Kohl's Tait defects* 1, 2, 3, 5 (mod 6).
7. **Controls as switches.**
   * Conjugating the one colour a by the class transposition h0 = τ_{2(6),4(6)} gives the 3n−1 group; this equals ν-conjugation followed by the relabelling b ↔ c.
   * The four mod-4 sign strategies of the 3n±1 game are four one-colour switches of G_C: Collatz (OPEN), 3n−1 (PROVED intransitive), "always one halving" (PROVED divergent) and "always ≥ 2 halvings" (PROVED transitive). **Collatz is the only open corner of this square.**
   * No rcwa switch gives DEFECT.
8. **Excluded minors (Theorem B, §4).**
   * SHEET (3n−1) and DRIFT (5n+1) are genuine excluded minors (incomparable atoms) of the rcwa order for connectivity.
   * DEFECT lies outside RCWA.
   * The three controls are pairwise independent; this is PROVED by a diagonal table of invariants: order, size with the multiplier, and finite description.
   * The rcwa order is not a well-quasi-order, so there is no finite Kuratowski list. For "every orbit reaches a finite set" the conjectural obstruction set {(q, ±1) : q ≥ 5} is an infinite antichain cut out by the *continuous* drift inequality log q > 2 log 2.
9. **Triples (§3).**
   * The exact Kuratowski–Tutte shape the session keeps meeting is "dual pair + self-dual", as in Tutte's regular matroids {F7, F7*, U2,4}: the sheets, the means and the 4-tournaments are REAL instances.
   * The "twins + container" shape (K5, K3,3 ⊂ Petersen) is REAL only for SHEET/DRIFT as twin atoms, with no container, and literally inside U. Everything else is ANALOGY or NUMEROLOGY.

## 1. Kohl's group as a Tait-coloured Schreier graph

### 1.1 Source

Kohl, "The Collatz conjecture in a group theoretic context", Prop. 1.2(a) (proved as Prop. 2.1). I re-read the statement and the proof table in the author-PDF text downloaded by the atlas lane (`scratch/procgen_atlas/lit_rcwa/dl/kohl_collatzgroups.txt`, git-ignored).

> G_C = ⟨τ_{1(2),4(6)}, τ_{1(3),2(6)}, τ_{2(3),4(6)}⟩ acts transitively on N \ 0(6) iff the Collatz conjecture holds.

* **Kohl's proof.** For n ∉ 0(3): {n, n^a, n^b, n^c} = {n} ∪ {C(n)} ∪ C^{−1}(n). Also C^{−1}(0(3)) = 0(6), and a = C on 3(6).
* **The journal reference.** "J. Group Theory 20 (2017)" is as recorded by the atlas lane. The text file carries no journal header, so volume and year are UNVERIFIED here.
* **Prop. 1.2(b), G_T.** Also read; see §1.9.

### 1.2 The graph (Theorem K1, PROVED; checked on |n| ≤ 60000)

The generators:
* a: n ↔ 3n+1 for n odd;
* b: m ↔ 2m for m ≡ 1 (3);
* c: m ↔ 2m for m ≡ 2 (3).

Each fixed point carries a semi-edge, so the Schreier graph Σ is a cubic graph with semi-edges in which every colour class is a perfect matching: a Tait colouring. The vertex types:

| n mod 6 | moved by | degree | semi-edges | Collatz meaning |
|---|---|---|---|---|
| 0 | — | 0 | a b c | deleted (the D-ray above an odd multiple of 3) |
| 1 | a, b | 2 | c | odd, one C-preimage |
| 2 | b, c | 2 | a | even, one C-preimage |
| 3 | a | 1 | b c | leaf (odd multiple of 3) |
| 4 | a, b, c | 3 | — | **branch point**: two C-preimages 2n and (n−1)/3 |
| 5 | a, c | 2 | b | odd, one C-preimage |

**Theorem K1.**
* The non-loop edges of Σ are exactly the undirected C-graph (C(n) = n/2 or 3n+1) induced on Z \ 0(6).
* The only multiple edge is the C-2-cycle {−1, −2}, coloured a and c.
* C maps N \ 0(6) into itself, and C^{−1}(0(3)) = 0(6).

*Proof.* The edges of Σ are {n, 3n+1} for n odd and {m, 2m} for 3 ∤ m. An edge {m, 2m} meets 0(6) iff 3 | m. ∎

Why 0(6) must go: along a doubling ray x·2^j with 3 ∤ x, the mod-3 residue alternates, since 2 ≡ −1. That is what lets two class transpositions 2-colour the doubling edges. A ray of multiples of 3 has no mod-3 sign; colouring it would need the parity of v₂, which is not a residue condition. Deleting 0(6) is harmless because multiples of 3 are transient.

**Corollary K1′ (PROVED).** Every component of an undirected functional graph has at most one cycle. Hence Σ is planar, and has no K3,3, K5 or Petersen minor, *whether or not Collatz holds*.
* **On N \ 0(6).** The unique cycle is the **rainbow triangle** 1 –a– 4 –c– 2 –b– 1 iff Collatz holds. All n ≤ 2·10⁵ reach it.
* **On Z_{<0} \ 0(6).** Starts down to −2·10⁵ meet exactly three cycles:

  | cycle | colour word |
  |---|---|
  | digon {−1, −2} | `ac` |
  | 5-cycle through −5 | `acacb` |
  | 18-cycle through −17 | `acacacacbacacacbcb` (7 a, 8 c, 3 b) |

### 1.3 Kempe chains (Theorem K2, PROVED; every vertex with |v| ≤ 60000 checked against the closed forms)

* **⟨b,c⟩.** For v = 2^j o with o odd:
  * if 3 ∤ o, the chain is the ray {o·2^i : i ≥ 0}; its colours alternate, starting with b if o ≡ 1 (3) and with c if o ≡ 2 (3);
  * if 3 | o, it is the singleton {o}.
* **⟨a,b⟩.** P_n = (2n –b–) n –a– 3n+1 –b– 6n+2 for n odd, with 2n present iff n ≡ 1 (6).
  * It has 3 or 2 edges and exactly one a-edge, and the P_n partition Z \ 0(6).
  * The ends are a-fixed (they lie in 2(6)).
* **⟨a,c⟩.** For y odd with y ≡ 1, 3 (6) and k = v₂(y+1), the chain is the odd run y –a– 3y+1 –c– T(y) –a– … –c– T^k(y).
  * It has 2k edges, and T^j(y) = 3^j(y+1)/2^j − 1.
  * The interior odd points lie in 5(6), and the end lies in 2(6).
  * The only other ⟨a,c⟩-orbit is {−1, −2}.
* **Unique closed chain.** On Z \ 0(6) the only closed Kempe chain is the ⟨a,c⟩-digon {−1, −2}. On N \ 0(6) every Kempe chain is a path or a ray.

*Proof.*
* **⟨b,c⟩.** b ∪ c consists of the doubling edges {m, 2m} with 3 ∤ m, and 2m ≡ −m (mod 3).
* **⟨a,b⟩.** a sends odd n to 3n+1 ∈ 4(6) ⊂ 1(3). Its b-partner 6n+2 lies in 2(6), which a fixes.
* **⟨a,c⟩, forwards.** The steps a then c give x → 3x+1 → (3x+1)/2 = T(x). The chain continues while T(x) is odd.
* **⟨a,c⟩, backwards.** The steps c then a give E(x) = (2x−1)/3. The chain stops at the first value ≡ 1 or 3 (mod 6), which c fixes.
* **Finiteness.** Forward runs are finite unless y = −1. Backward runs strictly decrease |x| unless x = −1. ∎

**The Banach reading (PROVED).**
* **Two contractions.** Walking a ⟨b,c⟩-chain iterates D(x) = 2x. Two backward steps along an ⟨a,c⟩-chain apply E(x) = (2x−1)/3. D and E are the inverse branches of T, and both are 2-adic contractions with factor 1/2 (THM-4471 §4).
* **Their fixed points.** D fixes 0, which lies in 0(6) and is deleted by Kohl. E fixes −1, which is the unique closed Kempe chain.
* **Chain length is 2-adic distance.** The ⟨a,c⟩-chain from y has length 2k, where 2^−k = |y+1|₂. So Kempe length measures 2-adic proximity to the Banach fixed point −1.
* **Mersenne starts.** y = 2^k − 1 opens the chain of length 2k ending at 3^k − 1 (checked k ≤ 11; THM-4473's T^k(2^k−1) = 3^k − 1).
* **One neighbourhood, four views.** The following all live in the same 2-adic neighbourhood of −1, now read as long Kempe chains:
  * THM-4471(D): stretching by 3/2 at Mersenne starts;
  * THM-4470(5): the pair-0 obstruction to bounded lookahead;
  * Applegate–Lagarias: "−1 resists".

### 1.4 Kempe quotients (Theorem K3, PROVED)

Contracting every chain of one type (each chain is connected) preserves and reflects connectivity. Since the chains are trees, apart from the digon, it also preserves the cyclomatic number. So **Collatz ⟺ each quotient restricted to N is connected**, and under Collatz each quotient is a tree plus one loop.

* **Σ/⟨b,c⟩ = the undirected Syracuse graph** on odd integers: edges {n, S(n)} with S(n) = oddpart(3n+1), and a loop at 1. Syracuse acceleration *is* this Kempe contraction.
* **Σ/⟨a,b⟩ = a subcubic "sibling" graph** on odd n (checked for |n| ≤ 7500). Its edges are:
  * {n, 4n+1} for every odd n;
  * {n, T(n)} for n ≡ 3 (4);
  * {n, T²(n)} for n ≡ 1 (8).

  The degree is 3 if 3 ∤ n and 2 if n ≡ 3 (6), and P₁ carries the loop {2, 4}. The owner's sibling ladder p ↦ 4p+1 (E D² = S E, mod-192 note) is an edge family of this quotient.
* **Σ/⟨a,c⟩ = the run graph.** Its vertices are the maximal odd runs, the Terras blocks 1^k 0. The degree is 1 + [y ≡ 1 (6)] + v₂(y+1).

### 1.5 Tutte: the colouring as a Z2×Z2 chain (Theorem K4, PROVED)

Put a = (1,0), b = (0,1), c = (1,1). The colouring φ is a nowhere-zero Z2²-chain. Its boundary ∂φ(v) = Σ_{e∋v} φ(e) is the sum of the missing colours at v:

| v mod 6 | 1 | 2 | 3 | 4 | 5 |
|---|---|---|---|---|---|
| ∂φ(v) | c | a | a | 0 | b |

So φ is a flow exactly at the branch points 4(6).
* **Where a flow survives.** A Z2-projection of φ is a flow on a subgraph iff that subgraph is a union of closed Kempe chains. By K2, the only such subgraph is the digon at −1.
* **Tutte's hypothesis fails everywhere else.** Tutte's theory needs a bridgeless graph: no finite graph with a bridge has a nowhere-zero flow (CITED). Under Collatz, Σ on N \ 0(6) is a tree plus a triangle, so every edge except three is a bridge. The Tait colouring is therefore "all boundary". This is a precise statement, and it gives no leverage.

### 1.6 What the Petersen graph obstructs here (d)

* **The general fact.** A finite connected cubic graph is the Schreier graph of a transitive action of Z2∗Z2∗Z2 (generators = colours) iff it is 3-edge-colourable. So the Petersen graph, the smallest snark, can never be an orbit of any group generated by three involutions.
* **The snark theorem.** It says every snark has a Petersen minor (CITED: announced by Robertson–Sanders–Seymour–Thomas; steps published 2016 and 2019; a 2026 apex preprint, arXiv 2608.22870, claims to complete it; primary texts unread).
* **Why this does not bite.** Σ is such an orbit graph *by construction*, and by K1′ it has no Kuratowski or Petersen minor at all. Collatz asks for connectivity, which is neither a colouring nor a minor-monotone property.
* **Verdict: Petersen obstructs nothing in a single sheet (PROVED).** Where it does appear (§2), it still obstructs nothing: the kernels there are 3-edge-colourable.

### 1.7 Switches (e)

* **Kempe switches are invisible (PROVED).** Swapping x ↔ y on a union of ⟨x,y⟩-chains does not change the graph, so transitivity is unchanged. For example, swapping a ↔ b on all P_n with n ≡ 1 (4) gives new rcwa involutions a′, b′; a′ has modulus 24 and differs from a at 50000 window points. The Schreier graph is identical (139999 edges checked).
* **The sheet swap is a one-colour switch (PROVED).** νaν = τ_{1(2),2(6)} (n ↔ 3n−1) = a^{h0} with h0 = τ_{2(6),4(6)}, while νbν = c and νcν = b. Hence ν G_C ν = ⟨a^{h0}, c, b⟩ = G_C⁻, Kohl's group for 3n−1.
  * The sheet swap = one-colour conjugation of a by the "±2" transposition h0, followed by the global relabelling b ↔ c, which fixes the uncoloured graph.
  * G_C⁻ has at least 3 orbits on N \ 0(6) (PROVED: the three cycles through 1, 5, 17).
* **The one-colour switch census (FINITE-EXACT).** Take h = τ_{c1,c2} with c1, c2 ⊆ 2(6) ∪ 4(6). Then the Schreier graph of ⟨a^h, b, c⟩ is exactly the functional graph of F_h (odd n ↦ h(3n+1), even n ↦ n/2) on Z \ 0(6). For moduli 6 and 12 (10 switches):

  | h | odd step | type |
  |---|---|---|
  | τ_{2(6),4(6)} | 3n−1 | **SHEET** (cycles 1, 5, 17; PROVED intransitive) |
  | τ_{2(12),10(12)}, τ_{4(12),8(12)}, τ_{4(12),10(12)} | mixed offsets (3n−7 / 3n+1, …) | SHEET-type (3–4 cycles; PROVED intransitive) |
  | τ_{4(6),2(12)} | 6n−4 | **DRIFT**, PROVED divergent: one-step ascent, the next odd value 3n−2 > n |
  | τ_{2(12),4(12)} | 3n−1 (n ≡ 1 mod 4), 3n+1 (n ≡ 3 mod 4) | **DRIFT**, PROVED divergent: exactly one halving per odd step, next odd value (3n∓1)/2 > n |
  | τ_{8(12),10(12)}, τ_{2(6),10(12)} | 3n∓1 with 4 ∣ 3n∓1; resp. (3n−5)/2, 3n+1 | **PROVED transitive**: odd-step descent, next odd value ≤ (3n+1)/4 < n |
  | τ_{2(6),4(12)}, τ_{4(6),8(12)} | (3n+1)/2 or 3n+1; 6n+2 | Collatz-equivalent re-timings of C (F ∈ {C, C²} up to inserted points) |

  The certificates are exact. Intransitivity comes from two cycles, or from a divergent orbit together with the cycle at 1. Transitivity comes from descent below n0 = 1 plus a direct check.

  For moduli ≤ 24 (114 switches):

  | type | count |
  |---|---|
  | SHEET-type | 48 |
  | DRIFT, PROVED divergent | 8 |
  | DRIFT, heuristic | 16 |
  | PROVED transitive | 3 |
  | open (Collatz-like) | 37 |
  | mixed | 2 |
* **The strategy square (PROVED).** Choose the sign in 3n±1 by n mod 4, written (sign used at n ≡ 1, sign used at n ≡ 3): (+,+), (−,−), (−,+), (+,−). These four residue strategies of Althöfer's 3n±1 game are exactly four one-colour switches of G_C:

  | strategy | switch | status |
  |---|---|---|
  | (+,+) | identity | Collatz, OPEN |
  | (−,−) | τ_{2(6),4(6)} | 3n−1, PROVED intransitive |
  | (−,+) | τ_{2(12),4(12)} | PROVED divergent, hence intransitive |
  | (+,−) | τ_{8(12),10(12)} | PROVED transitive |

  In the Althöfer lane's notation ([game note](collatz_procgen_20260922_althofer_game.md), §3), (−,+) is the all-u strategy, whose orbits are that note's ascending rays (L1), and (+,−) is the all-d strategy.

  **Collatz is one corner of a square whose other three corners are proved.** The difficulty is not in the Tait/switch structure.
  * The two pure corners are trivial: (−,+) has one halving per odd step, and (+,−) has at least two.
  * The two mixed corners (+,+) and (−,−) have a halving count per odd step that depends on n mod 4, with the same Haar log-drift log(3/4) per odd step on both.
  * Of the mixed corners, (−,−) is decided only by its explicit extra cycles, which is the SHEET phenomenon; (+,+) is Collatz.
* **DEFECT is unreachable (PROVED).** A nonidentity rcwa permutation has support of positive density.
* **Reading.** The Tait/Kempe/switch framework is not where the difficulty sits: *provably* transitive Kohl-shape groups lie one class-transposition switch away from G_C.

### 1.8 Why q = 3, and the mod-p generalisation (PROVED)

Take Kohl's mod-3 shape ⟨a_q = τ_{1(2),(q+1)(2q)}, b, c⟩ for the map qn+1. It deletes the doubling edges at multiples of 3. This is harmless iff multiples of 3 are transient, i.e. iff qn+1 ≢ 0 (mod 3) for every odd n, i.e. iff **3 | q**.
* **q = 5 fails.** 36 = 5·7+1 becomes a leaf and 18 is isolated, so the C₅-path 7 → 36 → 18 → 9 → 46 is cut.
* **Two ingredients.** Kohl's presentation combines two facts: 2 ≡ −1 (mod 3), and 3 | q.
* **General criterion.** Colour a doubling edge {m, 2m} by the parity of m's position on its ⟨2⟩-coset cycle in (Z/p)^×, and delete the multiples of p. This gives a proper Tait colouring iff **−1 is a power of 2 mod p** (every doubling cycle has even length). The shape tracks C_q iff, in addition, **p | q**, which makes multiples of p transient.
  * The first condition holds for p = 3, 5, 11, 13, 17, 19, 29, 37, 41, 43 and fails for p = 7, 23, 31, 47 (checked).
  * 5n+1 has a mod-5 Tait group ⟨τ_{1(2),6(10)}, b₅, c₅⟩, verified to have Schreier graph equal to the C₅-graph on Z \ 0(10). It is intransitive (cycles through 1, 13, 17).
* **Consequence.** The Tait/Kohl framework is DRIFT-blind as well as SHEET-blind. The mod-3 sign of §0 is the p = 3 case: "2 ≡ −1" generalises to "−1 ∈ ⟨2⟩".

### 1.9 Kohl's second group G_T (FINITE-EXACT)

G_T = ⟨τ_{0(2),1(2)}, τ_{1(2),2(4)}, τ_{1(4),2(6)}⟩ acts on N₀. On [0, 20000] its Schreier graph has cyclomatic number 1: the only cycle is the b/c digon {1, 2}. Colour a is the consecutive pairing {2k, 2k+1}, which is THM-4470's 3n−1 pairing. Kohl's conjugating map r(n) ∈ {2n−2, 2n−1} lands in the a-pair of index n−1. So G_T is another Tait-coloured tree plus a digon, and is equally minor-poor.

## 2. The sign fold and Althöfer's union graph

U is the simple graph on N with the halving edges {m, 2m} and, for odd n, both up-edges {n, 3n+1} and {n, 3n−1}: the move graph of the 3n±1 game. U_N is its restriction to [1, N].

### 2.1 Exactness (Theorem F, PROVED; checked on |x| ≤ 30000)

Let Γ₊ be the undirected 3n+1 graph on Z \ {0}, and π(x) = |x|.
* **(i) Homomorphism.** π is a graph homomorphism Γ₊ → U. Each halving edge has two preimages, one per sign; each up+ edge has one (in Z_{>0}); each up− edge has one (in Z_{<0}). By the transport theorem, ν C₊ ν = C₋.
* **(ii) A fold, not a covering.** π is locally injective but not locally surjective: deg_U(|x|) − deg_Γ₊(x) = 1 exactly for x ≡ 1, 2, 3, 5 (mod 6), and 0 for x ≡ 0, 4 (mod 6), with the only exceptions x = 1, 2 (the digon). So π is an immersion (a fold). **The deficient vertices are exactly Kohl's Tait defects**: the fold attaches the other sheet's up-edge precisely where the Tait colouring has a boundary.
* **(iii) The natural cover is trivial.** W = Γ₊ ∪ Γ₋ → U is a 2-fold covering with deck involution ν. But no edge of W joins x > 0 to y < 0, because every map x/2, 3x±1 preserves sign. So W = U ⊔ ν(U) is the **trivial** double cover.
* **The control.** The dodecahedron → dodecahedron/antipode = Petersen is a *connected* double covering: a planar cover of a projective-planar, nonplanar graph. It is checked in the script.
* **Verdict: ANALOGY.** The Collatz involution ν does not act on one connected graph; it exchanges the two sheets. The exact residue is (ii).

### 2.2 The union graph (FINITE-EXACT)

* **Degrees.** U is subcubic. Odd n has neighbours 2n, 3n+1 and 3n−1, all larger than n. Even n has neighbours n/2, 2n and whichever of (n∓1)/3 is odd. The degree is 3 except at 0 (mod 6) and at the simple-graph digon points 1 and 2.
* **Short cycles.** The girth is 3. The only triangle is {1, 2, 4} (checked in U₂₀₀₀). The cycles of length 3–8 in U₆₀₀ number 1, 0, 2, 0, 2, 33.

All the properties below are monotone in N, because U_N ⊂ U_{N+1}:

| property of U_N | holds iff | certificate |
|---|---|---|
| planar | N ≤ 51 | rotation system of U₅₁ verified by face tracing (V−E+F = 2 on all 15 components); explicit K3,3 subdivision in U₅₂ (below) |
| no K5 minor | N ≤ 67 | K5 branch sets in U₆₈ verified; kernel(U₆₇) is cubic on 10 vertices, where a K5 model needs a perfect matching M with kernel/M = K5, and all 6 perfect matchings fail; CP-SAT agrees |
| planar double cover (⟺ projective-planar, Negami 1986, CITED) | N ≤ 75 | a planar signing for each kernel N ≤ 75; all 256 signings of kernel(U₇₆) nonplanar, so U₇₆ is not projective-planar (easy direction) |
| no Petersen-family minor (linkless; Colin de Verdière μ ≤ 4) | N ≤ 91 | minors of the 6 non-Petersen family members verified in U₉₂; at N = 91, K6 is excluded by a count (a branch set of s vertices in a cubic graph has ≤ s+2 exits) and the others by CP-SAT (solver-certified) |
| no Petersen minor | N ≤ 103 | explicit Petersen subdivision in U₁₀₄ (below); kernel(U₁₀₃) (18 vertices) has none, by exhaustive subdivision search and CP-SAT |

**The first K3,3 (U₅₂).** The shores are {5, 7, 26} and {11, 14, 20}, and every edge is an arithmetic move:
```
5-10-20 [halvings]      5-14 [3*5-1]           5-16-32-11 [3*5+1; 32/2; 3*11-1]
7-14 [14/2]             7-20 [3*7-1]           7-22-11 [3*7+1; 22/2]
11-34-17-52-26 [3*11+1; 34/2; 3*17+1; 52/2]    14-28-9-26 [28/2; 3*9+1; 3*9-1]
20-40-13-26 [40/2; 3*13+1; 26/2]
```
* **What completes it.** Vertex 52 closes the obstruction through the ear 17 –(3·17+1)– 52 –(52/2)– 26.
* **Kernels.** kernel(U₅₁) is the prism, i.e. K4 with a truncated vertex. kernel(U₅₂) is K3,3 with a truncated vertex.
* **Sheet interaction.** Both kernels carry the 3n−1 5-cycle 5–14–7–20–10. So the first nonplanarity is sheet interaction: a 3n−1 cycle with 3n+1 ears.

**The first Petersen (U₁₀₄).** The branch vertices are {10, 13, 19, 20, 23, 26, 28, 34, 38, 44}. The 15 paths are listed with their moves in the `.out` file (section U6); the model is deterministic, with CP-SAT on 1 worker. For example, 34–11–32–16–8–3–10 [3·11+1; 3·11−1; halvings; 3·3−1; 3·3+1]. Vertex 104 enters through 35–104–52 [3·35−1; 104/2].

**Colin de Verdière (CITED equivalences + the table).** μ(U_N) ≤ 3 iff N ≤ 51; μ(U_N) = 4 for 52 ≤ N ≤ 91; and μ(U_N) ≥ 5 for N ≥ 92. The equivalences used: μ ≤ 3 ⟺ planar; μ ≤ 4 ⟺ linklessly embeddable ⟺ no Petersen-family minor.

**Kernels are Tait-colourable (FINITE-EXACT, SAT).** For every N with 52 ≤ N ≤ 1000, and for N = 2000 and 4000, the kernel of U_N is connected, bridgeless, cubic and 3-edge-colourable. So the Petersen minors of U obstruct nothing.

**Planar double covers are real but not arithmetic.** The Petersen-type phenomenon, a nonplanar graph with a planar double cover, occurs exactly for 52 ≤ N ≤ 75. But the signing by sheet and every edge-type signing (up− negative, up+ negative, halving negative) give nonplanar covers from N = 52 or earlier. The planar covers are not the negation cover.

**Reading.** The Kuratowski–Tutte graphs measure *sheet interaction*. They are absent from each sheet (K1′) and appear only in the superposition U. They arrive in the minor-forced order K3,3 (52) < K5 (68) < Petersen (104), so the "container" comes last, as Kuratowski–Tutte requires; this order is automatic.

## 3. The triple dictionary

The Kuratowski–Tutte (KT) pattern has two exact shapes:
* **KT-a, "twins + container".** X, Y are incomparable minimal obstructions, and Z ≥ X, Y is the minimal obstruction for a weaker property. Example: K5, K3,3 < Petersen for planarity versus 3-edge-colouring or 4-flows; "no Petersen minor ⟹ 4-flow" extends Tait's form of the four-colour theorem.
* **KT-b, "dual pair + self-dual".** A duality D satisfies D(X) = Y and D(Z) = Z. Example: Tutte's regular matroids, regular ⟺ no U2,4, F7, F7* minor, with F7 ↔ F7* dual and U2,4 self-dual. Graphic ⟺ additionally no M*(K5) or M*(K3,3).

Neither shape covers K5 versus K3,3 exactly: no duality swaps them, and they are one ΔY move apart up to an edge (ΔY(K5) = K3,3 + e, checked). ΔY is Kennelly's 1899 star–triangle transform, which preserves effective resistance: an electrical, continuous move. The ΔY/YΔ class of K6 has exactly 7 graphs and contains the Petersen graph, i.e. it is the Petersen family (computed).

| triple | twins | exact duality swapping them? | third object and order | shape | verdict |
|---|---|---|---|---|---|
| K5, K3,3 / Petersen | K5, K3,3 | none (ΔY up to an edge) | Petersen ≥ both (minors) | KT-a | REAL (reference) |
| F7, F7* / U2,4 | F7, F7* | matroid duality | U2,4 self-dual, not ≥ F7 | KT-b | REAL (reference) |
| μ ≤ 3 / μ ≤ 4 | K5, K3,3 (μ = 4) | — | Petersen family (μ ≥ 5): the next level of one spectral invariant | KT-a as consecutive levels | REAL (CITED) |
| sheets b ∈ {−1, 0, +1} | ±1 | **ν, exact** (transport theorem) | b = 0 (Mahler's 3x/2) is ν-fixed; it is the common scaling limit, *below* both | KT-b | **REAL** (KT-b); KT-a reversed |
| means GM, AM, QM | GM², QM² | **s ↦ 2AM² − s, exact** (QM² + GM² = 2AM², the parallelogram law) | AM fixed, between the twins | KT-b | **REAL** (KT-b); ANALOGY to the sheets via the cross term ±2ab (THM-4471) |
| 4-tournaments | the two diamonds | **converse, exact** (THM-4472: time reversal) | TT and strong are both self-converse; H = 1, 3, 3, 5 | KT-b with two fixed points | **REAL** (KT-b) |
| Kohl colours a, b, c | b, c | ν swaps b ↔ c (the mod-3 sign) but sends a to a^{h0} on the other sheet | a (up), no order | KT-b across {G_C, G_C⁻} only | REAL duality, ANALOGY as a triple |
| Kempe types | ⟨b,c⟩ (D, → 0), ⟨a,c⟩ (E, → −1) | none | ⟨a,b⟩ bounded | "two contractions + one bounded" | ANALOGY (the Banach reading is REAL) |
| controls SHEET / DRIFT / DEFECT | SHEET, DRIFT: incomparable **atoms** of the rcwa order, both non-Conn (PROVED) | none: order data versus size data | DEFECT lies outside RCWA; it realises both failure modes (a flip gives a cycle, density-zero flips give divergence), but is ≥ neither | KT-a twins without a container | ANALOGY (REAL atoms and independence, §4) |
| inputs: sign law, gate integrality, null-set avoidance | — | none | a chain of necessities | chain | NUMEROLOGY for KT |
| places {∞, 2, 3} | 2, 3 | no automorphism swaps places; time reversal swaps only the reading roles (parity forwards, mod 3 backwards) | ∞ (order and sign live only there; the product formula fixes \|·\|∞ on {2,3}-units) | none | ANALOGY; REAL sub-fact: the Tait colours are the 3-adic sign, flipped by 2 ≡ −1 and by ν |
| swap words Sturmian / square / cube | Sturmian (Theorem S, approximation), square (Theorem Y, Padé) | none | the cube Y3 (HYP-9127) escapes both methods | "minimal object outside the union of two methods" | ANALOGY: the "snark of current methods" is an apt metaphor with no order |
| negative cycles −1, −5, −17 | — | none | clocks 1/1, 3/2, 11/7 lie on the Stern–Brocot path to log₂3 at depths 0, 2, 5 (a chain); 11/7 = 3/2 ⊕ 8/5 involves the upper node 8/5; the −17 word 11110111000 is unbalanced, not the Christoffel word 01011011011 | chain | NUMEROLOGY |
| Banach, Brouwer, Sharkovskii | — | none | Caristi ⟹ Banach is real; Sharkovskii was refuted as a sheet separator | none | NUMEROLOGY (but D, E ↔ Kempe chains is REAL, §1.3) |
| U-thresholds K3,3 / K5 / Petersen | K3,3 (52), K5 (68) | — | Petersen (104) ≥ both | KT-a literally | REAL but forced by minor monotonicity |

**Summary.**
* **KT-b is REAL and recurrent.** "A Z/2 duality with a fixed point" appears in the sheets, the means and the 4-tournaments, as in Tutte's regular matroids. In all three the duality is a sign: the sign of b, the sign of the cross term ±2ab, and the direction of time.
* **KT-a is REAL only in part.** It holds only for the twin atoms SHEET/DRIFT, which have no container, and literally inside the union graph U.

## 4. An excluded-minor calculus for Collatz proof methods

### 4.1 Systems, minors, controls

* **Systems.** A system is (f, X) with f : X → X and X ⊆ Z \ {0}.
  * **Conn(f, X):** the undirected functional graph is connected; Collatz = Conn(T, N).
  * **P_fin:** every orbit enters a fixed finite set.
* **Minor operations.**
  * (m1) restriction to an f-invariant subset;
  * (m2) conjugation by an affine bijection, including x ↦ x/d on dZ and ν;
  * (m3) acceleration by contracting connected orbit pieces. The Syracuse acceleration is the ⟨b,c⟩ Kempe contraction (K3).
* **Lemma.** Conn and P_fin pass to invariant restrictions and to conjugates. Conn is preserved and reflected by (m3). *Proof.* Two points of an invariant Y whose orbits meet in X meet inside Y. ∎
* **Control moves** (not minor operations):
  * σ: change the side (the same map on −X; equivalently ν on the map);
  * δ: change the multiplier inside the 2-adic conjugacy class (every T_q with q odd is conjugate to the shift on Z₂ by its parity-vector map; Terras/Lagarias, and Akin 2004 as recorded by the atlas lane [R]);
  * ε: modification on a density-zero set.

### 4.2 Theorem B (Kuratowski-type barrier, PROVED)

Let (Φ, A) be a sound method: Φ is any observation of systems, and Φ(S) ∈ A ⟹ Conn(S). Suppose it proves C = (T₃,₁, N). Then Φ separates C from each of the three controls:
* σC = (T₃,₋₁, N);
* δC = (T₅,₁, N);
* εC = THM-4470's pair-flip maps.

In particular, a method that proves Collatz is invariant under none of σ, δ, ε.

*Proof.* Each control is non-Conn:
* σC has cycles through 1, 5 and 17;
* δC has cycles through 1, 13 and 17;
* εC with the single flip of pair 4 has the new cycle (3 5 8 12 6), and with the density-zero flips of THM-4470(4) the orbit of 3 diverges (300 steps re-verified).

Soundness plus Φ(control) = Φ(C) ∈ A would force Conn(control). ∎

This is the Kuratowski *necessity* direction. The converse, "separating all three controls suffices", is not a theorem. The HARD words (cube-swap Y3, HYP-9127) resist the current methods even though those methods do see order, size and pointwise data (Theorems S, D, Y). Sufficiency is the open problem itself, and Y3 is not a minor-type obstruction.

### 4.3 Minimality and independence (PROVED)

* **Atoms.** In the two-branch rcwa order (§4.4), (3, −1) on N and (5, 1) have no proper minors and are non-Conn. So SHEET and DRIFT are literal excluded minors of Conn. They are incomparable: the sheet sign b·side differs, and the multiplier is a minor invariant.
* **DEFECT is outside RCWA.** An rcwa map that agrees with T on a density-one set agrees with T on every residue class of a common modulus, hence everywhere. So εC is not a minor of any rcwa system. This is also why no class-transposition switch reaches it (§1.7).
* **Independence (diagonal table).** Each invariant separates C from exactly one control:

  | invariant | separates σ | separates δ | separates ε |
  |---|---|---|---|
  | drift sign (log q − 2 log 2) | no | **yes** | no |
  | sheet sign b·side | **yes** | no | no |
  | rcwa (finite description) | no | no | **yes** |

  So separating one control never forces separating another.
* **What a proof must see.** It must use, simultaneously:
  * the archimedean **order** (for σ);
  * **size with the true multiplier** (for δ; 2-adic data cannot help, since all T_q are 2-adically conjugate);
  * **pointwise, finite-description** data (for ε; density data cannot help).
* **How the classical families fare.**

  | method family | blind to |
  |---|---|
  | residue, Haar and density arguments | σ and ε |
  | 2-adic dynamics | σ, δ and ε |
  | contraction in \|x−y\| (THM-4471) | σ |
  | pairing statistics (THM-4470) | ε |
  | Kohl's wildness | σ and δ |
  | Kohl's Tait/Schreier framework (this note) | σ and δ (3n−1 through ν, 5n+1 through its mod-5 Tait group); it separates ε only trivially, because DEFECT maps have no rcwa presentation |

### 4.4 Is there a finite list of excluded minors? (PROVED / FINITE-EXACT / OPEN)

* **The order.** Two-branch maps T_{q,r}(x) = x/2 (x even), (qx + r)/2 (x odd), with q, r odd. The rcwa order is (q, r′) ≼ (q, r) iff r′ | r, since T_{q,r}(dx) = d·T_{q,r/d}(x). This is verified, and 123 scaled cycles are checked.
* **Not a well-quasi-order.** The multiplier q is a minor invariant, so the q-slices are pairwise incomparable. {(q, 1)} is an **infinite antichain**, so no Robertson–Seymour-type finiteness is available, and it fails.
* **Conn has infinitely many excluded minors (PROVED).** For every prime p ∤ 2q, both pN and N \ pN are invariant, so (q, p) is never Conn. Every q-slice therefore contains an excluded minor.
  * In the slice q = 3: if Collatz holds, the excluded minors are {(3, p) : p ≥ 5 prime}; if Collatz fails, (3, 1) alone is. So "finite list" ⟺ "Collatz false". This is a restatement, not a reduction.
* **P_fin (census: q ≤ 15, r ≤ 25, starts ≤ 3000, escape bound 10⁴⁰).**
  * q = 1: every orbit enters a cycle (PROVED: (x+r)/2 < x for x > r).
  * q = 3: no escape (Lagarias's 3x+k conjectures; OPEN).
  * q ≥ 5: mass escape.
  * Conjecturally the excluded minors of P_fin are exactly {(q, ±1) : q ≥ 5}. This infinite antichain is cut out not by finitely many minors but by the **continuous** inequality log q > 2 log 2. No member is PROVED: divergence is open for every q ≥ 5.

## 5. Bridges between discrete and continuous (the owner's second request)

1. **Kempe length is 2-adic distance (REAL).**
   * ℓ_{ac}(y) = −2 log₂|y+1|₂.
   * The unique closed Kempe chain is the Banach fixed point −1 of E.
   * The deleted vertex 0 is the Banach fixed point of D.
   * So the combinatorial Kempe decomposition of the Tait-coloured Schreier graph *is* the orbit decomposition of the two contractions of Z₂ in THM-4471.
2. **A spectral invariant sees the Kuratowski levels (REAL for KT, CITED; FINITE-EXACT for U).** Colin de Verdière's μ is the maximal corank of a generalized Laplacian (Schrödinger operator) with one negative eigenvalue and the Strong Arnold Property.
   * μ ≤ 3 ⟺ planar (no K5, K3,3 minor).
   * μ ≤ 4 ⟺ linkless (no Petersen-family minor).
   * The twins and the container are *consecutive levels* of one continuous-origin invariant.
   * For Althöfer's graph: μ(U_N) = 4 exactly for 52 ≤ N ≤ 91, and ≥ 5 from 92.
3. **Drift replaces finite obstruction sets (PROVED antichain; ANALOGY).**
   * In minor-closed graph classes, obstruction sets are finite.
   * In the rcwa order they are infinite antichains, and the natural substitute is a real-valued invariant, the drift log q − 2 log 2.
   * Its sign changes between q = 3 and q = 5 (zero at q = 4). THM-4470 singles out q = 3 as the unique AM-fair odd multiplier, whose drift is the AM–GM gap log(√3/2) < 0.
4. **Star–triangle is electrical (REAL for KT).**
   * ΔY(K5) = K3,3 + e, and the ΔY/YΔ class of K6 is the Petersen family.
   * Nothing analogous was found on the Collatz side; the only triangle of U is {1, 2, 4}.
5. **The fold fills the Tait defects (REAL, Theorem F(ii)).** The boundary of Kohl's colouring is exactly where the other sheet attaches.

## 6. Negative results, stated plainly

* The Petersen graph obstructs nothing:
  * no Kuratowski graph occurs in a single sheet (PROVED);
  * where they occur (U), the kernels are Tait-colourable (FINITE-EXACT: every N ≤ 1000, and N = 2000, 4000).
* Kempe switches cannot touch Collatz (the graph is unchanged).
* The natural sign cover is trivial. The projective window 52 ≤ N ≤ 75 exists, but its planar double covers are not arithmetic.
* The mediant "11 = 8+3" holds for clocks only. There is no word-level mediant, and the −17 cycle word is unbalanced.
* Nothing here reduces Collatz. The new content is structural:
  * the Banach/Kempe identification;
  * the defect-filling fold;
  * the PROVED independence of the controls;
  * PROVABLY transitive Kohl-shape groups one switch from G_C;
  * the non-wqo verdict with a continuous substitute.

**Next obligations (OPEN).**
* (i) The 37 "Collatz-like" switches of modulus ≤ 24 are each a new Collatz-type transitivity question. Are any of them equivalent to Collatz, or provable by a two-step descent?
* (ii) An arithmetic (non-SAT) 3-edge-colouring of the kernels of U, extending Kohl's colouring by the up− colour.
* (iii) A second, independent certificate of linklessness for N ≤ 91 (a linkless embedding or a second solver).

## 7. Reproduction and sources

```bash
bash 04-computation/experiments/procgen_kuratowski_20260925_run.sh > 05-knowledge/results/procgen_kuratowski_20260925.out
```

* **Runtime and resources.** One process at a time. Peak memory is 288 MB. The final run took 216 s wall-clock (311 s CPU) on a loaded machine, dominated by CP-SAT infeasibility proofs; CP-SAT runs with 1 worker when it prints models (deterministic: the U6 model was identical across two runs) and 2 workers for infeasibility.
* **What is reproducible.** Output lines are deterministic except timings.
* **Checks.** Every claim in the scripts is a `check(...)` that raises on failure. Solver claims are cross-checked:
  * the K5 absence by the perfect-matching argument;
  * the Petersen absence by an exhaustive subdivision search, with a positive control;
  * the Petersen-family absence at N = 91 is solver-certified only (K6 by a count).
* **Sources read this session.**
  * Kohl's paper (the local author-PDF text, see §1.1).
  * Wikipedia raw pages fetched 2026-09-25 with a generic user agent: *Planar cover*, *Colin de Verdière graph invariant*, *Petersen family*, *Snark (graph theory)*, *Nowhere-zero flow*, *Regular matroid*, *Graphic matroid*, *Kuratowski's theorem*, *Petersen graph*, *Y-Δ transform*.
* **CITED through those summaries only (primary texts UNVERIFIED):**
  * Kuratowski 1930 and Wagner 1937;
  * Tutte 1958 and 1965;
  * Colin de Verdière 1990 and Lovász–Schrijver 1998;
  * Robertson–Seymour–Thomas (linkless embeddings);
  * Negami 1986 (planar double cover ⟹ projective-planar, for connected graphs);
  * Archdeacon 1981 (35 projective obstructions);
  * the snark theorem (announced by Robertson–Sanders–Seymour–Thomas);
  * Inoue et al. 2026 (apex preprint).
* **Inputs used.**
  * THM-4470, THM-4471, THM-4472;
  * the mod-192 note (transport theorem);
  * the atlas §5;
  * codex's [decoder minors](decoder_minors_20260925.md) (the +2/doubling graph G_N, planar iff N ≤ 15) and [seam threshold](seam_threshold_20260925.md) notes, which are not redone here: their graph is a different operation graph;
  * THM-261 and THM-4113 (half-Kempe chains end at degree-2 vertices, exactly as the Kohl chains end at the Tait defects here).

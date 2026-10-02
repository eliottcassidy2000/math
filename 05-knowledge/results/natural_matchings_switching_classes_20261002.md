# OPEN-Q-060: natural matchings of tournament switching classes with untwisted Euler graphs exist only for n <= 4

**Provenance:** Claude thread session (project thread, branch `claude/project-thread-96889h`), 2026-10-01/02;
temporary identity, no `.machine-id`. This is the "tbij" lane named in the wave-24 letter of
collatz-procgen-20260922 (OPEN-Q-060's bijective form), continued here; no tbij result of that session was on
`origin/main` at 05791bbd. It continues the OPEN-Q-060 item of
[procgen_selfie_20261001_selfie_tournaments.md](procgen_selfie_20261001_selfie_tournaments.md) §5 and
[THM-4524](../../01-canon/theorems/THM-4524-selfie-tournaments-loop-gauge-arc-parity-and-the-odd-mallows-sloane-partner.md),
which records "the bijective form remains open".
**Status:** COMPLETE for this pass; independent audit OWED. **Promoted:** nothing; candidates are listed in §11
(no THM/HYP IDs reserved).
**Labels used:** PROVED, CITED, FINITE-EXACT, VERIFIED, CONDITIONAL, OPEN, EMPIRICAL.
**Code:** `04-computation/experiments/natural_matchings_switching_classes_20261002_run.py` (runner, ends ALL CHECKS PASSED) with helpers
`04-computation/experiments/natural_matchings_switching_classes_20261002_lib.py`. Needs numpy, sympy and nauty (`geng`, `gentourng`,
`labelg`, `dreadnaut`, with or without the `nauty-` prefix).
**Output:** `05-knowledge/results/natural_matchings_switching_classes_20261002.out` (runner stdout, `--full`).
**Parallel work.** [THM-4531](../../01-canon/theorems/THM-4531-twisted-mallows-sloane-principle-and-no-natural-bijection-for-switching-classes.md)
(PROVED + INDEPENDENTLY AUDITED) and its note [procgen_tbij_20261001_natural_bijections.md](procgen_tbij_20261001_natural_bijections.md)
landed on `origin/main` (80cbda5c) while this note was being finished; the two were written independently.
Overlap: Lemma 0 = THM-4531 N1; Theorem 1 = P1; Theorem 3 = O4 (one direction); the remark at the end of §1 = O1;
Theorem 11 = P3. The census deficits of Theorems 4 and 9 agree with THM-4531's exhaustive matchings, which go one
size further (n <= 10 for classes, n <= 9 for tournaments). THM-4531 also has Theorem G, P2, Theorem B and S,
which are not here.
New here: Theorems 5-6 extend O2 from 5 <= n <= 100 to 5 <= n <= 10^6. Rigid blocks are closed under
lexicographic products (Lemma 5.3(d)), so every order in Sigma is available, a superset of THM-4531's twist-rigid
sizes (prime powers = 3 mod 4 or 5 mod 8, and 9). The forced classes are told apart by hand, and A(n) was computed
up to 10^6 (A(n) >= 2 except at n = 8, 16, 40, which have computer certificates). This reduces the switching half of
HYP-9172(a) to binary additive problems over Sigma (density about x/(log x)^(1/4)). Theorem 10 adds two infinite
families to O3 (n = 21, 30, 35, 39, 42, 45, ... are new). Theorem 2 (the level lemma), the twisted Euler graph
counts and the correction of §7.4 are not in THM-4531.

Notation. V = {0, ..., n-1}. A tournament T has scores s_T(v) (out-degrees). Switching T at U reverses every
arc with exactly one end in U; [T] is the switching class (2^(n-1) labelled members). An Euler graph has all
degrees even. For a graph F and g in Aut(F), eps_F(g) = (-1)^(number of edges whose orientation g reverses),
for any fixed reference orientation of F; this does not depend on the orientation and is a homomorphism on
Aut(F) (THM-4524; Royle et al. Lemma 4.1). F is **untwisted** (in the language of Royle et al.: **even**) if
eps_F = 1 on Aut(F), and **twisted** (odd) otherwise. Cycle form (THM-4524): eps_F(g) = (-1)^k, k = number of
even-length cycles of g whose antipodal pairs {x, g^(m/2) x} are edges of F. A permutation is **level** if all
its cycle lengths have the same 2-adic valuation (THM-479, Babai-Cameron).

## 0. Headlines

| # | Result | Status |
|---|--------|--------|
| 1 | "Natural" = S_n-equivariant. A natural bijection on isomorphism classes exists iff a perfect matching exists in which each class is matched to an object whose stabilizer contains (a conjugate of) the class stabilizer (Lemma 0). | PROVED |
| 2 | **Odd n:** every switching class contains exactly one tournament whose scores are all congruent to (n-1)/2 mod 2. So for every odd n, A049313(n) = number of tournaments with all scores = (n-1)/2 mod 2, and "class -> that member" is a natural bijection (Theorem 1). | PROVED. n = 1 mod 4 is THM-1470; n = 3 mod 4 (all scores odd) was found independently of THM-4531 P1. Also the tournament case of Higashitani-Ueyama Thm 4.4 (CITED) |
| 3 | **Level lemma:** if g is level, every g-invariant Euler graph has eps_F(g) = +1; if g is not level, exactly half of them have eps_F(g) = -1 (Theorem 2). Twisted Euler graphs: 0,0,1,1,5,10,42,164,1246,13544,295716,12353336,1030186824 (n = 1..13). | PROVED (new direct proof; equivalent to THM-4524 + THM-479). Counts FINITE-EXACT |
| 4 | **Graphs:** for even n >= 4 no natural map from two-graphs to Euler graphs is injective on isomorphism classes; for odd n Seidel's map is a natural bijection (Theorem 3). | PROVED (elementary) |
| 5 | **The bijective form of OPEN-Q-060:** a natural bijection (switching classes of tournaments) <-> (untwisted Euler graphs) exists for n <= 4 and does not exist for n = 5, ..., 9 (Theorem 4). | FINITE-EXACT; n = 5 also by hand |
| 6 | It does not exist for any 5 <= n <= 10^6 (Theorems 5, 6): two non-isomorphic "forced" classes, built from rigid blocks (Paley and quartic-residue tournaments and their lexicographic products), must both go to the empty graph. | PROVED, with a computer check of an additive condition. For all n >= 5: CONDITIONAL on that condition |
| 7 | **Royle-Praeger-Glasby-Freedman-Devillers setting** (tournaments vs even graphs; they ask for a natural bijection): a natural bijection exists for n <= 4 and not for n = 5, 6, 7, 8 (Theorems 8, 9), nor for infinitely many further n, e.g. 15, 21, 30, 33, 35, 39, 42 (Theorem 10). | FINITE-EXACT (n <= 8); n = 5 by hand; families PROVED |
| 8 | Side result: for every l >= 2 and n >= 1, switching classes in Alt_n(Z/l) and modular Eulerian matrices are equinumerous up to isomorphism (Theorem 11). Higashitani-Ueyama prove this for l prime or gcd(n, l) = 1, and for l = 4 say they expect it without proof. | PROVED (= THM-4531 P3, found independently); VERIFIED for l in {2, 3, 4, 6, 8}, small n |

**What this settles in OPEN-Q-060.** The question's "still open" part had two items. (a) The bijective
form: unlabelled bijections exist trivially (the counts agree by THM-4524), so the content is in naturality,
and natural bijections exist only for n <= 4 (all n <= 10^6 proved; beyond, conditional). For odd n the
natural second incarnation is a class of tournaments, not of graphs (Theorem 1). (b) "How the twist sees
THM-479's level split": exactly through Theorem 2. The same obstruction answers, in the equivariant sense
and for n = 5..8 and infinite families, the open problem of Royle et al. Two earlier repo statements are
refined: THM-1470's n = 3 mod 4 clause (§2), and the claim in HYP-3799 and `07-reflections/one-word-two-even-graphs.md`
that the odd/even |Aut| mismatch is what forbids naturality (§7.4).

## 1. Natural maps and the matching criterion

Let X, Y be finite S_n-sets with the same number of orbits. A **natural matching** is an S_n-equivariant map
f: X -> Y that induces a bijection X/S_n -> Y/S_n. Any construction that uses only the structure of its
input, never the vertex names, is equivariant (in Joyal's language, f is a morphism of species in degree n).

**Lemma 0 (PROVED).** A natural matching exists iff the bipartite graph on X/S_n and Y/S_n, joining O to O'
when Stab(x) <= Stab(y) for some x in O, y in O', has a perfect matching.
*Proof.* If f is equivariant then Stab(x) <= Stab(f(x)). Conversely, given a perfect matching, choose for each
O a representative x_O and y_O in the matched orbit with Stab(x_O) <= Stab(y_O) and set f(g x_O) = g y_O.
This is well defined because g x_O = h x_O implies h^(-1) g in Stab(x_O) <= Stab(y_O). QED

We use X = labelled switching classes of tournaments (stabilizer Stab(C) = {g : g(C) = C}) and Y = labelled
untwisted Euler graphs (stabilizer Aut(F)); in §7, X = labelled tournaments and Y = labelled even graphs.
Only this direction is possible: for n >= 3 the empty graph is S_n-fixed, but no switching class is (a
transposition is not level, so it fixes no class, THM-479 / Corollary 2.2), and no tournament is.

## 2. Odd n: the canonical member

**Theorem 1 (PROVED).** Let n be odd and c = (n-1)/2 mod 2. Every labelled switching class of tournaments
contains exactly one tournament T* with s_{T*}(v) = c (mod 2) for all v.
*Proof.* Switching at U changes s(u) by |V \ U| (mod 2) for u in U and by |U| (mod 2) for u not in U, since u
gains or loses arcs only across the cut. As n is odd, exactly one of U, V \ U has even size; call it E. The
two rules say that switching at U flips the score parity exactly on E. The 2^(n-1) members of [T] are the
switchings at the unordered pairs {U, V \ U}, so they correspond bijectively to the even subsets E of V.
Let D = {v : s_T(v) != c (mod 2)}. Summing scores, C(n,2) = c + |D| (mod 2), and C(n,2) = c (mod 2) because
n is odd, so |D| is even. The member with E = D has all scores = c (mod 2), and it is the only one because E
is the set of flipped vertices. QED
**Corollary 1.1.** C -> T*(C) is S_n-equivariant, so Stab(C) = Aut(T*(C)) has odd order, and for odd n
A049313(n) = number of tournaments with all scores = (n-1)/2 (mod 2). *Verified:* 1, 2, 12, 792 for
n = 3, 5, 7, 9 (gentourng), and the labelled statement for every tournament on 5 vertices and a sample on 7.

*Relation to repo work.* For n = 1 mod 4 this is THM-1470 (even tournaments). For n = 3 mod 4 THM-1470 records
only that no member has all scores even and exactly n members have a single odd score; Theorem 1 adds the
all-odd member, and the n single-odd members are its switchings at one vertex. *Literature.* Write a
tournament as the matrix M_T over Z/4 with entries +1 (i -> j) and -1 (j -> i). The row sum of row v is
2 s(v) - (n-1), so "all scores = (n-1)/2 mod 2" means "all row sums = 0 mod 4", which Higashitani-Ueyama
(arXiv:2409.10904, Def 4.1) call modular Eulerian; Z/4-switching by a vector with constant parity is
tournament switching. Their Theorem 4.4 (gcd(n, l) = 1 gives exactly one modular Eulerian member per class,
up to isomorphism) therefore implies Theorem 1. They do not mention tournaments; the specialisation and the
direct proof above are ours.

## 3. The level lemma: how the twist sees the level split

**Theorem 2 (PROVED).** Let g in S_n and let W_g be the F_2-space of g-invariant Euler graphs.
(a) If g is level, eps_F(g) = +1 for every F in W_g.
(b) If g is not level, F -> eps_F(g) is a nontrivial character of W_g, so exactly half of W_g has eps_F(g) = -1.
*Proof.* For F in W_g and a cycle c of g of length m_c, let k_c(F) = 1 if m_c is even and the antipodal orbit
{x, g^(m_c/2) x} (x in c) lies in E(F), else 0; then eps_F(g) = (-1)^(sum_c k_c(F)) and F -> sum_c k_c(F) mod 2
is linear on W_g (orbit membership is additive under symmetric difference).
(a) If all cycles are odd there is nothing to prove. Otherwise all m_c = 2^a u_c with a >= 1 and u_c odd. Fix
x in c. Inside c, x has two neighbours in each non-antipodal pair orbit of F and one in the antipodal orbit,
so it has k_c(F) neighbours there mod 2. Let d_{cc'} = number of F-neighbours of x in another cycle c' (the
same for all x in c). Counting F-edges between c and c' gives m_c d_{cc'} = m_{c'} d_{c'c}, so
u_c d_{cc'} = u_{c'} d_{c'c} and d_{cc'} = d_{c'c} (mod 2). Every degree is even, so k_c(F) + sum_{c' != c} d_{cc'} = 0,
and summing over c gives sum_c k_c(F) = sum over unordered pairs {c, c'} of (d_{cc'} + d_{c'c}) = 0 (mod 2).
(b) Let c' be a cycle of maximal 2-adic valuation a' (>= 1, as g is not level) and c a cycle of smaller
valuation a. Let F0 = (antipodal orbit of c') + (one g-orbit O of pairs between c and c'). In F0 a vertex of
c' has degree 1 + |O|/m_{c'} = 1 + lcm(u_c, u_{c'})/u_{c'} (even), a vertex of c has degree
|O|/m_c = 2^(a'-a) lcm(u_c, u_{c'})/u_c (even), and every other vertex has degree 0. So F0 is in W_g and
sum k(F0) = 1. QED

**Corollary 2.1 (counts; FINITE-EXACT).** dim W_g = orb_2(g) - orb(g) + [g has an odd cycle] (orb_2 = orbits on
pairs; checked by GF(2) rank for every cycle type, n <= 13). By THM-4524's per-permutation identity
#Fix_g(classes) = sum_{F in W_g} eps_F(g), Theorem 2 gives Fix_g = |W_g| for level g and 0 otherwise (the
Babai-Cameron level formula, THM-479), and
#twisted Euler graphs = A002854(n) - A049313(n) = (1/n!) sum_{g not level} |W_g|
= 0, 0, 1, 1, 5, 10, 42, 164, 1246, 13544, 295716, 12353336, 1030186824 (n = 1..13),
matching the A002854 / A049313 values recorded in `05-knowledge/results/switching_classes_level_burnside_cbx2.out`.
This is the switching analogue of Royle et al.'s #graphs = #tournaments + #odd graphs, with "g has odd order"
replaced by "g is level".
*Relation to repo work.* Given THM-4524's identity and THM-479, Theorem 2 follows (a sum of signs that equals
|W_g| or 0); the converse also holds. The new content is the direct degree-parity proof and the explicit
twisting witness F0.
**Corollary 2.2.** g fixes some switching class iff g is level (Babai-Cameron Lemma 7.1, THM-479), recovered
from Theorem 2 and THM-4524. In particular no transposition fixes a class when n >= 3.

## 4. Graphs: the even-n wall

**Theorem 3 (PROVED, elementary).** For even n >= 4 there is no S_n-equivariant map from labelled two-graphs
(switching classes of graphs) to labelled Euler graphs that is injective on isomorphism classes. For odd n,
"class -> its unique Euler member" (Seidel 1974) is a natural bijection.
*Proof.* The classes [empty graph] = {K_{U, V\U}} and [K_n] = {K_U + K_{V\U}} are S_n-fixed and distinct for
n >= 3, so an equivariant f sends both to S_n-fixed Euler graphs. The only S_n-invariant graphs are the empty
graph and K_n, and K_n is not Euler when n is even. So f sends both classes to the empty graph. QED
Mallows-Sloane (1975) prove the equality for all n by counting; the OPEN-Q-060 record already noted their
remark that it is not bijective for even n. Theorem 3 makes the obstruction precise. The same argument says
nothing for tournaments, because no tournament class is S_n-fixed; there the obstruction appears only at
n = 5 and needs the forced classes of §6.

## 5. The census: n <= 9

**Theorem 4 (FINITE-EXACT).** A natural bijection between switching classes of tournaments and untwisted
Euler graphs exists for n <= 4 and does not exist for 5 <= n <= 9. Maximum natural matchings (Lemma 0, with
class stabilizers computed exactly and every invariant untwisted Euler graph enumerated):

| n | classes = untwisted Euler graphs | Euler graphs | maximum natural matching | deficiency |
|---|---|---|---|---|
| 3 | 1 | 2 | 1 | 0 |
| 4 | 2 | 3 | 2 | 0 |
| 5 | 2 | 7 | 1 | 1 |
| 6 | 6 | 16 | 5 | 1 |
| 7 | 12 | 54 | 10 | 2 |
| 8 | 79 | 243 | 76 | 3 |
| 9 | 792 | 2038 | 788 | 4 |

At n = 9 a Hall violator has 11 classes (stabilizer orders 81, 21, 15, 15, 9, 9, 9, 9, 9, 5, 5) and 7 allowed
partners. *n = 5 by hand.* By Theorem 1 the two classes are those of the rotational tournament R_5
(Aut = Z_5) and of source => 3-cycle => sink (scores 0, 2, 2, 2, 4; Aut = Z_3). The untwisted Euler graphs on
5 vertices are the empty graph and C4 + K1 (C5, K5, K3, the bowtie and K5 - K3 are twisted). Aut(C4 + K1) has
order 8, so it contains no element of order 3 or 5, and both classes must go to the empty graph. QED

## 6. All n from 5 to 10^6: rigid blocks and forced classes

**Definition.** A *rigid block* is a tournament B on w vertices with a subgroup H_B <= Aut(B) that is
transitive on V(B) and such that every nonempty H_B-invariant graph on V(B) is twisted. (H_B has odd order,
as Aut(B) does, so w is odd.) A switching class C is *forced* if some H <= Stab(C) has the empty graph as its
only invariant untwisted Euler graph. By Lemma 0 a natural matching sends every forced class to the empty
graph, so **two non-isomorphic forced classes rule out a natural bijection.**

**Lemma 5.1 (module lemma, PROVED).** Let M be a module of a graph F (each vertex outside M is joined to all
or none of M) and g in Aut(F[M]) with eps_{F[M]}(g) = -1. Then g extended by the identity is in Aut(F) and
has eps_F = -1. *Proof.* Orient every edge leaving M away from M; the extension maps these edges to edges of
the same kind without reversing them. QED

**Lemma 5.2 (blow-ups, PROVED).** Let F = Gamma[W_1, ..., W_m] (vertex i of Gamma replaced by an independent
set W_i, an edge ij by all W_i-W_j edges) with every |W_i| odd. Then F is Euler iff Gamma is Euler, and an
automorphism gamma of Gamma with |W_gamma(i)| = |W_i| lifts to an automorphism of F with the same eps.
*Proof.* deg_F(x) = sum_{j ~ i} |W_j| = deg_Gamma(i) (mod 2) for x in W_i. With the orientation induced from
Gamma, each Gamma-edge carries |W_i||W_j| (odd) edges, all reversed or all not. QED

**Lemma 5.3 (rigid blocks exist for every w in Sigma, PROVED).** Sigma = odd numbers with no prime factor = 1
(mod 8).
(a) One vertex (vacuous).
(b) Paley P_q, q = 3 mod 4 a prime power, H = {x -> ax + b : a a nonzero square}: H is transitive on unordered
pairs (-1 is a non-square), so the invariant graphs are empty and K_q, and K_q is twisted (a transposition
reverses one edge).
(c) R_q, q = 5 mod 8 a prime power: arcs x -> y iff y - x in D = M u gM, M the fourth powers, g a primitive
element. Since -1 lies in g^2 M, D and -D partition the nonzero elements, so R_q is a tournament. H = {x -> ax
+ b : a in M} has odd order q(q-1)/4 and two orbits on pairs (difference square or not), so the invariant
graphs are empty, the Paley graph, its complement and K_q. The Paley graph and its complement are twisted by
x -> -x: it has (q-1)/2 two-cycles {x, -x}, of which (q-1)/4 (odd) are edges.
(d) Lexicographic products A[B] of rigid blocks, with H = H_B wr H_A: an invariant graph is either nonempty
inside a copy of B (twisted by Lemma 5.1) or a blow-up of a nonempty H_A-invariant graph on A by independent
sets of odd size |B| (twisted by Lemma 5.2), or empty.
Every w in Sigma is a product of primes = 3, 5, 7 (mod 8), so (b) with q = p = 3 mod 4, (c) with q = p = 5 mod
8, and (d) give a rigid block of order w. QED

**Theorem 5 (forced classes, PROVED).** Let T = Q[B_1, ..., B_m] with rigid blocks (B_i, H_i) and H = H_1 x ...
x H_m <= Aut(T). If m <= 2, or m = 3 and two blocks have the same order, then [T] is forced.
*Proof.* Let F be H-invariant, Euler and untwisted. H_i x H_j is transitive on B_i x B_j, so between two blocks
F has all edges or none, and each block is a module of F. A nonempty F[B_i] is H_i-invariant, hence twisted,
hence F is twisted (Lemma 5.1); so F = Gamma[B_1, ..., B_m] with Gamma Euler (Lemma 5.2). On at most two
vertices the only Euler graph is empty. On three vertices the other one is the triangle, giving the complete
3-partite graph; if |B_1| = |B_2| = w, swapping B_1 and B_2 reverses the w^2 (odd) edges between them and no
others, so it is twisted. QED
With four or more blocks no class of this shape is forced: Gamma = C4 gives a complete bipartite graph with
even parts, which is Euler and untwisted.

**Theorem 6 (PROVED for 5 <= n <= 10^6; CONDITIONAL beyond).** For n >= 5 let A(n) be
* n odd: #{(w, w') : 2w + w' = n, w, w' in Sigma} + [n is a prime in Sigma];
* n = 2 (mod 4): #{w1 <= w2 in Sigma : w1 + w2 = n};
* n = 0 (mod 4): #{primes q = 3 (mod 4) : n/2 < q <= n - 1, n - q in Sigma}.
If A(n) >= 2, there is no natural bijection at n. A(n) >= 2 holds for all 5 <= n <= 10^6 except n = 8, 16, 40;
for those (and for every 5 <= n <= 100) two non-isomorphic forced classes were exhibited and separated by
computer (the multiset of canonical forms of the n tournaments obtained by switching each vertex to a source
and deleting it, an isomorphism invariant of the class). **Hence no natural bijection exists for any
5 <= n <= 10^6.**
*Proof.* Forced classes (Theorem 5): for odd n, TT3[B_w, B_w, B_w'] for each representation, plus the rigid
block itself if n is a prime in Sigma; for even n, B_{w1} => B_{w2} (for n = 0 mod 4 with B_{w2} = P_q).
Non-isomorphism:
* *Odd n.* Rigid blocks are regular, so score parities are constant on blocks and the canonical member T* of
  Theorem 1 is again Q*[B, B, B'] with Q* on 3 vertices. If Q* is transitive, the strong components of T* are
  the blocks; if Q* is the 3-cycle, the blocks are the maximal proper modules (the quotient is prime). Either
  way the block-size multiset {w, w, w'} is an invariant, and it determines (w, w'). The single block of prime
  order n is a prime tournament (a circulant tournament of prime order has no nontrivial module: translates of
  a minimal module of size >= 3 would be disjoint, and a 2-element module {x, x+d} forces, for every vertex
  c, c -> c+d iff c -> c-d, which is impossible), while the three-block T* are decomposable.
* *n = 2 mod 4.* For even n, switching at U flips no parity (|U| even) or all parities (|U| odd), so the
  unordered partition of V by score parity is a class invariant. For B_{w1} => B_{w2}, vertices of B_1 have
  score (w1-1)/2 + w2 and those of B_2 have (w2-1)/2; as w1 = w2 (mod 4), these parities differ, so the
  partition is {B_1, B_2} and its part sizes {w1, w2} are an invariant.
* *n = 0 mod 4.* Now the parity partition is trivial. For a pair {u, v}, restrict the class to V \ {u, v}
  (n - 2 = 2 mod 4 vertices) and let beta(u, v) be the smaller part of its parity partition. The parity of
  x after deleting u, v changes by [x -> u] + [x -> v], so the parts are separated by S(u, v) = {x : x beats
  exactly one of u, v}. Call {u, v} balanced if |S(u, v)| = (n-2)/2. For B_{w1} => P_q: mixed pairs give
  |S| = (w1-1)/2 + (q-1)/2 = (n-2)/2 (balanced); pairs inside B_1 give |S| <= w1 - 2 < (n-2)/2; pairs inside
  P_q give |S| = (q-1)/2 (P_q is doubly regular), balanced iff q = n - 1. So the number of balanced pairs is
  w1 (n - w1) for w1 >= 5 and C(n, 2) for w1 = 1, which determines w1. QED
*Computation.* `arith`: A(n) for all n <= 10^6 by FFT convolution; the least value of A(n) over n >= 1000 is
21, 169, 156, 146 in the classes n = 0, 1, 2, 3 (mod 4) (EMPIRICAL growth about n / sqrt(log n)). Allowing
prime powers q = 3 mod 4 in the n = 0 mod 4 case (Lemma 5.3(b) covers them) removes n = 40 from the
exception list but needs Paley tournaments over F_q, so the runner uses primes.
*What is missing for all n.* A(n) >= 2 for all n > 10^6 is a binary additive problem: Sigma has density about
x / (log x)^(1/4) (primes = 1 mod 8 sifted out, a sieve of dimension 1/4, two variables give dimension 1/2),
and the n = 0 mod 4 case also involves primes. We expect it to be provable for large n with sieve methods
but did not attempt it, so the statement for n > 10^6 stays CONDITIONAL.

## 7. The Royle-Praeger-Glasby-Freedman-Devillers setting

Royle, Praeger, Glasby, Freedman and Devillers (J. Algebraic Combin. 57 (2023) 515-524, arXiv:2204.01947)
prove by Cauchy-Frobenius that tournaments and even graphs (graphs with no automorphism reversing an odd number
of edges; our "untwisted" for arbitrary graphs) are equinumerous, and write: "It is an open problem to find a
natural bijection between the sets of unlabelled even graphs and tournaments on n vertices." With "natural"
read as in §1 (the only direction possible is tournaments -> even graphs):

**Theorem 8 (PROVED, n = 5 by hand).** There is no natural bijection for n = 5. The five tournaments with a
nontrivial automorphism (R_5 with Aut = Z_5, and four with Aut = Z_3: scores 0,1,3,3,3; 0,2,2,2,4; 1,1,1,3,4;
1,2,2,2,3) can only go to even graphs invariant under an element of order 5 or 3. For order 5 the invariant
graphs are empty, C5 (twice) and K5, and only the empty graph is even (a reflection of C5 reverses 3 edges).
For order 3 the even ones are the empty graph, K_{1,3} + K1, K_{1,4} and K_{2,3}. Five tournaments, four
possible partners. QED

**Theorem 9 (FINITE-EXACT).** A natural bijection exists for n = 3, 4 and does not exist for n = 5, 6, 7, 8.
Maximum natural matchings 11/12, 54/56, 451/456, 6868/6880 (deficiencies 1, 2, 5, 12); Hall violators 5/4,
13/11, 49/44, 292/280 (tournaments / allowed even graphs), all violators having |Aut| in {3, 5, 9, 15, 21}.
At n = 3 the bijection is C3 -> empty graph, TT3 -> path P3.

**Theorem 10 (PROVED).** Call a tournament *Royle-forced* if some H <= Aut(T) has the empty graph as its only
invariant even graph. Rigid tournaments (Lemma 5.3) are Royle-forced, and so is B => B' for rigid B, B' of the
same order w: inside blocks Lemma 5.1 applies, and the only other candidate, K_{w,w}, is odd (swapping the
sides reverses w^2 edges). Two non-isomorphic Royle-forced tournaments rule out a natural bijection. Hence none
exists for n in Sigma with at least two distinct prime factors (B_p[B_{n/p}] for two different primes p | n
have different prime quotients in their modular decomposition) and for n = 2w with w such a number
(B => B and B' => B' differ in their strong components). Examples: 15, 21, 30, 33, 35, 39, 42, 45, 55, 57, 63, 65, 66, 69, 70, ...
With three or more blocks no class of this kind is Royle-forced (a path quotient gives K_{w, w'+w''}, which is
even), so this method does not reach all n here; the switching setting has the extra Euler constraint, which
is why §6 covers much more.

**7.4 Correction to the repo's explanation.** HYP-3799, `07-reflections/one-word-two-even-graphs.md` and the
boxeph S159 broadcast say the natural bijection is blocked because tournament automorphism groups have odd order
while even graphs (n <= 5) have even-order groups. That mismatch rules out bijections preserving the
automorphism group, but not natural ones: equivariance only needs Aut(T) <= Aut(G), and at n = 3, 4 natural
bijections exist although the |Aut| profiles are disjoint. The actual obstruction is a Hall violation, first at
n = 5 (Theorem 8).

## 8. Side result: the l = 4 question of Higashitani-Ueyama

Alt_n(Z/l) = skew-symmetric n x n matrices over Z/l with zero diagonal; switching at v adds X_v (-1 in row v,
+1 in column v); D = <X_1, ..., X_n>; E = {M : every row sum is 0} (modular Eulerian); S_n acts by simultaneous
permutation. s_{l,n} = |(Alt/D)/S_n|, t_{l,n} = |E/S_n|. Higashitani-Ueyama (arXiv:2409.10904v4) prove s = t
for gcd(n, l) = 1 (Thm 4.4, bijectively) and for l prime (Thm 4.7, Burnside), and in Example 4.10 write that for
l = 4 they expect s_{4,n} = t_{4,n} for every n without a proof (read through a web summariser in two separate
fetches; the paper's own text was not machine-readable here).

**Theorem 11 (PROVED).** s_{l,n} = t_{l,n} for all l >= 2 and n >= 1.
*Proof.* <M, N> = sum_{i<j} m_ij n_ij is a well-defined (the product m_ij n_ij is unchanged by swapping i, j),
S_n-invariant, perfect Z/l-valued pairing on Alt_n(Z/l). The coefficient of s_v in <delta(s), N> is minus the
v-th row sum of N, so D^perp = E and E is the Pontryagin dual of Alt/D as an S_n-module. For a finite abelian
group X with automorphism g, |X^g| = |ker(g-1)| = |coker(g-1)| = |(X^)^g|. Burnside's lemma gives equal orbit
counts. QED (This is Brauer's permutation lemma for the abelian group Alt/D; Mallows-Sloane is l = 2.)
*Verified* by exact Burnside over both sides: l = 2: 1, 2, 3, 7, 16 (n = 2..6, = A002854); l = 3: 1, 2, 4, 14,
120; l = 4: 1, 3, 8, 62, 1760; l = 6: 1, 4, 17, 461; l = 8: 1, 5, 32, 2332 (n from 2). The tournament count
A049313 is the analogous count on the all-odd fibre of Alt_n(Z/4), where S_n acts affinely and Theorem 2's
sign appears.

## 9. Literature (searched 2026-10-01; 31 sources, mostly through a web summariser)

* Mallows-Sloane, SIAM J. Appl. Math. 28 (1975): two-graphs = Euler graphs, by counting (not read directly;
  secondary sources agree). Seidel 1974: the odd-n bijection.
* Babai-Cameron, Electron. J. Combin. 7 (2000) R38: Lemma 3.1, Lemma 7.1 and Thm 7.2 (level permutations,
  A049313). No Euler graphs, bijection or score-parity member.
* Royle-Praeger-Glasby-Freedman-Devillers (2023): tournaments = even graphs by counting; natural bijection
  posed as open; no switching classes or Euler graphs.
* Higashitani-Ueyama, arXiv:2409.10904 (2024-25): §2 and §8 above.
* Cameron, Math. Z. 157 (1977) (cohomology of two-graphs), Harries-Liebeck (1978, a Klein four-group stabilising
  a class without fixing a member), Liebeck (1982), Cheng-Wells, JCTB 40 (1986) (the l = 3 count), Moorhouse
  (1995, skew two-graphs): background, no equivariance statements found.
We found no source stating an equivariance obstruction for either equinumerosity, the odd-n score-parity
member for n = 3 mod 4, or the twisted Euler graph counts. Mallows-Sloane and Cameron 1977 could not be read
directly, so priority is not claimed beyond this search. Blocked during the search: arxiv.org and oeis.org
from the shell (HTTP 403), the OEIS search endpoint for the fetch tool (robots.txt).

## 10. OPEN

* A(n) >= 2 for all n > 10^6, or another source of two forced classes, to make Theorem 6 unconditional.
* Royle setting: THM-4531 O3 settles 5 <= n <= 15 and 16 larger n, and Theorem 10 adds two infinite families.
  The smallest n settled by neither is 16.
* The deficiency of the best natural matching: 1, 1, 2, 3, 4 for n = 5..9 (switching) and 1, 2, 5, 12 for
  n = 5..8 (Royle); THM-4531 adds 10 at n = 10 and 44 at n = 9 (Royle). Is it unbounded? A(n) gives a lower
  bound of A(n) - 1 for n <= 10^6.
* Weaker notions of naturality (equivariance under a subgroup, or natural correspondences that are bijective
  only on isomorphism classes of rigid objects) remain unexplored.

## 11. Promotion candidates

* Theorems 1 and 11 are already canon as THM-4531 P1 and P3 (found independently); nothing to promote.
* Theorems 5-6 as an extension of THM-4531 O2 to 5 <= n <= 10^6, and as progress on HYP-9172(a) (the reduction
  to A(n) >= 2).
* Theorem 10 as an extension of THM-4531 O3 (two infinite families).
* Theorem 2 (the level lemma) as the direct proof of how the twist sees THM-479's split, with the twisted counts.
* The §7.4 correction, against HYP-3799 and `07-reflections/one-word-two-even-graphs.md`.

## 12. Reproduction

`nice python3 -u 04-computation/experiments/natural_matchings_switching_classes_20261002_run.py --full > 05-knowledge/results/natural_matchings_switching_classes_20261002.out`
(35 min single-threaded on a shared 4-core machine: 2106 s, of which S4 with the n = 9 census took 683 s and
S12 took 827 s; the deposited run ended with 4573 checks and ALL CHECKS PASSED).
Without `--full` it takes 49 s and makes 4161 checks. Sections: S1 Theorem 1; S2 Theorem 2 and
the counts; S3 Theorem 3; S4 census (n <= 8, `--full` adds 9); S5 Royle census (n <= 7, `--full` adds 8) and
the n = 5 partners; S6 rigid blocks; S7 forced classes and class keys (5 <= n <= 60, `--full` to 100); S8 the
arithmetic condition to 10^6; S9 beta counts; S10 prime and doubly regular blocks; S11 Royle-forced families;
S12 Theorem 11.

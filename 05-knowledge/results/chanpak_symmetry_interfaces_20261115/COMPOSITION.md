# Exact composition of path covers and two notions of tournament symmetry

The point of the composition formula is to keep the physical continuation state. A scalar count does not, in general, determine what a quotient sees. This note gives an all-size recursion, an arithmetic symmetry test for exact targets, and a physical family showing that this test is stricter than testing the ordinary automorphism group. It does not prove spectrum coverage.

## 1. The substitution compiler

For a tournament T on b vertices let F_T(k) count **ordered** partitions of its vertices into k **nonempty directed path words**. In particular F_T(1)=H(T), and F_T(k)=0 for k>b. Let E_T(r) be the optional r-slot count, where slots are ordered and may be empty. Precisely,
```
E_T(r) = sum_k binom(r,k) F_T(k).
```
Finite binomial inversion recovers F_T from E_T. In an ordinal join, E factors at every r: each slot restricts to a possibly empty path word in every component; its blocks concatenate in the component order, with all connecting arcs present.

For an arbitrary quotient tournament Q and nonempty modules T_i, form Q[T_i] by orienting every arc between T_i and T_j as i→j in Q. A **schedule** for k slots is an ordered k-tuple of nonempty words over V(Q), with consecutive letters distinct and arcs in Q, and every quotient label appearing somewhere. Each occurrence of label i is one run in its module; let m_i be its total multiplicity across slots. Then
```
F_{Q[T_i]}(k) = sum over schedules of product_i F_{T_i}(m_i).
```
Indeed compress each of the k paths into maximal consecutive runs in its module. For each i, the ordered runs, in slot order then left-to-right within the slot, give a physical member of F_{T_i}(m_i). Conversely inserting those paths in their scheduled occurrences gives legal paths with the specified runs. The two constructions are inverse. Only schedules with m_i≤|T_i| need be enumerated. This is an exact recursion for arbitrary-size compositions, not a factorization of H under cyclic substitution.

For the canonical strong-module decomposition, the quotient at each nonleaf is a linear order or prime; see the statement of Gallai's tournament theorem in [Boussaïri–Chergui–Ille–Zaï](https://cdm.ucalgary.ca/article/download/67935/54027). The slot recursion and the linear-order factor rule apply at every such node.

## 2. Profile symmetry forces a divisor, even in a rigid tournament

Let d_i be the gcd of the nonzero F_{T_i}(m), and let the **normalized profile color** of i be the whole sequence F_{T_i}/d_i (including its length). Set D=product_i d_i, and let G be the subgroup of Aut(Q) preserving these colors. Then, simultaneously for every k,
```
D |G| divides F_{Q[T_i]}(k).
```
Factor D out of each summand of the schedule formula; every module occurs. The action of G on **ordered** schedules is free: to fix a schedule an element must fix each letter of each word, and every label is used. The normalized integer weight is constant along an orbit by the color hypothesis. Thus each orbit contributes an integer multiple of |G| after factoring D. The resulting divisor propagates through every ancestor because each ancestor reads some nonzero slot count from this child.

The true group obeys a different rule. For a prime quotient Q with at least three vertices, the children are exactly the maximal proper modules, so any automorphism induces a quotient automorphism preserving **isomorphism classes** of children; its kernel is the product of the child automorphism groups. For completeness, a module X meeting more than one child has a set J of occupied labels which is a module in Q: an outsider block has uniform arrows toward X. If J is proper, primeness forbids |J|≥2; if J is all, an outsider vertex in a child cannot have uniform directions toward X because its quotient label has both an incoming and an outgoing neighbor. (A source or sink in Q would give the complementary proper nontrivial module.) Thus every proper module X lies within a child. Every color-stabilizing quotient permutation lifts by chosen isomorphisms, and representatives of child colors can be used to choose these lifts coherently. In a linear SCC order, the automorphism group simply preserves each SCC and is the product of their automorphism groups.

These two quotient colors can disagree. The schedule divisibility is not a claim that there is an actual action on the vertices or that it detects every arithmetic divisor.

## 3. A size-by-size physical consequence

Take the rigid nearly transitive tournament R of order four, with rows (6,12,8,1), using outgoing bitmasks in zero-based label order. Its profile is (5,22,36,24). It is rigid by direct checking (also, its automorphism group has odd order and its two degree classes each have size two). For t≥2 define three ordinal children
```
A_j = TT_j → R → TT_(t−j),   j=0,1,2,
```
omitting TT_0. They are rigid, pairwise nonisomorphic by the position of their unique nontrivial SCC, and have identical complete slot profiles: E_{A_j}(r)=r^t E_R(r). Set
```
T_rig = C3[A_0,A_1,A_2],    T_sym = C3[A_0,A_0,A_0].
```
Both are strong; the first is rigid and the second has automorphism group C3 by §2. Their entire ordered cover profiles are identical, not merely H. Every entry is divisible by 3, even in the rigid case. The basic size is 3(4+t).

For clarity the single-slot schedule sum for cyclic Q=C3 is explicit. When the three run multiplicities are equal, there are three possible start labels. If they differ by exactly one, the advanced labels must form the first one or two labels on a cyclic word, which uniquely fixes the start. Other multiplicity patterns give zero. In particular if all child profiles equal a sequence f, with f(0)=0, then
```
H/3 = sum_{j≥1} [ f(j)^3 + f(j)^2 f(j+1) + f(j) f(j+1)^2 ].
```
Since j! divides F(j) by permuting the nonempty slots, only j=1,2 contribute to H/3 mod 3, giving (f(1)+f(2))(f(1)^2+f(2)^2). For t odd, f(1)=5 and f(2)=32·2^t−10 is divisible by 3. Hence H is divisible by exactly one factor of 3 in both examples (at t=3 the common H is 72465879938757). Thus the ghost divisor survives even when every physical vertex symmetry has been removed.

This gives a simple construction at **every n≥20**, not just infinitely many matched orders. Start with t=2, order 18, and apply either C3[X,TT_1,TT_1] or C3[X,TT_2,TT_1] to each of the two members, adding two or three vertices. Since the large child is distinguished, no nonidentity quotient rotation survives; the two groups remain respectively trivial and C3. The common slot profile remains identical by the schedule formula, and its common divisibility by 3 propagates. Every n−18≥2 is a sum of twos and threes. This is profile-equivalent physical expansion, not a density theorem for H.

### Every prescribed odd cyclic symmetry

The same construction is available with any odd q=2m+1≥3, not only three. Use the carousel quotient on Z/qZ, with i→j iff j−i is one of 1,...,m. Its automorphism group is C_q: its out- and in-neighborhoods are both transitively ordered, so fixing a vertex fixes every other vertex; the q translations are automorphisms. It is also prime. If a module contains the oriented pair i,i+d (1≤d≤m), it must contain i+m+1,...,i+m+d because each of them distinguishes the pair. Thus it contains a consecutive pair (immediately if d=1, or among those forced vertices otherwise). Rotate this pair to 0,1. The pair forces m+1; then 1,m+1 forces m+2,...,2m,0; then m+1,0 forces 1,...,m. Hence every non-singleton module is the whole quotient.

Take t=q−1 and all q distinct rigid children TT_j→R→TT_(t−j) for 0≤j<q; contrast with q copies of the first child. The two strong tournaments of order q(q+3) have identical full profile, divisible coordinatewise by q, and ordinary automorphism groups respectively trivial and C_q. The same +2 and +3 wrappers give this at every n≥q(q+3)+2 without changing either group or losing divisibility. This makes a prescribed divisor and exact cyclic symmetry size independent of each other; it does not prescribe the full Hamiltonian count.

### Finite modular dynamics of the size wrappers

The +2 and +3 wrappers above are exact linear operators on the whole optional-slot sequence (E_r=E_X(r)). Writing every (E_0=0) for a nonempty core, their outputs respectively are
```
E'_r = r [ E_(r+1) + r E_r − E_(r−1) ],
E''_r = r [ (2r+1)E_(r+1) + (r²−2)E_r − (2r−1)E_(r−1) ].
```
Here is a direct derivation. Use the cyclic orientation (X\to B\to C\to X), and locate the distinguished vertices of (B,C) in the ordered slots. In each scaffold list, an (X)-segment is optional before a first (B) or after a last (C), and mandatory between consecutive (C,B) letters; no other internal (X) segment can fit. An empty scaffold slot allows optional (X). If one segment is mandatory, subtract the optional count without it.

For the two singletons, separate slots contribute (r(r−1)E_r); a shared slot in order (BC) contributes (rE_{r+1}), and in order (CB) contributes (r(E_r−E_{r−1})). For (B=TT_2) and (C) singleton the three occupancy types give
```
r(r−1)(r−2) E_r;
r(r−1) [2E_(r+1)+3E_r−2E_(r−1)];
r [3E_(r+1)−E_r−E_(r−1)].
```
For the middle line, either the two (B)'s share a slot (unique increasing order), or a (B) shares with (C) (both orders, one with a mandatory segment). For the last, of the six permutations of (B_1,B_2,C), two have the invalid adjacent (B_2 B_1); (B_1 B_2 C) permits two optional segments, two (BCB)'s require an internal segment and permit a prefix, and (C B_1 B_2) requires an internal segment only. Adding gives the claimed operator.

For any fixed (M≥2), M divides E_M(X) by the free action of cyclic relabeling on the (M) slots (at least one slot has a vertex). Thus these operators close exactly on the vector of optional-slot counts of indices 1 through M−1 modulo M, with size transitions +2 and +3. This gives a finite modular inverse test for this particular grammar, not a coefficient-positivity test for integer counts. The first coordinates are (H+E_2) and (3E_2−H), respectively. In particular the [odd-cycle packing formula](https://arxiv.org/abs/2307.05569) gives E_2(X) ≡ 2^|X|·H(X) (mod 3): a packing of (j) odd cycles covering (d) vertices can allocate each whole cycle to either slot and distribute its remaining vertices freely, multiplying its ordinary weight by 2^(|X|−d+j) ≡ 2^|X| (mod 3). So every odd-order output of the +2 wrapper is divisible by 3, whereas the +3 wrapper reads −H(X) modulo 3 without that order restriction. Wrapping strong witnesses to hit a prime target must respect this distinction.

The finite-state description yields certified **all-residue coverage at every sufficiently large order inside this specific grammar**. Starting at the singleton, let R_n be all its reachable modular vectors at order n; give absent sizes the empty set. Then R_n=A2 R_(n−2) ∪ A3 R_(n−3). A repeat of three consecutive layers forces an all-future repeat by this exact recurrence. The modular replay gives:

| M | first certified repeated triple endpoint pair | all n at least | exact coverage |
|---:|---:|---:|---|
| 3 | 14, 16 | 12 | all odd residues mod 6 |
| 5 | 27, 29 | 12 | all odd residues mod 10 |
| 7 | 42, 44 | 13 | all odd residues mod 14 |
| 9 | 52, 54 | 26 | all odd residues mod 18 |

The first coordinate is H modulo M. For odd M the count's odd parity makes these equivalent to the asserted residues mod 2M. Every nonempty wrapper word gives a strong tournament because its quotient is the directed triangle with nonempty children. Backtracking through the finite transitions gives literal witnesses. These are separate modulus statements, not simultaneous coverage modulo their least common multiple unless a corresponding state certificate is provided. In fact, synchronizing the two state vectors for moduli 3 and 5 under the **same** sequence of +2/+3 choices gives a period-two repeated triple at orders 146 and 148. There are exactly fourteen of the fifteen residue classes modulo 15 at every order n≥146. The missing class is 10 at even order and 5 at odd order (respectively the odd residues 25 and 5 modulo 30). This is an exact obstruction for the singleton-seeded wrapper grammar, even though both marginals are full. It is not an obstruction for tournaments outside that grammar.

#### The infinite-palette invariant and an explicit bypass

The same scaffold computation works for arbitrary B in the one-singleton quotient C3[X,B,TT1]. Write b_r=E_B(r), x_r=E_X(r), so b_0=x_0=0. The exact output optional-slot operator is
```
E'_r = r [(b_(r+1)−b_r)x_(r+1)
         +(3b_r−b_(r+1)−b_(r−1))x_r
         −(b_r−b_(r−1))x_(r−1)].
```
Choose the slot of the distinguished singleton (r choices). Outside it, each of the r−1 slots has an optional B segment preceded by an optional X segment. In its own slot, the singleton has two possible B regions, before and after. The four occupancy patterns of the B regions have respective B counts b_(r−1), b_r−b_(r−1), b_r−b_(r−1), b_(r+1)−2b_r+b_(r−1). Use occupancy order (neither, before only, after only, both). The X choices are x_r, x_(r+1), x_r−x_(r−1), x_(r+1)−x_r respectively: a B-region before the singleton permits an X-prefix, one after it requires an X-transition from the singleton to B, and with no B-region after it an X-suffix is optional. This proves the operator, which is symmetric in X and B. For B=TT_a it specializes to b_r=r^a with no restriction on a.

For any tournament of order n, the optional-slot polynomial is the sum of 2^j r^(n−d+j) over all packings of disjoint nontrivial odd cycles, where j cycles cover d vertices. This follows by letting each entire cycle choose its slot and each uncovered vertex choose its slot independently. Thus E(−r)=(−1)^n E(r) and evaluation is periodic in r modulo any M.

Write h=E1, e=E2, Q(h,e)=h²−2e² modulo 5, and let χ be the quadratic character modulo 5. The following invariant holds for the **entire binary composition grammar** built from transitive seeds by C3[X,B,TT1]:
```
If h ≠ 0 (mod 3), then Q(h,e) ≠ 0 (mod 5) and
h mod 3 as a sign = (−1)^n χ(Q(h,e)).
```
A transitive seed of size n has h=1, e=2^n; its Q is 4 at even n and 3 at odd n, so the signs agree. For a composition, modulo 3, e=(−1)^n h, and the first output coordinate is (3h_B−e_B)h_X+(e_B−h_B)e_X. If both orders are odd this is zero. Otherwise it is −h_X h_B; thus nonzero output requires two nonzero inputs and at least one even order. In that case E3=(−1)^n E2 modulo 5 for each side, and inserting these identities in the operator at r=1,2 gives
```
(h',e') = [[A,L],[3L,A]] (h_X,e_X),  A=3h_B−e_B, L=e_B−h_B   (mod 5).
```
In the norm Q, this matrix scales by A²−3L²=Q(h_B,e_B) modulo 5. The quadratic characters accordingly multiply, and n'=n_X+n_B+1 supplies the minus sign. This proves the invariant. If h=0 modulo 5 and h nonzero modulo 3, Q=−2e² is nonsquare, forcing sign −1 at even order and +1 at odd order. So the missing classes 10 (even order) and 5 (odd order) modulo 15 are permanent for this grammar regardless of transitive reservoirs, size, or parenthesization. The finite replay above proves that all other classes occur eventually in the smaller singleton-seeded two-wrapper grammar.

**Bypass by padding one nontrivial module ordinally.** Take B=R→TT1 where R is the same rigid order-four tournament, giving (h_B,e_B)=(5,64) at order five. It fails the sign invariant: its nonzero h modulo 3 is −1 but the formula would require +1. Take X of order six obtained from a singleton by +3 then +2; (h_X,e_X)=(37,352). The strong output C3[X,B,TT1] has order 12 and
```
H = (15−64)·37 + (64−5)·352 = 18955 ≡ 10 (mod 15).
```
This is exactly the former forbidden class at even order. Replacing the last +2 inside X by TT_(1+4j) with singleton preserves the mod-3 and mod-5 state, because the relevant powers 2^a and 3^a modulo 5 have period four, and the mod-3 first-coordinate formula depends only on parity. It gives the same bypass at every order 12+4j. Then a TT4-with-singleton wrapper adds five vertices, multiplies (h,e) by 2 modulo 5, and reverses the nonzero sign of h modulo 3. Thus k=0,1,2,3 such wrappers produce the appropriate missing class at orders 12+4j+5k. These cover every order n≥24 with the missing class for its parity. Combining with the certified fourteen other classes, **every odd residue modulo 30 has a strong witness at each order n≥146**. This is a synchronized modular construction, not exact absolute coverage; its significance is that the obstruction is broken by changing the core, not by adding more transitive leaves to the same grammar.

A smaller all-order cutoff is certified without a growing sleeve: add just the one fixed B=R→TT1 operation (increment +6) to the original +2/+3 palette, synchronizing the mod-3 and mod-5 vectors as above. The resulting six-layer transition state repeats at orders 83 and 85; every odd residue modulo 30 is attained at every n≥42. If the TT4 sleeve (+5) is also allowed, that six-layer state repeats at 68 and 70 and its all-residue cutoff is n≥30. The [joint replay](cutoff_bypass.py) checks these exact layer equalities and all intermediate projections. The single external sleeve therefore does not merely give two isolated witnesses.


## 4. Exact-target implications and a selective switch

For an attained prime p≥5, SCC products leave only one nontrivial SCC, of count p; its size is less than p by the Busch minimum. At every node of that SCC's canonical substitution tree, a nontrivial gcd of the slot profile would propagate to p by §2. But the gcd divides the factorial of the node size, since F(b)=b!, so no such gcd is possible. Likewise the normalized profile-color group at every quotient must be trivial: its order would divide p and also the factorial of its quotient size. In particular actual rigidity is not a sufficient diagnostic for the prime lane. This is only a necessary condition, not a method of realizing every prime. More generally the **product** of all normalized color stabilizer sizes over the canonical tree divides every slot count: induct upward using the gcd at a child to retain each factor already found below it. The product is also a multiple of the actual automorphism-group order, since at every prime node the isomorphism-color stabilizer is a subgroup of the profile-color stabilizer; at a linear node both are trivial.

For the other hard lane 7p with prime p≥11, the known absence of count 7 and the SCC product law again leave one nontrivial strong component. It has order <p since f(p)>7p (base orders 11,12,13, then the three-step floor recurrence). The tree product just described can therefore only be 1 or 7: it divides 7p, and each local quotient group has permutation degree <p so contributes no factor p. At most one node can have any normalized profile symmetry at all, and in that event its order is 7. This is a scope-specific filter for candidate constructions on the 7p lane, not an exclusion of those counts.

There are exact **one-visit switches** which isolate one coefficient of the continuation. The transitive triangle and directed triangle have profiles (1,6,6) and (3,6,6). Also the order-five masks (24,1,3,6,14) and (18,8,3,5,14) have profiles (11,78,198,240,120) and (13,78,198,240,120). Replacing either pair as a module in any fixed physical quotient alters H by exactly twice the total receiver weight of its one-visit schedules; all repeat schedules are unchanged. This follows directly from §1. The multiplier is odd modulo two: schedules with any repeated label in another module carry even weight by k! divisibility, and the schedules with all modules visited once give H(Q) times a product of odd counts. Thus a single such switch always changes H by 2 mod 4, and two opposed same-size-of-step switches cannot give an exact difference of 2 in the same fixed receiver.

The missing construction is a receiver that puts these selective corrections at the prescribed absolute target. Ordinary quotient symmetry cannot be used as a substitute for that step, and merely breaking it with distinct actual module isomorphism classes does not eliminate the arithmetic divisor.

## Replay

`python3 verify_composition.py` independently enumerates constituent cover words using cuts at all invalid adjacent arcs, then checks the cyclic formula against a literal tournament subset DP in small cases. The all-size assertions use the schedule bijection and quotient arguments, not extrapolation. Run `python3 verify_modular.py` for the exact modular layer repetition and `python3 cutoff_bypass.py` for the enriched joint palette.

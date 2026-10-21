# Every fixed symmetry-slack layer has a rational generating function

Date: 2026-10-21. Status: PROVED in the stated scopes, with independent adversarial audit PASS and FINITE-EXACT controls. No historical-priority claim. All counts below are unlabelled tournament isomorphism classes, not Hamilton-path presentations or converse-merged classes unless stated.

## 0. Inheritance, advance, and connection contract

Inherited: `05-knowledge/results/continuation_width_symmetry_20261020/README.md`, Theorems 3–6 and Corollary 5C. For a tournament T, let f(T) be its minimum feedback-arc count and let s(T) be the number of vertices moved by at least one element of Aut(T). The nonnegative symmetry slack is δ(T)=3f(T)−s(T). The inherited bounded-exceptional-skeleton theorem gives an Aut(T)-invariant set B with |B|≤2δ, outside which every vertex is either fixed by Aut(T) or belongs to a homogeneous directed-triangle module that is an actual Aut(T)-orbit. Its proof also gives a representation by at most |B| singleton exceptional vertices and at most |B|+1 arbitrary ordered words of singleton/triangle atoms. That representation was not injective. Exact generating functions were proved only for δ=0,3,4; δ=1,2,5 are empty. The same packet explicitly leaves all-fixed-slack rationality open.

New proved advance: replace the nonunique exceptional presentation by the canonical partition into maximal modules belonging to the word class W=SEQ(K1,C3), then separate bounded short words from independent long words. This gives an exact finite compiler, proves bivariate rationality for every fixed slack, and proves the corresponding self-converse and converse-merged channels. It does not compute every finite configuration or produce an efficient general compiler.

Connection contract:
- Source: the inherited bounded-exceptional symmetry grammar and actual automorphism action.
- Target: a finite positive sum of rational word-language series.
- Map: take maximal W-modules, pass to their quotient, and retain the exact short words and long-slot set.
- Preserved: tournament isomorphism class, actual feedback, actual moved support, and converse.
- Lost only after scalar summation: module addresses and internal triangle phases; these can remain attached as sidecars.
- Decisive test: check maximality, weighted-feedback straightening, and actual quotient stabilizer. The diagonal C3 synchronization example at δ=6 is mandatory.

## 1. Notation and theorem

A module M of a tournament T is a nonempty vertex subset such that every vertex outside M either beats all of M or loses to all of M. Let W be the class of ordered joins of singleton tournaments and directed triangles, including the empty word. Its isomorphism classes are precisely words over the alphabet {1,3}; 1 has vertex weight one and triangle count zero, while 3 has vertex weight three and triangle count one. Write n(w),t(w),u(w)=n(w)−3t(w) for a word's vertex count, triangle count, and number of singleton letters.

Set

A(z,y)=1/(1−z−yz^3),
D_d(z,y)=sum over W words w with n(w)≤d of z^{n(w)}y^{t(w)},
L_d(z,y)=A(z,y)−D_d(z,y).

Thus L_d counts words of vertex length strictly greater than d. The empty word is included in D_d.

For d≥0 define

F_d(z,y)=sum over [T] with δ(T)=d of z^{|T|}y^{f(T)},
C_d(z,y)=the same sum restricted to self-converse classes,
M_d(z,y)=the same sum over converse-merged classes.

The empty tournament belongs to d=0. Then

F_0=A(z,y),
C_0=(1+z+yz^3)A(z^2,y^2),
M_0=(F_0+C_0)/2.

For every fixed d≥1 there is a finite, explicitly decidable set of configurations Γ_d and nonnegative integers a_c,b_c,ℓ_c with ℓ_c≤2d+1 such that

F_d(z,y)=sum_{c in Γ_d} z^{a_c}y^{b_c} L_d(z,y)^{ℓ_c}.   (1)

In particular F_d is rational and its denominator divides (1−z−yz^3)^{2d+1}. Empty layers have the zero series. Configurations are counted once up to the finite equivalence defined in §5; there are no unknown input tournament counts.

A finite subset Γ_d^sc consists of the configurations admitting an anti-isomorphism. Define the palindrome series and its long tail by

P(z,y)=(1+z+yz^3)/(1−z^2−y^2z^6),
P_{>d}(z,y)=P(z,y)−[vertex degree≤d]P(z,y).

Then

C_d(z,y)=sum_{c in Γ_d^sc} z^{a_c}y^{b_c}
 L_d(z^2,y^2)^{floor(ℓ_c/2)} P_{>d}(z,y)^{ℓ_c mod 2}.   (2)

Its denominator divides (1−z^2−y^2z^6)^{d+1}, and M_d=(F_d+C_d)/2. These are ordinary class series; Burnside's final factor 1/2 does not identify a Hamilton-path or marked carrier.

## 2. Canonical maximal W-modules

### Lemma 2.1: W is hereditary and closed under ordered join

This is immediate from its atom words: an induced subset of a C3 is empty, K1, the transitive two-vertex tournament, or C3. A transitive two-vertex tournament is the word 11.

### Lemma 2.2: overlapping W-modules have W union

Let A,B be modules with nonempty intersection. Nested cases are immediate. Otherwise write X=A\\B,Y=A∩B,Z=B\\A, all nonempty. The union is a module, because every outside vertex sees A and B uniformly and the common vertex makes those two directions agree. The module conditions force the three nonempty pieces into a transitive order, either X→Y→Z with X→Z or its reverse: for x in X,y in Y,z in Z, the fact that x lies outside B gives the same direction from x to y and z, and z outside A gives the same direction from x and y to z. These directions are uniform across the three sets by varying x,y,z. Each of T[X],T[Y],T[Z] belongs to W by heredity, so their ordered join T[A∪B] belongs to W.

Consequently the maximal W-modules are pairwise disjoint, cover V(T), and are intrinsic. They yield a unique quotient tournament Q and nonempty words (w_i) at its vertices, up to isomorphism of this decorated quotient. Conversely the decorated quotient reconstructs T exactly by substitution Q[W(w_i)].

### Lemma 2.3: a finite exact maximality test

For a proposed quotient Q and nonempty W words w_i, the displayed blocks are maximal W-modules if and only if there is no module I of Q with |I|≥2 whose inflated tournament lies in W. To prove the only nontrivial direction, a larger W-module intersects a displayed block and its union with that block remains a W-module by Lemma 2.2. Iterating over the finitely many intersected blocks produces a W-module that is a union of whole blocks. A union of whole blocks is a module exactly when its index set is a module of Q.

The condition that the inflation on I belongs to W depends only on Q[I] and the flags [w_i=1]. Its SCCs must all be singletons or directed triangles; every vertex in a directed-triangle SCC must be expanded by the one-vertex word 1. Indeed substitution into a strong quotient on at least two vertices is strong, and the only nontrivial strong member of W is the three-vertex triangle. Singleton SCCs may be expanded by arbitrary W words.

Thus exact maximality does not need the internal contents of any word longer than one vertex.

## 3. Feedback straightening and the actual support formula

### Lemma 3.1: arbitrary homogeneous modules straighten

For any substitution T=Q[T_i], with n_i=|T_i|>0,

f(T)=sum_i f(T_i)+C_Q(n),
C_Q(n)=min over linear orders π of V(Q) of sum_{i→j backward in π} n_i n_j.  (3)

For the lower bound, start with any vertex order of T. The cost inside each module is at least f(T_i). Choose one uniformly random vertex of each module and order the modules according to these representatives. For each pair i,j, the expected contribution n_i n_j times the indicator of a backward quotient arc equals exactly the number of backward cross arcs in the original order. Some representative choice therefore has weighted quotient cost no greater than the original cross cost. Conversely an optimal quotient order, with each module internally placed in an optimal order, attains the right side. This proves (3), including cyclic quotient hostiles.

Apply this to the canonical W partition. Let K=Aut(Q;(w_i)) be the actual stabilizer of the full word assignment, and let J be the set of quotient vertices moved by at least one element of K. Each W block has only its independent C3 rotations internally. Hence, writing t_i=t(w_i),u_i=u(w_i),

f(T)=sum_i t_i+C_Q(n),
s(T)=3sum_i t_i+sum_{i in J}u_i,
δ(T)=3C_Q(n)−U,   U=sum_{i in J}u_i.                 (4)

Crucially U is moved SINGLETON mass, not the entire mass of moved blocks. The latter is V=sum_{i in J}n_i and only supplies a bound. For a concrete hostile, C3[C3,C3,C3] has f=12, s=9 and slack27: its canonical words are (3,3,3), C=9, U=0 and V=9. The correct expression is3C−U=27; subtracting whole moved-block mass would incorrectly give18. This example is included in the independent audit.

### Lemma 3.2: moved canonical blocks are bounded

For δ(T)=d, use the inherited invariant exceptional set B with |B|≤2d. If a canonical W block is moved by some automorphism, every one of its vertices is moved. None can lie in a quiet actual triangle orbit outside B: that orbit is itself a W-module; it meets the block, and maximality plus Lemma 2.2 forces it to lie entirely in that block. But the orbit of a vertex of a moved block also meets a different, disjoint block. Contradiction. Thus every moved block is contained in B. Consequently

0≤U≤V≤2d,    C_Q(n)=(d+U)/3≤d.                      (5)

This is the only inherited bounded-skeleton input required after (4). It retains actual group orbits, not mere score classes, and allows synchronized exceptional actions.

### Lemma 3.3: canonical quotient size and free slots

Choose an optimal quotient order with weighted cost C≤d. It has at most d backward arcs. Let S be their endpoints, so |S|≤2d. Every vertex outside S is incident to no backward arc. In each gap of S in that order, the consecutive outside-S vertices form a transitive homogeneous module of Q. If a gap contained two or more vertices, the union of their W blocks would be a W-module, contrary to maximality. Therefore at most one outside-S vertex lies in each gap:

|Q|≤2|S|+1≤4d+1.                                    (6)

Any word of length >d cannot lie at a backward-arc endpoint because even one incident term n_i n_j exceeds d. Such long slots lie outside S, so their number ℓ is at most 2d+1. In particular the subtournament induced by the long slots is transitive.

## 4. Removing the last unbounded group ambiguity

For a fixed finite Q, let L be the long slots. Prescribe every word w_i at i outside L exactly, with 1≤n(w_i)≤d. The words at i in L may be arbitrary words with n(w_i)>d.

Call a quotient order admissible if every arc incident with a vertex of L is forward. Let C_0 be the minimum cost over admissible orders, counting only the short-short arcs with their prescribed masses; set C_0=∞ if there is no admissible order. Keep only configurations with C_0≤d. For every choice of long words, inadmissible orders have cost>d, while the minimizing admissible order costs C_0. Thus the actual weighted feedback cost is exactly C_0, independent of long-word contents and lengths.

The configuration stabilizer K_0 consists of automorphisms of Q preserving L setwise and preserving each prescribed short word. Since Q[L] is transitive, K_0 fixes L pointwise. It is therefore exactly the stabilizer of EVERY filling of the long slots: it cannot compare two long words or permute their positions. The moved set J lies among the short slots and is determined before any long words are supplied.

The maximality test of Lemma 2.3 likewise is independent of long contents: each long word has more than d≥1 vertices, so is not the one-vertex word. Actual slack is consequently the constant 3C_0−sum_{i in J}u(w_i).

This is the key finite-state step. It excludes fake configurations using the exact slack and maximality tests, rather than treating the inherited surjective grammar as an injective census.

## 5. The exact finite compiler and proof of (1)

For fixed d≥1:
1. Enumerate each tournament Q up to isomorphism with 1≤|Q|≤4d+1.
2. Choose a subset L of its vertices and a nonempty W word of vertex length≤d at every vertex outside L.
3. Apply Lemma 2.3's exact maximality test with long slots flagged non-singleton.
4. Compute C_0 from all admissible orders and require C_0≤d.
5. Compute the actual finite stabilizer K_0, its moved short slots J, and require 3C_0−sum_{i in J}u(w_i)=d.
6. Identify configurations under Aut(Q). Keep exactly one representative of each orbit.

All steps are finite and decidable. Theorem 3.3 proves completeness; §§2–4 prove soundness for EVERY long filling. Distinct retained configurations or distinct long fillings cannot give the same tournament: the canonical maximal-W quotient would identify them, and a stabilizer fixes every long slot. More explicitly, two transports to the chosen representative of a configuration orbit differ by a configuration stabilizer; because that stabilizer fixes every long slot individually, both transports give the same filling at the representative. Thus the correspondence is a literal bijection.

The term of a configuration is z^a y^b L_d(z,y)^ℓ where

a=sum_{i outside L} n(w_i),
b=C_0+sum_{i outside L}t(w_i),
ℓ=|L|.

Summing proves (1). No atom census beyond the stated finite explicit enumeration is presumed known, and the procedure is not asserted practical at large d.

## 6. Converse and palindrome word channels

Converse sends Q to Q^op and each word to its reversal; K1 and C3 are individually self-converse. This operation preserves vertex size, feedback and moved support. It induces an involution of the finite configuration orbits.

A tournament can be self-converse only if its configuration orbit is fixed. For such a configuration choose one anti-isomorphism τ of its quotient carrying prescribed short words to the reversals of their corresponding words and carrying L to itself. On the transitive tournament Q[L], τ must reverse the unique linear order. Every other such anti-isomorphism differs by a configuration automorphism, which fixes L pointwise. Thus the action on long slots is uniquely determined: floor(ℓ/2) pairs and, if ℓ is odd, one fixed central slot.

A paired pair of long slots is filled by w and reverse(w), contributing L_d(z^2,y^2). The possible middle long slot must carry a palindrome. Palindromes are v reverse(v), v 1 reverse(v), or v 3 reverse(v), giving

P(z,y)=(1+z+yz^3)A(z^2,y^2).

Removing those of vertex length≤d gives P_{>d}. These conditions are necessary and sufficient for an actual anti-isomorphism of the expanded tournament; they do not require τ itself to be an involution on the short blocks. They give (2). Since ceil(ℓ/2)≤d+1, its stated denominator follows. Burnside for converse then gives M_d=(F_d+C_d)/2.

## 7. Consequences and boundaries

At y=1 write p(z)=1−z−z^3 and let ρ be its positive zero. The ordinary exact-slack sequence has an eventual constant-coefficient recurrence annihilated by p(E^{-1})^{2d+1}; equivalently its reduced generating denominator divides p(z)^{2d+1}. The self-converse channel has denominator dividing p(z^2)^{d+1}.

For a nonempty d≥1 layer, let m_d be the maximum ℓ occurring in the finite compiler. Then the positive representation (1) gives a pole at ρ of exact order m_d: every top term has positive coefficient at ρ after multiplying by p^{m_d}. Therefore

[z^n]F_d(z,1) ~ c_d n^{m_d−1}ρ^{-n},  c_d>0,

with 1≤m_d≤2d+1. No explicit value of m_d or c_d for uncompiled layers is claimed. Nonempty layers are infinite because adding any outside W word preserves slack. The empty layer is excluded from this asymptotic statement. The unique dominant root follows as in the inherited W-word growth argument.

The results concern symmetry slack, not fundamental black degree. The two parameters are not interchangeable. They do not determine metagraph adjacency, bridges, connectivity, or Hamilton-path port allocations. The synchronized two-triangle δ=6 example remains a permitted finite quotient with diagonal C3, not two independent triangle atoms. The exact stabilizer test is what preserves it.

## 8. Other recovered connections and stopping reasons

1. Ordered-backward-graph trijection → fixed feedback rationality. A k-backward-edge order has ≤2k exceptional endpoints and ≤2k+1 transitive clone gaps; representative sampling straightens the weighted core. This was recovered independently here but is ALREADY INHERITED, explicitly credited in continuation_width_symmetry_20261020/incoming.md to defect_boundary_20261020/sparse_theory.md and commit602286dc6f. Its old linked file is absent in this checkout; no novelty claim is made. Preserved invariant: exact feedback via weighted quotient. Lost sidecar: optimal observer order and clone word positions. Decisive test: canonical clone quotient and actual weights, not the count of backward edges in one display.
2. Order-join metric product → path-colour rational sums. THM-2183 gives exact l1 cancellation, while THM-3166 and THM-3181 give multiplicative Q_T and fixed-depth recurrences. A distance- or path-profile sum over a finite SCC alphabet is a finite transfer-matrix series. This alone is a direct specialization of proved mechanisms, not enough novelty. The symmetry-slack result above differs because it proves a canonical bounded quotient and eliminates a genuinely open word-identification obstruction.
3. Marked cycle-plus-solo trijection → converse-reflection fixed loci. This remains a separate candidate, not needed above. Rotation-equivariant contraction alone does not justify converse equivariance: the spoke-to-cycle-edge replacement shifts the reflection gauge. An explicit reflection action on the original free-edge carrier is the decisive test. No speculative arithmetic or local-system analogy is promoted.

## 9. Explicit low-slack converse controls and the strong-family contrast

The inherited classifications of slack three and four immediately refine to

F_3(z,y)=y z^2 (A(z,y)−1−z) A(z,y)^2,
C_3(z,y)=y z^2 (P(z,y)−1−z) A(z^2,y^2),
F_4(z,y)=y^3 z^5 A(z,y)^2,
C_4(z,y)=y^3 z^5 A(z^2,y^2).

Here subscripts denote slack, not a directed triangle. The slack-three strong component is x→W→y→x with n(W)≥2. Its endpoints are intrinsic by the inherited score argument; feedback is t(W)+1 and moved support is 3t(W). Thus its slack is three for arbitrarily large W, and it stays strong. At the excluded boundary W=1 it becomes the triangle, gains quotient rotations, and has slack zero; W empty is transitive and also has slack zero. This is a concrete reason the new theorem is not the fixed-H finite-strong-atom proof: fixed symmetry slack contains unbounded strong components.

The converse of the slack-three central component replaces W by its reversed atom word. Self-converse outer words must be paired reversals, and the middle word must be a palindrome, proving its formula. The regular five-tournament is self-converse with feedback three, proving the slack-four formula.

Finite independent checks use `python3 verify_symmetry_slack.py 7` (omit 7 for the smaller order-six run). The full run passes 21,749 exact checks on all532 nonempty unlabelled tournaments through order seven; normal and optimized-mode transcripts byte-match. The checker constructs every unlabelled tournament through order seven by complete one-vertex extensions and exact score-compatible canonicalization (class totals1,1,2,4,12,56,456), then computes feedback by subset dynamic programming, actual automorphisms directly, every W-module, and the weighted quotient by all quotient orders. It checks every claimed structural equality, all 67 generated long-word replacement controls through order seven, the inherited exact slack0/3/4 counts, the zero-slack palindrome channel, the exact configuration-level converse criterion, and the unique long-slot reflection action. This is a finite audit, not a substitute for the all-order proof. It also retains the phase-locked two-triangle hostile with actual automorphism group C3, feedback4, moved support6, and slack6. In one orientation its canonical W blocks are three transitive pairs, so the canonical decorated quotient is C3[11,11,11]. Another order-six symmetric slack-six core has six singleton maximal W-modules; neither possibility is discarded by the compiler.


Independent review: `session23-slack-review/AUDIT.md` checks the inherited bounded-skeleton dependency and every main proof obligation. Its separate implementation passes4,576 exact assertions on139 actual tournaments,548 proposed W partitions including deliberately nonmaximal ones, and27 cyclic quotient weight profiles; named tests include long-word/converse controls through order15. This review is a scoped mathematical audit, not a formal verification or an unbounded census. No mathematical repair was required. Final changes after its first PASS only promote status and add the explicit representative-transport explanation and moved-singleton-mass hostile above.

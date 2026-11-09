# The stabilizer sidecar: identical scalar counts, different quotient gaps

**Status: PROVED; independent audit PASS.** All-order statements
below are elementary consequences of the cited OCF identity and action laws.
No completeness theorem for Chan–Pak is claimed.

## 1. Three different quotient statistics

Let T be a finite tournament, G <= Aut(T), P(T) its directed Hamilton paths
(regarded as ordered spanning words), and X(T) its two-coloured vertex-disjoint
directed odd-cycle packings, including the empty packing. A directed cycle is
identified under cyclic rotation, not reversal. Its colour is 0 or 1. G relabels
vertices and leaves colours unchanged.

The OCF identity (01-canon/theorems/THM-002-ocf.md, Grinberg–Stanley) says |P|=|X|=H(T).
The classical free action on spanning words (01-canon/theorems/LEM-003-aut-acts-freely-on-ham-paths.md) gives

    p_G(T) = |P/G| = H(T)/|G|.

For a finite action groupoid, define its mass and its coarse orbit count by

    m_G(X) = sum_[x] 1/|G_x| = |X|/|G|,
    c_G(X) = |X/G|.

Orbit–stabilizer proves the equality in the definition of mass, so

    m_G(X(T)) = p_G(T),
    c_G(X(T)) - p_G(T) = sum_[x] (1 - 1/|G_x|).              (1)

The right side is nonnegative. It vanishes precisely when the G-action on X
is free. Since the empty packing has stabilizer G, equality holds precisely
when G is trivial. In particular for every nontrivial G there is no
G-equivariant bijection P(T) -> X(T), despite the cardinality identity.
This is a direct obstruction, not a claim that no non-equivariant bijection
exists.

For G=Aut(T), inherited THM-4812 gives m_G(X) not in {7,21}; no new proof of
that inherited finite-certificate theorem is supplied. The session's separately
audited subgroup theorem, if promoted, extends the same conclusion to all G.
Equation (1) does not imply that c_G(X) excludes either value.

## 2. Exact reappearance of 7 and 21

Let R5 have vertex set Z/5 with i->j exactly when j-i is 1 or 2 mod5.
Its automorphism group is the translation group C5. Indeed an automorphism
fixing 0 must preserve its out-neighbourhood {1,2}; the arc 1->2 fixes these
individually, and likewise fixes the in-neighbourhood {3,4}. Translations
are automorphisms, so there are exactly five.

There are five directed triangles. Translation acts freely on them. There
are exactly two directed Hamilton cycles: all steps1, and all steps2. For
a proof of completeness, a length5 directed closed walk has steps1 or2;
its sum is divisible by5 and between5 and10, hence all steps are identical.
There are no two disjoint odd cycles on five vertices. Therefore X(R5) has

- one empty packing, fixed by C5;
- ten coloured triangles, in two free orbits;
- four coloured Hamilton cycles, each fixed by C5.

Thus H(R5)=15, p_C5(R5)=3, but c_C5(X(R5))=7. Its stabilizer histogram is
five orbits of stabilizer5 and two of stabilizer1, with mass5/5+2=3.

For C3, all three packings (empty and the two colours of its unique directed
cycle) are fixed by Aut(C3)=C3. Its path orbit count and packing mass are1,
but its coarse packing orbit count is3.

If A=>B is the order join, every directed cycle lies inside one factor and
every Hamilton path concatenates a path of A and a path of B. The ordered SCC
condensation has no nontrivial automorphism. In particular for strong A,B,
Aut(A=>B)=Aut(A)xAut(B). Packing actions and path actions are product actions,
so coarse counts, masses and path quotient counts all multiply.

Consequently C3=>R5 has H=45 and full automorphism group C3 x C5 of order15.
Its path quotient count and packing mass are3, while its coarse packing orbit
count is21. This is a genuine action-theoretic reason for both numbers to
reappear after changing the statistic, not an exception to Chan–Pak or THM-4812.

## 3. An all-order stabilization / nonstabilization contrast

For k>=0 let T_k=R5 => C3 => ... => C3 with k ordered C3 factors. Equal-sized
SCCs cannot be permuted because their condensation is a linear order. Thus

    H(T_k)=15*3^k,     Aut(T_k)=C5 x C3^k,
    |P(T_k)/Aut(T_k)|=3,
    m_Aut(T_k)(X(T_k))=3,
    c_Aut(T_k)(X(T_k))=7*3^k.                              (2)

Deleting the final C3 factor defines restriction maps on these particular
Hamilton paths (not on arbitrary induced subsets). On Hamilton-path orbit
sets these maps are bijections: the removed factor has exactly one orbit.
They give a constant three-element inverse system, up to canonical factor
identifications.

On coarse packing orbit sets the same restriction map has exactly three
points in every fibre. No transition or composite transition is injective;
the inverse limit is a seven-element set times {0,1,2}^N, hence uncountable.
The images at each fixed level stabilize immediately because the maps are
surjective. This emphasizes the distinction between stabilization of images
(the Mittag–Leffler condition) and eventual constancy. No claim about derived
inverse limits is needed here.

Even equal mass does not make the packing action groupoids equivalent as k
changes: the empty packing has stabilizer of order5*3^k, and every isotropy
group at level k-1 has order at most5*3^(k-1). Equivalences preserve isotropy
orders. Thus (2) is mass stabilization, not groupoid stabilization.

## 4. Connection contract and inheritance

Source: directed Hamilton words and the OCF packing set of the same T.
Target: free path quotient, coarse packing quotient, and packing action groupoid.
Map: G relabelling on each carrier; orbit projections. The OCF identity is an
equality of cardinalities, not an equivariant map between carriers.
Preserved: H before quotient, and H/|G| after weighted mass.
Lost: stabilizers under passage to coarse orbits or scalar mass.
Sidecar: stabilizer group/order and actual transition map.
Cheapest hostile: C3 already rules out an equivariant path/packing bijection;
R5 produces the missing value7 in the coarse statistic; an order join produces21.

This extends the quotient warnings in THM-4892 (the global colour-swap quotient
has relative-colour torsors) to the distinct vertex-automorphism quotient.
It reuses the explicit correction in
LEM-003-aut-acts-freely-on-ham-paths.md of the failed cycle-freeness step
preserved in the historical THM-212-paley-divisibility.md: the two arithmetic R5 Hamilton cycles
are fixed, not free. THM-212's final path-action proof and divisibility result
remain correct. No claim of external novelty is made.

Shao's paper motivates separating image stabilization from constancy; Zhang's
paper motivates retaining genuine gluing/automorphism data before taking a
coarse scalar. These are analogies of method only. Equation (1), the examples,
and the tower are proved here in finite combinatorics, not imported from
p-adic geometry. See the session's source pin memo for the exact paper scope.

Reproduce: python3 check.py and python3 -O check.py; outputs must agree.

## Full inherited routes

- 01-canon/theorems/THM-4812-presentation-fibres-7-and-21-are-absent-at-every-order-orbit-product-lemma-vertex-transitive-lemma-closed-window.md
- 01-canon/theorems/THM-4892-tournament-sites-plus-construction.md
- 01-canon/theorems/THM-212-paley-divisibility.md (historical failed sketch, not a dependency).
- 01-canon/theorems/LEM-003-aut-acts-freely-on-ham-paths.md

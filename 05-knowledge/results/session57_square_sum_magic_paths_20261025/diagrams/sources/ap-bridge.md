# AP biclique root bridge: exact certificate and controls

## Source and scope

Ajai Choudhry, *Two sets of integers such that all elements of the sumset of the two sets are perfect squares*, arXiv:2508.07806v1 (11 August 2025), Lemma 2 and equation (2.9): https://arxiv.org/html/2508.07806v1#S2 . The primary-source result identifies square-valued additive 3-by-3 matrices with semimagic square matrices by an explicit permutation; its converse supplies integer vertex potentials. Equation (2.9), reproduced in the notation below for checking, is the source of the semimagic permutation. The AP characterization, graph certificate, gauge specialization, and enumeration below are deductions checked here, not claims of priority. None establishes an integer magic square of nine distinct squares.

## Exact integer equivalence

Let a magic square of squares mean a 3-by-3 square with nine distinct positive integer square entries, equal row/column/diagonal sums. Its existence is equivalent to existence of two integer three-term APs A and B whose nine cross sums are distinct positive squares. Potential vertex labels need not initially be positive or globally distinct; a universal square scaling and gauge supplies those extra properties.

Every 3-by-3 magic array over the integers has its entry multiset

    { C+i*x+j*y : i,j in {-1,0,1} }.

For example, a valid magic placement is

    [ C+x,   C-x-y, C+y   ]
    [ C-x+y, C,     C+x-y ]
    [ C-y,   C+x+y, C-x   ].

The usual row/column/diagonal equations imply opposite entries add to twice the center and give this placement with integer x,y. Distinctness permits sign changes and interchange to select x>y>0; it additionally requires x != 2*y. Positivity is exactly C>x+y. Conversely, this placement has all eight line sums 3*C. Therefore a distinct-square AP sumset supplies the magic array simply by rearranging its entries.

### Universal positive, disjoint, induced realization

Given such C,x,y, set

    A = {1, 1+4*x, 1+8*x},
    B = {4*(C-x-y)-1, 4*(C-x-y)+4*y-1,
         4*(C-x-y)+8*y-1}.

Then A_i+B_j = 4*(C+(i-1)*x+(j-1)*y), for i,j=0,1,2. These are squares, with roots twice the original roots. All vertices are positive. The A vertices are 1 modulo 4 and the B vertices 3 modulo 4, so the parts are disjoint. Within either part, sums are 2 modulo 4 and are not squares. Thus the square-sum graph induced on these six positive integers is exactly K3,3. This is stronger than merely exhibiting nine edges in some larger square-sum graph.

Conversely, a six-vertex square-sum K3,3 with both parts APs and distinct cross-edge labels gives a magic square of distinct positive squares by the displayed rearrangement. The induced condition is unnecessary for this direction. Square scaling preserves the existence question but does not establish existence.

## Choudhry permutation and the precise two extra guards

For arbitrary potentials a1,a2,a3 and b1,b2,b3, the primary-source permutation is

    E = [ a3+b3, a2+b2, a1+b1 ]
        [ a1+b2, a3+b1, a2+b3 ]
        [ a2+b1, a1+b3, a3+b2 ].

All row and column sums equal T=sum(A)+sum(B). The main diagonal minus T is 2*a3-a1-a2. The other diagonal minus T is 2*b1-b2-b3. Hence E is fully magic precisely when a3 is the AP midpoint in A and b1 is the AP midpoint in B. Arbitrary ordering of AP parts does not suffice: the midpoint roles must be selected correctly.

For increasing indexed parts A0,A1,A2 and B0,B1,B2, use (a1,a2,a3)=(A0,A2,A1) and (b1,b2,b3)=(B1,B0,B2).

## Six independent edge-label constraints

Write L_ij=r_ij^2 for nine proposed square edge labels, i,j=0,1,2. Additive realizability is equivalent to the four relations

    L_ij-L_i0-L_0j+L_00 = 0,  i,j in {1,2}.

Once these hold, the two AP guards are

    L_00 - 2*L_10 + L_20 = 0,
    L_00 - 2*L_01 + L_02 = 0.

The 9-by-6 unsigned edge/vertex incidence matrix has rank 5 and one-dimensional gauge kernel (1,1,1,-1,-1,-1). The four cycle constraints have rank 4; all six displayed constraints together have rank 6. The square-entry problem is therefore the intersection of the coordinate-square locus with this three-dimensional linear magic subspace after permutation. No square-root sign restrictions are needed, because all constraints concern squared roots.

A certificate may supply roots r_ij, the four cycle residuals, the two AP residuals, distinctness/positivity checks, and one gauge t. Reconstruct

    A_i = t + L_i0 - L_00,
    B_j = L_0j - t.

All labels are recovered, and this is the general potential solution. After multiplying all L by 4 and arranging parts increasingly, choose t=1: positive, disjoint potentials follow exactly as above. A complete certificate should also record the actual six potentials, the nine labeled edges, the magic permutation, the eight line sums, both path orders below, and the shared edge. These are independently checkable without an existence search.

## Two Hamilton paths, optimal within the stated model

Use vertex orders

    P = A0,B0,A1,B1,A2,B2,
    Q = B1,A0,B2,A1,B0,A2.

The two five-edge paths cover all nine edges and share exactly A1--B0. Since a six-vertex Hamilton path contains only five edges, one path cannot cover nine distinct labels; two are optimal. Exact enumeration gives 36 unoriented Hamilton paths, six unoriented Hamilton cycles, and 36 unordered pairs of paths whose union covers all nine edges. Each path has exactly two covering partners. For a fixed path, its complement has components that are a length-three path and a single edge; precisely two endpoint-connecting edges from the original path yield a partner. The partner graph on the 36 paths is a disjoint union of six 6-cycles.

Six vertices are minimal ONLY among simple bipartite graphs carrying nine distinct edge labels: a bipartite simple graph on n<=5 vertices has at most floor(n^2/4)<=6 edges. At n=6, nine edges force K3,3. This does not prove minimality among general graphs or among other certificate models. Likewise, two AP guards are necessary in this fixed additive edge-coordinate presentation, not a universal lower bound on all possible formulations.

## Finite-field supplement, conditional and separate

Over an odd finite field F_q, the additive and AP identities remain valid, including characteristic three. Write A=a+{-x,0,x} and B=b+{-y,0,y}, and suppose the nine sums are distinct. Under gauge A->A+t, B->B-t, a cross-part vertex collision occurs exactly at t=(B_j-A_i)/2. Since {-x,0,x} is symmetric, the difference set B-A is a translate of the same nine-point grid as A+B. There are therefore exactly nine forbidden gauges and exactly q-9 gauges with six distinct vertex potentials. This count is conditional on nine distinct sums, not on nine distinct square roots.

In characteristic three, nine distinct labels mean x and y are linearly independent over F_3. The union of the two AP directions is then a two-dimensional affine F_3-plane. Merely describing two APs does not guarantee this independence. Positivity, ordering, and the integer mod-4 induced-graph argument have no finite-field analogue here. No claim that the finite-field graph is induced is made.

An independently enumerated nonzero-square witness over F_59 is C=1, x=11, y=14. Centered parts with gauge t=0 are A={48,0,11} and B={46,1,15}. The nine labels are {1,4,12,15,26,35,46,49,57}; all are nonzero quadratic residues and distinct. Exactly 50 gauges avoid cross collisions. This modular witness is not an integer solution or evidence of an integer lift.

## Reproduction

- verify.py: exact SymPy symbolic identities/ranks and exhaustive graph enumeration.
- results.json: all 36 covering pairs and summary values.
- finite_fields.py: exhaustive prime-field search with exact residue and gauge checks.
- finite-field-results.json: searched primes and the F_59 witness.

No repository or remote was modified.

## General induced-biclique normalization

For any positive-integer square-sum biclique with parts L and R, replace every left vertex u by 4*u+1 and every right vertex v by 4*v-1. Positivity is retained, cross sums become 4*(u+v), parts become disjoint modulo 4, and all within-part sums are 2 modulo 4. Thus the square-sum graph induced on the new vertices is precisely the prescribed complete bipartite graph. This proves equivalence between bare and induced positive-integer square-sum biclique existence after a permitted square scaling and gauge. It is not a height-preserving transformation and says nothing about finite fields. If the original parts are APs, the transformed parts remain APs.

## Exactly eight layouts for a real distinct magic multiset

For any nine distinct real numbers admitting a 3-by-3 magic arrangement, the center C is their mean. The eight other numbers occur in opposite pairs C±d, with four distinct positive gaps. After replacing the two parameters by their absolute values and interchanging them, write x>y>0. The four gaps are {x,y,x-y,x+y}; distinctness excludes x=2*y. The largest gap is x+y and the second largest is x. Therefore x and y are uniquely recovered from the unordered entry set: x is the second-largest gap and y is the largest minus the second largest.

The general magic-array parametrization then leaves exactly eight choices for the ordered signed parameter pair: signs on x and y and interchange. These give eight different layouts, all obtained from one another by the dihedral symmetries of the square. No other layout is possible. This is conditional on the multiset admitting a magic layout; an arbitrary nine-element real set generally admits none.

The real ordering argument does not transfer to characteristic three. For an affine two-dimensional F_3-plane of nine labels, the same parametrization has nine possible centers and 48 ordered independent direction pairs. It therefore has 432 distinct magic layouts, or 54 free D4 orbits. These assertions concern layouts of a fixed affine plane, not whether that plane lies in the square locus.

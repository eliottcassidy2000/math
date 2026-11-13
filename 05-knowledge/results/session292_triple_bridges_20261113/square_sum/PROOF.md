# Guarded reflection words and the birth of three-colour cycles

Status: **PROVED below; computational controls in `verify.py`; independent review pending.** No claim of Chan–Pak coverage or a new general square-sum Hamiltonicity theorem.

## Inheritance and the object

Let Q_N have vertices 1,...,N, distinct x,y joined exactly when x+y is a square. More generally allow targets in a set S. A target A is a colour and acts by the partial reflection R_A(x)=A−x: both endpoints must lie in [1,N] and differ.

Reflection composition and the two-colour forest are **INHERITED**, not new: [THM-4510](../../../../01-canon/theorems/THM-4510-sum-graphs-as-reflection-orbits-pythagorean-zigzags.md) and [its full note](../../procgen_sumgraph_20260926_reflection_orbits.md) already prove them, as well as tight three-target Hamiltonicity and primitive-Pythagorean spanning chains. [Crossroads](../../crossroads_family_20260926_squares.md) already proves first four-cycle birth at 46 and records first triangle at 30 by finite computation. [THM-4505](../../../../01-canon/theorems/THM-4505-small-graphs-encoding-arithmetic-square-sum-zigzag-fences-collatz-alphabet.md) supplies the first Hamiltonian window 15,16,17. We extend the guarded word calculus, not originate reflection algebra.

Closest mechanism: partial reflections. Hostile: a formal identity can force a repeated vertex or leave the positive interval. Corrected near miss: a closed colour word is not a simple cycle. Least-used coordinate: exact starting-vertex interval and collision set. The scope of novelty is the explicit three-colour cycle formula and its uniform consecutive-square birth law; no literature priority claim is made.

## 1. All-word affine guard theorem

For a word w=(a_1,...,a_l), put c_0=0 and c_i=a_i−c_{i−1}. Then its successive vertices are

v_i(x)=(-1)^i x+c_i, 0≤i≤l.

For a path, require 1≤v_i≤N for every i. Even i give 1−c_i≤x≤N−c_i; odd i give c_i−N≤x≤c_i−1. Intersect these integer intervals. Two indices of equal parity collide identically if their c_i agree, and never otherwise. Two indices of opposite parity collide only at the explicitly determined value x satisfying 2x=±(c_j−c_i). Delete those integral values. The result is exactly the set of starts realizing a simple path with this colour word.

For a simple cycle require l≥3, use i<l for pairwise-distinctness and impose v_l=v_0. The length restriction excludes a loop and traversal of one undirected edge in both directions. If l is even this requires c_l=0; if l is odd it requires 2x=c_l. The range guards still apply. These conditions are necessary and sufficient, since every edge sum is the prescribed a_i by construction and distinctness excludes loops. This is a finite exact compiler for each word, valid for every N; it is not a finite enumeration of the infinite square-sum graph.

The unguarded affine action discards the interval and collision conditions. Retaining those sidecars makes word concatenation exact: the second word's guard is pulled back by the first word's affine map, and cross-word vertex collisions must also be excluded. This last condition matters when asking for a simple path rather than a walk. No bounded-state claim is made for all words.

## 2. Exact three-colour hexagon theorem

Fix distinct positive integer targets A<B<C. Define p=B−A, q=C−B and K=A+q. The word ABCABC produces

(x, A−x, x+B−A, A+C−B−x, x+B−C, C−x),

then returns to x. Its six vertices lie in [1,N] exactly when

max(C−B+1,C−N) ≤ x ≤ min(A−1,N−B+A).

Set y=2x−K and

R=min(A+B−C−2, 2N−B−C+A).

The interval becomes |y|≤R with y≡K mod2. Vertices of the same parity position never collide because p,q>0. Comparing the nine cross-parity pairs gives exactly

y∈{0, ±p, ±q, ±(p+q)}.

Therefore a simple cycle of this colour pattern exists precisely for

0<|y|≤R, y≡K mod2, y∉{±p,±q,±(p+q)}.                 (1)

Replacing y by −y changes the starting position by three places along the same cycle. These are the only two starts with the prescribed oriented word ABCABC on that unoriented cycle: the cycle has exactly two A edges, and reversing direction changes ABC to ACB. Thus the exact number of distinct unoriented cycles with this pattern, ignoring any extra chords, is

0, if R<0; otherwise

floor((R+(K mod2))/2)
 − #{d∈{p,q,p+q}: d≤R and d≡K mod2},                (2)

where the brace is a set, so p=q is counted once. A cycle supported by these three distinct colours has a unique target triple, so different triples do not double count.

This is a piecewise-affine growth law with an eventual constant value. The breakpoints are precisely the positivity cap and the three collision distances. It classifies the entire infinite-N evolution of this word family from the three targets.

### Exact birth polynomial, not just a counting formula

Put D=A+B−C−2 and L=B+C−A. The admissible positive parameters are

Y={y:1≤y≤D, y≡K mod2, y∉{p,q,p+q}}.

Since L≡K mod2, each cycle has exact integral birth order b_y=(L+y)/2. Therefore, writing h_N for (2),

sum_{N≥0} h_N z^N = (1−z)^(-1) sum_{y∈Y} z^((L+y)/2).  (4)

The numerator is an interval of consecutive powers of z with at most three powers removed. In particular a fixed three-colour family introduces at most one new cycle per size, its missing births are exactly the eligible collision distances, and after its last birth it persists unchanged. This is a complete recurrence/generating-function description for this infinite-N family. Different triples can have simultaneous births; (4) is not a claim that the full square-sum graph has a rational global generating function.

## 3. A uniform consecutive-square family and the central triangle

Take A=t², B=(t+1)², C=(t+2)², t≥4. Then p=2t+1, q=2t+3,

A+B−C−2=t²−2t−5≥3,
B+C−A=t²+6t+5.

If t is even, K is odd, and the smallest positive permissible y is 1. If t is odd, K is even, and it is 2. Neither equals p,q,p+q. Equation (1) therefore proves that the first simple hexagon occurs exactly at

N_t=floor((t²+6t+5)/2)+1.                            (3)

This is an infinite constructive family with explicit vertex list from y=1 or2. At the first size there is exactly one such cycle. Examples: N_4=23, N_5=31, N_6=39, N_7=49, N_8=59. This is not a Hamiltonian cycle on [1,N_t].

At y=0 the hexagon repeats a triangle rather than being simple. For arbitrary A<B<C the triangle exists exactly when A+B>C and A+B+C is even; its three vertices are

(A+B−C)/2, (A+C−B)/2, (B+C−A)/2.

Its first size is (B+C−A)/2. Hence for consecutive squares with odd t≥5, a central triangle is born one size before the first hexagon; for even t the central parameter is a half-integer and no triangle exists for that triple. This is a literal parity-controlled bifurcation, not an analogy inferred from matching cardinalities.

## 4. First three-colour hexagon and first triangle in the square graph

If a square-colour hexagon exists with N≤22, then C≤2N−1≤43, hence C≤36. Positivity requires A+B≥C+2. For C≤25, the sum of the two largest smaller square colours is at most C (at C25 equality9+16=25), impossible. For C36 the only possibility is A16,B25, whose lower bound N≥ceil((C+B−A)/2)=23 is impossible. At N23, y=−1 or1 gives the single cycle

13,3,22,14,2,23,13.

So **23 is the first birth for the ABCABC three-colour family**, not asserted to be the first even cycle of any length. At N24 a second cycle is present, 1,24,12,4,21,15,1 (from x12); at N25 there are still exactly two. Targets49 are possible at25, but their two largest smaller squares25+36 give positivity while the size lower bound is30; all other triples are worse or fail positivity.

For any square-sum triangle, C≤2N−1. Below N30, C≤58, so C∈{4,9,16,25,36,49}. Targets≤25 fail A+B>C. At C36 only16,25,36 satisfy positivity, but their sum77 is odd. At C49 the positivity possibilities are16,36,49 and25,36,49; the first has odd sum101 and the second has first size30. Thus the finite triangle record in Crossroads upgrades to an elementary first-birth proof: the unique first triangle is {6,19,30} with colours25,36,49.

## 5. What transports toward tournaments, and what does not

The useful connection is an exact *guarded operation calculus*: source object → affine operation word → algebraic identity → admissible interface/positivity conditions → collision-free realization. It parallels Session291's distinction between a cluster exchange identity and a legal marked-endpoint construction. An identity alone has lost both domain and repeated-state data. A closed word at y0 demonstrates this loss using a triangle traversed twice.

A candidate tournament core move should similarly carry its source class, named ports, endpoint responses and positivity guards, and compose those guards rather than infer realizability from scalar recurrences. This is a methodological transfer with a precise hostile control, **not an isomorphism** between square-sum graphs and tournaments or a compiler of cycles into specified Hamiltonian-path counts. No Chan–Pak target is newly claimed attained.

The historical set {4,2,5,3,7,13} is separately identified in `borwein_ladder_owner_numbers_vsc_twelve_20261028.md` as {d+1:d divides12}; it is not a square-sum cycle in the displayed order (4+2=6). The actual first-three-colour hexagon's data above are not replaced by that label coincidence.

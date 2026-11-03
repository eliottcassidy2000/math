# Actual tournament hexagon covers and reversible local computation

**Status.** PROVED elementary constructions, supported by self-contained finite controls. The six seed profiles and class identifications are FINITE-EXACT. The arbitrary-width cover and compiler statements follow from the proofs. These are pinned labelled tournaments and their single-arc-flip class graphs, **not** the blue/black Hamilton-path-complement metagraph. No universality, spectrum-conjecture, or literature-priority claim is made.

## 1. Inheritance and the missing coordinate

The S167 ternary compiler supplies strong six-vertex classes A_0,A_1,A_2 with common complete path-cover profile and pairwise single-flip class adjacency. Its cellular clock was synchronous; one CA step could require several physical flips, and every triangle necessarily introduced relabelling holonomy. The S160 ordered engine correctly retained an external rule order and intermediate phase: forgetting those data did not yield an autonomous unmarked machine.

The new sidecar is much smaller than a general permutation frame. A **single bit per ternary cell** suffices, and it is carried by actual arc orientations. This gives a concrete regular tournament cover, in contrast to a cover whose added labels are only formal observation states.

Live concepts: full-profile invisibility; literal versus class moves; covers and holonomy; reversible local gates; finite rings and clock phases. The hostile example is a logical three-cycle: a labelled tournament cube is bipartite, so three literal arc flips cannot close. The strongest survivor is a six-cycle with a nontrivial deck bit.

## 2. An actual six-state switch

Use the lexicographic pair order 01,02,...,45, with a 1 bit meaning i→j for i<j. Put

    B_0=153, B_1=8345, B_2=8337,
    B_3=8369, B_4=177, B_5=185.

Indices on B are modulo six. Their successive differences reverse the arcs

    35, 04, 12, 35, 04, 12.

**Theorem 1 (literal hexagon).** The six literal tournaments induce C_6 in the labelled single-arc-flip cube. Every B_h is strong and has

    P(z)=29z+131z²+162z³+77z⁴+15z⁵+z⁶,
    F(r)=r⁶+12r⁴+16r².

There are exactly three isomorphism classes among them, with [B_h]=A_(h mod 3). The permutation

    p=(4,2,1,5,0,3)

maps B_h to B_(h+3). It is an involution, the product (0 4)(1 2)(3 5). Thus the quotient map h↦h mod3 is the ordinary connected double cover C_6→K_3.

**Verification.** `check.py` explicitly checks every arc difference, strong connectivity, every relabelling orbit, the action of p, and two independent profile counts. The two omitted corners of the same three-arc cube are masks145 and8377; both have H=33, not29. Thus precisely the six displayed corners have this profile within that cube. This is a finite seed certificate, not an exhaustive claim about all order-six tournaments.

The double cover is necessary for the three-cycle in this pinned model: a literal one-flip graph is bipartite, whereas K_3 is not. A fixed representative of each of three logical states cannot preserve all three transitions.

## 3. An exponential-sheet cover with an elementary deck group

For h=(h_1,...,h_N)∈(Z/6Z)^N take the ordinal sum

    U_h=B_(h_1)⊕...⊕B_(h_N),

on fixed, disjoint labelled blocks. All cross arcs point from the earlier block to the later block. Let x=h mod3.

**Theorem 2 (physical product cover).** The displayed labelled tournaments induce C_6^□N. Passing to their tournament isomorphism classes gives exactly the S167 graph K_3^□N, and

    q:C_6^□N→K_3^□N

is a connected regular 2^N-sheeted covering. Its full deck group is (C_2)^N, acting by h_i↦h_i+3 independently. Every displayed labelled tournament has order6N and

    F_U(r)=(r⁶+12r⁴+16r²)^N,    H(U)=29^N.

**Proof.** Two literal words are one arc apart exactly when one block traverses an edge of its induced hexagon. Each vertex has exactly two neighbors in each coordinate, mapping bijectively to the two logical neighbors, so q is a covering. All six blocks are strong, hence the SCCs of U_h are precisely its six-vertex blocks, in their uniquely ordered condensation. Isomorphisms cannot exchange positions; within each position the class is exactly h_i mod3. This proves the quotient description and fiber size. The displayed deck translations act freely and transitively on each fiber. A connected cover has at most its fiber size many deck transformations (one vertex image determines a deck map), so these are the full deck group. The path-slot product F_(T⊕V)=F_T F_V proves the profile identity. ∎

The cover is of this displayed class subgraph. It is not claimed to cover the entire tournament-class graph, nor to prohibit shortcuts through undisplayed classes. The S167 strong wrapper also transports the class dynamics to order30N+10 and common full profile. With fixed identical guard blocks, the same literal core flips give the pinned cover there; its Hamilton count is common but is not asserted to remain29^N.

## 4. Every controlled trit permutation has a literal reversible lift

A controlled single-cell gate leaves its controls fixed and replaces a target trit a by π_c(a), where π_c is a permutation of F_3 depending only on the logical control values c. The controls may include any designated neighboring cells, but not the target itself as an additional independent control.

For physical target state h∈Z/6Z define

    d(h,c)∈{−1,0,1},    d(h,c)≡π_c(h mod3)−h (mod3),
    h↦h+d(h,c) (mod6).

**Theorem 3 (gate lifting).** This is a reversible lift of the logical gate. It commutes with every deck-bit translation. If d is nonzero it reverses exactly one literal tournament arc; if d=0 it does nothing. A local gate stays local in the designated cell geometry.

**Proof.** The projection is the required logical permutation. For each target logical state, its two physical lifts map bijectively to the two lifts of its image, since replacing h by h+3 preserves d. Distinct logical inputs have distinct logical outputs. Hence the physical map is bijective. Fixed logical controls and unchanged physical controls give the controlled assertion. The hexagon proves the arc statement. ∎

For a layer whose targets are distinct and whose controls are never changed in that layer, the gates commute. It can be applied in parallel or in any sequential order. Every sequential nontrivial gate is a real edge; a parallel layer is generally **not one edge**. Locality means proximity of the designated cell indices, not a sparse neighborhood in the complete underlying tournament.

**Corollary 4 (exact orbit lifting).** Fix a repeated finite circuit of these gates, with unchanged cell positions. If a logical configuration has least period p under the circuit, each physical lift has least period p or2p. It has period p exactly when every cell undergoes an even number of nontrivial updates during those p rounds.

**Proof.** At a logical return, every physical displacement is0 or3 modulo6. The sum of the signed unit increments at cell i has parity equal to its number of nontrivial updates; since it is divisible by3, this parity identifies whether its displacement is0 or3 modulo6. The second logical lap repeats the increments, so its deck translation squares to the identity. No earlier physical return can project to an earlier logical return. ∎

The bound is2p, not2^N p: the number of sheets is exponential, but every deck element has exponent two. The logical gate a↦a+1 has period3 and lifts to the literal hexagon clock of period6, proving the doubled case occurs.

## 5. A constant-depth compiler for reversible cellular automata

Let F:F_3^m→F_3^m be a bijection on a periodic row. Suppose F_i depends only on cells within radiusr of i and its inverse G_i depends only on radiuss. Translation-equivariance may be imposed when F is called a CA. Use two tracks x,y, so N=2m.

Apply these five layers, evaluating each right side before changing its layer's targets:

    1. y_i ← y_i+F_i(x)
    2. x_i ← x_i−G_i(y)
    3. x_i ← x_i+y_i
    4. y_i ← x_i−y_i
    5. x_i ← x_i−y_i                       (all arithmetic in F_3).

The last three layers swap x and y. Layers1–2 are translations of each target by a function of the other track; layers3–5 are affine permutations of one target with the other fixed.

**Theorem 5 (clean reversible compiler).** The five-layer circuit is a permutation of all3^(2m) logical states. On blank workspace it acts exactly as

    (x,0)↦(F(x),0).

Every layer consists of commuting local single-cell permutations, with radius at most max(r,s). Lifting the gates as in Theorem3 yields a reversible circuit on actual labelled tournaments of order12m with common H=29^(2m), and common complete profile. Each round uses at most5m nontrivial literal arc flips, or five simultaneous layers. At round boundaries the *logical* workspace is blank; its physical deck bits need not be reset. A period-p F-orbit has physical period p or2p under this lifted circuit.

**Proof.** The first two layers give (x,F(x)) and then (0,F(x)), and the swap finishes the displayed formula. Every layer is reversible on the whole state space, independently of the blank condition. The general whole-space formula is

    C(x,y)=(y+F(x), x−G(y+F(x))).

Its inverse reconstructs x=d+G(c), y=c−F(x) from (c,d). Gate locality and all tournament statements follow above. ∎

The *five-layer* locality radius is max(r,s); the composed CA can have radiusr+s. A local inverse is a hypothesis, not a consequence claimed here for arbitrary finite-width bijections. This is a compiler for already specified reversible local dynamics, not a proof of computational universality.

There is a sharp irreversibility boundary. Let E be injective and assign exactly one physical codeword to each logical state. An injective physical map Q cannot satisfy Q(E(x))=E(F(x)) for a noninjective F: two different x with the same F(x) would collide. This obstruction assumes that the output is the prescribed codeword, including every workspace and deck bit. Allowing different deck bits for the same logical output can retain erased information and therefore lies outside that clean-encoding requirement. A reversible simulation of an irreversible CA must retain such information, allow a nonclosed representation, or weaken its simulation requirement. A profile-invisible encoding alone does not remove the obstruction.

## 6. An interacting radius-one clock with no physical period penalty

A particularly explicit family avoids any need to compute a local inverse for F. On two tracks of a length-m ring let S be the cyclic shift (Sb)_i=b_(i+1), and put

    R(a,b)=(b, a+b+Sb).

Its inverse is

    R^−1(c,d)=(d−c−Sc,c),

so both maps have radius one for **every** m. Implement R by the local shear

    a_i←a_i+b_i+b_(i+1),

followed by the three swap layers3–5 above. This is four local layers and at most4m nontrivial literal flips after lifting.

**Theorem 6 (exact physical clocks).** If m=3^e, e≥0, then R has exact order8m and

    R^(4m)=−I.

The state (δ_0,0) has least period8m. The four-layer lifted circuit on the literal tournament family has the **same exact order8m**, and every physical lift of (δ_0,0) has period8m. These are real pinned tournament states of order12m, all with the same complete profile and H=29^(2m).

**Proof of logical order.** Work over

    A=F_3[z]/(z^m−1)=F_3[n]/(n^m),  n=z−1,
    M(z)=[[0,1],[1,1+z]].

At z=1, M(1)^4=−I and M(1) has order8. Consequently M(z)^4=−I+nD(z) for a matrixD overA. The two summands commute, so taking the3^e-th power gives M^(4m)=−I, and then M^(8m)=I.

For exactness, M(z)^8−I=nB(z), where a direct differentiation at z=1 gives

    B(1)=[[1,2],[2,2]],    det B(1)=1 in F_3.

Thus B(z) is invertible over the local ringA. If m≥3 and q=m/3, then

    (M^8)^q=I+n^q B(z)^q≠I,

because n^q≠0 and B(z)^q is invertible. Therefore M^8 has exact orderm; together with the order8 reduction modulo n, M has order8m. The same argument works on v=(1,0): its reduction has period8, and B(z)^q v has a nonzero reduction, so n^q B(z)^q v≠0. The vector1 corresponds to δ_0 under the shift-module identification. For m=1 the claim is the explicit2×2 calculation.

**Proof of physical order.** Each gate in the four-layer circuit is linear overF_3. The logical state after4m rounds is the negative of its original state. Therefore the second half of the8m-round logical microtrajectory is the termwise negative of the first half, including intermediate gate states. A signed physical increment is the representative in {−1,0,1} of target minus source. Negating both states negates that integer increment. The second half's physical increments cancel the first half at every cell. Every literal lift consequently returns after8m rounds. The logical basis orbit proves that neither its lifts nor the global physical permutation have smaller order. ∎

The calculation B(1) is checked by exact polynomial matrix arithmetic in `check.py`, not by floating-point fitting. Controls run through m=81, including all intermediate physical arc flips. The proof, not those controls, gives arbitrary e.

This is a genuinely interacting local rule: the shear reads its neighboring b-cell. It is not the earlier single-track I+S clock. Its inverse remains radius one at every ring size, although the special exact-period formula is only asserted at powers of three.

**FINITE-EXACT off-family controls.** For the same basis state and four-layer circuit, the logical/physical periods at widths2,4,5,7 are respectively8/16,40/80,328/328,2920/2920. The first two are hostile controls against extending the no-doubling conclusion to every width. The last two give substantially longer concrete clocks (width7 uses84 tournament vertices), but these four finite observations are not an infinite exponential-period theorem. `check.py` checks first return, not merely return at the advertised time, and checks the deck-bit parity criterion at the first logical return.

## 7. Scope, hostile controls, and the useful structural conclusion

The preserved object is the complete path-cover profile. The lost data include the logical CA state, physical deck bits, and clock phase. The sidecar consists of labelled cell locations, their row-neighbor relation, and the finite layer schedule. Every recorded micro-move reverses a real arc, but one synchronous round changes many arcs and is not a single flip-class edge or a blue/black metagraph line.

In particular:

- The triangles cannot be realized by three fixed labelled representatives; the explicit double cover repairs this.
- The two excluded three-cube corners have H=33 and are hostile profile controls.
- Irreversible rules do not become reversible merely because all profiles agree.
- Equal full profiles do not determine the phase, period, or gate history of even these explicitly realized reversible dynamics.
- The general lift can double a period; the linear anti-period mechanism proves why the interacting8m clock avoids doubling.
- A fixed finite ring is a finite machine. Arbitrarily large rings and growing periods do not, by themselves, prove Turing universality or unbounded memory at fixed tournament order.

The exact new bridge is therefore

    reversible local gates → C_6^□N literal arc dynamics
                           → K_3^□N class dynamics
                           → one constant complete-profile fiber.

## 8. Reproduction and sources

Run from this directory:

    python3 check.py
    python3 -O check.py

Both write byte-identical `certificate.json`. The program uses only the Python standard library and explicit checks that remain enabled under optimization. It independently verifies the seed profiles by subset/set-partition DP and permutation/cut enumeration; all720 relabellings; every controlled single-trit permutation; finite product covers; exhaustive small compiler states; actual physical gate transitions; the polynomial order calculation; and basis clocks through width81.

Inherited, exact sources:

- [S167 ternary profile switch](https://github.com/eliottcassidy2000/math/blob/596ae1ea47f640689304af59237e5f737c75c0b8/05-knowledge/results/session167_recurring_structures_20261103/ternary/TERNARY_CONSTRAINTS.md), Theorems1–2 and5: the seed class triangle, profile product, strong wrapper, and the earlier synchronous clock.
- [S160 ordered metagraph engine](https://github.com/eliottcassidy2000/math/blob/596ae1ea47f640689304af59237e5f737c75c0b8/05-knowledge/results/session160_ordered_programs_20261103/metagraph/NOTE.md): actual complement-line moves with a retained controller and intermediate phase. Its native lines are different operations from the flips in the present packet.
- [S167 octonion transport](https://github.com/eliottcassidy2000/math/blob/596ae1ea47f640689304af59237e5f737c75c0b8/05-knowledge/results/session167_recurring_structures_20261103/octonion/NOTE.md): the warning that an abstract cover need not be physically realized by tournament states. The present hexagon cover supplies such a realization for a different, explicitly specified base.

Open directions: find smaller strong wrappers that keep the physical cover; classify profile-preserving local switches beyond this hexagon; determine which controllers can themselves be encoded using only intrinsic unmarked tournament data. None is settled here.

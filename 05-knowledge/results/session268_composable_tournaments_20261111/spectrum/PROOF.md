# Physical scalar tomography and exact cyclic constructors

**PROVED (self-contained):** scalar tomography, arithmetic inversion, the affine lift, and the sparse multiplicative-closure bound below. The ordinary path-cover substitution congruence is inherited from this session's `../modular/PROOF.md` section 1. Computational results in `results.json` are FINITE-EXACT controls. Chan–Pak remains OPEN. No literature-priority claim is made for the elementary cyclic blowup formula.

## 1. Definitions and recovered mechanism

Let q_r(T) be the number of ordered lists of r nonempty vertex-disjoint directed paths covering T, with arbitrary endpoints. Thus q_1=H, q_n=n!, and q_r=0 for r>n. Put q_0(empty)=1, and q_0(T)=0 for nonempty T. For a transitive tournament TT_a, q_r(TT_a)=F(a,r)=r! S(a,r), because distributing labelled vertices among ordered nonempty boxes determines each increasing path uniquely.

Write C(A,B,C) for the tournament A→B→C→A with the internal tournaments unchanged. All three modules are nonempty. A Hamiltonian path decomposes uniquely into maximal runs in modules. Its run labels follow the directed cyclic order; the chosen ordered covers inside the modules determine the path uniquely. Consecutive runs cannot have the same label. This is the preserved predicate and exact composition map; a scalar H alone loses higher cover coefficients.

For arbitrary nonempty A,B,C, their run counts in a cyclic word are either (k,k,k), a cyclic permutation of (k+1,k,k), or a cyclic permutation of (k+1,k+1,k). Therefore

    H(C(A,B,C)) = sum over k>=0 of
      3 q_k(A)q_k(B)q_k(C)
      + q_(k+1)(A)q_k(B)q_k(C)
      + q_k(A)q_(k+1)(B)q_k(C)
      + q_k(A)q_k(B)q_(k+1)(C)
      + q_(k+1)(A)q_(k+1)(B)q_k(C)
      + q_(k+1)(A)q_k(B)q_(k+1)(C)
      + q_k(A)q_(k+1)(B)q_(k+1)(C).                 (1)

There are three starts for the equal run-count case, and exactly one start for each of the other six count patterns. This proves (1), including its multiplicities, without assuming Hamiltonian paths remain contiguous in a module. Naively replacing each module by H would be false.

## 2. Main theorem: physical tomography of the entire cover response

For d>=1 define the physical scalar experiment

    P_d(T)=H(C(TT_d,TT_d,T)).

Set F(d,r)=0 outside 1<=r<=d. Formula (1) gives

    P_d(T)=sum_{r=1}^{d+1} K(d,r) q_r(T),
    K(d,r)=3F(d,r)^2+2F(d,r)F(d,r+1)+F(d,r+1)^2
           +F(d,r-1)^2+2F(d,r-1)F(d,r).           (2)

In particular K(d,d+1)=(d!)^2. All coefficients are nonnegative integers and depend only on d. The highest coefficient comes from a path starting and ending in T, with d runs in each transitive module and d+1 runs in T. No higher cover degree is physically possible.

**Tomography theorem.** For every n-vertex nonempty tournament T, the n ordinary integers H(T),P_1(T),...,P_(n-1)(T) uniquely determine its entire ordered-cover vector. Explicitly,

    q_1=H(T),
    q_(d+1)=[P_d(T)-sum_{r=1}^d K(d,r)q_r]/(d!)^2.  (3)

The divisions are exact on physical inputs. The largest probe has 3n−2 vertices. Computing the coefficient triangle and decoding uses O(n^2) arithmetic operations; bit costs depend on the output integers. This is a physical readout, not a claim that counting each probe from a raw adjacency matrix is polynomial time.

**Known-order optimization.** For n>=3 the last two cover coordinates are universal: q_n=n! and q_(n−1)=(n−1)! binom(n,2). The second identity chooses the unique two-vertex path and orders it with the n−2 singleton paths. Therefore H(T) plus only P_1,...,P_(n−3) reconstructs all remaining coordinates. There are n−2 measured scalars including H; no enlarged probe is needed when n=3, and the largest required tournament has3n−6 vertices when n>=4. At n=1,2 the whole response is universal. The unoptimized n-scalar version above remains useful when testing normalization or when treating arbitrary formal cover vectors.

**First-difference theorem.** If two tournaments have equal q_1,...,q_(r−1) and different q_r, then C(TT_(r−1),TT_(r−1),T) separates their H values by exactly ((r−1)!)^2 times the q_r difference. For r=1 the unmodified H already separates. At most 2r−2 auxiliary vertices suffice. The assertion does not require the input orders to agree.

**Optimal observation-budget theorem.** The cost2r−2 is minimal among *all* one-hole pure tournament-substitution contexts, not just cyclic probes. Indeed T remains a module in the completed context: every exterior vertex either beats every vertex of T or loses to every vertex of T. Between two successive maximal T-runs in a Hamiltonian path there must therefore be at least two exterior vertices. A single exterior vertex cannot both receive an arc from T and send an arc back to T. If m vertices are outside T, there are at most floor(m/2)+1 T-runs. Grouping paths by their external run pattern expresses the output H as a nonnegative integer linear combination of q_1,...,q_(floor(m/2)+1). Coefficients depend only on the context, not on the internal arcs of T.

Consequently equal q_1,...,q_k imply equal outputs for every context with at most2k−1 exterior vertices. A first difference at q_r is invisible below2r−2 exterior vertices, and the cyclic sensor above detects it at exactly2r−2. The context-budget filtration is therefore sharp. For comparison, arbitrary pointed arc-seam interactions are not pure substitution contexts and may break module homogeneity; this lower bound does not apply to them. Hole duplication is also excluded.

**Full abstraction for substitution.** Two nonempty tournaments have the same full cover vector iff no one-hole pure tournament-substitution context can distinguish them by H. The forward implication is the exact substitution congruence in `../modular/PROOF.md` section 1, iterated through the context tree. For the reverse implication choose the least differing q_r and use the explicit context above. Different orders are covered by zero-extending q; q_n=n! ensures a first difference exists. Thus no information in q can be discarded universally if every physical substitution context must be predicted. This is behavioral minimality, not a statement that q cannot be compressed by a nonlinear encoding or that each coordinate is algebraically independent on physical tournaments. Graph-isomorphism information is still lost by q.

A sibling's n=6 exact witness has qA=(9,150,756,1704,1800,720) and qB=(9,126,684,1656,1800,720). Their first scalar probe is 177 versus153, because P_1=3q_1+q_2. Higher probes distinguish them by the explicit coefficient formula; scalar equality is not a valid merge rule for compositional enumeration.

## 3. Exact affine families from an arbitrary seed

Let Q(T)=q_2(T). Specializing (1) with a singleton gives

    H(C(TT_1,TT_b,T))=(2^b+1)H(T)+(2^b−1)Q(T)
                         =(H+Q)2^b+(H−Q).          (4)

Every b>=1 gives an actual strong tournament on |T|+b+1 vertices, regardless of whether T itself is strong. Strongness follows by travelling through the other two nonempty modules. Given a seed vector and target N, membership in this one-parameter family reduces to checking whether (N−H+Q)/(H+Q) is an integral power of two >=2. Construction emits the literal cyclic adjacency blocks, with no enumeration of candidate tournaments.

If two seeds have the same H and different Q, the scalar gap after this lift is (2^b−1)ΔQ. The six-vertex pair above therefore gives infinitely many equal-order strong pairs with counts 159*2^b−141 and 135*2^b−117, separated by24(2^b−1), using b+7 vertices. This is a physical exponential amplifier of previously hidden state. It supplies actual new outputs but does not imply dense range or universal arithmetic digits.

## 4. Singleton transitive atoms: exact inverse and chirality

Taking T=TT_c in (4) gives the basic infinite strong family

    h(b,c)=2+(2^b−1)(2^c−1),  b,c>=1.              (5)

The vertex cost is b+c+1. Among this family alone, the unordered pair {b,c} is uniquely determined by h. Here is a direct proof that does not rely on factoring arbitrary integers. Put M=h−2 and suppose b<=c. If b<c, v_2(M−1)=b. If b=c>=2, v_2(M−1)=b+1. Square presentations are unique by monotonicity. Two nonsquare presentations have the same minimum exponent from the valuation, then the same maximum. A square of 2^b−1 cannot equal an unequal product: its smaller exponent would be b+1, making both factors larger than 2^b−1. The exceptional b=c=1 is M=1. This proves injectivity.

An implementation may simply try b<=c with b<=log2(M+1), test divisibility by 2^b−1, and test whether the quotient plus1 is a power of two. It uses O(log h) integer tests. The valuation proof permits a constant-size candidate list but is not needed for the safe decoder.

The three maximal proper strong modules of this tournament are its three transitive blocks, the usual unique strong modular decomposition with quotient C3. Isomorphisms act on the block-size triple by cyclic rotation. Consequently the count h(b,c) gives one self-converse core class when b=1 or b=c, and two converse core classes when 1<b<c. These are exact multiplicities within this family, not all tournament classes at that h. The only nontrivial automorphism is the order-three rotation at b=c=1; every other core is rigid. The class GF for arbitrary singleton padding and ordered SCC products is thus exact after assigning each scalar atom weight w(h)z^(b+c+1), with w=1 or2 as above. Scalar atoms may be composite integers; they are not claimed to be irreducible in the full realizable-count monoid.

## 5. Quantitative obstruction: even arbitrary products are very sparse

Let M be the multiplicative monoid generated by every h(b,c) in (5). Then

    #{m in M : m<=X} <= exp(C (log X)^(2/3))         (6)

for an absolute constant C, for X>=3. In particular this entire infinite atom bank, even closed under arbitrary SCC composition, has X^o(1) distinct scalar counts and cannot prove Chan–Pak.

Proof. Put j=b+c. For j>=4, h>=2^(j−2)>=2^(j/2); j=2,3 satisfy the same bound directly. A product<=X therefore admits a generator multiset with total j-weight at most L=2log_2 X. There are at most j generator choices of weight j. Counting such multisets overcounts distinct products, which is safe. For every t>0 their number is at most

    exp(tL) product_{j>=1}(1−exp(−tj))^(−j).

Its logarithm is at most tL+ζ(3)/t^2, since
sum_{j>=1}j exp(−uj)=1/[4sinh²(u/2)]<=1/u². Choose t=(2ζ(3)/L)^(1/3). One may take exponent 3*ζ(3)^(1/3)*L^(2/3)/2^(2/3). This proves (6), even counting chiral versions separately. A sharper constant is not needed.

This diagnoses the compositional debt precisely: large explicit atom inventories and multiplication do not supply value-changing arithmetic interactions. The theorem concerns only this bank, not all SCC compositions of arbitrary strong tournaments. More generally O(j^(d−1)) generators in logarithmic weight j imply an exp(O((log X)^(d/(d+1)))) upper bound by the same partition argument.

## 6. Scope and open construction debt

The closest inherited mechanism is the full path-cover substitution algebra. The hostile is H=13, absent from (5) and its product closure: 13 is prime and 13−2 is not a product of Mersenne factors. It already has a two-apex/crossing-arc realization in inherited THM-4878, so this is a grammar limitation, not a forbidden count. The corrected near miss is treating the scalar count as a closed module state. The least-used useful sidecar here is physical scalar tomography of higher ordered covers.

The live concept board is cyclic runs, Stirling surjections, physical readout, arithmetic decoding, multiplicative partition growth, and hidden-state amplification. The source-to-target map preserves every path-cover coefficient precisely when the entire vector is kept; scalar and isomorphism quotients lose different data. A finite alphabet providing target-directed H-changing interactions is still required for Chan–Pak. No universal compiler or convergence of a target search is established here.

Reproduce: `python3 05-knowledge/results/session268_composable_tournaments_20261111/spectrum/replay.py`. No inherited conjecture is used in these proofs. The standard substitution congruence is a session sibling dependency, and tomography itself is proved directly above.

## 7. Rooted seam interaction creates two decodable arithmetic surfaces

Use the sibling's proved one-back-arc seam in `../automata/PROOF.md` section 1. All arcs run from A to B except b→a. Its scalar law is

    H(A seam B)=H(A)H(B)−e_A(a)s_B(b)+L_A(a)R_B(b),

where L counts an arbitrary path followed by a path starting at a, and R counts a path ending at b followed by an arbitrary path. Empty arbitrary paths are allowed. For A=C(TT_1,TT_b,TT_c) rooted at its singleton,

    e_A=s_A=1,   L_A=1+(2^b−1)2^c=h+(2^b−2).

Indeed a root-start path either is the singleton, or uses a nonempty subset of the B module and any subset of the C module; its complement is transitive. This is L, not R: orientation matters.

The sibling's literal five-vertex targets B0 and B1 have masks8 and10, respectively, in lexicographic pair order01,02,03,04,12,13,14,23,24,34, with bit1 meaning i→j. At root0 both have H=9,s=1, whereas R=15 and11. Consequently the two exact outputs are

    N0(b,c)=24h+15(2^b−2)−1
           =24(2^b−1)2^c−9*2^b+41,
    N1(b,c)=20h+11(2^b−2)−1
           =20(2^b−1)2^c−9*2^b+37.                (7)

Their difference is4L_A, so replacing one fixed scalar-invisible module yields exponentially growing scalar differences without changing the order b+c+6. These are genuine nonmultiplicative interactions, not products of scalar counts. For example b=1 gives the two rays24*2^c+23 and20*2^c+19. For general b, an exact membership compiler tries b, divides N+9*2^b−41 by24(2^b−1) (or N+9*2^b−37 by20(2^b−1)), and tests for a power of two2^c with c>=1. It needs O(log N) integer tests: from (7), N_i>=2^b for b,c>=1, so b<=floor(log2 N) is safe. Emitting the adjacency matrix takes O((b+c)^2) bits. Values outside these surfaces require another constructor; no completeness is asserted.

The surfaces remain sparse in isolation. Their role is to exhibit an explicit controlled interaction between two retained port states and a decoded exponential atom family. The sibling supplies further two-sided mixed-switch interference; we do not claim (7) preserves the input port slice under iteration. That missing state-return property is exactly what prevents reading these equations as a full radix compiler.

## 8. Target-directed minimum-cost compiler and actual size gains

`compile_target.py` defines an explicit infinite generator language, with three optional layers:

1. the singleton transitive atoms (5), both parameters positive;
2. the two five-vertex seam surfaces (7), both parameters positive;
3. the two arbitrary-seed lift rays from (4), using literal six-vertex seeds
   A=TT_1→C(TT_1,TT_1,TT_3) and B=C3→C3.

The seeds' vectors are exactly the two vectors in section2, and the executable derives them from the literal graphs. It does not merely insert a claimed vector into a witness. The generated tournament language consists of all ordered SCC sums of these generator tournaments. Every generator is strong: cyclic lifts are strong, and a single reversed inter-block arc joins two strong blocks into a strong tournament. Both five-vertex seam targets contain the directed cycle0→4→3→2→1→0, so their strongness is literal. Different generator presentations can have the same scalar; the compiler retains a least-size one. This is safe only because subsequent composition in this compiler is SCC concatenation, not an arbitrary context.

For target N, enumerate all generator parameter pairs up to bitlength(N) and discard values>N. Every relevant parameter is in this range by the explicit formulas. There are O(log²N) candidate generators. Define D(1)=0 and

    D(m)=min_{generator g dividing m} [size(g)+D(m/g)],

with infinity for an empty minimum. All generators exceed1, so the recurrence terminates. Memoization visits only divisors of N and yields the exact minimum vertex count in the declared language. No integer factorization oracle is needed: only divisibility tests against the explicit generator bank. There are O(τ(N)log²N) candidate divisibility tests and comparisons, excluding big-integer bit costs and writing the output. No polynomial-in-logN worst-case assertion is made for τ(N). Target1 is emitted as a singleton, separately from the empty product used internally.

The literal witness is the ordinal sum of the selected generator graphs. Each primitive has at most4log_2(g) vertices: for basic atoms this follows from h>=2^((b+c)/2); the seam counts are at least their base h and at least59, absorbing the six-vertex overhead; the two lift rays have counts at least2^(b+6), absorbing their seven-vertex overhead. Summing over factors gives at most4log_2(N) vertices for every N>1 this compiler succeeds on. This is a convenient uniform upper bound, not a sharp bound on all tournament realizers.

Executed exact optimization results (`compiler_results.json`) are:

- H153: minimum10 vertices in the atom-product language, still10 with seam generators, but8 with hidden-seed lifts.
- H177: no atom-product witness;11 vertices after adding seams,8 after adding lifts.
- H495:14 vertices using atom products,13 with seam generators,9 with the lift generators.
- H27081=153*177:20 vertices in the seam-extended language versus16 after adding lifts.
- H59 and H71: new8-vertex seam generators, absent from the original atom monoid.
- H393 and H1000001: no witness in this language. This is not a nonrealizability claim.

The H153,H177,H495,H27081 full target searches visit6,4,12,16 memoized divisor states respectively. Their minimum sizes are relative to these three explicitly nested languages, not global tournament minima. Every returned witness of order<=16 in the measurement set was independently recounted from its adjacency matrix by the subset Hamilton-path dynamic program. At10^6 the banks contain94,307,330 distinct scalar generators in the three layers.

This is actual constructive efficiency and a measurable size benefit from retaining hidden cover state. It does not solve the remaining arithmetic problem: the enlarged finite collection of exponential surfaces still has O(j) logarithmic-shell parameters, and the partition argument continues to give exp(O((log X)^(2/3))) scalar counts after arbitrary products. Dense coverage requires a richer recursive interaction with a controlled returning state. Repeatedly using the compiler's scalar quotient as if it retained arbitrary root ports would be invalid.

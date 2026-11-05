# Complete unrestricted equality classification at amplitude 72

This is a separate extension of `../minimum-positive-clock/PROOF.md`. That frozen packet proves that 72 is the globally least positive amplitude for an eventual transitive-root-clone identity H_m=A*3^m+1. Here no symmetry assumption is made, and the equality sources are completely classified.

## Theorem

A rooted finite tournament (T,r) satisfies H_m=72*3^m+1 eventually, equivalently for every positive integer m, if and only if it is an ordered join

    Tr_l -> (S,a) -> Tr_k,  l,k>=0,

where the unique rooted strong core has rows

    S=[8,5,9,50,39,327,159,63,223], a=4.

Empty padding blocks are allowed. This core is rooted-isomorphic to the earlier attaining source

    [414,260,8,50,294,327,415,318,12], root4.

It is rigid and has exactly one root with this clone identity. Consequently, for each n>=9 there are exactly n-8 unlabelled rooted sources and n-8 unlabelled unmarked sources, (n-8)n! labelled rooted pairs and labelled unmarked sources, and (n-8)(n-1)! sources when the root label is prescribed. There are none at smaller orders.

These statements concern one fixed root replaced by a transitive module. They do not classify other graph-changing clocks or solve the unrestricted Hamilton-path spectrum.

## 1. The all-order cutoff, including the order-eleven boundary

The run expansion in the preceding minimum packet makes any eventual identity exact at every positive m. The first two values are 217 and 649. The marked root cannot be a singleton SCC, since cloning such a vertex preserves the count. The Hamilton-path product of every SCC outside its strong component divides

    gcd(217,649)=1.

Thus every exterior SCC is a singleton. Removing them leaves a strong root core of H=217 and preserves all clone counts. This reduction does not require an anti-involution, so the numbers of singletons before and after the core are independent.

Busch's Theorem1 gives the safe bound 5^(n-1)<=H^3 and the exact strong-order-eleven minimum 9*5^2=225. Because

    5^11 > 217^3,  and 225>217,

the core has order at most ten. The exact order-eleven bound is essential: 5^10<217^3, so the uniform bound alone would leave order eleven unresolved. The published source is Arthur H. Busch, *A Note on the Number of Hamiltonian Paths in Strong Tournaments*, EJC13(2006)N3, Theorem1, DOI10.37236/1141, https://www.combinatorics.org/ojs/index.php/eljc/article/download/v13i1n3/pdf.

## 2. Complete finite equality certificate

The standalone C++ verifier regenerates every strong tournament of order at most ten and H<=217, up to exact isomorphism. The strong-deletion lemma, induced Hamilton-path monotonicity, and exact degree-class permutation canonicalizer are proved in the preceding minimum packet and are unchanged here. Thus the cap does not omit any relevant intermediate graph or final core.

The retained class counts at orders three through ten are

    1,1,6,35,353,1144,726,190.

Every root of matching rank two in every retained graph with H=1 modulo6 is tested against

    H_2=3H_1-2,  H_3=9H_1-8.

At matching rank two the root-run degree is exactly three, so these equations are necessary and sufficient for the entire pure positive-unit ternary identity. The complete search tests 2,817 roots and finds exactly the displayed order-nine root, with H=217. In particular, no order-ten equality source exists.

The Python verifier independently recounts every retained graph by uncapped reverse-endpoint dynamic programming, checks strongness and exact canonical keys, recomputes root matching by augmenting paths, and directly constructs and counts the clone graphs. It independently compares the complete labelled strong universe through order six with the generated bank. It also enumerates every degree-compatible bijection of the equality core: only the identity is an automorphism. Checking all nine roots, without a balance or symmetry restriction, gives root4 as the sole equality root. A separate exact bijection identifies it with the earlier attaining source.

The actual clone counts at multiplicities one through five are

    217,649,1945,5833,17497.

They give positive run coefficients (217,216,72), while the rank-two bound excludes higher runs. Hence the displayed core realizes the identity for every m>=1.

## 3. Padding and counting

Ordered singleton SCCs form a transitive prefix and a transitive suffix. They contribute one to every Hamilton-path product, so every pair (l,k) yields an equality source. Conversely Section1 proves that every equality source has this form, and Section2 identifies its rooted strong core uniquely.

An isomorphism preserves the ordered SCC list. Therefore it preserves l and k, fixes every singleton SCC, and restricts to an automorphism of S. The full graph is consequently rigid. Its root is unique, because a singleton SCC cannot be a nonconstant clone root and the core has exactly one equality root. There are n-8 choices of l in 0,...,n-9, which proves all unlabelled and labelled formulas.

The independently reviewed marked self-converse equality theorem is recovered by additionally requiring l=k. This gives only odd orders and one unlabelled source per such order. It is strictly narrower than the present unrestricted theorem. No claim about equality of the two padding lengths is used here.

## 4. Reproduction and provenance

Run `python verify.py` and `python -O verify.py`; their deterministic JSON outputs agree. The adjacent `verify.cpp` regenerates its own bank, with no imported discovery data. The Python file checks ten physical prefix/suffix choices, including unequal padding, at clone multiplicities one, two and three.

The smaller-amplitude exclusion and its proof remain frozen in the preceding directory. The present new work is the complete H217 equality bank, the unique unrestricted core, and its arbitrary-padding enumeration. The earlier symmetry-restricted equality packet is `session197-positive-clock-audit/EQUALITY.md`; it is preserved unchanged. No literature-priority claim is made.

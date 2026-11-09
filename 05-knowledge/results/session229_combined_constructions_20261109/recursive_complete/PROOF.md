# Three-block escape and bounded-template recursive compilation

Status: PROVED candidates pending independent audit; computations FINITE-EXACT.
Chan–Pak completeness remains OPEN. No unrestricted coverage interval or new-to-literature scalar is claimed.

## Inheritance and live board

The source is the session229-work checkout. Inherited proofs are THM-3121-path-cover-walk-content-substitution-kernel.md, THM-3134-tournament-endpoint-jet-and-c3-newton-profile-transform.md and THM-3166-falling-factorial-order-join-path-colour-transform.md under 01-canon/theorems. The immediate predecessor is 05-knowledge/results/session228_recursive_constructions_20261108/recursive/PROOF.md and cyclic_inverse/PROOF.md. Session227's multiroot/PROOF.md supplies the positive-cover viewpoint. These sources, not historical session prose, authorize the identities below. The directed-triangle kernel itself is inherited, not new.

The board is physical substitution; full cover state; arity rather than depth; strong-node size; target-directed compiler. The closest mechanism is the session228 singleton-C3 compiler. The hostile is the scalar 19, still unavailable without a new seed/template. The corrected near miss is believing arbitrary binary depth includes arbitrary three-block substitutions. The underused sidecar is the third ordered-cover coordinate. The source-to-target map is literal tournament substitution; it preserves physical Hamilton counts, but scalar-only projection loses the future profile. Full Q restores this loss.

## 1. A sharp size bound survives unrestricted three-factor cyclic synthesis

For nonempty T, write h=H(T), F(r)=r! pc_T(r), and Q(t)=sum pc_T(r)(t)_r. Put n=|T|. Then h>=1 and

    F(2)=sum_(empty proper S) H(T[S]) H(T-S) >= 2^n-2.

For C=C3[A,B,D], the inherited diagonal law is

    H(C)=sum_r [3 F_A(r)F_B(r)F_D(r)
      +sum_i F_i(r+1)prod_(j!=i)F_j(r)
      +sum_i F_i(r)prod_(j!=i)F_j(r+1)].

Let n_i=|A_i| and a_i=2^(n_i)-2. Keeping only all terms in which each run count is 1 or 2 gives

    H(C) >= 3 + sum_i a_i + sum_(i<j) a_i a_j + 3 a_1 a_2 a_3
         >= 2^(n_1+n_2+n_3-2)+1.                     (1)

The difference between the two lower bounds is

    (1/2)sum_(i<j)a_i a_j + (11/4)a_1 a_2 a_3 >=0.

All discarded terms are nonnegative. The bound is sharp at every total order n>=3: C3[TT_(n-2),K1,K1] has H=2^(n-2)+1. Thus every arbitrary C3 node of scalar at most X has at most 2+floor(log2(X-1)) vertices, exactly the same cap as in the prior singleton-C3 grammar. In particular adding the third moving block does not incur a worse order cap. Equality in the second inequality requires at least two blocks to be singletons; no claim about equality in the first inequality beyond the displayed sharp witnesses is needed.

Consequently the session228 exact target/minimum-depth compiler extends to arbitrary C3[A,B,D]: enumerate all child profile triples by total order under this cap, using THM-3121 for the whole resulting profile; include binary ordered joins and finite physical seeds. Flatten exterior joins and take multiplicative scalar closure as before. Full-profile equality is a context congruence, and h is nondecreasing from child to parent. All children are smaller than the parent. These facts prove termination, completeness, and exact minimum nesting depth by the same induction as the predecessor. Neither profile truncation nor a fixed depth bound is legitimate.

## 2. A literal escape from the old complete grammar

Set E(A)=C3[TT2,TT2,A]. Since F_TT2(1)=1, F_TT2(2)=2 and higher terms vanish, the inherited kernel gives

    H(E(A))=11 H(A)+17 F_A(2)+4 F_A(3)
           =4 Q_A(3)+5 Q_A(2)-11 H(A).                (2)

This is positive in the cover coordinates; the negative coefficient after changing basis does not authorize subtracting physical graphs. For A=TT_m,

    H(E(TT_m))=4*3^m+5*2^m-11, m>=1.                (3)

At m=3 this gives 137 on seven vertices. The complete session228 compiler proves 137 is not attainable from all order<=5 seeds by arbitrary joins and C3[A,B,K1] nodes, even with unrestricted nesting. Equation (3) therefore gives an exact grammar escape, not just a new syntactic expression or higher nesting depth. For m=1,2 it gives 11,45, already in the old grammar. No assertion is made that infinitely many values in (3) escape the old grammar.

The finite companion computes all old profiles up to the proven cap for X=4095 and every single new C3[A,B,D] node with its three children any old-grammar profiles of orders <=5. Exactly 17 distinct scalar values <=4095 from these nodes are absent from the entire old scalar language. These are not the full image of the expanded recursive grammar. The entire old bank has 577 values, matching the inherited receipt. The new list, physical witnesses and reproduction checks are in CHECK.json.

The scalar 19 remains a hostile for the full C3 extension with these seeds. A three-factor node with no singleton has n_i>=2, hence already H>=45 by the first bound in (1). All nodes below 45 must have a singleton and are old operations. Since all ancestors have scalar at least each descendant, any expression of scalar19 cannot use a genuinely new node. The old complete compiler excludes19. Thus merely removing the singleton restriction still cannot establish all-odd coverage.

## 3. General finite-template theorem

Let R be any fixed quotient tournament on q>=3 vertices with a specified directed Hamilton cycle. Inflate its vertices by nonempty tournaments T_i of orders n_i. Define h_i=H(T_i), f_i=F_T_i(2), r_i=h_i+f_i. One traversal of the specified cycle, cut at each of its q arcs, gives q distinct quotient Hamilton path words. Starting at vertex i and traversing that cycle back to i gives a quotient word containing i twice and all other vertices once. Their physical decorations are pairwise disjoint by their maximal block-run words. Therefore

    H(R[T_i]) >= q prod_i h_i + sum_i f_i prod_(j!=i)h_j
               = sum_i r_i prod_(j!=i)h_j
               >= sum_i (2^(n_i)-1).                 (4)

The last inequality follows from h_i>=1 and r_i=h_i+f_i>=2^(n_i)-1. Convexity gives H>=q(2^(n/q)-1), but no floating point is needed: if n=qa+b, 0<=b<q, the minimum of the last sum over positive integer n_i with fixed n is

    L_q(n)=(q-b)2^a+b2^(a+1)-q.                      (5)

Pairwise balancing of two parts differing by at least2 proves (5). Hence H<=X implies n<=max{n>=q:L_q(n)<=X}. This is a finite logarithmic bound for fixed q. If X<q the set is empty: discard that template, since it has no output below the cutoff. If all templates are discarded, use profile cap0 and handle seeds entirely through exterior scalar multiplication.

Fix a finite physical seed bank and a finite bank of such cycle-certified templates. Close under all these substitutions and ordered joins. For each X, take the maximum cap supplied by (5) over the templates. The same target-bounded full-profile compiler is complete for this grammar, with all finite seeds of scalar<=X included as possible exterior factors even when their orders exceed the cap. Store any template cycle as an easily checked certificate; no unverified connectivity theorem is required. The q=3 cap from (1) is substantially sharper and should replace (5) for triangles.

This combines old quotient mechanisms in one terminating recursive language. It does not assert efficient runtime, an all-target predecessor, or that any one finite template bank is sufficient. Exterior joins still invalidate any order cap for whole outputs: C3=>TT_m has constant scalar3 and arbitrarily large order. Only template-rooted nodes are bounded.

## Reproduction and audit scope

Run `python3 -I -B verify.py > CHECK.json`, and similarly with `-O` for optimized replay. The standalone script uses standard-library code only and has explicit checks active under optimization. It includes the inherited complete binary compiler (credited in its source), direct subset Hamilton DP for every new witness, all 216 transitive C3 size triples with block orders1..6 for the bound (via exact inherited kernel), direct graph checks for size triples1..3, and independent subset counts for every reported new scalar. It checks (4) on all strongly connected labelled quotient tournaments of orders3 and4 and all transitive block orders1 or2, using explicitly found Hamilton cycles. The all-order conclusions are proved above, not inferred from those tests.

## 4. Complete recursive ternary image and the inherited unsaturated seed

The supplemental standalone-in-directory compiler `ternary.py` implements the full THM-3134 Gregory–Newton transform, not just its scalar slice. It loads only the accompanying `verify.py` by explicit local path and supports isolated Python. For every triple it first applies the positive exact scalar gate, then retains the whole polynomial and minimum nesting depth. The sharp cap (1) makes the image exhaustive even at arbitrary depth.

At X=4095 the complete counts are:

- Prior singleton-C3 grammar, all labelled seeds through order5: 577.
- Arbitrary three-child C3 grammar, same seeds: 620, with depth counts 79,431,109,1 at depths0,1,2,3. Exactly43 values are added; the first is137.
- Arbitrary three-child C3 grammar, same seeds plus the inherited unsaturated six-vertex H19/Q256 module with outgoing masks (40,1,3,6,15,30): 660, with depth counts109,445,105,1. Exactly83 values are added over the original577; the new primitive seed removes19, but29 is the first remaining admissible missing scalar.

The extra seed is inherited from session228's unsaturated_module/CHECK.json and is directly verified from its graph. It is not a newly discovered tournament. Each of these is an exact result for its specified grammar, not a claim of minimal order over all tournaments or unrestricted spectrum coverage. In particular finite computational exclusion of29 does not mean no tournament has29 paths.

The stronger calculation supersedes the limited17-value one-node comparison as the complete grammar result; the17-value calculation remains an independently reproducible control. TERNARY.json stores full added-value lists, first misses, profile counts by order, minimum-depth counts, sorted scalar/depth digests, recipe examples and literal-H verification counts for new recipes of order<=14. Both normal and optimized runs byte-match. Reproduce with `python3 -I -B ternary.py > TERNARY.json` and `python3 -I -B -O ternary.py > TERNARY.optimized.json`.

A structural explanation of the remaining small barrier is available without interpreting finite search as all-tournament exclusion: a node genuinely using three nonsingleton blocks has count at least45, while the first cyclic ancestor of an H19 seed, if it is a singleton-C3 node, has count at least H19+Q2=275: the intervening join-only child has h>=19 and r=Q2-h>=237, and the other nonempty child has h>=1,r>=1. Exterior multiplication using19 first changes a nonseed value at57. Thus adding the H19 seed cannot fix29 at any recursive depth. A further physical primitive or a different operation is necessary at that target.

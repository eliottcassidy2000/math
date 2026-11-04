# Cyclic two-singleton exterior: exact lift and finite C3-nine certificate

Producer result, awaiting independent review. This packet imports the 27-name C3 frontier's census, filter, and exact modulo-eight Hamilton response bank. It is not an arbitrary-core scalar-lift theorem.

## Exact ambient formula

Let M be a module and let the two exterior vertices have arcs M -> u -> v -> M. At any actual two-port core query, give both physical ports x,y exactly M's incidences to u,v: x,y -> u and v -> x,y. Thus S=M+x+y is a module of the augmented ambient. Write H=q1 and q2 for the complete ordered two-path-cover count. Every ambient Hamilton path has either one S-run (the three block orders S,u,v; u,v,S; v,S,u), or two S-runs (S,u,v,S). Consequently, exactly over the integers,

    F_ambient = 3 H(S) + q2(S).

No cubic or quartic contribution has been dropped: this is a bijection on complete physical paths. For a transitive three-block quotient, select ports with M's incidences; all three blocks occur in their forced order, and F_ambient=H(S).

## A general modulo-eight two-cover identity

For every tournament S with n>=3 vertices, the imported odd-cycle expansion gives H(S)=1+2 O(S) modulo four, where O counts all directed odd cycles of length at least three. This is the exact same input recorded in ../../session189-incoming-anchor-audit/INCOMING_PROOF.md Section 1 (Grinberg--Stanley, Theorem 1.39).

Expand q2(S) as twice the sum of H(A)H(B) over unordered nonempty bipartitions. Modulo eight its constant term is 2^n-2=-2, and its remaining term is four times the total number of odd cycles contained entirely in one bipartition side. A fixed nonspanning odd cycle on r vertices occurs 2^(n-r)-1 times, an odd number. A spanning cycle never occurs. If n is odd and C(S) counts directed Hamilton cycles, this proves

    q2(S) = 2 H(S)+4-4 C(S) modulo eight,
    F_ambient = 5 H(S)+4-4 C(S) modulo eight.

For even n, the same statement holds without C(S). Thus at odd order an exact sufficient condition is 5 Delta H - 4 Delta C != 0 modulo eight. In particular a core gap four lifts if and only if the Hamilton-cycle parity agrees. This also explains why arbitrary scalar lifting cannot be inferred from the scalar gap alone.

## Cycle parity contracts the nine triangles

Take M=Q[C3^9], child-uniform ports, so n=29. Anchor every Hamilton cycle of S at the singleton port x. Decompose its unique traversal into maximal child runs and the two singleton ports. If a child has h runs, its physical fillings number q_h(C3), namely (3,6,6). Modulo two every word with a repeated child vanishes. All remaining words use each of the nine children and both ports once, and each has odd filling multiplicity 3^9. Therefore C(S) modulo two equals the Hamilton-cycle parity of the 11-vertex augmented quotient.

This reasoning is unchanged by any original-module changes inside the physical triangles: their complete cover profiles stay (3,6,6). Consequently both H and C at these child-uniform queries, and hence the exact ambient counts, are independent of those interior changes.

## Retained finite computation

cycle_parity.cpp counts Hamilton cycles of each augmented quotient modulo two by an x-anchored subset DP. It consumes the inherited bank_input.txt, and cycle_parity.tsv retains source parity, target parity, and their XOR for all twelve inherited queries at all 2348 endpoint-zero quotient representatives. The inherited residues.tsv supplies their exact physical H residues.

For each entry evaluate 5(H_source-H_target)-4(C_source-C_target) modulo eight. The twelve inherited queries separate 2347 of the 2348 quotient classes. The sole miss is node 56186, with outgoing quotient rows

    8,13,49,100,139,147,55,335,127.

The additional query (a,c)=(2,4) has physical H residues (3,7) and equal Hamilton-cycle parity, hence ambient gap four modulo eight. Thus the twelve inherited queries together with (2,4) give a 13-query fixed separating battery for the cyclic two-singleton exterior in the entire endpoint-zero class. extra_input.txt and its two output TSVs retain an exploratory 96-query bank on the sole miss; index 20 is the added (2,4) witness. witnesses.tsv retains one ambient witness per quotient.

There are 7624 inherited quotient/query pairs with odd cycle-parity gap, so Delta q2 is not universally zero modulo eight. No unrestricted arbitrary-core response-injectivity or arbitrary-exterior theorem is asserted.

The endpoint-nonzero class has an inherited child-uniform core query with Delta H=2 modulo four. The identity above gives q2(S)=2H(S) modulo four, so G(S)=3H(S)+q2(S)=H(S) modulo four. Consequently its core contrast two lifts directly to the cyclic ambient, without importing an arbitrary-ambient theorem. For a transitive exterior quotient, the exact equality G=H already gives the same conclusion. A simultaneous isomorphism/converse transport of the quotient and both masks, leaving physical interiors paired on the same names, transports every witness. Converse exchanges source/target. The opposite cyclic exterior orientation is obtained by swapping exterior names.

## Reproduction and physical convention

Run python3 verify.py and python3 -O verify.py. The verifier checks pinned imports, rebuilds the entire cycle-parity bank from source, rebuilds the extra core-response bank from the inherited exact recurrence, and checks the modular identities on 1296 literal small tournaments (all n=3,4,5 and 100 fixed-seed examples each of n=6,7). It validates all physical graphs and writes every chosen witness. Imported full twelve-query H residues are hash-pinned rather than re-proved here.

physical_witnesses.jsonl.gz has old core names 0 through 26, u=27, v=28, x=29, y=30. Thus there are 29 old names and 31 augmented names. Each physical mask expands its nine coarse bits to triples and additionally includes bit 28, since v points to each port; bit 27 is absent. u points to v, core and ports point to u, and v points to core and ports. Each record retains both complete physical augmented graphs and both modulo-eight ambient residues.

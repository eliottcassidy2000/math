# Module-indecomposable unit ternary clocks and explicit arithmetic witnesses

Session 197. **Two explicit physical sources, with complete finite structural certificates.** Unit ternary clone clocks need not have any nontrivial homogeneous module, even away from the root. This is a scope counterexample to a possible structural extrapolation from small clocks, not a classification or a new spectrum-coverage theorem.

## 1. The two sources

The bit in position v of row u means u→v. The marked root is the last vertex.

Negative unit source, order10:

    [46,320,2,1014,999,966,901,775,517,7].

Positive unit source, order11:

    [46,640,2,2038,2023,1990,1927,1797,1543,1029,7].

Both are strong, rigid tournaments with **no nontrivial homogeneous module**. In particular neither has a transitive twin pair anywhere. Here a nontrivial module means a vertex subset of size between2 and n−1 whose members have identical incidences with every outside vertex. This property is distinct from arithmetic primality of H.

Their root-run vectors (a0,a1,a2,a3) are respectively

    (259,779,780,260),
    (521,1561,1560,520).

Consequently replacing the root by a transitive chain of length m gives, for every integer m>=1,

    H^-_m=260*3^m−1,
    H^+_m=520*3^m+1.                              (1)

The deletion values at m=0 also agree with (1). Nontrivial clone chains at m>=2 are themselves modules, so the module-indecomposability assertion is deliberately about the two original m=1 sources only.

## 2. Physical construction and the exact all-order identity

Use O={a,b,u} in order a→u→b, with a→b. The I side is a transitive chain whose cross-neighborhood word is

    1,0,1,0^j,2,0,2,

where j is0 for the negative source and1 for the positive source. Type1 vertices lose to a and beat b,u. Type2 vertices lose to b and beat a,u. Type0 vertices beat all three O vertices. Every I vertex beats the marked root, which beats every O vertex. Labels are a=0,b=1,u=2, followed by I in chain order and then the root. This gives precisely the two displayed row arrays. Replacing the last vertex by a transitive m-chain gives the graphs in (1).

The crossing graph O→I has matching rank exactly two. All crossing edges meet {a,b}, and a type1 edge and a type2 edge form a matching of size two. In a collapsed clone-run Hamilton word, each internal root-to-root fragment needs a disjoint O→I crossing edge. Thus there are at most three root runs. The verifier independently counts the complete run vector by the exterior-word/gap recurrence. Its coefficient a4 is also checked to be zero, while the cover argument rules out all higher coefficients.

For m positive labelled transitive clones, distributing them among j nonempty ordered runs gives j! S(m,j) choices: within each run their order is forced, but the sets of labels assigned to distinct runs are arbitrary. Therefore

    H_m=a1+(2^m−2)a2+(3^m−3*2^m+3)a3.

The displayed vectors have a1=3A+beta, a2=3A, a3=A, so the binary terms cancel and the constant becomes beta. This proves (1) at every positive length. At length zero the separate deletion count is a0=A+beta.

These sources are the i=k=1 specializations of the already independently reviewed four-seed family, published at:

https://github.com/eliottcassidy2000/math/blob/ef10dc937a51f36154189a26965e8624e4895eb8/session197/clock-gap-families/README.md

The local named predecessor is `session197-prime-atoms/two-module-unit-clocks/PROOF.md`, SHA256 `eaf139d56e65acfbea8e5b7550770b2ea21bddbcbbe0ad88764c7f3a5f7445ea`. That all-parameter family theorem is corroborating context; the explicit run counts and matching argument above suffice for these two rays.

## 3. Complete structural certificates

For each original source the verifier visits every vertex subset with size between2 and n−1 and checks the defining module incidences against every outside vertex. None passes. This is a complete finite proof for the explicit graphs, not a random search or a claim about larger families.

Every automorphism preserves outdegree. The verifier enumerates every within-degree-group permutation and checks all adjacency entries. There are only8 such permutations for the order10 source and16 for the order11 source. Only the identity survives. Strongness is independently checked by directed reachability from every vertex.

This refutes the universal claims that every unit ternary clock has a twin pair, has a nontrivial module somewhere, or reduces by a nonroot module contraction. It does not conflict with the unit-boundary/twin theorem for **rank-one binary** clocks: these sources have rank two.

No least-order or least-amplitude claim for module-indecomposable unit clocks is made.

## 4. Actual prime and 7p members

The receipt includes full adjacency rows and independently counted H for all of these graphs:

- Negative ray, m=2: order11, H=2339, prime.
- Negative ray, m=3: order12, H=7019, prime.
- Positive ray, m=1: order11, H=1561=7*223, with223 prime. This arithmetic witness is itself module-indecomposable and rigid.
- Negative ray, m=6: order15, H=189539=7*27077, with27077 prime.

Primality is certified by trial division through the exact integer square root:48,83,14 and164 respectively. The source H=779 factors as19*41 and is not claimed prime. These are finite arithmetic witnesses on all-order graph families. Nothing here asserts infinitely many primes or 7p values, covers all prime atoms, or proves an inverse representation theorem.

## 5. Reproduction

Run `python3 verify.py` and `python3 -O verify.py`. The checker is standalone and uses only the standard library. It constructs the graphs directly, computes Hamilton counts by a subset endpoint DP, computes the root-run vectors separately, exhausts modules and degree-preserving automorphisms of the original sources, and verifies every displayed primality claim by exact division. It writes no files; redirect its output to preserve a receipt. Normal and optimized outputs must agree.

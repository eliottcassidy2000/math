# Independent review of Q/Z commutation and sparse compiler outputs

Status: **APPROVED / PASS.** The producer applied the explicit nontrivial-source qualifier to both transfer laws and stable closed forms, and separately handled singleton inputs. Graph normal form includes the singleton. The finalnewfreeze/ directory freezes the final reviewed inputs including the approved sharpness addendum; approved/ and earlier snapshots preserve the review lineage.

## Findings

1. **Graph theorem verified.** After canonical relabeling, deletion is TT_q => (G-root) => TT_z, with root beating the Q block and losing to the Z block. Q vertices are ordered newest first; Z vertices oldest first. Every Q beats every Z under either insertion order. Thus the full rooted graph, not merely its seven coordinates, depends only on (q,z). There are at most k+1 isomorphism types at word length k.

2. **Transfers and reset verified for nontrivial inputs.** L is the root-end cover statistic, R the root-start statistic. Both displayed transfers and their common composite have the correct orientation and endpoint subtraction. The QZ reset sets s=e=d; deletion count stays d. The coordinates need not have been stable before reset.

3. **Stable formulas verified for nontrivial stable inputs.** Both exponent-zero axes are included. Direct substitution gives the displayed h,a,rho,b formulas, expanded X formula, and constant h-a-rho+b+2j. The conserved-defect rewriting is the inherited identity, not an independent assumption about signs. Potentially negative expanded coefficients do not support the counting bound; the positive unexpanded h expression does.

4. **The explicit logarithmic-square bound is valid.** For nonempty deletion, singleton-root covers and full root-ending paths with the other bin empty are disjoint, giving L>=d+e=2d in a stable source. Likewise R>=d+s=2d; h>=d. Therefore h(q,z)>=d2^q and >=d2^z. With compiler r in {1,2}, j in [0,8], X=h+a+2j>=h and H_compiled=3X+2r, each output <=N satisfies the claimed exponent bounds. Invalid residue/budget states contribute no outputs. The bound counts exponent pairs, not distinctness or multiplicity.

5. **Transient axes verified at their minimum exponent1.** On Q, the coefficient in X is 2(L-d)>=0, with zero coefficient allowed. On Z it is R-d+b>0, since R>=d and b>0. Each positive-coefficient affine dyadic sequence has O_G(log N) outputs below N despite any negative constant term; a zero coefficient gives one value. Empty word gives at most one further output. Words with both letters are covered by one stable reset. The finite-source-bank O(log²N) conclusion follows; infinite banks and other moves remain excluded.

6. **Singleton distinction is necessary and repaired.** Its profile is (1,1,1,1,1,1,1). The nontrivial Q transfer incorrectly preserves e=1 there, whereas the literal two-vertex Q graph has e=0. The proposed common composite would also give h=2 instead of the actual QZ triangle's h=3. Graph commutation itself is unaffected. The separate singleton argument is correct: pure axes are transitive, Q-only nonempty compiler base is constant, Z-only has start0 and is inadmissible, and QZ is the stable directed triangle with d=1. Consequently the same asymptotic conclusion includes singleton seeds after handling them separately. Section2's closed-form hypothesis now explicitly excludes the singleton too.

## Independent reproduction

Run `python3 independent.py` and `python3 -O independent.py` in this audit directory. The script imports no producer code or numerical output. It uses direct permutation enumeration on every induced subset for endpoint and cover profiles, rather than the producer's subset dynamic programming. Its separate graph-normal-form construction compares full canonical adjacency matrices.

Results, with matching normal and optimized outputs:

* 9,525 canonical graph checks: every rooted labeled tournament of orders1..4, root0, and every Q/Z word of lengths0..6, including the empty word and singleton.
* 80 literal permutation-profile checks of transfers and stable descendants.
* 250 stable exponent pairs q,z in0..4, with compiler-base expressions checked at every j in0..8.
* 120 transient-axis checks at exponents1..6; three source profiles have the degenerate Q coefficient0.
* The singleton transfer is a deliberate hostile control and fails as expected.

The frozen revised producer verifier was also rerun and its output matches the retained producer receipt. Finite checks supplement the all-length proof; they do not establish its quantifiers by themselves.

## Consequence and novelty boundary

The inherited paired-source synchronization criterion remains valid. What fails is the interpretation of two move names as independent binary digits: their graph-level commutation collapses words to two counters. This gives a rigorous sparse-family limitation even if a fully paired seed is later found. It does not rule out all-length coverage by different operations or an unbounded changing source bank, and it does not settle the unrestricted Hamiltonian-path spectrum.

## Approved sharpness addendum

**APPROVED / PASS for final Section4.** The exact amended theorem and both producer verifiers are frozen in `finalnewfreeze/`, with `INPUT_SHA256SUMS`. This is the final reviewed version; `approved/` preserves the pre-sharpness approval.

The directed-triangle profile is correct. Root-start remains d=1 throughout, so the width-eight compiler always has r=1,j=0. Substitution in the stable formula gives exactly the displayed X and H=3X+2.

For fixed n=q+z>=2, the subtraction 4·2^q+2^(n-q+1) is convex as a function of q. Its maximum on integer1<=q<=n-1 is attained at an endpoint; the q=n-1 endpoint dominates the other by2^n-4>=0. Thus the stated lower layer bound is attained, including the boundary n=2. The upper bound is strict. Successive layer intervals are disjoint for every n>=2.

A sum of two powers of2 uniquely determines the unordered pair of exponents: distinct exponents give two set bits; equal exponents give a single set bit and cannot coincide with a distinct-exponent sum. Accordingly each scalar value has at most two exponent-pair preimages. This justifies the retained lower bound (n-1)/2 even when that number is nonintegral. The affine H map preserves separation and multiplicities.

For K=floor(log2((N-11)/24))>=2, one has24·2^K+11<=N. Every included layer n<=K satisfies

    H < 24·2^n+11 <= 24·2^K+11 <= N.

Thus all counted outputs are strictly below N, hence also at most N, and summing the layer lower bounds gives K(K-1)/4. Since K grows as log2 N, the constructive lower bound and previously audited upper bound prove Theta((log N)^2) for this single physical source. No coverage or positive-density conclusion is claimed.

`independent_sharpness.py` uses the independent audit's permutation-profile routines, not producer code. Normal and optimized outputs match. It verifies ten literal compiler-base graphs (q,z>=1, q+z<=5), 3160 exponent pairs through layer80, and79 exact integer N-threshold checks. It also checks the sharper optional arithmetic observation that layer n has exactly ceil(n/2) distinct values, because its only collision involution is (q,z)↔(z-1,q+1). This observation is an audit control, not needed for the published sufficient bound. The amended producer sharpness verifier was rerun with output matching its retained receipt.

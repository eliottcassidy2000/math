# A composite-circulant hostile to the constant-sum path test

Status: PROVED counterexample to a test implication, not to the strong pure-blue classification. FINITE-EXACT retest of the incoming displayed paths through order 23. No all-order repair of those particular path formulas is claimed.

## 1. The precise failed implication

A Hamiltonian path P=(v_0,...,v_(n-1)) is palindromic here when its positional reversal r_P(v_i)=v_(n-1-i) is an anti-automorphism of the tournament. In a circulant on Z/nZ, constant paired sums v_i+v_(n-1-i)=c imply r_P(x)=c-x, which is indeed an anti-automorphism.

The converse is false for general circulants: additional anti-automorphisms need not have this affine form. Thus nonconstant paired sums alone do NOT prove that a path is nonpalindromic.

[THM-4730](../../../01-canon/theorems/THM-4730-strong-pure-blue-alphabet-complete-points-of-symmetry-are-k1-c3-t5.md), section2, uses this converse when it rules out its three-interval and descending-pair paths by unequal sums. Its `atom_alphabet_proofcheck.py` originally implemented the same weaker condition in `palindromic(n,path)`. The original finite script checked directedness correctly but did not test every arc under the reversal in its odd-order lane. Commit eb6c3e9818e independently identified this gap, repaired the checker to test arcs, corrected the universe totals and added a proof-gap override. Those incoming repairs are preserved; the discussion here records the original failed implication and our independently produced hostile.

## 2. Exact order-nine hostile, with a proof

Use n=9 and the skew connection set S={1,3,4,7}, so x->y iff y-x lies in S modulo 9. Consider

    P=(3,1,2,6,4,5,0,7,8).

The increments are (7,1,4,7,1,4,7,1), all in S. Thus P is a directed Hamiltonian path. Its paired sums are

    (2,8,2,2,8,2,2,8,2) modulo 9,

which are not constant. Nevertheless its reversal is an anti-automorphism.

For a conceptual proof, let q add 3 to vertices congruent to 0 modulo 3 and fix the other vertices. Between different residue classes, orientation depends only on their residues, so q preserves it. Within the residue-zero class, q is the rotation of its directed C3; the other two classes are fixed pointwise. Hence q is an automorphism. The ordinary cyclic Hamiltonian path P0=(0,1,...,8) has anti-automorphism r(x)=8-x as its positional reversal. Since P=q(P0), its reversal is q r q^(-1), again an anti-automorphism.

The portable verifier independently tests every unordered vertex pair under r_P, in addition to the increment and unequal-sum checks. This verifies the actual predicate rather than another affine proxy.

## 3. What survives the incoming explicit-path computation

For every odd n=7,9,...,23, the verifier enumerates every skew connection set containing 1. It checks each three-interval path associated with an unordered pair a<b in S\{1}, a+b<=n+1; for the exceptional no-small-element set it also checks the stated descending-pair path. Every tested path is Hamiltonian and its reversal genuinely fails the arc condition. The numbers of tested paths are:

    n:      7   9  11  13  15   17   19    21    23
    tests:  4  17  61 193 561 1537 4033 10241 25345

Total: 41,992 actual arc tests with zero invalid witnesses. In particular all 2,044 incoming selected normalized-set witnesses are valid over this range. The original incoming summary's 2,047 was an arithmetic transcription error, now corrected upstream: the retained output has 4+8+...+1024=2,044 sets, consisting of 2,035 selected three-interval cases and 9 descending-pair exceptions.

This finite result does not show that the proposed formulas always yield a nonpalindromic path. Nor does the hostile above belong to the incoming displayed family: it invalidates the test implication, not those formulas. A future proof could validate them by a true anti-automorphism argument, but none is supplied here.

## 4. The classification remains proved by a separate route

The [uniform prefix-rerouting proof](../pure_blue_atom_classification_20261020/UNIFORM_PREFIX_REROUTING_CLASSIFICATION.md) rules out all odd circulants of order at least7 using explicit violated arcs and unit-step reindexing. It does not assume that all anti-automorphisms are affine. Combined with its audited strong-cycle reduction, it proves the alphabet `{1,C3,regular cyclic C5}` independently of THM-4730's incomplete odd proof.

The incoming rotation lemma and even-order mixing-cycle proof are different matters. Both are valid on inspection: the latter constructs a Hamiltonian cycle whose shift-by-two sends an even reference vertex to an odd one, making the automorphism group transitive and contradicting even-order regularity. This is an elementary alternative to the independent proof's cited Lu-crossing route. Preserve this credit rather than treating the entire incoming result as false or replacing the earlier independent proof history.

## 5. Incoming independent correction preserved

Before this packet was frozen, commit eb6c3e9818e published an independent
[weighted-core inheritance audit](../weighted_core_transport_20261020/inheritance.md#10-incoming-pure-blue-alphabet-claim-a-proof-test-gap-not-a-refutation).
It gives a different path in the same order-nine circulant,
`(0,1,2,3,7,5,6,4,8)`, with nonconstant partner sums and a true positional
anti-automorphism. It also reports a broader special-family finite probe
through order 25. We verified that new hostile directly on all 36 arcs;
we do not present its broader finite run as our own rerun.
The earlier public warning and checker repair are retained verbatim in the
proposed exact-prior update. Our contribution to the current theorem route
is the separately audited all-order prefix proof, not priority over this
independently discovered correction or a repair of the old special paths.

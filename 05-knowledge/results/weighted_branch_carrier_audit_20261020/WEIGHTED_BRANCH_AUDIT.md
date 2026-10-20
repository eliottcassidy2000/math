# Weighted level-branch carrier: independent audit and corrected flip census

Date: 20 October 2026. Audited source: commit `127698803818637f07ae89c32ac5d04c59161292`, THM-4726, `01-canon/theorems/THM-4726-level-branch-as-a-weighted-count-over-two-element-decorations-and-the-decorated-flip.md`.

**Verdict: the general centralizer-weight identity and the branch totals pass; two material corrections are required.** The normal-complement formula has the wrong pointwise denominator, and the claimed decorated-flip count is an unmarked count. The correct number is **7 fixed points among 31 involution-decorated types at n=8**, not 9. The 49-pair weight histogram and total 15 are correct. This audit does not review or endorse THM-4727 or the unrelated Jacobian material in the source commit.

## 1. Exact carrier and proof

Let C run over isomorphism types of tournament switching classes, set H=Aut(C), and let rho(H) be the fraction of odd-order elements. The level branch has the exact identity

    N_lev(n) = sum_[C] (1-rho(H))
             = sum_([C],[x]_H), x!=1 of 2-power order, rho(C_H(x)).

The second index is an H-conjugacy class of **2-elements**, not merely an involution. Equivalently the indices are S_n-orbits of marked pairs (C,x). Every even-order element g has a unique commuting decomposition g=x y, where x is a nonidentity 2-element and y has odd order; both are powers of g. For fixed x, y ranges over the odd-order elements of C_H(x). A conjugacy class of x has |H|/|C_H(x)| elements. Dividing its contribution by |H| gives exactly rho(C_H(x)). No normal-complement assumption is needed for this identity.

This is a weighted identity, not a claim that these indices already form a set of cardinality N_lev. Forgetting x, its conjugacy class, or its centralizer loses information needed by the sum. At n=8 there are 27 unmarked even-Aut types, 31 involution-marked types, and 49 nontrivial-2-element-marked types.

### Correct normal-2-complement specialization

Suppose H=K semidirect P, where K is a normal odd-order 2-complement and P is a Sylow 2-subgroup. For x in P,

    rho(C_H(x)) = 1/|C_P(x)|,

not 1/|P| in general. The projection C_H(x) -> C_P(x) is onto with kernel C_K(x); that kernel is its normal odd-order complement. H-conjugacy between elements of P is exactly P-conjugacy: projecting a conjugacy equation to H/K gives the corresponding conjugacy equation in P. Consequently

    sum_[x]_P, x!=1 1/|C_P(x)|
      = sum_[x]_P, x!=1 |[x]_P|/|P|
      = (|P|-1)/|P|.

Equivalently one may sum 1/|P| over every individual nonidentity element of P. These two indexing conventions must not be mixed. If P is abelian, each conjugacy-class weight does equal 1/|P|.

The smallest relevant hostile is H=P=D8, of order 8. Its nonidentity 2-element conjugacy classes have weights 1/4, 1/4, 1/4, 1/8, summing to 7/8. A reflection has centralizer V4 and weight 1/4, not 1/8. The independent controls also check C3 x D8, so the correction is not confined to a trivial odd complement.

For n=2 mod4 with n>=6, the inherited semiregularity theorem forces every even-Aut Sylow 2-subgroup to have order 2. Thus the single involution decoration contributes 1/2 and THM-4531's Theorem B pairs the terms. The n=2 exception must remain explicit: its matching is a cut, its flip is fixed, and N_lev(2)=1/2.

For an arbitrary finite group, rho(C_H(x)) need not be the reciprocal of a 2-power: for the central involution in C2 x A5 it is 3/8. This example is a boundary for the unrestricted finite-group identity, **not** a claim that C2 x A5 occurs as a switching-class automorphism group. A reciprocal-centralizer formula for a special class of automorphism groups requires its own normal-complement argument; the exact rho formula needs no such replacement.

## 2. Independent exhaustive finite census

The C++ program enumerates all 2^((n-1)(n-2)/2) labelled switching states after vertex 0 has been switched to a source. It applies every permutation in S_n, partitions states into orbits, explicitly computes their stabilizers and 2-element conjugacy classes, and uses rational arithmetic throughout. It imports no repository mathematics, nauty code, or producer orbit data.

| n | labelled quotient states | class types | even-Aut types | 2-element decorations | involution decorations | true flip-fixed | N_lev |
|---|---:|---:|---:|---:|---:|---:|---:|
| 2 | 1 | 1 | 1 | 1 | 1 | 1 | 1/2 |
| 4 | 8 | 2 | 2 | 4 | 2 | 2 | 1 |
| 6 | 1024 | 6 | 4 | 4 | 4 | 0 | 2 |
| 8 | 2097152 | 79 | 27 | 49 | 31 | 7 | 15 |

At n=8 the weight counts are exactly

    19 at 1/2, 14 at 1/4, 16 at 1/8,
    19/2 + 14/4 + 16/8 = 15.

The Sylow-type masses also pass: C2 contributes 19/2, C4 contributes 3/2, V4 contributes 1, C8 contributes 7/4, and D8 contributes 5/4. The order-8 Sylows are certified from explicit cyclic elements or an order-8 centralizer with five involutions and two elements of order 4. The order-168 group is not identified by name solely from its order; its enumerated centralizers and weights suffice for this audit.

A second Python implementation uses directed adjacency matrices rather than affine bit transformations. It constructs every permutation commuting with each fixed-point-free involution by permuting and flipping its vertex pairs. It performs 12,114 centralizer-constrained isomorphism tests over the n=2,4,6,8 decorations and verifies the marked flip is an involution and preserves the weights. It separately enumerates all of S8 for the decisive D8 stabilizer. Both Python ordinary and -O runs agree byte-for-byte. The C++ optimized and address/undefined-behavior-sanitized runs agree byte-for-byte (leak detection disabled because the host's ptrace prevents LeakSanitizer).

An independent partition-formula calculation reproduces N_odd and N_lev through n=16, including the n=2 exception and totals 1,2,15 at n=4,6,8.

## 3. Why the 9-of-31 claim is false

For the marked involution flip

    tau(C,g) = (C+M_g,g),

fixedness requires a relabelling h with hC=C+M_g **and hg=gh**. The producer script comments that it checks this, but its executed test compares only the two unmarked canonical forms. Thus its output counts nine decorations whose underlying types agree, not nine fixed marked types.

The discrepancy occurs precisely in the D8 automorphism type. Under the independent low-bit-first encoding of edges (1,2),(1,3),...,(6,7), with vertex 0 a source, take state 4722 and

    g1 = [1,0,3,2,6,7,4,5],
    g2 = [2,7,0,4,3,6,5,1].

These are representatives of the two reflection conjugacy classes. Flipping their matchings gives states 594445 and 269715 respectively. Unmarked isomorphisms exist:

    h1 = [0,5,1,2,3,6,4,7],
    h2 = [0,2,3,4,6,1,5,7].

But neither commutes with its marked involution. In fact h1^-1 g1 h1=g2 and h2^-1 g2 h2=g1, and the two involutions are not conjugate inside Aut(C). All 384 commuting relabellings fail in each case. The two decorations are exchanged by tau, rather than fixed.

The seven true fixed decorations are three from C2-type classes, two from C4-type classes, and two from the V4 class. No D8 decoration is fixed. Therefore the true flip has **12 transpositions**, not 11.

There is a useful exact structural refinement: C_Aut(C)(g)=C_Aut(C+M_g)(g), as literal subgroups of S_n, because any h commuting with g preserves M_g. Thus tau preserves its weight. At n=8 the involution decorations have total mass 12, of which 5/2 is fixed mass and 19/2 lies in transpositions. The noninvolution 2-element decorations supply the remaining mass 3. Neither the number of transpositions nor its mass equals N_lev=15.

The source's prescribed unions C2 union C4 and C8 union D8 union V4 are arithmetic groupings, not consequences specifying a mandatory pairing. Nonintegral Sylow-stratum masses show that an integral-packet interpretation cannot keep every Sylow stratum separate. They do not force those exact unions or define any canonical carrier.

## 4. Optional finite continuation: the full marked carrier

There is a natural extension to all nonidentity 2-elements x. Put j(x)=x^(ord(x)/2), its terminal involution, and define

    tau_full(C,x)=(C+M_j(x),x).

This is well defined on marked isomorphism types, is an involution, and preserves C_Aut(C)(x), hence the exact weight: every relabelling commuting with x preserves j(x) and its matching. The marked element x still fixes the flipped class.

A separate matrix-action probe gives fixed counts 2/4 at n=4, 0/4 at n=6, and **9/49 at n=8**, with seven fixed order-2 decorations and two fixed order-4 decorations at n=8. The two additional fixed decorations have state 291 and weight 1/4 each; explicit commuting witnesses are in the receipt. The n=8 fixed mass is 3. This is a different carrier from the source's 31 involution decorations. The coincidental numeral nine does not repair the original claim.

This probe remains FINITE-EXACT. It does not show that fixed mass or transposition mass is integral in general and does not yield the open integral-packet construction. Even a fixed-point-free, weight-preserving involution would not by itself prove integer total: a 2-cycle of weight 1/4 per element contributes 1/2, and one of weight 1/8 contributes 1/4. Packets must clear the full 2-primary denominator, or be replaced by a separately justified integral combinatorial species. The marked x and its full centralizer must remain in the carrier; replacing x by j(x) alone loses its order and centralizer data.

## 5. Precise relation to the prime-local branch/code/lift packet

Source: `05-knowledge/results/prime_local_branch_code_lifts_20261020/PRIME_LOCAL_BRANCHES_CODES_AND_LIFTS.md`.

There is an actual mathematical bridge for the finite-group-action lane. For any prime p and any finite G-set X, with rho_p(H) the p-regular proportion,

    B_p(X)=sum_[a in X/G] rho_p(G_a),
    |X/G|-B_p(X)=sum_([a],[x]), x a nontrivial p-element,
                         rho_p(C_(G_a)(x)).

The proof is the same unique commuting p-part/p'-part decomposition. The map retains the **same action**, the object stabilizer, the marked p-element conjugacy class and its centralizer, and exactly preserves the complementary restricted Burnside contribution. Passing to the number loses the object, marked class, automorphism embedding, and possible pairing data. Restoring these requires the full decorated action/groupoid, not only an integer or prime-support list.

At p=2 on the affine F2 gain carrier, this is exactly the tournament switching carrier and the complement is N_lev. On the linear graph-switching carrier it is not N_lev. The decisive n=3 check is: both regular branches equal 1, but the linear total is 2 with complement 1, while the affine tournament total is 1 with N_lev=0. Equal regular sums do not identify the singular carriers. At odd p, the predecessor's affine-to-linear translation by 1/2 applies; it does not identify different prime actions or construct simultaneous mixed-prime witnesses.

The predecessor's numerical integrality proof is unchanged: Frobenius places the restricted Burnside sum in Z[1/p], and the explicit p-regular fixed-space formula independently places it in Z_(p). Their intersection is Z, with its stated exception. The new weighted refinement does not supply the missing canonical partition/pairing.

There is **no established map here** from the decorated switching carrier to Pascal codes or Jordan lift clocks. The code identification preserves its own coefficient/evaluation weights, and Jordan reduction preserves a factor trajectory but can lose exact period. Their concrete analogical lesson is to retain a lost coordinate; that is not a shared carrier or transferred theorem. The faithful CRT map in the predecessor retains the whole local-state tuple and period data. It does not license combining separate p-marked class choices into a single object.

## 6. Precise relation to the Collatz local-prime-power packet

Source: `05-knowledge/results/collatz_local_prime_power_gates_20261020/LOCAL_CARRY_AND_PRIME_POWER_GATES.md`.

Its exact source map is a fixed chronological word w -> (B(w),D(w)), with D=2^q-3^p, and then to the valuation deficits of B against every prime power in D. Integer-cycle recognition is D|B for the **same word**; positivity additionally requires D>0. This is not a finite-group stabilizer or conjugacy-class map, and no map from those words to the weighted carrier is supplied by either packet.

The correct connection is an obstruction-retention lesson, not a Collatz theorem. A quotient is safe only after specifying the predicate and discarded data:

- In the branch problem, forgetting the marked involution changes the fixed-point predicate. The D8 example falsely turns a transposition into two fixed points.
- In Collatz, forgetting prime-power depth changes the integer predicate: `11101000100` has B=475, D=1805=5*19^2, gcd=95 and reduced denominator 19. Prime support passes but integrality fails.
- Forgetting the common word also changes quantifiers. At (q,p)=(15,8), D=73*359, the independently reproduced necklace contents are 1:418, 73:10, 359:1, and D:0. Each prime has a witness, but no common word satisfies both.
- Forgetting primitive repetition can inflate content without changing the reduced denominator. The chronological order, clock, common necklace and full valuations remain indispensable.

The comparison controls reproduce the first two arithmetic hostiles and the n=3 linear/affine boundary directly. They do not rerun or reclassify the predecessor's full q<=400 edit census. No Collatz convergence result, all-cycle exclusion, mixed-prime symmetry correspondence, or transfer from branch integrality to Collatz descent follows.

## 7. What remains open

Numerical branch integrality is already proved. THM-4531's n=2 mod4, n>=6 pairing remains valid. At n divisible by four, what remains open is an explicit natural decorated or weighted combinatorial construction explaining the level-branch integer. Neither the 27 unmarked even-Aut types, the 31 involution-marked types, nor the natural full 49-marked terminal-involution flip at n=8 gives the desired fixed-point-free unit-mass pairing. A different enrichment, grouping rule, or construction may exist; these controls do not rule out every carrier or force an involution as the only possible proof architecture.

## 8. Reproduction and review scope

From this packet directory:

    g++ -std=c++17 -O2 -Wall -Wextra independent_carrier.cpp -o independent_carrier
    ./independent_carrier > independent_carrier.replay.json
    cmp independent_carrier.json independent_carrier.replay.json
    python3 independent_marked_check.py > independent_marked_check.replay.json
    cmp independent_marked_check.json independent_marked_check.replay.json
    python3 -O independent_marked_check.py > independent_marked_check.optimized.replay.json
    cmp independent_marked_check.json independent_marked_check.optimized.replay.json
    python3 comparison_controls.py > comparison_controls.replay.json
    cmp comparison_controls.json comparison_controls.replay.json
    python3 probe_full_decorated_flip.py > probe_full_decorated_flip.replay.json
    cmp probe_full_decorated_flip.json probe_full_decorated_flip.replay.json

Source provenance is pinned to the exact commit above. No claim about unrelated THM-4725 or THM-4727 is promoted by this review.

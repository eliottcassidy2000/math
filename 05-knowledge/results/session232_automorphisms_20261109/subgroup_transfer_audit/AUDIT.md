# Session 232 independent subgroup and transfer audit

**Verdict: PASS.** The subgroup gap theorem, strict-window prime cost and equality characterization, threefold separation, sharp 5/7/25 strata, optional module-interaction inequality, and pointwise-stabilizer path-orbit recurrence are mathematically valid with their stated hypotheses. No substantive theorem repair was necessary. Review prompted a simpler three-run proof of the sharp strata, corrected a temporarily stale supplemental checker-coverage sentence, and clarified cyclic-extension notation. These refinements are incorporated in the reviewed producer supplement.

## Scope and dependency boundary

Reviewed producer inputs:

- `session232-aut-spectrum/PROOF.md` and `SHARP_STRATA.md`, plus both verifiers.
- `session232-paper-bridge/PROOF.md` and its verifier, specifically path-prefix orbit transfer, terminal stabilization, and the elementary group-action explanations for seven-coordinate cut states.
- The exact THM-4812 file named in those proofs, including its appended audit repairs and the relevant MISTAKE-681 correction. Parts (b), (c), and the established full-group gap from (d) are inherited. Part (g) is not used. This audit does not rerun the historical all-order certificates or claim an independent proof of their classical dependencies.
- The supplement's prior F2 credits exist at pinned commit `4b23190e9c39`: `session196/symmetry/PROOF.md` and `session197/orbit-wildcard/PRIME_CUBES.md`. Their stated constants and conventions match the present use.

External preprint interpretation and the independent graph-encoding theorem are outside this audit's mathematical scope. The subgroup producer's complete verifier also replays its graph controls, but that does not replace the graph worker's separate audit. No claim of complete odd spectrum, resolution of Chan–Pak, or new priority for general group-action methods is accepted or needed.

## 1. Subgroup prime cost and exact equality

Write G=Aut(T), take its actual vertex orbits O_i, and identify G with its faithful restriction subgroup of P=product Aut(T[O_i]). OCF gives H(T)>=product H(T[O_i]); free Hamilton-word action gives q_A=[G:A]q_G. Thus

    q_A >= [G:A][P:G] product q(T[O_i]).

For q_A<23, every vertex-transitive orbit tournament lies in {K1,C3,R5,P7} by the inherited strict-window lemma. With c,r,p counting nontrivial orbit types, |P|=3^c 5^r 21^p and product q=3^(r+2p). Therefore prime divisors of |A| belong to {3,5,7}, and q_A>=3^(v5(|A|)+2v7(|A|)). The group inclusions, valuations, strict threshold, and direction of inequalities are correctly typed.

Equality forces both indices to be one, r=v5(|A|), p=v7(|A|), and equality in the OCF partition bound. In particular there is no directed triangle meeting more than one orbit block. Every orbit tournament is strong or a singleton: an automorphism preserves the intrinsic linear order of SCCs, so a vertex-transitive tournament cannot have multiple SCCs.

The no-cross-triangle step is valid. If an outside vertex x sees a strong block B both ways, a path in B from an outneighbor of x to an inneighbor has a first orientation change, producing a cross triangle. Hence B is a module. Applying this to both blocks gives uniform incidence between blocks; the quotient has no directed triangle and is transitive. The orbit blocks are exactly the SCCs. Conversely an ordinal sum of the stated strong blocks has product automorphism group and multiplicative H. Thus the iff characterization, including arbitrarily many K1/C3 factors, is exact.

For the subgroup gaps: q_A in {7,21} makes q_G a divisor of that value. The inherited full-group exclusion leaves q_G=1, or q_G=3 in the 21 case, while forcing index divisible by7. Hence 7 divides |G|, contradicting the prime-cost lower bound q_G>=9. This does not incorrectly infer gap preservation from arbitrary multiplication. Nonfaithful external group actions must use their image, as the text says.

## 2. Threefold separation and sharp strata

Let b=3^(v5(|A|)+2v7(|A|)) and retain q_A<23. If q_A<3b, oddness forces both subgroup indices to be one and the orbit exponents to equal the acting-group valuations. Thus A=G=P and the orbit blocks are independently acting modules.

A directed triangle in the quotient, anywhere in that quotient, has three disjoint whole-block Hamilton-word families. Its expansion U consequently satisfies H(U)>=3H(X)H(Y)H(Z). OCF applied to U and all remaining orbit blocks gives q_A>=3b, a contradiction. The quotient is transitive and q_A=b. This proves the sharpened alternative q_A=b or q_A>=3b without any F2 constants or additional classification.

The 7, 5, and 25 strata follow with prices 9,3,9 and respective bounds q_A<23,<9,<23. The equality structures then follow from the main theorem. The optional stronger triangle-through-module inequality is also sound: the three whole-block families and the two-run M-B-C-M family are disjoint, with cardinalities 3H(M)H(B)H(C) and F2(M)H(B)H(C). OCF with outside blocks is legitimate. Here F2 is ordered and both paths are nonempty. Independently reproduced F2(R5)=90 and F2(P7)=1638, and cyclic expansions by two singleton modules have (H,|Aut|,q)=(135,5,27) and (2205,21,105). The unrestricted bound q>=105 for all non-equality seven-symmetry tournaments would be false; the producer correctly excludes it and retains carousel7 with q25 as hostile.

## 3. Orbit transfer: orientation and stabilizer typing

For a source endpoint-state orbit b represented by (S,u), N(b,a) counts legal *outgoing* appended vertices v whose new endpoint state lies in target orbit a. Let K_b=G_(S) be the pointwise stabilizer, not the setwise stabilizer of (S,u). Every ordered path with this state has stabilizer exactly K_b. The available extensions form a K_b-set; an appended vertex has stabilizer of order k_a. Hence every orbit has size k_b/k_a and

    W(b,a)=N(b,a) k_a/k_b

is a nonnegative integer. Deleting the last vertex is equivariant and gives each child path orbit a unique parent orbit. Consequently c_a=sum_b c_b W(b,a), initialized with one per singleton state orbit, is exact. Different chosen representatives only conjugate stabilizers and preserve their orders. Impossible endpoint states correctly retain zero path counts.

The groupoid-mass version m_b=c_b/k_b obeys m_a=sum_b N(b,a)m_b with initial mass1/k_b. The raw N matrix is therefore valid for this different measure. The incorrect hybrid uses initial coarse count 1 with N. The source=>C3 hostile gives3 rather than the correct1.

**Additional independent typing hostile:** on C3=>K1, the terminal endpoint-state stabilizer has order 3 but its pointwise/path stabilizer has order 1. Substituting the former into W gives terminal count 3 rather than1. This is a concrete reason not to replace pointwise Fix(S) by setwise endpoint-state symmetry, even when the latter sounds natural. The producer's formula already uses the correct quantity; this is an added regression control, not a correction to its theorem.

Terminal stabilization is valid and sharp: an automorphism fixing n-2 vertices could only transpose the remaining pair, which reverses their arc. Thus all relevant pointwise stabilizers are trivial once |S|>=n-2. A chain of n-3 vertices followed by C3 has residual rotation at the preceding nonempty layer when n>=4. For n=3 the corresponding nontrivial layer is empty, as stated.

The seven-coordinate divisibility explanation also follows from free action on ordered spanning paths and ordered spanning path pairs, with cut-invariant endpoint restrictions. The stated bilinear ordinal product is compatible with dividing by the specified product action. Rigidity of the marked core is necessary, not sufficient, for a unit gain; C3 with singleton cut supplies the stated converse hostile. This does not enlarge the inherited moving-core construction or provide a new unit-gain sufficiency theorem.

## 4. Finite exact checks and reproducibility

`independent.py` has no producer or repository imports. It builds adjacency matrices, explicitly enumerates automorphism groups and their literal subgroups, uses a pull subset recurrence for H, constructs state orbits by orbit removal, compares with DFS enumeration of actual ordered prefix paths, and checks both coarse-orbit and rational-mass recurrences. It also checks the deliberately wrong endpoint-state-stabilizer recurrence on the hostile above.

Independent universe and checks:

- All 1,099 labeled tournaments on orders 1..5.
- All 1,301 literal subgroup/tournament pairs; 315 equality cases.
- All 74,521 endpoint-state orbits; 328 transfer edges needing correction from N to W.
- Every set partition of those tournaments, filtered for strong/singleton blocks, tested for the no-cross-triangle implication; exactly 1,099 admissible no-cross-triangle partitions, the unique SCC partition per tournament.
- 52 low-window arithmetic cases for exceptional product-group orders with c<=6,r<=2,p<=1.
- All subgroups for selected R5, P7, carousel7, C3=>C3, R5=>R5, and P7=>C3, plus the diagonal C3 action on C3=>C3. The latter has H9, acting group3, and quotient3, distinctly from the full-group quotient1.
- Independent two-path-cover counts and cyclic singleton-extension controls described above.

The separate producer replay covers all 532 isomorphism classes through order 7, all 640 literal class-subgroup pairs, 168 strict-low-window prime-cost tests, and 33 equality cases. Its normalized-Hamilton-word census and exact canonicalization were read and checked; mass H/|Aut| is independently compared with literal words. The supplemental replay applies 168 threefold gates and all applicable small sharp-stratum cases, plus the two-R5 order 10 control.

Every final verifier was run normally and with `-O`; each output pair is byte-identical. Explicit failure sentinels call each verifier's exception gate with a false condition in both modes and require nonzero exit with RuntimeError. AST inspection finds no assert statements in these verifiers. These are finite controls, not extrapolation to an all-order theorem.

Reproduction from this directory:

    python3 -B independent.py > independent.normal.json
    python3 -B -O independent.py > independent.optimized.json
    cmp independent.normal.json independent.optimized.json
    python3 -B producer_subgroup_replay.py > producer_subgroup.normal.json
    python3 -B -O producer_subgroup_replay.py > producer_subgroup.optimized.json
    python3 -B producer_transfer_replay.py > producer_transfer.normal.json
    python3 -B -O producer_transfer_replay.py > producer_transfer.optimized.json
    python3 -B strata_replay/verify_strata.py > producer_strata.normal.json
    python3 -B -O strata_replay/verify_strata.py > producer_strata.optimized.json

`failure_sentinels.json` records all explicit negative-gate results; `MANIFEST.json` pins reviewed inputs and replay outputs. Producer replay copies are verbatim; only the independent verifier is independent code. All audit writes stayed inside this audit directory. No producer files or repository files were edited.

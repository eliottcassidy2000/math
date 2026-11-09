# The exact small seven- and five-symmetry strata

Session232 supplement to `PROOF.md`, Section3. **PROVED with the same inherited dependencies as Theorem A.** The sharp strata require no extra constants. The optional stronger module-interaction bounds use the explicit ordered two-path-cover constants independently verified below. This supplement does not modify the frozen primary verifier or substitute for its independent audit.

## 1. Types and constants

For a nonempty tournament M, define F2(M) to be the number of **ordered pairs of nonempty directed paths whose disjoint vertex sets partition V(M)**. Neither path is allowed to be empty. In particular this is not the conventions-with-empty-parts quantity H⊛H; that quantity is F2+2H.

For R5=Cay(Z/5,{1,2}) and P7=Cay(Z/7,{1,2,4}),

    H(R5)=15,      F2(R5)=90,
    H(P7)=189,     F2(P7)=1638.                  (1)

These constants are inherited, with explicit prior credit to [S196 symmetry, Section4](https://github.com/eliottcassidy2000/math/blob/4b23190e9c39/session196/symmetry/PROOF.md) and [S197 prime cubes, Section1](https://github.com/eliottcassidy2000/math/blob/4b23190e9c39/session197/orbit-wildcard/PRIME_CUBES.md). The supplement checker recomputes F2 in two ways: subset products, and every vertex permutation with every nontrivial cut, checking only the two path interiors, not the arc across the cut.

## 2. A normalized module-interaction lower bound

Suppose G=Aut(T) is the full direct product of the automorphism groups of its vertex-orbit tournaments T_i. Each O_i is then a module: an element supported on O_i is available independently, and transitivity on O_i makes every exterior vertex have uniform incidence with it. Collapse those modules to a quotient tournament Q.

Suppose three quotient vertices form a directed triangle M→B→C→M, and let U be its expansion. Write h_X=H(X). Selected Hamilton words in U give

    H(U) ≥ (3h_M+F2(M)) h_B h_C.                (2)

There are three whole-block orders MBC, BCM, CMB; each yields h_M h_B h_C words. Separately, choose an ordered nonempty two-path cover (P1,P2) of M, and a Hamilton word of each of B,C, and concatenate P1,B,C,P2. The three cross-module arc types make this a directed Hamilton word. Its maximal block-run pattern is MBCM, so different input covers give different words and none overlap the previous three families. This proves (2). Other run patterns only add nonnegative contributions.

Let all other orbit blocks be D_j. OCF positivity applied to the partition U,D_1,… gives

    H(T) ≥ H(U) ∏ h_Dj.

Dividing by the full product group order therefore yields

    q(T) ≥ (3+F2(M)/H(M)) ∏_i q(T_i).           (3)

This is a lower bound for a graph containing the triangle expansion, not an assumption that it has no other vertices or that every Hamilton path has one of the four chosen run patterns.

Every vertex of a strongly connected tournament lies on a directed triangle. To see this, take a shortest directed cycle through that vertex. If its length exceeds3, its chord from the first to the third vertex either shortens it or produces a directed triangle through the starting vertex. Equivalently, for a named vertex x partition the other vertices into N+(x),N−(x). Strong connectivity forces an arc from N+(x) to N−(x); otherwise no directed path can return to x from N+(x). Such an arc completes a triangle through x. The second argument directly covers every named vertex and is the one used here.

## 3. A threefold separation after the prime cost

**Theorem C (no intermediate cost).** Suppose q_A(T)<23 and set

    b=3^(v_5(|A|)+2v_7(|A|)).

Then either q_A(T)=b or q_A(T)≥3b. In the first case the complete equality characterization in Theorem A applies.

**Proof.** It suffices to consider q_A<3b. In the notation of Theorem A, let α=v_5(|A|), β=v_7(|A|), and let r,p count R5,P7 full-group vertex orbits. We have

    q_A ≥ [G:A] d 3^(r+2p),
    r≥α, p≥β.

Both indices are odd. If their product were greater than1, it would be at least3, contradicting q_A<3b. Likewise r>α or p>β would increase the power of3 by at least one. Hence A=G=∏Aut(T_i), r=α and p=β. The full independent actions make every orbit block a module, as in Section2.

If the quotient Q contained any directed triangle X→Y→Z→X, its expansion would contain the three distinct one-run Hamilton-word families XYZ,YZX,ZXY. Therefore H(U)≥3H(X)H(Y)H(Z). Apply OCF positivity to the partition consisting of U and every other orbit block, then divide by the full product automorphism order. This gives

    q_A ≥ 3∏_i q(T_i)=3b,

again a contradiction. The quotient is thus triangle-free and transitive. Consequently T is an ordinal sum of the orbit blocks and q_A=∏q(T_i)=b. QED.

This proof needs no F2 evaluation, no classification of low strong fibres, and no directed triangle through a distinguished module. The three one-run families are enough. The stronger Section2 count remains a useful optional test for interaction involving a particular module.

## 4. Exact five-, seven-, and twenty-five-symmetry strata

**Corollary D.** For A≤Aut(T):

1. If7 divides |A| and q_A<23, then q_A=9, A=Aut(T), and T is an ordinal sum of K1,C3 and exactly one P7.
2. If5 divides |A| and q_A<9, then q_A=3, A=Aut(T), and T is an ordinal sum of K1,C3 and exactly one R5.
3. If25 divides |A| and q_A<23, then q_A=9, A=Aut(T), and T is an ordinal sum of K1,C3 and exactly two R5.

**Proof.** Theorem A bounds the price b by q_A. In the three cases the allowed exponents force respectively b=9,3,9: a higher exponent would give at least27 (or81). In each case q_A<3b, so Theorem C gives equality. Apply Theorem A's equality structure. QED.

In the range below23, the only possible path-orbit count with a seven-factor or a twenty-five-factor in the acting group is exactly9. In particular no arbitrary automorphism subgroup realizes11,13,15,17 or19 with seven-symmetry.

**Optional stronger module-interaction bound.** If the orbit modules are one P7 and otherwise C3/K1 and a quotient triangle passes through P7, Section2 gives q≥105. If the orbit modules are one R5 and otherwise C3/K1 and a quotient triangle passes through R5, it gives q≥27. These exact constants are stronger for those distinguished-module interactions but unnecessary for the sharp strata. They are conditional on the stated direct-product orbit-module structure. We do not claim every tournament with seven-symmetry and q≠9 has q≥105: the carousel Q7 has q=25.

## 5. Sharp and hostile controls

- R5: (H,|Aut|,q)=(15,5,3), and P7: (189,21,9), attain both conclusions.
- R5→R5: (225,25,9), so Corollary D(2)'s strict upper threshold9 cannot be relaxed to≤9 with the same conclusion.
- Q7=Cay(Z/7,{1,2,3}): (175,7,25), a nearby hostile to any unrestricted assertion that seven-symmetry forces q=9.
- The C3 substitutions with one R5 or P7 module and two singleton modules (the closing arrow returns to the same module, not a second copy) have (H,|Aut|,q)=(135,5,27) and(2205,21,105). Their quotient is literally a triangle; here (2) is an equality because the two external blocks are singletons. These attain the module-interaction constants used in the proof.
- P7→C3: (567,63,9), shows that arbitrary triangle SCC padding is compatible with the exact stratum.

Reproduce with `python3 verify_strata.py` and `python3 -O verify_strata.py`. The small full-group census applies the index identity to every literal subgroup, so it tests the threefold-separation theorem and both nonvacuous small subgroup strata, not only the A=G specialization. The order10 R5→R5 control tests the twenty-five stratum. It does not rerun the inherited THM-4812 all-order finite certificates.

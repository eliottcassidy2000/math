# Session185: larger-component escapes and the next closure barriers

2026-11-04. **PROVED:** explicit clock-uniform returns crossing the previous SCC≤4 boundary; two further children from each new carrier; an all-size ordered-context obstruction; two normalized SCC5 valuation families; a complete fixed-front/transitive-tail inverse; a persistent three-state unit edge gate; and an exact 7-adic classification of one inherited ray. **OPEN:** a recursively closed two-digit carrier, the necessary coverage theorem, native-operation certificates for the new replacements, and unrestricted Hamilton-path spectrum completeness.

## 1. What was inherited, and what now changes

[Session184](https://github.com/eliottcassidy2000/math/blob/602eac4f1e7687032cecdeb0f5a1d23a2129ca04/00-navigation/S184.md) proved that, in its root-attachment architecture, every suffix with SCCs of order at most four has the same complete normalized return language from the common all-marked TT3 entry:

    {empty,2,4,22,42,422,4222}.

In particular neither first-generation sibling had a clock-uniform +4 successor within that grammar. Some SCC of order at least five was necessary, but the theorem did not construct one or exclude larger components.

The physical seven-coordinate signature and ordered-sum law are inherited from [Session182's interface proof](https://github.com/eliottcassidy2000/math/blob/ad09a24c8fd25eca9ea320d1b240c04990ada96b/05-knowledge/results/session182_physical_returns_20261104/induction/PROOF.md). For a marked suffix(B,W), a nonempty transitive prefix TT_t and the specified root attachment give

    H_t=c+2^t F,       c=h−J.

The full tournament is strong exactly when the terminal suffix SCC meets W. A uniform digit d∈{2,4}, with fixed clock shift q, requires and is implied by

    2^q F_child=3F_parent,       c_child=3c_parent+d,

at legal clocks t≥1 and t+q≥1. These tools and the global exclusions7,21 are antecedents, not new theorems of this session.

## 2. Actual SCC7/8 escapes, with both next children

The [closure packet](../closure/PROOF.md) supplies the two missing +4 edges by explicit replacement suffixes, retaining their actual adjacency and port masks:

- Address24: strong order7, rows(2,20,9,3,45,79,31), port mask42, signature(35,676,252,267,70,13,14). Thus(F,c)=(252,−35), reached from address2=(21,−13) with shift−2.
- Address44: strong order8, rows(2,4,33,23,71,155,47,95), port mask52, signature(89,2144,504,627,118,17,24). Thus(F,c)=(504,−29), reached from address4=(21,−11) with shift−3.

Each new carrier has both requested next children, also literal strong suffixes:

- 24→242 has(F,c)=(756,−103), shift0; 24→244 has(1512,−101), shift−1.
- 44→442 has(F,c)=(756,−85), shift+1; 44→444 has(378,−83), shift+2.

All six new suffix graphs and full signatures are retained in [witnesses.json](../closure/witnesses.json). At original entry clock4 the simultaneously legal count paths are

    107→323→973→{2921,2923},
    107→325→979→{2939,2941}.

The exact coefficient identities prove these edges for every legal clock. Negative shifts consume prefix vertices: the displayed finite paths require original clock at least4. One cannot simply iterate clock-consuming steps indefinitely from a fixed initial graph. Nor does the packet provide all third-generation addresses or successors from its four terminal states.

The [independent audit](../closure_audit/AUDIT.md) recomputes all nine source/new signatures by partial permutations, counts literal full attachments by direct path enumeration, and independently checks strongness and shifts. The [carrier audit](../carriers_audit/AUDIT.md) separately recomputes the six new physical witnesses. These are constructions beyond the former exclusion, not a minimality theorem for SCC size, an old-name-preserving transformation, or a certified native Hamilton-complement program.

## 3. Why adding ordered context cannot continue these two carriers

The [mining lemma](../mining/PROOF.md), with [independent audit](../mining_audit/AUDIT.md), proves a stronger negative result than failure of a finite append menu. If a marked suffix B has

    c≤−3,       min(F,G)>−2c−2,

then no ordered extension X→B→Y, retaining B and its marking, has constant3c+2 or3c+4. X,Y may have arbitrary size and arbitrary markings.

The ordered-sum law gives c_(X→B→Y)≤h_X h_Y c. Any nontransitive context has Hamilton count at least3, overshooting either target in the negative direction. If both contexts are transitive, any left marking costs at least F and any right marking at least G; if both are unmarked, c stays unchanged. The displayed inequality closes every case.

Both new carriers satisfy the criterion: for24 the threshold is68 versus F252,G267; for44 it is56 versus F504,G627. Thus their successful next children genuinely require changing/replacing the internal carrier rather than merely adjoining ordered marked context. This is consistent with the positive replacement witnesses above. The order7 carrier has a unique proper nontrivial homogeneous module{3→0}, with3 marked and0 unmarked, and minimum feedback arc number3. Module inflation is a concrete next direction, not an established returning rule.

## 4. SCC5 defeats uniform clock normalization, but exact constants recover one finite inverse

Two complementary results prevent a misleading extrapolation of S184's finite normalization.

First, the [closure packet](../closure/PROOF.md) exhibits a marked-first SCC5 seed followed by repeated20-vertex transitive blocks, with

    F_k=2^(6+20k),       c_k≡1+k (mod3).

Hence normalized v2(F) is unbounded even on infinitely many constants coprime to3. The terminal port is real, and no initial unmarked prefix can be removed. Its [audit](../closure_audit/AUDIT.md) independently validates the seed and20-bit joint coefficient and checks the all-k geometric and modular proofs.

Second, the [growing-carrier packet](../carriers/PROOF.md) reaches the specific odd coefficient63 required by the blocked siblings. A marked strong order5 front with signature(11,100,12,27,0,1,4), followed by the transitive port word

    u_m=(63·2^(20m+8)−3)/25,

has F_m=63·2^(20m+10). Its actual constant is

    c_m=11−11j(u_m)−(u_m mod2)−27u_m,

which tends to negative infinity. High dyadic valuation on the correct odd ray is therefore not a missing +4 return at constant−35 or−29.

For a fixed marked front(h,p,F,G,J,e,s) and any positive transitive tail mask u, the exact inverse uses

    F_new=F+pu,
    c_new=h−J−h j(u)−e(u mod2)−Gu.

Prescribing C bounds0<u≤floor((h−J−C)/G). Together with the empty-tail case and the odd-part coefficient test, this covers every tail length and fixed integer clock shift in that subclass. The complete strong-front bank through order5 gives no return at either missing constant with odd coefficient63. The [independent audit](../carriers_audit/AUDIT.md) reconstructs the bank and867 literal candidate suffixes.

The exclusion is for one marked strong front of order at most5 followed by a transitive tail. It does not settle arbitrary SCC≤5 words. The two valuation examples disprove uniform scale bounds; they do not disprove target-dependent bounds.

## 5. A different construction: three adjacent odd counts with a three-run channel

The [alternative packet](../alternative/PROOF.md) extends the inherited Session176 occupied-gap technique to three strong six-vertex seeds differing by two successive single-edge reversals. Inflating their common vertex by TT_r gives

    H_j(r)=4·3^r+9·2^r−3+2j,       j=0,1,2, r≥1.

Thus every size has a three-point odd interval, initially{27,29,31},{69,71,73},{177,179,181}. Its run profiles are(27+2j,21,4): both higher-run channels balance, while the three-run channel is nonzero. This explicit family is the increment; the general occupied-gap law and nonuniform-feedback evaluator were already known.

For an arbitrary nonempty substituted tournament A, the exact law is

    H_j[A]=(27+2j)H(A)+21q_2(A)+4q_3(A).

The edge steps are2H(A), so they remain unit steps only for transitive A. Also no fixed shift between the three states implements a uniform3H+d rule: the coefficients would require both3^q=3 and2^q=3. The [independent audit](../alternative_audit/AUDIT.md) checks the full grouped skeleton table, all225 small arbitrary-module substitutions, and the all-order proof.

The upper state has exact7-adic valuation1 for r≡4 or6 modulo42, starting with469=7·67 and3493=7·499. These are infinite strong-count progressions, not infinitely many proven prime cofactors. Only O(log X) family outputs lie below X.

## 6. The unique 7-divisible ray in the inherited completed tree

The [prime-lane packet](../prime_lanes/PROOF.md) extracts new arithmetic information from the old address422 word S1 S0 C1 C1. Its inherited(F,c)=(189,−91) gives

    H_t=7(27·2^t−13),       t≥1,       order=t+9.

Among S184's seven normalized rays, every F is divisible by7 and only this c is. Thus422 is the only possible7p ray in that completed tree; the other six rays never give a multiple of7. This does not classify larger-SCC replacements or other entry states.

The exact clock densities are

    density(v_7(H_t)=1)=2/3,
    density(v_7(H_t)=e)=2/7^(e−1),       e≥2.

A unique lifted root modulo3·7^(j−1) describes divisibility of27·2^t−13 by7^j, proving every valuation occurs infinitely often. Deterministic trial division certifies prime cofactors at exactly clocks1,4,8,13,19,28,32,35 through35. The last gives a44-vertex strong graph with H=6493990551461=7·927712935923, derived from the exact family formula rather than enumeration of that graph's paths. The [independent audit](../prime_lanes_audit/AUDIT.md) checks the literal base interface, lifting, densities, and every finite primality decision.

## 7. What is still needed

The next proof must control the physical carrier, both F and its joint-defect constant c, terminal ports, and a legal clock simultaneously. It must return both children to a class on which the same construction works again, then prove entry and coverage of the intended target counts. If native operations are required, their guards need separate certificates.

This session supplies actual escapes and sharper constraints on continuing them. It does not prove all odd counts except7,21 realizable, either complete prime lane, infinitely many prime cofactors on an exponential ray, or unlimited branching. All five component packets have independent audits and normal/optimized verification. Targeted source checks support the stated comparison with prior work; historical novelty is not established.

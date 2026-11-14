# Odd coordinates: a true 3+5-bit return, and a physically complete four-branch slice

## Two origins, three products

Write an odd number H as 2i+1 or 2m−1, so m=i+1. The first chart sends index 0 to H=1, and the second sends index 1 to H=1. The forbidden pair {7,21} is therefore labelled {3,10} or {4,11}, respectively. Transporting the ring structure of the integers gives the following operations on the odd set:

| Origin used | Addition | Multiplication | Additive zero | Multiplicative one |
| --- | --- | --- | --- | --- |
| 0 mapped to 1 | H+K−1 | (H−1)(K−1)/2+1 | 1 | 3 |
| 1 mapped to 1 | H+K+1 | (H+1)(K+1)/2−1 | −1 | 1 |

Neither is the ordinary odd-number product HK supplied by ordinal joins of tournaments. The odd set under native multiplication cannot be made into a ring by merely choosing a new addition: such a ring would need a multiplicatively absorbing odd z, but 3z=z forces z=0. A change of coordinate is real, but it does not secretly preserve each operation in the old chart.

## The pair is literally an arithmetic two-cycle

The [older signed inverse Collatz work](../collatz_mod6_20260917_extended_collatz_scc.md) identified the positive inverse-greedy minus-sheet cycle {4,11}. Its map on m coprime to 3 chooses the least k≥0 with 2^k m in {2,5} modulo 9 and sends m to (2^k m+1)/3. The six unit residues modulo 9 give k in {0,1,2,3}. Under H=2m−1 this becomes

    C(H) = [2^k(H+1)−1]/3.

The genuine arithmetic two-cycle is now 7→21→7 with k-word (3,0). The affine relation does not turn this into ordinary Collatz or give a tournament-preserving map.

For a desired image M=2u−1 with 3 not dividing u, put X=3u−1 and e=v₂(X). The full predecessor menu in the m-chart is X/2^k for 0≤k≤e and k<d, where d=2 for u≡1 mod 3 and d=4 for u≡2 mod 3. Indeed the previous accepted residue backwards in the modulo-9 cycle of multiplication by 2 is at distance two from 2 (2→1→5) and four from 5 (5→7→8→4→2). If k is at least that distance then a smaller exponent would have been legal. In the odd chart the resulting predecessor is (3M+1)/2^k−1.

This gives exactly C⁻¹(7)={21} and C⁻¹(21)={7,15,31,63}. The last three are all attained Mersenne counts; this is a concrete refutation of any proposed forward closure of the spectrum under C. For k=0 the reverse operation M→3M is physical by ordinal multiplication by a triangle, but the example 63→21 shows why it cannot be inverted on an arbitrary realization. Divisibility is not factor extraction.

## What the hidden 3+5 bits can mean exactly

On the frozen two-letter (3,0) word the return map is m→(8m+4)/9, or m−4→(8/9)(m−4). Modulo 2⁸ the map depends exactly on m modulo 2⁵, because 9 is invertible; different low five bits give different outputs, and all eight high-bit choices over a fixed lower five are forgotten. Its image is the single congruence class m≡4 mod 8. After three complete returns, all eight low bits of the input are overwritten; the stable low-eight output is 4. This is the genuine three-plus-five splitting of a register on one known branch, not a splitting of the additive group.

The exact selector stability has a second, ternary coordinate. For r≥1, a positive start runs the word (3,0) for r consecutive repetitions if and only if

    m = 4 + 9^r t, and 8^(r−1)t is 0 or 2 modulo 3.

For a single repetition start at 4+9t. After branch 3 the value is 11+24t; its branch is 0 exactly for t≡0 or 2 mod 3. After branch 0 the value is 4+8t. Each required next occurrence of branch 3 forces another factor 9 in m−4. With m=4+9^r t all but the last repetition have their next return forced; the last one sees the indicated remaining ternary digit. Thus the density for r repeats is exactly 2/(3·9^r). Only m=4 follows this word forever among integer starts; this does not classify basins or other cycles.

For any ℓ≥1, the complete ordered first ℓ greedy exponents are determined by m modulo 3^(ℓ+1), and that precision is sharp at the 3-adic unit level. Equal starts modulo 3^(ℓ+1) have the same first selector mod 9; each identical branch divides the difference by 3 and multiplies by a 3-adic unit. Conversely, varying m by multiples of 3^ℓ leaves ℓ−1 branches fixed and makes the last input cover the three lifts of its unit residue from modulo 3 to modulo 9; the k table is nonconstant on either unit class. Varying by 2⁸·3^ℓ has the same effect while conserving all eight binary bits. The fixed word therefore contracts binary register data but is fragile at the greedy ternary selector.

There are two different finite algebra shadows. The underlying additive group of Z/256 cannot split as Z/8 × Z/32: its element 1 has order 256 but the latter group's orders divide 32. The 3+5 short exact sequence with kernel 32Z/256 is not a direct product; carrying is real. The finite field F₂₅₆ also has no subfields F₈ or F₃₂, since subfield degrees divide 8. The split ring F₈ × F₃₂ does have 256 points but a different Galois dynamics: componentwise Frobenius has four fixed points, four orbits of length 3, twelve of length 5, and twelve of length 15, from component orbit counts (2 fixed, 2 triples) and (2 fixed, 6 quintuples). In contrast F₂₅₆ Frobenius has orbit lengths dividing eight (with 2,1,3,30 orbits of lengths 1,2,4,8 respectively by subfield subtraction). The useful precedent is the repository's [Paley/Frobenius work](../collatz_paley_bridge_20261001.md), rather than an isomorphism of finite fields. A split observer has lost the native binary carry.

### PRIMEGAME: what a fraction can and cannot see

The frozen loop is literally the fraction 8/9 on positive d=m−4. A native FRACTRAN program consisting of [8/9] runs floor(v₃(d)/2) steps; the actual greedy loop can stop one step sooner because of the final unit digit. More generally any fixed ordinary FRACTRAN program on d uses a finite prime alphabet (the primes in its numerators or denominators). If d=9^r u with u coprime to this alphabet and to 3, choose u≡1 and u≡2 modulo 3 by CRT. These have identical tracked valuation vectors and hence exactly the same program choices for their entire orbit, but different r-repeat eligibility for the affine greedy loop. This forbids a native fraction-only detector of that condition without a changed data representation. It does not forbid the [prior encoded FRACTRAN counterexample search](../session323_computability_and_native_descent_20261113/FRACTRAN_ORBIT.md), or the use of separate valuation registers as in PRIMEGAME.

## From an arithmetic predecessor to a genuine tournament

Let A be a finite tournament with marked vertex v, and replace v by a transitive block of size a≥1; every external vertex sees the entire block as it saw v. A Hamilton path compresses to one of F_j legal words using every external vertex once and j nonconsecutive occurrences of a block marker. To inflate a word, choose a surjection from the a labelled vertices to the ordered j visits, then read each visit in its forced transitive order. Thus

    H(a) = Σ_j F_j · j! S(a,j)

where S(a,j) is the Stirling number of the second kind. If all visits have j≤2, this is H(a)=C+D·2^a, with D=F₂ and C=F₁−2F₂. Conversely a ray affine in 2^a must have no legal j≥3: positive F_j would force growth at least of order j^a. A sole exit at v suffices for at most two visits, since each departure other than the final possible end requires the same outside vertex. If A is strong and has order more than one, so is each inflation, by projecting the paths through the unchanged exterior.

For large enough a and a predecessor exponent k, the affine inverse transforms a physical profile as an arithmetic pair:

    C + D·2^a  →  ((3C+1)/2^k−1) + (3D)·2^(a−k).

The new offset must be odd and the greedy phase must hold; finding a second graph with the transformed profile is a separate requirement. Only the odd part of D matters if a bounded shift of the inflation length is allowed.

### A one-apex obstruction and its exact physical bypass

For a single apex over a transitive order, code a vertex by 1 if it beats the apex, 0 otherwise. Hamiltonians may be specified by the increasing vertices before and after the apex. An inversion 0 at position j before a 1 at i>j contributes 2^(i−j−1), and the remaining cases are adjacent 10 gaps and two endpoint cases. For a word U 0^a V the coefficient of 2^a is (Z(U)+1)B(V), where Z(U) reads the zeros of U as high-to-low binary 1-bits and B(V) reads the ones of V as low-to-high weights. This follows by summing the inversion pairs crossing or leaving the repeated zero block.

If the coefficient is 3 times a power of two, its odd part must come from one factor. If Z(U)+1 has odd part 1, U has only a terminal run of zeros after harmless leading ones, and V must be a zero run followed by 11 and then harmless zeros. This reduces to 0^b11 with H=3·2^b−1. If it is the Z(U)+1 factor with odd part 3, the binary pattern for Z(U)=3·2^j−1 forces U, up to leading ones, to 01 followed by j zeros. Then V has precisely one 1 and the word reduces to 010^b1, with H=3·2^b+3. In both cases b varies with a. Leading ones and terminal zeros contribute only a neutral source/sink; repeated-one padding is the converse case. Therefore no fixed one-apex word with a single repeated constant-bit run gives 3·2^s+1 for an unbounded sequence of s: an unbounded meeting would force the same odd part of the slopes and identical intercepts.

Nevertheless that whole ray is physical from s≥2. Inflate root 0 in the strong frame with outgoing rows (2,12,9,17,7), the same H=13 source used for a 69 Fermat seam. Its visit numbers are F₁=13,F₂=6 and no third visit. Thus the inflated count is 1+6·2^a = 1+3·2^(a+1). For odd r≥3 the k=1 predecessor of M=2^r+1 is N=1+3·2^(r−1); choose a=r−2. This supplies an actual family above the scalar obstruction.

### A complete four-branch menu

For every odd a≥5, the attained strong target M=5·2^a+5 has all four predecessors allowed. The relevant 2-adic valuation is v₂(3M+1)=4, and (M+1)/2≡2 mod 3; thus k=0,1,2,3 exactly. Use the following literal two-visit strong frames (the indicated root is zero-based):

| Role | Frame bit rows | Root | (C,D) |
| --- | --- | --- | --- |
| M itself | (2,4,17,7,43,15) | 0 | (5,5) |
| k=1 | (2,12,25,97,11,23,55) | 0 | (7,30) |
| k=2, before ordinal ×3 | (2,4,9,51,39,71,31) | 6 | (1,20) |
| k=3 | (2,4,41,19,7,91,31) | 3 | (1,15) |

Their inflation exponents are a, a−2, a−4, a−3, respectively. Use a triangle ordinal join on M for k=0 and on the listed k=2 inflation; the two joins are not claimed strong. The odd predecessor values are exactly 3M, 7+30·2^(a−2), 3[1+20·2^(a−4)], and 1+15·2^(a−3). This is a complete four-way physical slice, not four residue witnesses.

The exact marked-frame atlas through order seven has only one active odd slope D (coprime to 3) in the intersection of the four needed normalized sets (5,D), (1,D), (7,3D), (1,3D): D=5. This is a statement about using these profiles and an ordinal multiplication for k=2, not exclusive of other graph grammars. A more attractive power-of-two source ray M=2^r+5 is physical with frame (2,12,17,37,43,71,31), root0, (C,D)=(5,32). For even r≥6, k=0,2,3 are physical by the elementary ordinal constructions and the H=13 frame; its k=1 ray 7+3·2^(r−1) has no two-visit strong frame through order seven with that normalized slope and offset. There may be an unrelated construction for that ray.

### General single-run classification and one further physical return

The constant-run obstruction above has a factorization form for any odd normalized slope D. Every factorization D=d₁d₂ into positive odd numbers gives a canonical shallow template U 0^b V. Take U to be the binary complement of d₁−1 at its ordinary bit length (empty for d₁=1), and V to be the reversed ordinary binary digits of d₂. The intercept is H(U0V)−2D. These are all intercepts in the single constant-run one-apex grammar at that normalized slope. Indeed Z(U)+1=2^u d₁ gives Z(U)=2^u(d₁−1)+(2^u−1): those u low bits extend the middle zero run, and extra leading ones are neutral. Likewise B(V)=2^v d₂ merely extends the middle zero run and extra trailing zeros are neutral. The one-bit run case is still covered by converse.

For the concrete class U=0·1^m and V=1·0^(q−1)·1, with m,b≥1 and q≥2, the same inversion calculation gives

    H = (2^m+1)(2^q+1)·2^b + (2^m−2^(q−1)+1).

Indeed internal inversions in U and V contribute 2^m−1 and 2^(q−1)−1, endpoint and adjacent-downstep terms give 4, and removing the −B(V) term accounts for the displayed intercept. Taking q=m+1 gives an infinite library of offset-one rays with slopes D_m=(2^m+1)(2^(m+1)+1). Adjacent slopes obey 3D_m=D_(m+1) only when m=1, by cancelling the common factor and solving 3·(2^m+1)=2^(m+2)+1. This gives D=15→45, useful after the k=3 predecessor of the complete menu above. Namely M₁=1+15·2^(a−3) has k=1 predecessor M₂=1+45·2^(a−4), which is strong via the word 011 followed by a−4 zeros and then 1001. That predecessor is exactly greedy, because the source m=(M₂+1)/2 is 1 modulo 9. The next coefficient 135 has no offset-one member in the whole single-run one-apex grammar: checking all ordered divisor splits gives intercepts −117, −27, −17, −5, 5, 17, 35, 85 (the set is also checked by the replay). Thus a further physical return must change the operation or frame, not merely continue this binary archetype. This is a structural obstruction only for the named one-apex grammar.

The standalone replay computes the compressed-word counts by enumerating all external permutations and admissible marker gaps, then independently counts literal Hamilton paths of each inflated control. It exhausts nonisomorphic frames through seven for the tiny atlas and checks four generations of the exact ternary return cylinders. These are partial templates; neither a global reverse closure nor the unrestricted Chan–Pak conjecture follows.

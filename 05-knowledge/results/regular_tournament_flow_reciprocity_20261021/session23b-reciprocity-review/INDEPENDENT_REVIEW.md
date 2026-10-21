# Independent audit: regular-tournament flow reciprocity

2026-10-21. Verdict: **PASS**, in the stated regular-tournament and rational-function scopes. The proof invokes classical cone reciprocity correctly; the tournament-specific fixed-dimension parity argument is valid. No geometric inversion identification follows or is claimed. All reviewer files are outside the repository; no Git operations were performed.

## 1. Proof audit

The reviewed candidate is `session23b-reciprocity/REGULAR_FLOW_RECIPROCITY.md`.

1. **Action convention.** With incidence columns `e_tail-e_head`, an ordinary map satisfies `B A=P B`, whereas an anti-map sends `(u,v)` to `(a(v),a(u))` and satisfies `B A=-P B`. This minus sign is essential and is present in the candidate.
2. **Relative interior and lattice.** For a regular tournament, the all-ones vector lies in every fixed circulation space and is strictly positive. Hence the fixed cone is full-dimensional in that fixed space. Each coordinate functional is nonzero there, since it evaluates to 1 on this vector, so a zero coordinate is on a proper supporting hyperplane. Thus the relative interior is exactly the strictly positive vectors, even if some inequalities are redundant. In the lattice `L_a=V_a intersect Z^m`, addition of the all-ones vector bijects closed-cone lattice points with interior lattice points. The weight change is exactly `m`, not the number of arc orbits.
3. **Classical input.** Beck and Sottile, *Irrational proofs for three theorems of Stanley*, Theorem 1, states the required rational-cone reciprocity with sign `(-1)^(dim K)`. I independently opened the primary paper at https://arxiv.org/pdf/math/0501359 . The theorem is applicable in a basis of `L_a`; specializing to total weight is legitimate because every nonzero ray has positive weight. This is an identity of rational functions, not a power-series substitution at zero.
4. **Ordinary parity.** A tournament automorphism group has odd order, because an involution cannot exchange the endpoints of an arc while preserving it. Nontrivial real eigenspaces of an odd-order operator come in conjugate pairs. Therefore ordinary fixed-space codimension is even.
5. **Anti-involution parity.** An anti-involution has exactly one fixed vertex and `r` transposed pairs. Its trace on the arc space is `r`; its induced operator on the vertex zero-sum space is `-P`, of trace zero. Hence its trace on the circulation space is `r`, its fixed dimension is `(d+r)/2=r^2`, and this has the same parity as `d=r(2r-1)`.
6. **General anti-map.** The candidate's vertex-cycle argument is correct. A shorter argument is available: if `a` is anti then `a^2` belongs to the odd-order automorphism group. If `a^2` has odd order `s`, then `a` has order `2s`, and `j=a^s` is an anti-involution. The commuting odd part `h=a^2` gives `Fix(a)=Fix(h|Fix(j))`, so the dimensions again have equal parity. This shorter version avoids needing the stronger vertex-cycle classification.
7. **Burnside transfer.** Ordinary and anti-coset averages have the same inversion sign and translation factor, so averaging is valid. The anti-coset average counts converse-fixed ordinary orbits, not all decorations over a self-converse underlying tournament. The order-one tagged-group convention is consistent with the companion packet.
8. **Denominator sign.** If `F=N/product(1-q^ell_i)` and `L=sum ell_i`, then `N(q^-1)=(-1)^(d-s)q^(m-L)N(q)`. The displayed parity condition `s=d mod 2` is necessary for calling the numerator palindromic. If `N(0)=1`, it forces degree `L-m` and leading coefficient 1. These claims concern the displayed representation and do not establish minimal recurrence order.

No correction is needed for validity. The simplification in item 6 and explicit clarification of volume normalization below are recommended presentation improvements.

## 2. Exact controls and limits

The independent checker uses explicit `require` failures, not Python assertions. It exhaustively constructs regular orientations, partitions them into full permutation orbits, and verifies every ordinary and anti action in the following census:

| Order | Labelled regular orientations | Isomorphism classes | Automorphism orders |
|---|---:|---:|---|
| 1 | 1 | 1 | 1 |
| 3 | 2 | 1 | 3 |
| 5 | 24 | 1 | 5 |
| 7 | 2640 | 3 | 7, 3, 21 |

All three regular order-seven classes are self-converse. The classes with automorphism orders 7, 3, 21 have respectively 720, 1680, 240 labelled realizations. The last is Paley7.

For every action the checker independently computes fixed dimension in two ways: by arc-orbit quotient incidence rank and by the rank of the full incidence-plus-fixed-equations matrix. A third verification uses the average exact character over the cyclic subgroup. It also checks the signed incidence intertwiner, anti-involution traces, odd-part decomposition, all-ones translation in the integral orbit-coordinate lattice, and the common parity. Small closed boxes and shifted strictly positive boxes are independently enumerated whenever there are at most 12 arc orbits.

For order seven the exact dimensions are:

- Ordinary identity: 15.
- Ordinary order 3: 5.
- Ordinary order 7: 3.
- Anti order 2: 9.
- Anti order 6, present in Paley7: 3.

The script also constructs all circulation vectors of total weight at most 8 by sums of simple directed cycles, deduplicating vectors rather than treating cycle decompositions as distinct flows. Completeness of this finite construction follows from integral cycle decomposition. It checks ordinary orbit counts and anti-coset Burnside against literal orbit sets. A separate direct enumeration of nonnegative edge vectors validates the complete slices through weight 8 for orders at most 5 and through weight 4 for order 7. No large-weight raw-flow enumeration is used.

Exact symbolic substitutions check the C3 series; the independently supplied R5 raw, rotation, ordinary-orbit, anti-fixed and extended-orbit formulas; and the three Paley7 fixed-cone formulas derived below. The R5 formulas are corroboration only, not premises in the all-order proof.

### Strong nonregular hostile

The strong order-four tournament has

`R(q)=1/((1-q^3)^2(1-q^4))`,

`K(q)=1/((1-q^6)(1-q^4))`.

Its raw series has inversion `R(q^-1)=-q^10 R(q)`, not shift 6. More strongly, `K(q^-1)=+q^10 K(q)`: the raw and anti-fixed dimensions are 3 and 2. Consequently `(R+K)/2` satisfies neither choice of common sign with shift 10. This independently shows why the parity argument, and not only a common interior shift, is needed before taking an extended-group average.

## 3. Three exact Paley7 fixed-cone series

These are additional small explicit generating functions, derived from fixed-lattice equations, not inferred from a coefficient fit. Take Paley7 on `Z/7`, with arcs `u -> v` when `v-u` belongs to `{1,2,4}`.

### Order-seven translation

Any nonidentity translation has three arc orbits, each of size seven and each itself balanced. Therefore

`U7(q)=1/(1-q^7)^3`.

It has fixed dimension 3 and inversion `U7(q^-1)=-q^21 U7(q)`.

### Order-three automorphism `i -> 2i`

In the seven arc-orbit coordinates `(a,b,c,d,e,f,g)`, each orbit has size three. Direct incidence reduction gives exactly

`e=a`, `c+d=a+f`.

For fixed nonnegative `a,f,b,g`, there are `a+f+1` choices of `(c,d)`. The total weight is `3(3a+2f+b+g)`. Using

`sum_(a,f>=0) (a+f+1)x^a y^f = (1-xy)/((1-x)^2(1-y)^2)`,

one obtains

`U3(q)=(1-q^15)/((1-q^3)^2(1-q^6)^2(1-q^9)^2)`

`     =(1+q^3+q^6+q^9+q^12)/((1-q^3)(1-q^6)^2(1-q^9)^2)`.

The fixed dimension is 5. The second presentation has a palindromic numerator of degree 12, matching `33-21`. The first presentation has an anti-palindromic numerator because it has six denominator factors, illustrating the denominator-parity warning concretely. Both give `U3(q^-1)=-q^21 U3(q)`.

### Order-six anti-map `i -> 3i`

The four arc orbits have lengths `(6,6,6,3)`. In their coordinates `(a,b,c,d)`, the sole balance equation is

`a+d=2c`.

Thus `a,d` have the same parity. There is a unique parametrization

`a=2u+epsilon`, `d=2v+epsilon`, `c=u+v+epsilon`,

where `u,v,b>=0` and `epsilon` is 0 or 1. The total weight is

`18u+12v+6b+15epsilon`.

Consequently

`V6(q)=(1+q^15)/((1-q^6)(1-q^12)(1-q^18))`.

This fixed cone has dimension 3 and is not unimodular in the fixed circulation lattice: the two parity classes supply its two fundamental-parallelepiped terms. Its inversion is `V6(q^-1)=-q^21 V6(q)`, with numerator degree `36-21=15`.

The checker reconstructs all arc orbits and incidence equations from the named maps. It also compares both `U3` and `V6` with direct bounded parameter sums through total use 80. The sums are low-dimensional and do not enumerate raw Paley7 flows to that weight.

## 4. Audited cumulative asymptotic consequence

Let `T` be fixed, self-converse and regular on `n=2r+1>=5` vertices. Let `F_T(M)` count ordinary flow-orbits of total use at most `M`, and let `K_T(M)` count those orbits fixed by converse. Then

`K_T(M)/F_T(M) = C_T M^(-r(r-1)) (1+O(1/M))`,

for a positive constant `C_T` depending on the tournament and cutoff. This statement is cumulative: it does not assert pointwise coefficient asymptotics.

### Precise constant and proof

For every ordinary or anti action `a`, write

`P_a=K_a intersect {sum_e x_e<=1}`,

`alpha_a=vol_(L_a)(P_a)`.

The measure is the lattice measure in which a fundamental parallelepiped of `L_a` has volume 1. It is **not** the alternative convention that multiplies this measure by `d_a!`.

The polytope `P_a` is bounded, rational, full-dimensional in its `d_a`-dimensional lattice space, and has positive volume. Lattice-point counting gives

`#(M P_a intersect L_a)=alpha_a M^(d_a)+O(M^(d_a-1))`.

This is the leading-volume part of Ehrhart theory. It also has an elementary proof here: choose a lattice basis and compare the integer count with unit cubes; only cubes meeting the finitely many dilated facets cause an error, and there are `O(M^(d_a-1))` of them. Thus no period-collapse hypothesis is required.

Every nonidentity tournament automorphism acts nontrivially on the full circulation space when `n>=4`. Indeed, for every three-vertex set, its underlying triangle has a signed circulation vector, with entries `+1` or `-1` according to the tournament's arc orientations. If the action were the identity on the circulation space, it would preserve the support of every such vector and hence every three-subset. A vertex permutation preserving all three-subsets of an `n>=4` set is the identity. The odd-order eigenvalue argument then gives fixed codimension at least two for every nonidentity ordinary action.

Writing `G=Aut(T)`, ordinary Burnside therefore yields

`F_T(M)=(alpha_id/|G|) M^d+O(M^(d-1))`.

For every anti-map `a`, its fixed space lies in the fixed space of its associated anti-involution, so `d_a<=r^2`. At least one anti-involution exists, and every such involution has dimension exactly `r^2`; hence the maximum is attained with positive volume. Anti-coset Burnside yields

`K_T(M)=(sum_(a anti, d_a=r^2) alpha_a/|G|) M^(r^2)+O(M^(r^2-1))`.

Dividing proves the displayed asymptotic, since `d-r^2=r(r-1)`, with

`C_T=(sum_(a anti, d_a=r^2) alpha_a)/alpha_id > 0`.

The same proof works for the capacity cutoff `0<=f_e<=M`, replacing `P_a` by `K_a intersect [0,1]^m`. These polytopes are full-dimensional because half the all-ones vector is an interior point. The exponent is unchanged; the constant generally changes.

### R5 check and order-seven prediction

The explicit R5 formulas give lattice leading coefficients

`alpha_raw=17/1555200`,

`alpha_ordinary=17/7776000`,

`alpha_converse=7/17280`.

Thus its cumulative ratio is asymptotic to `(3150/17) M^-2`. This is different from the separately proved pointwise R5 coefficient-ratio constant `2100/17`; the two cutoffs must not be confused. The exact constants are checked symbolically from the `q=1` leading pole constants divided by the respective factorials.

For every self-converse regular order-seven tournament, the cumulative exponent is `-6`. This conclusion does not require a closed raw numerator, a minimal recurrence, or control of all root-of-unity poles.

## Reproduction

Run `python independent_check.py > normal.json` and `python -O independent_check.py > optimized.json`, then compare the files. The receipt records the exhaustive small-tournament census, all fixed-dimension histograms, finite orbit counts, explicit Paley fixed-cone equations and coefficients, hostile controls, and R5 cumulative constants.

## Final packaging audit

The final main theorem changes only its status from a proof candidate to proved with independent audit and adds a pointer to `COROLLARIES.md`. Its mathematical argument is unchanged. The corollary file contains exactly the reviewed sections 3 and 4 above, preceded by a provenance-and-scope paragraph. Both supplied source hashes were verified against actual file bytes.

The local `session23b-reciprocity/check_corollaries.py` independently reconstructs the two named Paley arc-orbit incidence matrices using the opposite, equivalent incidence sign convention; verifies their ranks and reduced equations; derives the order-three weighted geometric sum; verifies the order-six parity parametrization expression; checks all three inversion equations; and checks the distinct R5 pointwise and cumulative constants. Its normal and optimized executions have been reproduced separately in this review directory and compared with the producer receipt.

For clarity, the corollary text's bounded-parameter check through total use 80 is performed by this review directory's `independent_check.py`. The smaller local `check_corollaries.py` supplies an additional exact equation-and-identity check, not that bounded enumeration. The two checks have different, explicitly distinguished roles.

`manifest.json` pins the final theorem, corollary, producer checker and receipts, plus every reviewer source and receipt. Reproduction passed in both normal and optimized Python; no Git operation was performed. Final packaging verdict: **PASS**.

# Explicit fixed cones and cumulative symmetry rarity

Session23b. Root-derived cumulative consequence, independently reviewed. The explicit Paley fixed-cone derivations arose in independent review and have additionally been checked by the root from actual arc orbits. No literature-priority claim is made. See REGULAR_FLOW_RECIPROCITY.md for the general theorem and definitions.

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


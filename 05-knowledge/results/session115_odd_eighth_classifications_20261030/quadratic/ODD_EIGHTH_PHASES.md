# Odd eighths are the complete cubic sector of the coherent Fano gauge family

**PROVED elementary classification and explicit transport.** This extends the inherited signed-XOR norm/Hamming-code interface with a character-measurement classification and a uniform single-defect formula. The Hamming code, norm-composition locus, Walsh transform and affine-code facts are not claimed as new. No major open conjecture is resolved, and no literature-priority claim is made.

## 1. The actual objects and inherited map

The [Session112 norm-sign theorem](../../session112_counting_interfaces_20261030/norms/NORM_SIGNS.md) fixes the seven line bits A,B,C,D,E,F,G and proves that the Euclidean norm-composing XOR tables are precisely

 E=1+B+C+D, F=1+A+C+D, G=1+A+B+D.

All bits are in F2. Each table is a sign gauge of the fixed base table A=B=C=D=0. A gauge function s:F2^3→F2 acts by

 F_s(x,y)=F_0(x,y)+s(x)+s(y)+s(x+y).

Two gauges induce the same table exactly when they differ by a linear character. Fixing s(0)=s(1)=s(2)=s(4)=0 picks the unique representative of each of the16 tables. Its values are

 s(3)=A, s(5)=B, s(6)=D, s(7)=C+D.

Write x=(x1,x2,x3), with integer labels x1+2x2+4x3. The normalized gauge polynomial is therefore

 s(x)=A x1x2+B x1x3+D x2x3+κ x1x2x3,
 κ=A+B+C.                                                    (1)

This explicit polynomial is the missing bridge: a coherent line table becomes a Boolean phase, not a density merely because both mention eight.

For a named character u define

 W_s(u)=sum_x (-1)^(s(x)+u·x),
 p_s(u)=Pr_uniform x[s(x)+u·x=1]=(8-W_s(u))/16.                (2)

Changing the gauge representative by a linear character shifts the labels u; hence the multiset of readings descends to the table. Forgetting the chosen base table or changing arbitrary basis signs changes the table and may change this profile. It is not an invariant of the abstract octonion algebra: all these tables are isomorphic.

## 2. Complete classification of all16 tables

There are exactly three sectors.

1. **Trivial normalized gauge:** A=B=C=D=0. Walsh spectrum is8 once and0 seven times; readings are0 once and1/2 seven times.
2. **Seven nonzero quadratic gauges:** κ=0, s≠0. Walsh spectrum has four zeros, three values+4 and one value−4. Readings are1/2 four times,1/4 three times and3/4 once.
3. **Eight cubic gauges:** κ=1. Walsh magnitudes are6 once and2 seven times. Every reading belongs to

     {1/8,3/8,5/8,7/8}.                                    (3)

   Seven tables have profile1/8 once,3/8 four times,5/8 three times. The remaining table, A=B=C=D=1, has3/8 seven times and7/8 once. Thus the four fractions occur across the family; they do not all occur in every table.

### Proof of the quadratic sector

For κ=0 the polar form B_s(x,y)=s(x+y)+s(x)+s(y) is alternating. A nonzero alternating form in dimension3 has rank2: its rank is even and it is not zero because a nonzero normalized quadratic has a nonzero quadratic coefficient. Its radical has dimension1. Squaring the Walsh sum gives

 W_s(u)^2=8 sum_(z in rad B_s) (-1)^(s(z)+u·z).

This is16 on the four characters agreeing with s on the radical, and0 on the others. Since s(0)=0, Fourier inversion gives sum_u W_s(u)=8. Consequently the four nonzero coefficients have three plus signs and one minus sign. This proves all claims without a numerical spectrum search.

### Proof of the cubic sector

For κ=1, choose

 v=(1+D,1+B,1+A).

The point indicator δ_v(x)=product_i(1+x_i+v_i) has the same cubic and quadratic terms as s. Thus s(x)=ε+r·x+δ_v(x) for unique ε and r. The normalizations force

 ε=1[v=0], r_i=1[v=e_i]+1[v=0].                            (4)

A direct Fourier expansion gives

 W_s(u)=(-1)^ε (8·1[u=r]−2(-1)^((u+r)·v)).                (5)

So u=r is the unique exceptional character. Its value is+6 unless v=0, when it is−6. Formula(4) shows v=0 exactly when A=B=D=1 and hence C=1. In that case every other coefficient is+2. In the other seven cases, Fourier inversion forces four+2 and three−2 coefficients. Formula(2) proves the complete profiles.

The characteristic eighths here are **uniform probabilities over eight labels**. They are not eighth roots of unity, odd residues modulo8, or Collatz content rates unless a separate predicate-preserving map is supplied.

## 2a. Complete relative-gauge classification under GL(3,2)

Let linear coordinate changes act on the quotient of zero-constant Boolean functions by linear characters: [s] maps to [s composed with g]. The16 classes have exactly four orbits, of sizes

 **1,7,7,1:** zero; nonzero quadratic; cubic defect v nonzero; cubic defect v zero.

For quadratics, every nonzero polar alternating form has rank2, and choosing a symplectic pair plus a radical basis identifies all seven classes. Any two quadratic functions with the same polar form differ by a linear function, already quotiented out. For cubics, the unique representative description modulo linear functions is

 [s]=[δ_v+δ_v(0)].

Linear changes carry v to another point (with the inverse according to the pullback convention). They fix zero and act transitively on the seven nonzero points. Cubic parity and the zero class separate the stated sectors, proving completeness. The two cubic profiles are therefore exactly the two relative-gauge orbits.

**Action sidecar:** this is an action on relative gauges, or equivalently on a pair consisting of a multiplication table and its chosen base table, with both transported together. It is not an assertion about bare index relabeling on16 tables while silently holding a non-invariant base table fixed. Abstract algebra isomorphism forgets all four of these sectors.

An independent short spectrum check is also available: in dimension3, odd truth-table weight makes each Walsh coefficient congruent to2 modulo4. Parseval says their squared sum is64. Thus exactly one has magnitude6 and the other seven magnitude2. This argument independently recovers the cubic magnitudes, but the point-defect representation retains the labelled signs and the uniform-dimensional extension.

## 3. Uniform classification in every binary dimension

Let V=F2^d, N=2^d, d≥3. Define a single-defect affine phase by

 f(x)=ε+r·x+δ_v(x), ε∈F2, r,v∈V.

Then

 W_f(u)=(-1)^ε (N·1[u=r]−2(-1)^((u+r)·v)).                (6)

Consequently its character-twisted negative-sign probabilities lie in

 {1/N,1−1/N,1/2−1/N,1/2+1/N}.                            (7)

Conversely, a Boolean phase has some coefficient |W_f(r)|=N−2 **iff** it differs at exactly one point from an affine phase ε+r·x. This follows immediately by counting the agreements and disagreements contributing±1 to that coefficient. Thus(6) is a complete classification of the family detected by that extremal Fourier test, not just a sufficient construction.

For d≥3, N−2>2, so the exceptional character is unique. Its sign determines ε. The signs of the other Fourier coefficients determine v, giving unique parameters and exactly2N² such functions. Distinct affine words differ on at leastN/2≥4 positions, another proof of uniqueness of their radius-one defect representations.

At d=3 there are2N²=128 such functions, exactly the number of odd-weight Boolean functions. Affine words have even weight, so their one-point defects have odd weight; counting proves they exhaust the odd-weight sector. This recovers the complete cubic case because the cubic ANF coefficient is the parity of the truth table.

**Dimension boundary:** this exhaustion is special to d=3. At d=4, the indicator of {0,e1,e2} has odd weight3 but Fourier magnitudes10 or at most6, never14. It is not a one-point affine defect. Also d=2 does not have a unique exceptional character: N−2=2, so the uniqueness assertion must not be extended there.

## 4. Exact arithmetic selection sidecar

Session112 identifies a coherent line table z with the squarefree product

 m(z)=19^A 2^B 59^C 29^D 5^E 89^F 7^G.

Its divisor sum has squareclass30 iff the table is norm-composing. Within those16 promised-valid bank selections, the odd-eighth sector is exactly

 v_19(m)+v_2(m)+v_59(m)=1 mod2,                            (8)

because this is κ=A+B+C. The other eight selections have quarter/half profiles. The unique7/8-profile selection is m=19·2·59·29=65018, with sigma(m)=108000=30·60².

This is a literal finite-bank map preserving the coherent-table predicate and the measured Fourier profile. It does not identify a divisor sum numerically with a Fourier coefficient. A syndrome-zero difference uses symmetric difference of prime selections; repeated prime multiplication is not this operation. Outside the stated bank, equation(8) has no asserted meaning for norm composition.

## 5. Loss boundary and use

- The three-bit norm syndrome certifies coherence but does not determine the profile: all16 tables have the same valid syndrome.
- The extra bit κ decides whether every character reading is an odd eighth. It does not determine which character has the exceptional coefficient or the point v.
- For a promised cubic table, the exceptional Fourier location/sign and the remaining signs reconstruct the point defect via(5). Reading only the unordered profile loses that location: seven tables have identical profiles.
- Every table composes the same Euclidean norm. Hence norm-composition alone cannot distinguish this phase structure.

This gives a concrete synthesis with the earlier continuation-aware theme: coarse validity, phase class, and labelled reconstruction require different retained data. No implication for Collatz convergence, unrestricted LRC, factoring complexity or higher-dimensional norm composition follows merely from the shared fractions.

## Reproduction

Run `python check_odd_eighths.py` and `python -O check_odd_eighths.py`. The checker uses explicit exception gates. It enumerates all16 normalized gauges and all128 gauges with s(0)=0, computes Walsh sums literally, independently forms line-table coboundaries, checks the arithmetic bank, exhausts all256 Boolean functions in dimension3, checks320 literal uniform-defect controls through dimension7, and verifies the dimension4 hostile. Proofs above supply the all-dimensional statement; finite tests are controls.

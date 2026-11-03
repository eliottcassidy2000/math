# A complete denominator-sixteen classification and its triangular generating function

**PROVED by the elementary argument below. FINITE-EXACT checks are separate.**
These are scoped deductions with no claim of publication priority. The known
eight-point graph at 29 is inherited from
[THM-4139, rational three-cycle and order-six lift](../../../01-canon/theorems/THM-4139-rational-three-cycle-order-six-lift-and-horizontal-carrier.md),
not a discovery here. Infinity is excluded throughout.

## D1. Every shifted triangular map

For every integer k, put

    G_k(b)=b(b+1)/2-k,  b in Z,
    T_r=r(r+1)/2,      r>=0,
    B(k)=#{b in Z : the G_k orbit of b is bounded}.

Then

    B(k)=4([k is triangular]+[k-2 is triangular])
         +2[k in {1,2,5,6}]+8[k=4].                    (D1)

Here triangular numbers are nonnegative and include zero. Consequently

    sum_(k>=0) B(k) z^k
      =4(1+z^2) sum_(r>=0) z^(r(r+1)/2)
        +2(z+z^2+z^5+z^6)+8z^4.                       (D2)

This is an exact formal generating function, not a guessed numerical fit.
The maximum is eight, attained exactly at k=3,4. Period three occurs only
at k=4. Every other cycle has length one or two. There are seven affine
functional-graph types including the empty graph: four points with two fixed
points; four with a two-cycle; their two six-point extensions; the eight-point
union of fixed and two-cycle components at k=3; and the eight-point three-cycle
graph at k=4.

For K>=0 let N(K)=1+floor((sqrt(8K+1)-1)/2), and set N(K)=0 for K<0.
The number of parameters 0<=k<=K having any bounded integer orbit is

    N(K)+N(K-2)-[K>=3]+[K>=4].                         (D3)

The only intersection of the two triangular ladders is k=3: solving
T_r-T_s=2 gives (r-s)(r+s+1)=4, hence r=2,s=1. The extra k=4 is outside
both ladders. Summing (D1) also gives a closed cumulative count with
multiplicity. Thus the supported parameters have density zero, although
there are infinitely many parameters in each of the two families.

## D2. Equivalent complete rational preperiodic classification

Let c be any odd integer and f_c(x)=x^2-c/16. Every rational preperiodic
point has denominator exactly four. If c is not 5 modulo 8 there are none.
If c=8k-3, the affine bijection x=(2b+1)/4 gives

    f_c((2b+1)/4)=(2G_k(b)+1)/4.                     (D4)

It therefore transfers the COMPLETE rational preperiodic set to D1's
bounded integer set. This does not classify bounded *real* orbits starting
at arbitrary rationals; p-adic escape and real escape are distinct.

In numerator coordinates a=4x, write F_c(a)=(a^2-c)/4. The exceptional
small cases are:

| c | Bounded odd numerators | Cycles, in numerator coordinates |
|---|---|---|
| -3 | +/-1,+/-3 | (1), (3) |
| 5 | +/-1,+/-3,+/-5 | (-1), (5) |
| 13 | +/-1,+/-3,+/-5 | (-3,-1) |
| 21 | +/-1,+/-3,+/-5,+/-7 | (-3), (7), (-5,1) |
| 29 | +/-1,+/-3,+/-5,+/-7 | (-7,5,-1) |
| 37 | +/-3,+/-5,+/-7 | (-7,3) |
| 45 | +/-3,+/-5,+/-9 | (-5), (9) |

For c<=-11, c=5 mod8, the set is empty. For every c>=53, c=5 mod8:

* If c+4=s^2, the set is {+/-(2-s),+/-(2+s)}. The two fixed
  numerators are 2-s,2+s; the other two map to them.
* If c-12=s^2, the set is {+/-(-2-s),+/-(-2+s)}. The two cycle
  numerators are -2-s,-2+s; the other two map to them.
* Otherwise the set is empty.

The two large-case alternatives are disjoint. Their simultaneous square
conditions would give a difference of odd squares of 16, hence c=21.
For all odd c the cardinality has the concise form

    4([c+4 is a nonnegative square]+[c-12 is a nonnegative square])
      +2[c in {5,13,37,45}]+8[c=29].                  (D5)

## Proof: an escaping exterior and six surviving candidates

**Denominators.** At an odd denominator prime p of x, the negative
valuation of x doubles under f_c, since c/16 is p-integral. At 2,
v_2(c/16)=-4. If v_2(x)<-2, its valuation doubles; if v_2(x)>-2,
the next valuation is -4 and subsequently doubles. Preperiodicity therefore
forces x=a/4 with a odd. For the next denominator also to be four,
v_2(a^2-c)=2, which forces c=5 mod8. Conversely that congruence preserves
the odd-numerator lattice. On this discrete lattice boundedness and
preperiodicity are equivalent.

For c<=-11, F_c(a)-a=((a-2)^2-(c+4))/4 is uniformly positive, so every
odd-integer orbit diverges. For c>=-3 let M be the largest positive odd
integer satisfying

    M^2-4M-c<=0.

If |a|>M, then F_c(a)>|a|, and the subsequent positive odd iterates grow
strictly without bound. All bounded orbits lie in {-M,-M+2,...,M}.
The cases M<=9 are exactly the ten parameters
-3,5,13,21,29,37,45,53,61,69; the last four have M=9. Direct
substitution gives the table above, the empty sets at 53 and 69, and
{+/-5,+/-9} at 61. These last three agree with the large-case rule.

For M>=11, write uniquely

    c=M^2-4-8j,   1<=j<=(M-1)/2.                     (D6)

These inequalities are exactly M's maximality and its defining inequality.
A bounded point must also satisfy a^2>=c-4M. But for |a|<=M-6,

    a^2 <= M^2-12M+36 < M^2-8M <= c-4M,

because M>9. Thus the entire bounded set is contained in

    {+/-M, +/-(M-2), +/-(M-4)}.                       (D7)

There is no approximation in this six-state restriction. Its three images are

    F_c(+/-M)     = 2j+1,
    F_c(+/-(M-2)) = 2j-M+2,
    F_c(+/-(M-4)) = 2j-2M+5.                         (D8)

The middle magnitude can land in D7 only when j=1, when its image is
-(M-4). That target then leaves D7, since M>=11. Consequently +/- (M-2)
cannot be bounded. Either of the other magnitudes can land in D7 only at

    j=(M-5)/2, (M-3)/2, (M-1)/2.

The middle choice sends both to the excluded magnitude M-2. At the first
choice, c=(M-2)^2+12 and the remaining four points contain the two-cycle
(M-4,-M). At the last, c=(M-2)^2-4 and the remaining four points contain
the fixed points M and 4-M. These are precisely D2's alternatives.
This proves completeness, not merely the construction of the displayed cycles.
Substituting c=8k-3 proves D1, D2 and D3. QED.

The entire candidate interval is forward invariant exactly at
c in {-3,5,13,21,29}: its lower image is F_c(1), so the additional
condition is c<=4M+1. Combining it with c>=M^2-4M forces M<=7;
checking the remaining congruent parameters gives the list.
For M>=11 only six candidates remain after one restriction, and two further
restriction passes suffice to find the surviving four or empty set. This is
an effective classification step, not an experimental stopping cutoff.

## Sources and scope

The [September trace/reversal work](../arithmetic_braids_20260917_geometry.md)
already identifies 29/16 as the unique fixed marked-cycle trace under affine
time reversal, and the [September continuation](../arithmetic_seams_20260921_dynamics.md)
already supplies its complete eight-point graph. The
[recent modulo-290 model](../session146_carrier_primegame_20261103/prime290/README.md)
already proves the 145-state modular census and its divergent-lift control.
Here the addition is the explicit all-parameter triangular classification,
proof by six surviving candidates, and its count formula, not a new discovery
of the period-three example.

**CITED context, not a dependency of this proof.** Eliahou--Fares,
[2019 preprint](https://arxiv.org/pdf/1903.08865), Proposition 2.1 proves the
denominator constraint; Theorem 3.16 excludes rational periods above three
when the parameter denominator has at most two distinct prime factors.
Their [published article](https://doi.org/10.4153/S0008439521000588) confirms
that scope. Thus the period bound for denominator sixteen was already known.
Their earlier [Ramsey paper](https://doi.org/10.1016/j.dam.2015.07.038)
already supplies a graph-theoretic approach; its abstract was checked, not
its full proof. No general Poonen conjecture or Collatz result follows here.

Reproduce `python -X utf8 verify.py` and `python -X utf8 audit.py`.
The first prunes exact numerator graphs for 12,514 parameters; the second
walks integer triangular orbits independently for k=-20..1500. Neither
finite scan substitutes for D6--D8's all-parameter argument.

# A smaller proof dependency before coefficient stopping

2026-10-04. **PROVED** elementary all-height families, exact macro recognition,
completed-route splicing and the 19-adic digit theorem below. **FINITE-EXACT**
for the declared controls. Universal Collatz coverage remains **OPEN**.
No historical novelty claim is made for common-future or inverse-fibre rules.

The constructive gain is a family of sound choices for the adaptive selector:
an input can acquire a strictly smaller proof dependency while every observed
actual iterate, and every observed coefficient, is still growing. The length
of that growing prefix is unbounded. A second construction supplies a complete
route for selected members and organizes them into recursive 19-adic fibres.

## 1. Inheritance and the point of the extra rule

The [guarded sibling grammar](creative_sibling_20260925.md) already proves
`U(4y+1)=U(y)`, its repeated version, and certificate transfer through a common
future. It also proves that retaining only the smallest sibling can discard
a certificate:241 and483 are the hostiles. The new selector must retain all
available certified alternatives, not silently identify them with one normal
form. The same note's increasing-source closure is an inherited finite
benchmark, not a count recomputed or enlarged without qualification here.

The [boundary compiler](collatz_boundary_compiler_20261004.md) handles a
contracting word by isolating one least positive source. Its first-crossing
atlas does not accept a prefix whose coefficient has never contracted. Our
new family reaches that precise missing operation: **a smaller common-future
dependency before any coefficient crossing**. It does not contradict the
first-stopping coincidence conjecture; the smaller number is a proof
dependency, not an actual iterate.

The [recursive entry note](entry_20260927_recursive.md) already uses a source
valuation to recognize an unbounded macro counter. The
[inverse-ray codec](inverse_ray_ternary_addresses_20261004.md) supplies
checked first-hit route objects and modular evaluators. The recent
[fixed-observation alias theorem](paley_geometry_observations_20261004.md)
is the hostile to dropping the actual source in favour of its residue.

Anchor: adaptive guarded selection. Niche: intermediate sibling contraction.
Wildcard: recursively address the completed subfamily modulo19. Live board:
**immutable source / actual prefix / smaller dependency / suffix certificate /
dyadic guard / 19-adic address**. The least-used sidecar is the exact position
where the supplied smaller certificate joins the actual orbit.

## 2. All-height virtual contraction theorem

Let `U(n)=oddpart(3n+1)` on positive odd integers. For each integer `s>=1`,
choose the unique integer j satisfying

\[
 (3/2)^j<4^s<(3/2)^{j+1}.                         \tag{1}
\]

The inequalities are strict by unique prime factorization, and `3s<=j<4s`.
Define

\[
 P=3^j,\quad Q=2^j,\quad R=Q4^s,\quad
 c_s=\frac{4^s-1}{3},\quad C=P-Q(1+c_s),
\]
\[
 r=(3R-C)P^{-1}\pmod{4R},\qquad0<r<4R.             \tag{2}
\]

Here the inverse is modulo the power of two4R. For every

\[
 n=r+4Rk,\quad k\ge1,\qquad y=\frac{Pn+C}{R},      \tag{3}
\]

the following assertions hold:

1. n and y are positive odd integers, `y=3 mod4`, and `y<n`.
2. The first j valuations of n are1, and the next is `2s+1`.
3. Every actual iterate through time `j+1` exceeds n, and every corresponding
   coefficient `3^i/2^A_i` exceeds1.
4. Nevertheless
   \[
    U^j(n)=4^s y+c_s,\qquad U^{j+1}(n)=U(y).       \tag{4}
   \]
   A supplied certificate for the smaller y therefore gives a certificate
   for n by an explicit common-future splice.

The least member `k=0` may also be used whenever the same exact test `0<y<n`
passes. No unproved universal assertion about those boundary members is
needed.

### Proof of guards, inequality, and actual chronology

From (1), `2R/3<P<R`. Also

\[
 C=P-(R+2Q)/3>0,
\]

because `4^s>2`, and `C<P<R`. Congruence (2) gives `y=3 mod4`.
The identity

\[
 P(n+1)=Q(4^s y+1+c_s)
\]

has right-hand bracket congruent to2 modulo4. Consequently
`v2(n+1)=j+1`. This proves that the first j odd steps have exact valuation1,
with

\[
 U^i(n)+1=3^i(n+1)/2^i,\qquad0\le i\le j.
\]

At time j, substitution gives the first equality in(4). Then
`3U^j(n)+1=4^s(3y+1)`, and `y=3 mod4` makes the next valuation exactly
`2s+1`, proving the second equality.

The first j coefficients are `(3/2)^i>1`. The next is `3P/(2R)>1`
by(1). Their positive affine carries show that the actual iterates exceed
their original source as well. On the other hand,

\[
 y<n\quad\Longleftrightarrow\quad n>C/(R-P).
\]

Since `C<P<R` and `R-P>=1`, the threshold is below R. Every positive lift
in(3) is above R. This proves all claims. The root certificate for y is a
separate input to the splice; the inequality alone does not supply it.

### The first explicit family

At `s=1,j=3`, `P=27,Q=8,R=32,C=11`, yielding

\[
 n=79+128k,\quad y=67+108k,\quad k\ge0.             \tag{5}
\]

The least member also passes the inequality. Its diagram is

```text
79 ->119 ->179 ->269 ->101
                       67 ->101
```

The lower line means that67 reaches the same101; there is no Collatz edge
from269 to67. The actual word of79 is `(1,1,1,3)` through the join, with
coefficient81/64>1. The supplied route for67 can be spliced after101.

## 3. Fixed-input recognition and comparison with a bounded atlas

For a supplied positive odd n, compute `j=v2(n+1)-1`. If `j<3`, this rule
does not apply. Otherwise set

\[
 s=\left\lceil\frac{\operatorname{bitlength}(3^j)-j}{2}\right\rceil.
\]

Check(1), congruence(2), and `0<y<n`. These are finite exact tests. At most
one s can succeed, since the values j in(1) strictly increase with s.
The macro compresses an unbounded actual prefix; it does not pretend that
the prefix duration or parameter bit length is bounded.

The source families are disjoint: each has a distinct value `v2(n+1)=j+1`.
Their applicability density among positive integers is

\[
 \sum_{s\ge1}2^{-(j_s+2s+2)}.                      \tag{6}
\]

Removing or restoring finitely many boundary members per s does not change
this density: all sufficiently high-s families lie in one shrinking
dyadic neighbourhood of-1. More explicitly, `j_s>=3s` bounds the tail
after s=L by `1/(124*2^(5L))`.

For every `s>=5`, `j_s+1>=18`. The entire family lies outside the inherited
first-coefficient-crossing atlas through16 steps, because no coefficient
crossing has yet occurred. This is a proved extension of **smaller-dependency
applicability**, relative to that named finite atlas. It is not new complete
root coverage, and it is not claimed disjoint from every earlier grammar.

## 4. A completed subfamily with a ternary address

The inherited sibling grammar supplies

\[
 Y_b=\frac{32\,64^b-5}{9},\qquad
 Y_b\xrightarrow{1}\frac{16\,64^b-1}{3}
       \xrightarrow{4+6b}1,\quad b\ge0.            \tag{7}
\]

Every Y_b is3 modulo16. Fix s and its chart. Put `y_0=(Pr+C)/R`.
There is a unique `b_0 modP`, `0<=b_0<P`, for which

\[
 Y_{b_0}=y_0\pmod P.                              \tag{8}
\]

Indeed `v3(64^d-1)=2+v3(d)` for positive d, so
`v3(Y_b-Y_a)=v3(b-a)`. Thus b maps bijectively onto Y_b modulo every
power of3. Equivalently, the next ternary digit is a nonzero linear
function of the next b digit. The script recovers b_0 digit by digit,
testing the three candidates independently at each stage.

For `h>=1`, take

\[
 b=b_0+Ph,\qquad N_h=\frac{R Y_b-C}{P}.             \tag{9}
\]

Equation(8), together with `Y_b=y_0=3 mod4`, puts N_h in the source
cylinder(2). It is above R: `b>=P>j+2s`, hence `Y_b>R`; then `C<P<R`
implies `N_h>=Y_b>R`. The theorem therefore applies, and(7) supplies
the previously missing certificate for the smaller dependency.

The exact first-hit word is

\[
 \boxed{1^j,\ 2s+1,\ 4+6b.}                       \tag{10}
\]

Its odd rank is `j+2`, and its ordinary rank is `2j+2s+6b+7`.
All preterminal odd iterates exceed N_h; the final odd iterate is1.
The constructor can also build completed certificates at `h=0`. Smaller-source
assertions there are restricted to the explicitly checked instances; the
constructor does not assert that inequality for every `h=0` chart. These
instances are not needed for the all-height tail statement. The first
addresses `(s,j,b_0)` are
`(1,3,26),(2,6,459),(3,10,10798)`.

This is a concrete family that contains its route home. The selector's
separate job on an independently supplied integer is still to find a valid
rule and certify its smaller obligation.

## 5. Nineteen-way recursive refinement of each completed family

In(9), changing h multiplies the normalized exponential by `H=64^P`.
Since `P=3^j` and

\[
 64^3=1+19\cdot13797,\quad13797=3\pmod{19},
\]

we obtain

\[
 H=1+19P\pmod{19^2},\qquad v_{19}(H-1)=1.          \tag{11}
\]

Consequently, for distinct nonnegative h,h',

\[
 \boxed{v_{19}(N_h-N_{h'})=1+v_{19}(h-h').}         \tag{12}
\]

All coefficients multiplying the exponential difference in(9) are units
at19, so lifting the exponent proves(12). The whole completed family has
one fixed first digit modulo19. Modulo `19^a`, its parameter h runs
bijectively through all `19^(a-1)` extensions of that digit.

The recursive digit rule is particularly simple. For `0<=d<19`,

\[
 N_{h+d19^t}-N_h
 =\gamma_s d19^{t+1}\pmod{19^{t+2}},\qquad
 \gamma_s=32R64^{b_0}/9\pmod{19}.                  \tag{13}
\]

Here gamma_s is nonzero. Divide the next target-digit error by gamma_s
to select the next parameter digit. The resulting parameter class has
infinitely many nonnegative members; adding its period chooses `h>=1`
without changing the requested residue.

This is an actual recursive certificate constructor, not a formal orbit:
every selected h produces(10). Its depth-one fibre depends on s, so this
one construction does not claim all19 initial residues. The companion
[four-step route bank](mod19_route_lifts_20261004.md) supplies a separate,
sharp bounded-rank construction covering all of them. Neither statement
identifies a supplied integer with a constructed integer merely because
their residues agree.

The larger virtual source cylinder(3), before requiring the completed suffix,
also meets every residue modulo every power of19. This follows by solving
one linear congruence in k, because `4R` is a unit modulo19. Such a source
has a smaller dependency; its completed certificate still requires the
suffix or another supplied child proof. These are deliberately different
coverage predicates.

## 6. Reproduction, controls and connection ledger

Run:

```text
python -X utf8 -B 04-computation/experiments/virtual_contraction_ladders_20261004.py
python -O -X utf8 -B 04-computation/experiments/virtual_contraction_ladders_20261004.py
```

The [script](../../04-computation/experiments/virtual_contraction_ladders_20261004.py)
and [output](virtual_contraction_ladders_20261004.out) declare these universes:
30 charts; eight positive lifts each, with every actual prefix replayed;
90 completed symbolic certificates and seven literal controls below100000
bits; all target residues in the selected first-digit fibres for eight
charts through19-adic depth3; and an unexpanded depth24 residue request
at s20, producing a certificate of odd first-hit rank70. The codec's independent
modular evaluator checks the source formulas, including mixed moduli.

| Source -> target | Preserved predicate | Lost datum / required repair | Cheap decisive control |
|---|---|---|---|
| Actual growing prefix -> smaller sibling dependency | A real common future and strict source rank | The smaller node is not an actual iterate; keep the join |79 and67 meet at101|
| Source -> macro address s | Exact dyadic guard | Source magnitude and boundary inequality remain necessary | `recognize` verifies `0<y<n` |
| Supplied rank-two suffix -> complete long-prefix route | First arrival at1, exact ordinary and odd ranks | A residue match alone is not exact source membership | Independently replay the source formula and AST |
| Completed family ->19-adic parameter tree | One compatible source fibre at every depth | Integer height and supplied-source identity | Equation(12), plus the fixed-observation hostile |

The next useful test is whether the adaptive selector can repeatedly choose
these or other strictly smaller dependencies from the sole root seed. Finite
closure and supplied-child compilation are separate measurements. A pending
source remains pending; neither a finite residue address nor an available
formal macro is itself a completed certificate.

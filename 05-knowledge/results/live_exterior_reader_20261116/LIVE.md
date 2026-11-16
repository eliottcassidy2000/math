# Paired exterior choices improve the absolute reader error

**Proved.** For every sufficiently large odd N there is a literal
binary reader tournament with
```
|H − N| = O(N^(2/5)).
```
The argument keeps both exterior coordinates live and chooses them
jointly. It does not claim equality or density one.

Use the synchronized reader with I=−11, Z=82, g=9. The exact count
at address b of n bits is
```
H = I + D(Z+gb/D)(2−rev_n(b)/D) + gE
D = 2^n
0 ≤ E ≤ O(n).
```
See the [binary reader proof](../variable_reader_coverage_20261116/PROOF.md)
and [digit-reversal defect](../harmonic_reader_residual_20261116/RESIDUAL.md).
Set n=floor(log₂(N/88)), m=floor(n/2), k=n−m, B=2^m, C=2^k.
Thus C/B is 1 or 2 and D=BC. Write
```
b = BU+rev_m(X),      rev_n(b) = CX+rev_k(U)
0 ≤ U<C,              0 ≤ X<B.
```
Take L=2^floor(m/5) and restrict U and X to multiples of L.
The reversed remainders are then <B/L and <C/L respectively.
Hence, uniformly on the entire sublattice,
```
H = I + (CZ+gU)(2B−X) + O(B/L + log B).
```
A fixed choice of both U and X would leave an error of order B.
Instead, solve the continuous equation for one in terms of the other:
```
X*(U) = 2B − (N−I)/(CZ+gU).
```
Requiring X to be a multiple of L becomes a fractional-part problem
for f(j)=X*(Lj)/L with j integer.

## Uniform macroscopic interval

Write u=U/C, x=X/B. For all large N the normalized target
(N−I)/D ranges in [88,177]. At the two opposite corners of
[1/50,49/50]^2 the expression (82+9u)(2−x) is respectively
less than 88 and greater than 177. The line segment between
these corners provides one interior solution by continuity.
As the graph x=x*(u) has slope g (N−I)/[D(Z+gu)^2] bounded
above by 1, an interval of u of fixed positive length centered
at that solution stays in [1/100,99/100]^2. We use its middle
half to leave margin for the choice below.

Let J=B/L. Derivatives in the corresponding j interval satisfy
```
f'(j) = g(N−I)/(CZ+gLj)^2,
f''(j) = −2g² L(N−I)/(CZ+gLj)^3,
|f''(j)| ≍ 1/J.
```
The interval has length ≍J, so the range of f' has length
bounded below by a positive absolute constant. Thus one can
choose a rational r/s in that range, with s uniformly bounded.
For instance a fixed denominator whose reciprocal is less than
half the uniform range works. Let j* solve f'(j*)=r/s. Let j0
be the smallest integer ≥j*. On j=j0+sh set
```
F(h)=f(j0+sh)−rh.
```
Derivative F' is nonpositive on h≥0, with |F'(0)|=O(1/J)
and |F''|≍1/J (including the fixed denominator in the constants).
On 0≤h≤K√J the decrease is at least 2 once K is a sufficiently
large constant; that interval stays inside the reserved margin
for large J. Adjacent integer evaluations differ by O(J^−1/2).
Consequently this decreasing sequence crosses an integer with
error O(J^−1/2). Since rh is integral, f(j) is within the same
distance of an integer. Choose X to be L times that integer,
automatically in the permissible range. This argument can target
any prescribed fractional phase; here that phase is 0.

The main count misses N by at most O(B L/√J). With L≈B^(1/5),
```
B L/√J = √B L^(3/2) = O(B^(4/5)),
B/L = O(B^(4/5)).
```
Since B≍√N, this proves the claimed bound. One can split the
macroscopic interval into any fixed number of smaller intervals
and run the proof in each, yielding independently located compatible
exteriors (constants then depend on the number of intervals).

## Why this is not yet exact

Price's reciprocal chart supplies the right curve; a rational slope lets
integer samples revisit all fractional phases in a short segment instead
of selecting one point. However the lower reversed bits and the gap
between a crossed phase and its adjacent samples are still real integer
errors, of the traded-off scales B/L and √B L^(3/2). The proof
has not retained the right phase and the actual residual reader state
simultaneously. To close the remaining gap, one could use the retained
curve to impose the sign-compatible interface as a further condition,
or replace general phase crossing by a physically exact descent
intertwining this error with the next scale.

## A growing number of alternatives

The same proof can request M surviving address pairs, where
1≤M≤B^(1−η) for any fixed η>0. Choose a dyadic
L≍(B/M)^(1/5). Run until the decrease in F reaches 2M+1;
this requires O(√(MJ)) steps and individual steps have size
O(√(M/J))=o(1). The reserved interval remains available because
M/J=O((M/B)^(4/5))=o(1). Cross M distinct integer levels of F,
each yielding its own pair; successive crossings are distinct
once individual steps are <1. The bounds become
```
|H − N| ≲ B/L + √B L^(3/2)√M + log B
        ≲ B^(4/5) M^(1/5).
```
This counts distinct exterior addresses, not necessarily different
isomorphism classes or different realized Hamiltonian totals.

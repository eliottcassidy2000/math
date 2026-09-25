# Primitive triples seed square-sum paths with an exact translation clock

Status: **PROVED** elementary identities, complete local path structure, and stopping bounds; **FINITE-EXACT** computational controls. Date: 2026-09-25. This is the independent square/Pythagorean lane of the Zenodo discussion; source interpretation is handled by the coordinating lane. No source-priority claim, Hamiltonicity implication, or Collatz convergence claim is made.

## 1. Recovered coordinates and the small-number chain

Closest proved mechanism: [THM-3333, Gaussian-square light cone](../../01-canon/theorems/THM-3333-gaussian-square-farey-pythagorean-triangular-light-cone.md) and [collatz_mod6_20260917_pythagorean_semicircle.md](collatz_mod6_20260917_pythagorean_semicircle.md), theorems 1, 3, 5, and 7, supply the Euclid, odd-root, and Thales charts. [ternary_bridge_20260925.md](ternary_bridge_20260925.md) explains why marked square gaps, their signs, and primitive content matter for dynamical transport. [decoder_prime_square_20260925.md](decoder_prime_square_20260925.md) gives the exact square-class toggle under halving and the repaired endpoint-sensitive `Q_24` obstruction.

Canonical hostile: one square-sum edge does not specify a rational right triangle, and a local path does not cover a whole graph. Corrected near miss: degree two forces two path edges only after an interior-vertex condition is established. Least-used sidecar here: the square index labeling each edge reflection. Live board: odd-root gaps; Thales projection; square-edge reflections; odd-chain translation; primitive scaling; collision and positivity guards.

Write a primitive Pythagorean triple as

\[
a^2+b^2=c^2,\quad a\text{ odd},\quad b\text{ even},\quad a,b,c>0.
\]

Its Euclid parameters are `a=m²-n²`, `b=2mn`, `c=m²+n²`, with coprime opposite-parity `m>n>0`. Its odd-root gaps are

\[
c+b=(m+n)^2,\qquad c-b=(m-n)^2.
\tag{1}
\]

For `(a,b,c)=(3,4,5)`, these are `9` and `1`. Hence the numbers in question have exact geometric meanings:

\[
9=5+4=3^2,\quad 14=9+5=2(3+4),\quad 25=9+16=5^2,
\quad 11=25-14.
\]

But `9+5=14` is not a square-sum edge. The valid local edge involving the last two numbers is `11+14=25`.

Two further root identities are rigid, rather than universal recursions. For any positive real right triangle,

\[
(c+b)+c=2(a+b)
\iff 2c=2a+b
\iff (a,b,c)=\lambda(3,4,5),\quad\lambda>0.
\tag{2}
\]

Indeed `c=a+b/2` in the Pythagorean equation gives `3b²/4=ab`, hence `a=3b/4`. Thus among primitive integer triples, the user's identity selects exactly `(3,4,5)`.

Likewise, with semiperimeter `p=(a+b+c)/2` and inradius `r=(a+b-c)/2`, the area identity is `ab/2=rp`. Therefore `ab=a+b+c` is equivalent to `r=1`. In primitive Euclid coordinates `r=n(m-n)`, so `r=1` forces `n=1,m=2`. A general integer right triangle is an integer multiple of a primitive one; its inradius scales by that same integer, so this remains the unique positive integer solution. The root has product and perimeter `3*4=12=3+4+5`, semiperimeter 6, and `14=2(p+r)`.

## 2. A lossless primitive-triple to square-path map

Let `Q_N` be the simple graph on `1,...,N`, with distinct vertices adjacent exactly when their sum is a positive square. Every primitive triple above determines the ordered list

\[
\mathcal P(a,b,c)=
\bigl(a^2,\ b^2,\ a^2+2c+1,\ b^2-2c-1,\ a^2+2,\ b^2-2\bigr).
\tag{3}
\]

**Theorem.** This is a six-vertex simple path in `Q_(c²-5)`, and its consecutive edge sums are

\[
c^2,\quad(c+1)^2,\quad c^2,\quad(c-1)^2,\quad c^2.
\tag{4}
\]

The map is injective when the ordered path is retained: recover `a` and `b` by taking positive square roots of its first two entries, and `c` from their sum.

**Proof.** The five sums follow by direct expansion. Since `a,c` are odd and `c>a`, one has `c-a>=2`, so

\[
b^2=(c-a)(c+a)\ge2(c+a),\qquad b^2-2c-1\ge2a-1\ge5.
\]

All entries are consequently positive; the same inequality bounds the third entry by `c²-5`. The second is at most `c²-9`; the remaining bounds follow immediately. Primitivity gives `b=0 mod4` and `c=1 mod4`, so the six entries have respective residues

\[
(1,0,4,5,3,6)\pmod8.
\tag{5}
\]

They are therefore distinct. This proves the theorem.

At the root, the actual path is

\[
9\;--\;16\;--\;20\;--\;5\;--\;11\;--\;14,
\]

with edge labels `25,36,25,16,25`. This connects the user's numbers through genuine square edges while preserving their different roles as vertices and square targets. In particular 25 is a square target here, not a vertex of this six-vertex path.

The induced subgraph is also classified exactly. Square residues modulo eight are only `0,1,4`. Among the ten nonconsecutive pairs in (5), only positions 1 and 5 can have a square sum. Thus the induced graph is precisely the displayed path unless

\[
2(a^2+1)\text{ is a square},
\quad\text{equivalently}\quad a^2-2h^2=-1\text{ for an integer }h.
\tag{6}
\]

In that case it has exactly one extra edge, joining the first and fifth vertices. The first hostile by hypotenuse is `(7,24,25)`, where `49+51=100`. There are infinitely many such controls: odd solutions of the negative Pell equation are generated, for example, by powers of `1+sqrt2` of odd exponent; each odd `a>=3` also occurs in the primitive triple `(a,(a²-1)/2,(a²+1)/2)`. No Pell completeness claim is needed for this existence argument.

## 3. Why four square reflections give n+2

For a positive integer square index `s`, define the affine reflection

\[
R_s(x)=s^2-x.
\]

It specifies an edge of a positive square graph only when both endpoints are in the graph and distinct. Algebraically,

\[
R_vR_u(x)=x+v^2-u^2.
\]

Consequently, for integers `c>d>=1`, the four successive square indices `(c,c+d,c,c-d)` give

\[
R_{c-d}R_cR_{c+d}R_c(x)=x+2d^2.
\tag{7}
\]

This is the second difference of the quadratic function `s²`. For `d=1`, it is an exact realization of the odd-chain operation `n -> n+2` by four square-edge reflections. Formula (3) is the first four-reflection block seeded at `a²`, followed by one more `R_c`.

The distinction between scaling and translation matters. Scaling a square-edge picture by `t²` obeys

\[
R_{ts}(t^2x)=t^2R_s(x).
\]

It scales the index gap from `d` to `td` and the translation from `2d²` to `2t²d²`. In particular `d=2` gives translation by 8, not multiplication by 8. Positivity does not survive arbitrary gap changes: the root with `d=2` gives the algebraic walk

`9 -> 16 -> 33 -> -8 -> 17`,

so it is not a path of positive integers. The unit index gap is a necessary part of the valid primitive-triple construction.

## 4. Maximal repetition: exact collision and positivity bounds

Fix a primitive triple and use the `d=1` macro repeatedly. Put `x=a²`, `y=b²`. Block `j>=0` consists of

\[
\begin{aligned}
A_j&=x+2j, &B_j&=y-2j,\\
C_j&=x+2c+1+2j, &D_j&=y-2c-1-2j,
\end{aligned}
\]

followed by `A_(j+1)`. Define

\[
M=\frac{b^2}{2}-c,\qquad
L=\frac{b^2-a^2-2c-1}{2},\qquad h=L/2.
\tag{8}
\]

Here `M>=3`, and `h` is an odd integer, possibly negative. A whole block is positive exactly for `j<M`, because `D_j=2(M-j)-1`; in that range every vertex also lies below `c²`.

Within each of the four sequences, vertices never repeat. Parity excludes intersections between the odd sequences `A,D` and the even sequences `B,C`. The remaining cross-layer equalities are exactly

\[
A_j=D_k\iff j+k=L,
\qquad B_j=C_k\iff j+k=L.
\tag{9}
\]

If `h<0`, there are no collisions at nonnegative indices. If `h>=0`, the earliest is in block `h`, because an earlier collision would require two indices below `h` summing to `2h`. At that block,

\[
B_h=C_h=(c+1)^2/2.
\]

The reflection would be a forbidden self-loop. Moreover

\[
M-h=(c-1)^2/4>0,
\]

so this collision precedes the positivity bound.

It follows that the maximum number of complete simple blocks is

\[
q=\begin{cases}M&h<0,\\h&h>0.\end{cases}
\tag{10}
\]

The maximal simple prefix inside `Q_(c²-1)` has exactly `4q+2` vertices: after `q` complete blocks, one more `B_q` is new and valid. If `h>0`, the next `C_q` repeats `B_q`. If `h<0`, the next two algebraic vertices are `C_M=c²+1` and `D_M=-1`; the first leaves the graph and the second loses positivity. If the upper graph bound is removed, `C_M` adds one further positive vertex, but the next step still fails.

For `(3,4,5)`, `(M,h)=(3,-1)`. The full simple prefix in `Q_24` is

`9,16,20,5,11,14,22,3,13,12,24,1,15,10`.

The next vertex is 26 and the following one is -1. For `(5,12,13)`, `(M,h)=(59,23)`; the first failure is the repeated vertex 98, since `98+98=196`. These are two different exact stopping mechanisms, not a generic termination theorem for another dynamical system.

## 5. Circle meaning, quotient loss, and the descent boundary

The first edge of (3) has an inherited Thales interpretation. In a right triangle with hypotenuse `c`, the altitude divides it into segments `a²/c` and `b²/c`. Multiplying these segment lengths by `c` gives the first two square-path vertices. The normalized right-angle vertex is

\[
(X,H)=\left(\frac{a^2}{c^2},\frac{ab}{c^2}\right),
\qquad (X-1/2)^2+H^2=1/4.
\]

The complementary projections are `X` and `1-X`. This is a map from the primitive triangle to a rational point on the circle with diameter `[0,1]`, and to an edge whose two endpoints are themselves squares.

Most square-sum edges lack that last predicate. For example `11+14=25` describes real right-triangle legs `sqrt11,sqrt14`, not integer legs. More decisively, the four-reflection macro sends the seed `a²` to `a²+2`, and

\[
a^2<a^2+2<(a+1)^2\qquad(a\ge1).
\]

Thus the macro immediately leaves the integer-square seed class. It is an exact arithmetic walk on square-graph vertices, not a closed iteration on primitive triples. A map transporting a primitive-triple descent or a Collatz root certificate through this walk would need an additional state coordinate and a proof of the required guard. The finite positivity bound for these reflections does not supply such a proof.

Similarly, (3) and (10) certify local simple paths, not Hamiltonian paths. The root's fourteen-vertex prefix lies in `Q_24`, whose lack of a Hamiltonian path is independently proved in the inherited endpoint-aware note. This gives a concrete control against replacing a coverage obligation by the existence of a structured local walk.

## 6. Reproduction and scope

An independent proof audit of sections2--4 passed: positivity, residue
distinctness, the sole extra chord, the parity of h, and the exact first
repeat/escape boundary were checked from the displayed formulas.

```text
python 04-computation/experiments/zenodo_square_20260925.py
python -O 04-computation/experiments/zenodo_square_20260925.py
```

The standard-library script uses exceptions for every check. Its universe is all 1593 primitive triples with `c<=10000`, independently cross-checked against direct square-root enumeration for `c<=300`; full reflection-prefix replays for every primitive triple with `c<=100`; and 3900 exact tests of the gap-`d` macro. It checks the six vertices, all fifteen induced edges, inverse seed recovery, positivity, circle coordinates, endpoint formulas, the two rigid root identities, and all first-failure locations in the bounded replay. Path formula evaluation and direct reflection replay are independent implementations. Controls include the root, the negative-Pell extra chord, the `98` self-loop, the `d=2` negative vertex, and the first nonsquare translated seed. All infinite claims above have elementary proofs; the finite census corroborates their scopes.

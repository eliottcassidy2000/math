# Guarded reset-two debt families with a smaller common future

**Status: PROVED exact identities, guarded dependency reductions, densities,
and a restricted ansatz obstruction; FINITE-EXACT declared controls;
OPEN universal reset-two closure and Collatz convergence.**
2026-10-04. No historical priority claim.

## 1. Inheritance, object types, and the retained obligation

The closest inherited mechanism is the unbounded reset-at-least-three switch
in [checked switch, sections 2–3](checked_switch_phase19_20261004.md).
It gives a smaller proof dependency even when the actual trajectory is still
growing. The [virtual contraction ladders](virtual_contraction_ladders_20261004.md)
retain different smaller dependencies for some of the same sources; neither
choice should suppress the other.

The canonical hostile is first reset two: substituting two into the earlier
word switch creates an exponent zero. The corrected survivor is a debt state
with its child certificate retained, not a completed join. Sources 7 and 27
show that its terminal types differ. The least-used relevant sidecar here is
the difference of two ordered affine carries. A short search of that coordinate
produced a new guarded portal through the odd boundary, and an exact rational
orbit explains why the same ansatz fails at the next debt exponent.

The live board is: actual valuation words; smaller source certificates; the
two debt boundary types; sibling joins; exact dyadic cylinders; signed terminal
basins. The anchor is proof reuse. The niche is a rational orbit that certifies
failure of a whole symbolic ansatz. The wildcard is which carry identities
allow further boundary portals. No tournament or finite phase is used as a
substitute for an integer certificate.

Write \(U_s(n)=\operatorname{oddpart}(3n+s)\), \(s=\pm1\), on positive odd
integers, with the plus sheet denoted \(U\). A word lists actual positive
two-adic valuations. Its metadata \((P,Q,B)\) mean

\[
 F_{s,w}(n)=\frac{Pn+sB}{Q},\qquad P=3^{|w|},\quad Q=2^{\sum w}.
\]

Appending exponent \(a\) changes the metadata to
\((3P,2^aQ,3B+Q)\). The affine expression is an actual word only with its
integrality and exact-valuation guards.

## 2. A reset-two switch

**PROVED.** Let \(r,a\ge1\). If a positive odd source \(n>1\) has actual word

\[
 w=(1^r,2,1,6,a),
\]

then \(m=(n-1)/2\) is positive odd and smaller, has actual word

\[
 v=(1^{r-1},3,1,2,1,a+2),
\]

and both words reach the same endpoint after \(r+4\) odd steps.

For \(r=1\), the source word has metadata
\((243,2^{a+10},1279)\), and the child word has
\((243,2^{a+9},761)\). The identity \(2(761)-1279=243\) proves that their
inverse sources at a common endpoint satisfy \(m=(n-1)/2\).
The same metadata relation holds after prepending \(r-1\) exponent-one
steps to both words: \(H(n)=(n-1)/2\) commutes with the affine
exponent-one map. Equivalently, for every \(r\),

\[
 P_w=P_v,\qquad Q_w=2Q_v,\qquad 2B_v-B_w=P_w.
\]

The first valuation one forces \(n=3\bmod4\), so the inverse child is a
positive odd integer. Backward induction from the odd endpoint gives every
intermediate inverse as an integer: reduction modulo three first checks the
last inverse edge, then the remaining powers of three check earlier edges.
Each inverse numerator is odd, so the prescribed powers of two are exact.
All intermediate values are positive. This proves word legality, rather than
just equality of two unguarded affine formulas.

The same calculation on the minus sheet replaces the child by
\((n+1)/2\). It transports that child's supplied basin; it does not identify
the distinct basins containing 1, 5 and 17.

### The debt calculation behind the switch

For the inherited first-reset-two state, put

\[
 n=2^{r+1}t-1,\quad
 b=1+v_2(3^rt-1),\quad
 M=\frac{3^rt-1}{2^{b-1}}.
\]

The child reaches \(M\) by \((1^{r-1},b)\). The source after
\((1^r,2)\) is \(3\cdot2^{b-2}M+1\). A state
\(Y=3^e2^hM+1\) follows forced exponent-two steps when \(h\ge3\);
the two terminal types are

\[
 h=2:\quad U(Y)=\operatorname{oddpart}(3^{e+1}M+1),
 \qquad
 h=1:\quad U(Y)=3^{e+1}M+2.
\]

The new switch takes \(e=h=1\), hence \(b=3\), and imposes the exact
guard \(M=59\bmod128\). For \(Z=9M+2\),

\[
 W=\frac{27M+7}{64}\ \text{is odd},\qquad
 U^3(M)=\frac{27M+23}{16}=4W+1
\]

with actual child word \((1,2,1)\). Siblings \(W\) and \(4W+1\) have
the same next odd image, and the latter valuation is larger by two.
This supplies the last exponent \(a\) versus \(a+2\).

For example,

    315 --1,2,1,6,2--> 19,
    157 --3,1,2,1,4--> 19.

Here the source already falls to 25 before the join. This is an illustration
of the identity, not evidence for new growing-prefix efficacy. Source 7 has
debt \((e,h,M)=(1,2,1)\); source 27 has \((1,1,5)\), which fails the portal
guard. Neither is erased by this rule.

## 3. Exact selector, density, and a genuinely growing subfamily

For any positive exponent word with metadata \((P,Q,B)\), its exact source
cylinder is

\[
 Pn+B=Q\pmod {2Q},
 \quad\text{or}\quad
 n=(Q-B)P^{-1}\pmod {2Q}.
\]

This congruence gives the displayed endpoint odd. At each earlier step,
reduction at the required two-adic precision forces the next division,
and the final oddness rules out an extra factor of two at any intermediate
step. Thus it specifies exactly the word, not just a necessary residue.
Every positive representative is an actual positive route.

For the switch, \(\sum w=r+a+9\), so the full cylinder modulus is
\(2^{r+a+10}\). One can avoid guessing \(a\): parse the finite run
\(r=v_2(n+1)-1\), check the exact prefix \((1^r,2,1,6)\) modulo
\(2^{r+10}\), calculate its odd endpoint \(W\), and set
\(a=v_2(3W+1)\). The script's applicability function does precisely this,
without a root search. It returns a smaller dependency and a checked join,
not a home certificate.

The union over all \(a\ge1\) is just that prefix cylinder for each \(r\).
Different \(r\) are disjoint because they specify the first reset position.
Its natural density is

\[
 \sum_{r\ge1}2^{-r-10}=2^{-10}
\]

among all positive integers, and \(2^{-9}=1/512\) among odds. This is
not an unjustified countable addition of densities: the omitted runs
\(r>R\) lie in \(n=-1\bmod2^{R+2}\), whose upper density tends to zero.

**PROVED growing subfamily.** Restrict to \(a=1\), \(r\ge8\).
Every source-side prefix through the common endpoint is strictly larger
than the original source. For the initial ones, and the steps following
the first two reset letters, the affine multiplier exceeds one. The
smallest late multiplier, after the exponent six, is

\[
 \frac{3^{r+3}}{2^{r+9}}=\frac{27}{512}\left(\frac32\right)^r>1
 \quad\Longleftrightarrow\quad r\ge8.
\]

The last exponent one multiplies it by \(3/2\), and every carry is positive.
Hence all displayed source values exceed \(n\), uniformly over every positive
lift of the exact cylinder. The claim fails at \(r=7\): the least source
91903 falls to 82807 before the join.

At \(r=8\), the whole cylinder and join are

\[
 n=236031+524288t,\quad m=118015+262144t,\quad
 z=478505+1062882t,\qquad t\ge0.
\]

The source word has eight ones followed by \(2,1,6,1\); the child word
has seven ones followed by \(3,1,2,1,3\). Every source-side step exceeds
its original source. This is new relative to the initial-reset-at-least-three
schema, since the first reset is two. No comparison against every possible
older common-future rule is claimed.

The disjoint growing cylinders have density
\(\sum_{r\ge8}2^{-r-11}=2^{-18}\), or \(2^{-17}\) among odds.
The same initial-run tail bound proves existence of this density.

### Least-counterexample consequence, not an unconditional pending count

If there is a positive odd source that never reaches 1, let \(n\) be the
least such source. It is not 1. A first valuation at least two gives
\(U(n)\le(3n+1)/4<n\), a contradiction. Otherwise its initial one-run is
finite: an infinite run would require every power of two to divide \(n+1\).
An eventual reset at least three is handled by the inherited smaller-child
join. Its child is home by minimality; the small child-root case \(n=3\)
is directly \(3\to5\to1\). Therefore the least counterexample must have
first reset two.

For a run \(r\), write \(n=2^{r+1}t-1\) with \(t\) odd. Reset two is
equivalent to \(t=(-1)^r\bmod4\). These disjoint cylinders have modulus
\(2^{r+3}\) and total density \(1/8\), or \(1/4\) among odds. Again the
initial-run tail controls the infinite union. The new all-terminal portal
family is inside that kernel and cannot contain a least counterexample,
because its smaller child would be home. Excluding just this family leaves
a necessary candidate domain of density \(127/1024\), or \(127/512\)
among odds. This is neither a nonconvergence density nor a count of inputs
that the selector can certify without their child premises. The two further
portals below are deliberately not included in this density figure.

## 4. General portal lemma and the first failed ansatz

Let \(e\ge1\), and let a positive word \(u\) have length \(e+2\), total \(S\),
and carry \(B_u=2^S+7\). Then the same argument gives the guarded switch

\[
 (1^r,2^e,1,S+2,a)
 \quad\longleftrightarrow\quad
 (1^{r-1},2e+1,u,a+2),
 \qquad m=(n-1)/2.
\]

Indeed the inherited debt reaches \(Z=3^{e+1}M+2\); the next prescribed
exponent \(S+2\) gives \(W=(3^{e+2}M+7)/2^{S+2}\), while \(u\) sends
\(M\) to \(4W+1\). The same positive-odd inverse guard proves the alternative
word legal. Three exact examples are

| e | u | S | B |
|---|---|---|---|
| 1 | (1,2,1) | 4 | 23 |
| 3 | (1,6,1,3,1) | 12 | 4103 |
| 4 | (2,2,1,3,3,1) | 12 | 4103 |

This is a lemma for every word satisfying its carry condition, not a theorem
that such a word exists for every \(e\).

There is a stronger, exact stopping reason at \(e=2\). Consider the specific
two-step ansatz: send \(Z=3^{e+1}M+2\) by \((a,b)\), and \(M\) by
\((u,b+\delta)\), with \(|u|=L=e+2\), total \(S\), and equal formal endpoints
for all \(M\). Matching denominators and carries requires

\[
 \delta=a-S,\qquad 3B_u+2^S=21+2^a.
\]

Since \(L\ge3\), the smallest possible carry is \(3^L-2^L\ge19\);
therefore \(a>S\). Reduction modulo three forces \(a-S=2s\), \(s\ge1\).
Equivalently,

\[
 B_u=c_s2^S+7,\qquad c_s=(4^s-1)/3.
\]

Thus the odd rational start \(-7/3^L\) would have to follow the exact
accelerated word \(u\) to \(c_s\). Oddness of the endpoint and backward
induction make every intermediate an odd rational with odd denominator,
so its valuation word is unique. This is an arithmetic obstruction, not
a finite exponent search.

For \(e=2\), \(L=4\), that unique rational orbit is

\[
 -7/81\xrightarrow{2}5/27\xrightarrow{1}7/9
 \xrightarrow{1}5/3\xrightarrow{1}3.
\]

Its endpoint 3 is strictly between \(c_1=1\) and \(c_2=5\). No positive
\(s\) works. The two-step ansatz therefore fails for \(e=2\) at every
exponent height. It does not rule out longer joins, another smaller
dependency, or a different boundary representation. For lengths 3, 5 and 6,
the same rational test gives the three words in the table and endpoint 1.

## 5. Certificate transport and information retained

The source object is a supplied positive odd integer together with an exact
word guard. The target is a smaller integer and a common endpoint. The map
retains both ordered words and uses \(m=(n-1)/2\). It preserves the statement
that source and child have the same eventual plus-sheet root, or the same
supplied signed basin. It does not preserve the actual trajectory, source
label, ordinary route length, or an ordering of all candidate dependencies.
Required sidecars are the sign, source identity, guards, common endpoint,
and the child's first-hit certificate.

For the three listed portals, the source cannot hit 1 before the last edge.
A hit before the displayed large exponent would force subsequent valuations
to be two, contradicting that exponent. A hit immediately after it would
require \(3^{e+2}M=2^{S+2}-7\), but the right side has three-adic valuation
one for \(S=4\) and \(S=12\), smaller than \(e+2\). The child cannot hit 1
early because its last exponent is \(a+2\ge3\), whereas the root self-edge
has exponent two. Thus a supplied first-hit child certificate contains its
whole displayed child prefix. Splice its suffix onto the checked source word.
The code verifies the supplied certificate's actual source before splicing.

The words have equal odd lengths and the source total valuation exceeds the
child's by one. Consequently transport keeps the full odd certificate rank
and increases the full ordinary rank by exactly one. A dependency reduction
is therefore a decrease of source label, not a decrease of route length.
A word identity alone is not a completed root certificate.

The minus cycles \(1\), \(5\to7\to5\), and
\(17\to25\to37\to55\to41\to61\to91\to17\) remain separate.
Their valuations never include six, so none triggers the first portal.

## 6. Reproduction and exact finite scope

Companion program:
[reset_two_debt_families_20261004.py](../../04-computation/experiments/reset_two_debt_families_20261004.py).
Saved stdout:
[reset_two_debt_families_20261004.out](reset_two_debt_families_20261004.out).

    python 04-computation/experiments/reset_two_debt_families_20261004.py
    python -O 04-computation/experiments/reset_two_debt_families_20261004.py

Both modes use explicit exception checks. The declared universe is:

- 4608 exact rewrites: three portals, \(r=1,\ldots,24\), \(a=1,\ldots,8\),
  both signs, and the first four positive lifts of each cylinder.
- 96 direct debt-boundary-to-sibling checks.
- 363 growing-prefix examples for \(r=8,\ldots,128\), three lifts each,
  plus an exact multiplier test for every prefix in those cases.
- All 10000 odd sources through 20000 for independent applicability checks:
  20 match the first portal. This is an applicability census, not a
  convergence census.
- 26 supplied-child AST transports in the growing family, with bounded
  child-route searches explicitly supplying the premises and canonical
  source comparisons performed only after construction.
- A separate first-hit endpoint-one control:
  \(n=(2^{168}-1279)/243\) has word \((1,2,1,6,158)\); its child has word
  \((3,1,2,1,160)\). A proposed smaller example 68579 fails the source
  cylinder: an integral candidate hub alone did not supply the root-to-hub
  ternary guard. It is retained as a hostile, not used as a certificate.
- Exact rational ansatz controls, the \(r=7\) and 315 efficacy hostiles,
  the retained debts at 7 and 27, signed-basin controls, and six malformed
  source/type/certificate rejections.

The general results above use the proofs, not extrapolation from these
finite ranges. Universal debt closure remains open.

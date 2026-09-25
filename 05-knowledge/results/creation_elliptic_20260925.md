# A coordinate-controlled creation certificate on the fruit elliptic curve

**Status: PROVED** on the explicitly generated subgroup below; **FINITE-EXACT** for the stated controls. The extension from that subgroup to every rational curve point is **CITED**, not newly proved by saturation. This construction proves a terminating elliptic membership decoder, not Collatz convergence. Date: 2026-09-25.

## 1. Inheritance, corrected boundaries, and the actual seed

Use

\[
E:y^2=x^3+109x^2+224x,\quad G=(-4,28),\quad T=(56,728),
\quad Q=2T=(4,52),\quad U=3T=(0,0).
\]

Work first on the subgroup `L=<G,T>`, where `G` has infinite order and `T` has order six. Thus every point of `L` has a unique expression `mG+kT`, with integer `m` and `k mod 6`. No full-group generator claim is needed for the proofs on `L`.

The [original fruit audit](catalan_elliptic_20260921_elliptic.md) proves the projective coordinate map, six torsion points, and the positivity test. The [existing subgroup-tree note](prime_shells_20260921_elliptic.md) already constructs a coefficient ternary tree and proves that the literal inverse-tripling tree from `9G` has seven subgroup vertices. It also records the failed preservation of fruit positivity. The [rank and positive-multiples audit](collatz_mod6_20260921_fruit_rank_positive_multiples.md) proves rank one by a PARI descent, verifies saturation at primes below `10^4`, and corrects the BSD index formula to `|Sha|/[E(Q):L]^2`. Those finite saturation and analytic data do not themselves prove index one.

The primary [Bremner–Macleod paper](https://ami.uni-eszterhazy.hu/uploads/papers/finalpdf/AMI_43_from29to41.pdf), printed pp. 30–32, gives the coordinate map, torsion, and in Remark 2.2 identifies `G` as a generator modulo torsion and `9G` as the first positive multiple. Its generator statement is the **CITED** route to extending the decoder to all of `E(Q)`. The decoder below is proved without relying on that extension.

The exact primitive triple of `9G` has coordinate lengths `(81,80,79)` and is the repaired triple in the inherited audit. The script verifies its three full integers and its fruit sum exactly. This does not assert that any differently transcribed or permuted new input has already passed the same test.

Closest mechanism: the real-component parity and rational division guards. Canonical hostile: `9G` is positive, while the known ternary children and the halved creation descendants need not be. Corrected near miss: a point-generation tree is not a tree of positive triples. Least-used sidecar: the rational square class of `x`, which retains both free parity and torsion parity without first computing the free coefficient.

The live board is: square classes; rational halving; torsion phases; integer ranks; positivity; actual Collatz arrows.

## 2. The visible square class carries two parity bits

**E1 (PROVED).** Define the square-class map

\[
\alpha:E(\mathbb Q)\longrightarrow\mathbb Q^*/\mathbb Q^{*2},
\qquad
\alpha(O)=1,\quad\alpha(U)=224,\quad\alpha(x,y)=x\quad(x\ne0).
\tag{1}
\]

It is a group homomorphism. For an ordinary chord `y=lambda*x+nu`, the three intersection abscissae have product `nu²`, by substituting the line into the cubic. The third intersection is the negative of the sum and has the same abscissa as the sum. This proves the square-class law, with tangencies counted with multiplicity. A vertical chord gives a square product. For a chord through `U`, the translation formula `x(P+U)=224/x(P)` gives exactly the declared exceptional value; `U+U=O` and the identity cases agree.

Now `alpha(G)=-1` and `alpha(T)=14`. The four classes `1,-1,14,-14` are distinct over the rationals. Consequently

\[
\alpha(mG+kT)=(-1)^m14^k\pmod{\mathbb Q^{*2}}.
\tag{2}
\]

Within `L`, this gives the exact equivalence

\[
P\in2L\quad\Longleftrightarrow\quad\alpha(P)=1.
\tag{3}
\]

Indeed (2) detects that both `m` and `k` are even, precisely the condition for solving `2n=m, 2j=k mod 6`.

The two bits are computable by testing which of `x, -x, x/14, -x/14` is a rational square, with the exceptions in (1). No factorization of the enormous numerator is necessary: a rational-square test uses integer square roots of numerator and denominator. For a point in `L`, exactly one test succeeds. The sign of `x` detects free parity, not the sign of the free coefficient; negating a point leaves `x` unchanged.

## 3. Halves can be recovered by quadratic square-root tests

This makes the construction an actual coordinate algorithm, beyond relabelling an integer tree.

Let `R=(X,Y)` be a nonzero point with a rational half and `X!=0`. The duplication formula is

\[
X=\frac{(u^2-224)^2}{4v^2},\qquad (u,v)\in E.
\tag{4}
\]

Thus `X` is a rational square. Choose its positive square root `w`. Eliminating `v` and setting `z=u+224/u` gives

\[
z=2X\pm2Y/w,\qquad
u=\frac{z\pm\sqrt{z^2-896}}2.
\tag{5}
\]

For each rational candidate `u`, test both rational square roots of `u³+109u²+224u` and verify doubling exactly. These finitely many tests recover every rational half. The case `R=O` has halves `O,U`. The case `R=U` has none: (4) would require `u²=224`. No irrational branch is silently accepted.

If `R` belongs to `2L`, its two rational halves lie in `L` and differ by `U`. There cannot be others because the only rational two-torsion is `O,U`: the remaining two roots would require `10985=65*13²` to be a rational square. Among these halves, exactly one has even torsion coordinate. By (2), it is the unique half with `alpha` in `{1,-1}`. This gives an intrinsic, checkable kernel-phase selection.

## 4. A complete terminating decoder with a six-point terminal set

**E2 (PROVED on `L`).** Given `P`, read `r,t in {0,1}` from

\[
\alpha(P)=(-1)^r14^t.
\]

Subtract the indicated offset and select the even-torsion half:

\[
P'=\operatorname{selectedHalf}(P-rG-tT).
\tag{6}
\]

The input coefficient `m` is not provided to the algorithm. Nevertheless the proof gives its exact action:

\[
m'=\lfloor m/2\rfloor,\qquad
k'=\text{the unique }j\in\{0,2,4\}\text{ with }2j+t=k\pmod6.
\tag{7}
\]

Retain the two output bits and the point `P'`. The certificate equation is

\[
P=2P'+rG+tT.
\tag{8}
\]

The terminal set is the six explicit points

\[
\mathcal B=\{jQ,-G+jQ: j=0,1,2\}.
\tag{9}
\]

Define nonnegative integer ranks

\[
\rho(m)=\begin{cases}m&m\ge0,\\-m-1&m<0,\end{cases}
\qquad D(m,k)=2\rho(m)+(k\bmod2).
\]

Equation (7) gives `rho(m')=floor(rho(m)/2)` and `k' mod 2=0`. Therefore `D` strictly decreases at every step outside `B`. Its zero set is exactly `B`. This proves termination on **all** of `L`, including negative coefficients and torsion inputs. If `rho>0`, the required number of steps is its binary length; if `rho=0` but the torsion parity is odd, one normalization step suffices.

A certificate is a terminal point and the ordered bit pairs in (8). Reading it backwards creates the initial point by exact additions and doublings. A finite accepted certificate proves membership in `L` without a global generator theorem; every point in `L` receives one from the decoder. The intended target is subgroup membership and reconstruction, not arrival at the single identity point.

For the large supplied seed, the certificate is particularly small:

```text
9G -> 4G -> 2G -> G -> O,
digits (r,t): (1,0),(0,0),(0,0),(1,0).
```

All intermediate point coordinates can be checked independently by (8). The coordinate signs of their fruit triples are not preserved.

## 5. Positivity and the actual Collatz arrow remain separate

**E3 (PROVED).** No positive fruit solution is the double of a real elliptic point. The inherited positive-fruit criterion requires

\[
x<-14/3,\qquad x^2+112x+784>0,
\]

while (4) makes every finite real double have nonnegative `x`. Points of order two double to `O`, which is also outside the admissible positive-fruit locus. Thus the offset subtraction in (6) is essential; raw halving of the positive point `9G` is impossible even over the reals.

The step (6) sends `9G` to `4G`, whose fruit triple is not positive; the following `2G` and `G` also fail positivity. The terminal point `-G` is fixed by the decoder. The terminal set and coefficient-side convention cannot be omitted in a claim of descent to `O`.

There is an exact guarded elliptic encoding of the **actual shortcut Collatz** maps on the pure subgroup `{mG}`:

\[
m\text{ even}: P\mapsto\operatorname{selectedHalf}(P),\qquad
m\text{ odd}: P\mapsto\operatorname{selectedHalf}(3P+\sigma G).
\tag{10}
\]

Its output is `T_sigma(m)G`, for `sigma=±1`. The same square-class and kernel-phase tests make this lift explicit. But it is a different map from (6): the plus step at `9G` gives `14G`, not `4G`.

The inherited canonical-height identity also exposes the remaining obstruction. On pure multiples, the odd plus step multiplies height by `((3m+1)/(2m))²`, exceeding one for every positive `m`. The odd minus step also raises height for `m>1`. Thus the useful binary creation rank does not become a decreasing Collatz rank simply by passing to the curve.

## 6. Division quotients: four, nine, thirty-six, and sixteen

The subgroup presentation gives, for each positive integer `b`,

\[
L/bL\cong\mathbb Z/b\mathbb Z\;\times\;
\mathbb Z/\gcd(b,6)\mathbb Z.
\tag{11}
\]

Therefore the quotient has 4 classes for division by 2, 9 for division by 3, 36 for division by 6, and **16**, not 8, for division by 8. The second factor is the torsion residue. Three successive binary refinements have `|L/8L|=16` classes: four choices at the first stage, then two at each later stage after fixing the even-torsion section. Treating each stage as an independent four-way choice would incorrectly count 64.

These counts concern quotient classes, not counts of rational preimages of a point. A divisible point has two rational halves or three rational third-parts, according to the rational kernel. The inherited seven-vertex inverse-tripling tree at `9G` remains intact. An actual 36 occurs in the mixed division quotient; it does not identify this object with the earlier 36-cycle, graph edges, or representative pairs.

## 7. Reproduction and remaining obligation

Run

```text
python 04-computation/experiments/creation_elliptic_20260925.py
python -O 04-computation/experiments/creation_elliptic_20260925.py
```

The [script](../../04-computation/experiments/creation_elliptic_20260925.py) uses exact fractions and explicit checks. The [output](creation_elliptic_20260925.out) records: all 294 points `mG+kT`, `-24<=m<=24`, `0<=k<6`; each intrinsic parity read, half selection, rank step, and complete reconstructed certificate; 2,916 independent square-class homomorphism pair checks; 72 duplication-positivity controls; all 48 lifted plus/minus shortcut steps at `m=1..24`; the exact large triple and its four-step certificate; and the finite division quotient counts. The algorithm receives point coordinates alone, while separately generated coefficients are used only by the hostile controls and independent verification.

The positive result is a short, coordinate-verifiable creation certificate with a proved integer descent rank. The remaining Collatz obligation is to replace (6) by its actual arrow (10) while obtaining a decreasing rank or another valid termination mechanism. The curve equation and its division quotients do not supply that missing step by themselves.

An independent read-only proof audit of sections 2–6 re-derived the halving quadratics, checked the `O/U` exceptions, unique kernel-phase selection, strict integer rank and terminal-set convention, and all four quotient counts: **PASS**. Normal and optimized script outputs are identical, and the note's local links resolve.

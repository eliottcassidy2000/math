# Exact analytic terminals for the unchanged curves E17 and E627

2026-10-13. **FINITE-EXACT integer receipts + CITED arithmetic theorems.**
The final full-BSD implication uses the supplied accepted low-corank theorem;
it is labelled separately from the classical conclusions.

## 1. The first Selmer obstruction is bypassed, not renamed

The earlier [finite Selmer-growth packet](bsd_selmer_growth_20261013a.md)
left both

\[
E_n: y^2=x^3-n^2x,\qquad n=17,627=3\cdot11\cdot19,
\]

with pure 2-Selmer dimension two. That one-level statistic does not decide
the full 2-infinity Selmer corank. The present calculation keeps exactly
these curves and supplies a different arithmetic certificate: a nonzero
central L-value. It does not compute a full 4-Selmer group or order.

For a positive odd squarefree integer n, set

\[
A_n=\#\{(x,y,z)\in\mathbb Z^3:2x^2+y^2+32z^2=n\},\quad
B_n=\#\{(x,y,z)\in\mathbb Z^3:2x^2+y^2+8z^2=n\}.
\]

Tunnell identifies \(a(n)=A_n-B_n/2\) and proves

\[
L(E_n,1)=\frac{a(n)^2\beta}{4\sqrt n},\qquad
\beta=\int_1^\infty\frac{dx}{\sqrt{x^3-x}}>0.                 \tag{1}
\]

**CITED:** J. B. Tunnell, *A Classical Diophantine Problem and Modular
Forms of Weight 3/2*, Inventiones 72 (1983), 323–334, printed pp. 328–329,
Theorem 3 and the preceding coefficient identity.
[Original scanned paper](https://sites.math.rutgers.edu/~zeilberg/EM22/JT1983.pdf),
[publisher record](https://doi.org/10.1007/BF01389327).
This is the nonvanishing direction; no conjectural converse is used.

The exact counts give:

| Source n | A_n | B_n | a(n) | Exact positive L-value |
|---:|---:|---:|---:|---|
| 17 | 4 | 16 | -4 | \(4\beta/\sqrt{17}\) |
| 627 | 40 | 48 | 16 | \(64\beta/\sqrt{627}\) |

Thus both analytic ranks are zero. These are exact formulas, not decimal
evaluations with an assumed numerical error bound.

## 2. The arithmetic implication and its precise dependency

Both curves have complex multiplication: over \(\mathbb Q(i)\),
\((x,y)\mapsto(-x,iy)\) is an endomorphism squaring to negation.
The cited Coates–Wiles and Rubin nonvanishing theorems give rational rank
zero and finite Tate–Shafarevich group for these CM curves when (1) is
nonzero. See Karl Rubin, *The work of Kolyvagin on the arithmetic of
elliptic curves*, printed p. 129, final paragraph before acknowledgments,
which explicitly records this CM implication and its attribution:
[author's proof exposition](https://www.wstein.org/edu/Fall2003/252/references/rubin-kolyvagin/Rubin-The_work_of_Kolyvagin_on_the_arithmetic_of_elliptic_curves.pdf).
CM alone is not the hypothesis: the exact nonvanishing certificate is
essential. No present-day general theorem is being inferred from the
historical analytic conjecture also discussed on that page.

For every prime q, Kummer exactness now gives

\[
0\longrightarrow E_n(\mathbb Q)\otimes\mathbb Q_q/\mathbb Z_q
\longrightarrow\operatorname{Sel}_{q^\infty}(E_n/\mathbb Q)
\longrightarrow\Sha(E_n/\mathbb Q)[q^\infty]\longrightarrow0.
\]

The left term is zero because the rational group is finite, and the right
term is finite. Hence the **full q-infinity Selmer corank is zero**.
This is stronger than the formerly inconclusive first-level upper bound.

**Conditional application of a supplied accepted theorem.** Theorem 1.1
of *Exact Birch–Swinnerton-Dyer Formula from Low Selmer Corank*, already
typed in the [earlier packet](bsd_selmer_growth_20261013a.md), applies to
each unchanged curve. Under that accepted input, the full BSD formula
holds for both. The present computation neither re-proves that theorem
nor evaluates the regulator, all Tamagawa factors, or the order of Sha.
In particular, no numerical or inferred rank bound is called a Sel4
receipt. The two-dimensional pure 2-Selmer group is consistent with
nontrivial finite Sha[2], rather than a divisible corank-two tail.

More precisely, finite Kummer exactness and rank zero now identify the
pure 2-Selmer quotient with Sha[2]. Thus for each of E17 and E627,
Sha[2] is exactly `(Z/2)^2`, of order four. In the eight dimension-four
new terminals below it is exactly `(Z/2)^4`, of order sixteen. This
determines the 2-torsion subgroup, not the full finite group or its order.

## 3. Exact arithmetic engines and controls

The [checker](../../04-computation/experiments/bsd_tunnell_terminals_20261013c.py)
implements two independent counting paths:

1. enumerate all signed x,z within the positive-definite bounds and test
   the remaining y-square with integer square root;
2. convolve the three finite theta polynomials, using the multiplicity
   one for zero and two for each positive square.

For each form, every solution obeys
\(|x|\le\lfloor\sqrt{n/2}\rfloor\) and
\(|z|\le\lfloor\sqrt{n/c}\rfloor\), so the first search is exhaustive.
Nonnegative exponents make degree truncation in the second path exact.
They agree on both forms for every integer 1 through 627.

A third calculation uses the independent product expression

\[
q\prod_{m\ge1}(1-q^{8m})(1-q^{16m})
\sum_{z\in\mathbb Z}q^{2z^2}.
\]

Its coefficient agrees with \(A_n-B_n/2\) at all 314 odd indices through
627. The product and its identification with a(n) are in Tunnell's
introduction and printed pp. 327–328. Only finitely many factors can
contribute at the declared degree.

The zero control n=5 has a(5)=0 and the exact point
\((25/4,75/8)\in E_5(\mathbb Q)\). `terminal(5)` returns `None`, meaning
this nonvanishing gate supplies no conclusion. It is not evidence of
corank zero. `audit` reauthenticates the integer source, all counts and
the coefficient; bool/float aliases and altered receipts are rejected.
The old Monsky API independently confirms the same two pure dimensions.

No PARI executable is available on this host: PowerShell's `gp` is an
unrelated alias. No software was installed and no point-search or
numerical L-function oracle was used. Source scans were inspected in
scratch outside the tracked repository.

## 4. Exact gain on the previously declared finite universe

The extension retains **all 4,056 positive odd squarefree n <= 10,000**
from the previous packet. It adds no residue or point-search filter.
The prior baseline consists of its 2,886 cases with pure 2-Selmer
dimension at most one. Apply the independent Tunnell gate to the 1,170
remaining sources, without replacing or twisting any marked curve.

| Pure 2-Selmer dimension | Previously unresolved | New nonzero terminals | Still unresolved by this union |
|---:|---:|---:|---:|
| 2 | 988 | 716 | 272 |
| 3 | 164 | 0 | 164 |
| 4 | 18 | 8 | 10 |
| Total | 1,170 | **724** | **446** |

Thus **3,610 of 4,056** sources meet a low-corank hypothesis through
these two gates. The 724 new sources have the stronger rank-zero and
finite-Sha conclusion proved from exact nonvanishing in section 2;
the original baseline may also include corank-one curves. Under the
supplied accepted theorem, all 3,610 satisfy its full BSD conclusion.
No claim is made that the 446 residual curves fail BSD or have corank
at least two. Their least source is 41. A zero central value leaves
possible higher leading terms and arithmetic information to be resolved.

`census()` returns every triple `(n, pure_dimension, a(n))` in increasing
source order. The exact residual list is in the saved output; the full
row table has SHA256
`791511bd55a2dd936c544e840df278a3816dab3235c3038e547626f96ad45ff2`
under ASCII JSON with separators `(',', ':')`. This records unchanged
source identities, not just histogram totals.

The bulk theta coefficients agree with the separate eta-product path
at **all 5,000 odd indices through 10,000**, including nonsquarefree
indices as coefficient controls. The arithmetic gate still requires
squarefreeness. A declared 59-source subset is additionally rechecked
by signed enumeration: the first and last twelve gains, every 53rd
gain, every 37th residual, all eight dimension-four gains, and
5, 17, 627, with duplicates removed. This does not describe a full
second signed enumeration of the 10,000 universe.

## 5. Reproduction

Reproduce in normal and optimized modes:

```text
python -B 04-computation/experiments/bsd_tunnell_terminals_20261013c.py
python -B -O 04-computation/experiments/bsd_tunnell_terminals_20261013c.py
```

The [saved output](bsd_tunnell_terminals_20261013c.out) records the exact
finite universe and retains the cited/conditional distinction.
Normal and optimized runs agree at **1,755 checks**, LF SHA256
`ea2f95ad73f6d21eab8582c5e17e1836d737333eb7ee74862410b092d8c1b681`.

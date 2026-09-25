# The actual 9 to 14 to 7 passage on a primitive-triangle spine

**Status: PROVED** elementary constructions and guarded descent; **INHERITED PROVED** spine and Gaussian-square mechanisms; **FINITE-EXACT** stated controls. This lane does not audit the separately supplied DOI construction, prove Collatz convergence, or equate the entire Collatz graph with the Berggren tree. Date: 2026-09-25.

## 1. Recover the actual state map before changing coordinates

The [primitive-triple filtration](ternary_triples_20260925.md) already identifies the integer boundary

\[
\Phi(n)=\left(n,\frac{n^2-1}{2},\frac{n^2+1}{2}\right),\qquad n>0\text{ odd},
\tag{1}
\]

including the degenerate terminal `(1,0,1)`. Its numerator/denominator and orientation information must be retained. The [elliptic creation decoder](creation_elliptic_20260925.md) distinguishes the creation step `9G -> 4G` from the actual shortcut step `9G -> 14G`.

The nearest geometric mechanism is already in [THM-3334 — Berggren parabolic spine and Gaussian collision torsor](../../01-canon/theorems/THM-3334-berggren-parabolic-spine-gaussian-collision-torsor.md), equations (9), (18), and (23): the spine is the `C-B=1` section and extends to the degenerate root. [THM-3341 — U-spine square-hypotenuse transplant](../../01-canon/theorems/THM-3341-u-spine-square-hypotenuse-transplant-and-triangular-plane-torsors.md) supplies the Gaussian-square relation producing `(7,24,25)` from `(3,4,5)`.

Canonical hostile: the genuine plus step `3->5` increases the spine depth. Corrected near miss: treating an even shortcut state as an odd spinor without changing its normalization. Least-used sidecar: the marked leg and the gap `C-B`, which remember which integer and clock were encoded.

Live board: primitive normalization; Berggren state spine; exact valuation guard; creation digits; actual elliptic group edges; local square identities. The new content below is their explicit compatibility for an actual Collatz descent cylinder, with an infinite family of completed root certificates.

## 2. Even 14 has a valid primitive triangle if its parity marker survives

**Z1 (PROVED).** For every positive integer `n`, put

\[
d(n)=\begin{cases}2&n\text{ odd},\\1&n\text{ even},\end{cases}
\qquad
\Psi(n)=\frac{1}{d(n)}(2n,n^2-1,n^2+1).
\tag{2}
\]

This is an integral primitive Pythagorean triple with marked first and second legs. At `n=1` it is the permitted degenerate terminal. The common divisor of the displayed raw entries is exactly two for odd `n` and one for even `n`: any odd common prime would divide both `n` and `n²-1`, which is impossible, and the power of two is read directly. Its inverse and parity marker are

\[
n=\frac{A}{C-B},\qquad C-B=2/d(n)\in\{1,2\}.
\tag{3}
\]

Unlike the odd-only chart, the marked first leg may be even. Sorting or swapping the legs without retaining the mark would change the inverse.

Consequently the **actual two-shortcut route** is

\[
9\longrightarrow14\longrightarrow7,
\qquad
(9,40,41)\longrightarrow(28,195,197)\longrightarrow(7,24,25).
\tag{4}
\]

All three triples are primitive. Passing to odd first-return time discards the middle state only after recording the exact two-step clock: `3*9+1=28=4*7`.

The quotient that simply replaces every integer by its odd part cannot carry the one-step dynamics: `14` and `7` have the same odd part, but their shortcut successors have odd parts `7` and `11`. Retain the parity/valuation phase, or explicitly use the induced first-return map on odd states.

Primitive hypotenuse alone is not a globally ordered replacement for the integer: `Psi(4)` has hypotenuse 17 while `Psi(5)` has hypotenuse 13. On the odd section (1), however, the hypotenuse is strictly increasing in `n`, so a genuine odd-return drop has a genuine geometric drop.

## 3. The Collatz return is an actual guarded jump on the inherited spine

Write the inherited Berggren spine matrix as

\[
V=\begin{pmatrix}1&-2&2\\2&-1&2\\2&-2&3\end{pmatrix}.
\]

It satisfies

\[
V\Phi(n)=\Phi(n+2),\qquad
\Phi(n)=V^{(n-1)/2}(1,0,1).
\tag{5}
\]

Thus `j=(n-1)/2` is its depth from the degenerate root; the usual positive Berggren root `(3,4,5)` has `j=1`. This parabolic spine is inherited, not a new parameterization claim.

Let

\[
U_+(n)=\frac{3n+1}{2^k},\qquad k=v_2(3n+1),
\]

and put `j'=(U_+(n)-1)/2`. Its state-triangle edge is exactly

\[
\Phi(n)\longrightarrow V^{j'-j}\Phi(n).
\tag{6}
\]

The exponent is a guarded, state-dependent integer; (6) is not a fixed Berggren branch word implementing every Collatz step. It also concerns **state triangles** `Phi(n)`, not the consecutive-pair triangles whose Berggren incomparability was proved in [ternary_berggren_20260925.md](ternary_berggren_20260925.md).

**Z2 (PROVED guarded contraction).** For `n=1 mod 4`, the halving count is at least two, and

\[
j'\le\frac34j.
\tag{7}
\]

This follows from `U_+(n)<=(3n+1)/4` and `n=2j+1`. It is a strict decrease for `n>1`. If instead `n=3 mod 4`, the halving count is exactly one and

\[
j'=(3j+1)/2>j.
\tag{8}
\]

Thus the same geometric coordinate exposes both descent and the remaining ascent obstruction.

More sharply, on the exact cylinder

\[
n=8a+1,\qquad a\ge1,
\]

one has

\[
3n+1=4(6a+1),\qquad U_+(n)=6a+1,
\quad j=4a\longrightarrow j'=3a.
\tag{9}
\]

Its ordinary shortcut route uses exactly two edges. Its geometric hypotenuse drop is

\[
C(n)-C(U_+(n))=2a(7a+1)>0.
\tag{10}
\]

At `a=1`, the edge is literally one Berggren-parent step: `V^(-1)(9,40,41)=(7,24,25)`. For general `a`, it is a jump back by exactly `a` spine edges. At `a=0`, the terminal `n=1` is fixed in odd-return time and the drop is zero.

## 4. The even shortcut appears in a genuine middle Berggren child

This is a second positive compatibility, retaining a different marked coordinate. Write `n=2p-1`, so `Phi(n)` has Euclid parameters `(p,p-1)`. The middle Berggren child acts by

\[
(p,q)\longmapsto(2p+q,p).
\]

On the spine its new parameters are

\[
(3p-1,p)=\left(\frac{3n+1}{2},\frac{n+1}{2}\right).
\tag{11}
\]

The first parameter is exactly the plus shortcut value. The determinant-one-in-magnitude Euclid transformation preserves coprimality and opposite parity. Both parameters are retained, so the child remains a primitive triangle and the map is reversible on its marked image, characterized by `m'=3r'-1`.

For the displayed example the source parameters `(5,4)` become `(14,5)`, giving the actual middle child

\[
(9,40,41)\longrightarrow(171,140,221).
\]

Its parameter 14 encodes the shortcut value, but the child itself is not `Psi(14)`. Extracting that marked parameter and applying (2) is an explicit change of representation. Forgetting the second parameter would lose the primitive-triangle reconstruction.

## 5. The creation parity certificate can verify the actual two-step macro

Use the fruit curve point `G` and the intrinsic creation decoder from [creation_elliptic_20260925.md](creation_elliptic_20260925.md). For `a>=1`, let `Q=aG`. Then the same guarded family has

\[
P=(8a+1)G=8Q+G.
\tag{12}
\]

Three creation-decoder steps read low free bits `(1,0,0)`, with torsion bits zero, and recover `Q`. This is a coordinate-checkable certificate of the cylinder, not a claim that those three creation arrows are the actual Collatz arrows.

The actual two-shortcut endpoint is instead

\[
R=\frac{3P+G}{4}=6Q+G=(6a+1)G,
\tag{13}
\]

where the fraction means the twice-selected rational half with the retained torsion phase. Thus the creation certificate supplies a guard and a reusable tail for a valid Collatz macro edge. The canonical height drop is exact:

\[
\widehat h(P)-\widehat h(R)=4a(7a+1)\widehat h(G)>0.
\tag{14}
\]

This repairs the local difficulty at `9G`: `9G->14G` grows, while the certified first return `9G->14G->7G` descends. It does not supply such a decreasing macro for every initial integer.

## 6. Why the 3,4,7,14,11,25 configuration is exact and special

The equality `7=3+4` gives Euclid parameters `(4,3)`, hence `(7,24,25)`. Equivalently this is the inherited Gaussian square of `(3,4,5)`:

\[
(|3^2-4^2|,2\cdot3\cdot4,5^2)=(7,24,25).
\]

The numbers 14 and 11 are respectively twice the state 7 and its next plus shortcut. The proposed sum has a precise uniqueness statement:

\[
2y+T_+(y)=C(\Phi(y))\quad(y>0\text{ odd})
\Longleftrightarrow\quad
\frac{7y+1}{2}=\frac{y^2+1}{2}
\Longleftrightarrow\quad y=7.
\tag{15}
\]

So `14+11=25` is a real local junction, but it cannot be imposed as a universal recurrence.

It nevertheless sits on an **infinite family of complete root certificates**. For every `r>=1`, let

\[
n_r=\frac{7\cdot4^r-1}{3}.
\tag{16}
\]

These are positive odd integers, beginning `9,37,149,597,...`, and `3n_r+1=7*2^(2r)`. Every one reaches 7 after exactly `2r` shortcut steps, with the final junction `14->7`. The finite accelerated tail

```text
7 -> 11 -> 17 -> 13 -> 5 -> 1
```

is checked directly. Appending it proves convergence for **all** of (16), and transports to their state triangles and pure elliptic multiples. This uses the inherited inverse-fibre mechanism with its exact valuation; the new geometric reading is that every such certified route meets the same marked square-hypotenuse junction. It makes no claim that all odd integers belong to this family.

## 7. Lorentz similarity and the limit of changing a linear height

For a fixed sign `sigma` and exponent `k`, the spinor operation `(s,t)->((3s+sigma*t)/2^k,t)` induces a rational linear operator on `(A,B,C)`. On the integer boundary it is

\[
M_{k,\sigma}=\frac{1}{2\,4^k}
\begin{pmatrix}
3\,2^{k+1}&-\sigma2^{k+1}&\sigma2^{k+1}\\
6\sigma&8+4^k&10-4^k\\
6\sigma&8-4^k&10+4^k
\end{pmatrix}.
\tag{17}
\]

With `J=diag(-1,-1,1)`, direct expansion gives

\[
M_{k,\sigma}^{t}JM_{k,\sigma}=\frac{9}{4^k}J.
\tag{18}
\]

The actual exact-valuation guard ensures that (17) lands on the proper primitive integer state triangle. Cone preservation alone does not give a height decrease: the factor is `9/4` on the ascending `k=1` branch.

More generally, any fixed linear functional on the triple restricts to a polynomial of degree at most two in `n` along (1). If this polynomial is a proper height tending to positive infinity, it is eventually strictly increasing. There are arbitrarily large `n=3 mod 4` with `U_+(n)>n`, so that height increases on infinitely many actual edges. This rules out a universal strictly decreasing **fixed linear height** on this spine, including one obtained by a fixed Lorentz coordinate change. It does not exclude adaptive ranks or certificates grouping several edges.

## 8. Exact controls

Run

```text
python 04-computation/experiments/zenodo_triples_20260925.py
python -O 04-computation/experiments/zenodo_triples_20260925.py
```

The [script](../../04-computation/experiments/zenodo_triples_20260925.py) and [output](zenodo_triples_20260925.out) check all marked states `1<=n<=10001`; all 5001 odd spine states in that interval; both signed Lorentz operators for `k=1..12`; the guarded family at `a=1..2500`; the completed root-family certificates at `r=1..24`; and exact curve-group creation/Collatz compatibility at `a=1..4`. Controls retain the primitive parity normalizer, the exact first-return exponent, both Euclid parameters, and the terminal exception.

An independent algebra audit checked (2)–(3), (9)–(11), and (15), including the strictness condition `a>=1` and the necessary marked parameters: **PASS**. The mechanism proves an actual infinite guarded descent family and complete certificates for (16); the remaining universal obligation is to control the ascending spine jumps (8).

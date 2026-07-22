---
id: THM-2091
title: "Centered triple energy is necessary for a rank-seven terminal obstruction"
status: >
  PROVED. Averaging THM-2086's exact relative-Hunter edge identity over the
  uniform Cayley spanning tree of K_7 gives an exact margin formula in the
  total global pair mass, seven mixed-fold charges, and one centered
  guard/danger triple energy. THM-1234 and an exact odd-guard charge census
  imply that every rank-seven guard containment has centered triple energy at
  least 2059/90090. On the live branch 7 does not divide h, every speed
  divisible by 7 has zero charge; one through four such speeds sharpen the
  necessary energy thresholds explicitly. Below the relevant threshold some
  spanning tree supplies safe mass outside the guard. The energy also has an
  exact multiplicity-moment form.
  A hostile control proves the average-tree test is sufficient only and cannot
  replace THM-2081's maximum weighted tree. This is not LRC(14).
source: codex-2026-07-22-LRC-centered-triple-energy
depends_on:
  - THM-1234
  - THM-2081
  - THM-2086
related:
  - THM-1122
  - THM-1166
  - THM-2087
script: 04-computation/lrc14_centered_triple_energy_codex_20260722.py
output: 05-knowledge/results/lrc14_centered_triple_energy_codex_20260722.out
script_sha256: 725cdb1f7f997bc3fb04ae23a9abac8bf26cd75f1f84c6c02aee30e854be5dff
output_sha256: 421d2e111061237c3dff6eaa0373d5e71394419927ce9419aa4a569ed302d9bc
hash_basis: repository blobs with LF line endings
---

# THM-2091 -- centered triple energy obstruction

Let `Q={q_1,...,q_7}` be seven distinct positive speeds and let `h` be odd.
Put

```text
X_i(t)=1_(D_(q_i))(t),             H(t)=1_(E_h)(t),
e_i=measure(D_(q_i) intersect E_h)-2/49,
rho_ij=measure(D_(q_i) intersect D_(q_j)),             (1)
```

and define the centered triple atom

```text
delta_ij
 =integral (X_i-1/7)(X_j-1/7)(H-2/7).                 (2)
```

The two aggregate quantities are

```text
R=sum_(i<j)rho_ij,              D=sum_(i<j)delta_ij.   (3)
```

THM-2086's absolutely convergent channel formula is, in this notation,

```text
w_ij=measure(D_(q_i) intersect D_(q_j) intersect E_h^c)
    =(5/7)rho_ij-(e_i+e_j)/7-delta_ij.                (4)
```

Indeed THM-2086's genuine complement channel `R_h(q_i,q_j)` is
`-delta_ij`: the nonzero Fourier coefficients of `1_(E_h^c)` are the
negatives of those of `1_(E_h)`.

## 1. Exact uniform-tree average

Let `T` be a uniformly random labelled spanning tree of `K_7`, and put

```text
A(Q,h)=E_T[sum_(ij in T)w_ij-Delta_h(Q)],
Delta_h(Q)=2/7-sum_i measure(D_(q_i) intersect E_h)
          =-sum_i e_i.                                (5)
```

Then

```text
A(Q,h)=(10/49)R+(37/49)sum_i e_i-(2/7)D.              (6)
```

### Proof

By Cayley's formula, a fixed edge belongs to a uniformly random tree on seven
labels with probability `2/7`. Therefore

```text
A=(2/7)sum_(i<j)w_ij+sum_i e_i.                       (7)
```

Insert (4). Each `e_i` appears in six unordered pairs, so

```text
(2/7)(-1/7)sum_(i<j)(e_i+e_j)+sum_i e_i
 =(-12/49+1)sum_i e_i
 =(37/49)sum_i e_i.                                   (8)
```

The other two terms give `10R/49` and `-2D/7`, proving (6). QED.

If `A(Q,h)>0`, at least one tree has positive margin. THM-2081 then gives

```text
measure(G_Q minus E_h)>0.                              (9)
```

Thus guard containment forces `A(Q,h)<=0`.

## 2. The sharp seven-charge lower ledger

For one speed `q`, write the reduced THM-2080 variables

```text
a=h/gcd(h,q),                     b=q/gcd(h,q).        (10)
```

Since `h` is odd, `a` is odd. THM-2080 gives

```text
e(q,h)>=-1/(4ab).                                      (11)
```

Among coprime positive pairs `(a,b)` with odd `a`, the seven largest possible
negative charges `-e` are exactly

```text
5/294, 8/539, 3/245, 3/245, 4/441, 4/441, 9/1078.     (12)
```

They occur respectively at

```text
(a,b)=(1,6),(11,1),(3,5),(1,5),(9,2),(9,1),(11,2).    (13)
```

Their sum is

```text
41/495.                                                (14)
```

Consequently seven distinct speeds satisfy

```text
sum_i e_i>=-41/495.                                    (15)
```

### Proof

For `ab>=30`, (11) gives

```text
-e<=1/120<9/1078,                                     (16)
```

so every member of the top-seven list occurs in the finite set `ab<=29`.
Direct substitution in THM-2080's exact fold formula on the coprime odd-`a`
rows in that finite set gives (12)--(13), including ties. For fixed `h`, the
reduced ratio `b/a=q/h` determines `q`; distinct speeds therefore use distinct
rows. The worst sum is the sum of the seven largest entries, proving
(14)--(15). QED.

This is a charge **spectrum**, not seven independent copies of the sharp
`1/42` floor. Distinctness is what prevents repeating `(a,b)=(1,6)` seven
times.

## 3. The centered-energy obstruction

THM-1234's five-subset average gives, for every seven-set,

```text
R>=22/65.                                               (17)
```

Indeed each of the `21` five-subsets has total pair mass at least `44/273`,
and every pair occurs in ten such subsets.

Combining (6), (15), and (17) gives

```text
A(Q,h)
 >=(10/49)(22/65)-(37/49)(41/495)-(2/7)D
 =2059/315315-(2/7)D.                                  (18)
```

Therefore

```text
D<2059/90090
  implies measure(G_Q minus E_h)>0,                    (19)
```

while every guard containment must satisfy

```text
D>=2059/90090=0.022854... .                            (20)
```

### Proof

The constant identity in (18) is exact. If (19)'s strict inequality holds,
the right side of (18) is positive, hence `A>0`; Section 1 supplies a positive
tree and THM-2081 proves the claimed safe mass. Taking the contrapositive gives
(20). QED.

The result complements THM-2085/2087. Those theorems say a survivor has a
height-57 cut of relation channels. Equation (20) says those channels must in
aggregate carry a fixed positive signed mass; mere existence of several short
relations is not enough.

### The live modulo-seven sharpening

Assume now that `7` does not divide `h`, and put

```text
k=#{q in Q:7 divides q}.                               (20a)
```

If `7|q`, then in the reduced variables (10), `7|b`. In the exact THM-2080
fold, `y={b/7}=0`, so `F(x,0)=0` and hence

```text
e(q,h)=0.                                              (20b)
```

Only the other `7-k` speeds can spend negative charge. Summing the first
`7-k` entries of (12) and repeating (18) gives the following exact necessary
energy thresholds:

```text
k   cap on -sum_i e_i       containment forces D >=
0        41/495                    2059/90090
1       3613/48510                 396587/8828820
2       3173/48510                 608227/8828820
3        911/16170                 273289/2942940
4        713/16170                 368527/2942940.      (20c)
```

THM-2086 already restricts its all-height residual to `1<=k<=4`, so the last
four rows, not merely the universal first row, apply there. This modular
ledger is still only necessary. In the bounded terminal bank it closes 29
rows before the exact average margin is evaluated; all remaining bounded rows
are closed by that exact margin.

## 4. Multiplicity-moment form

Let

```text
N(t)=sum_i X_i(t)                                      (21)
```

be the danger multiplicity. Then

```text
D=(1/2) integral (H-2/7)
       [N^2-(19/7)N+6/7].                              (22)
```

### Proof

Pointwise,

```text
sum_(i<j)(X_i-1/7)(X_j-1/7)
 =(1/2)[(N-1)^2-sum_i(X_i-1/7)^2]
 =(1/2)[N^2-(19/7)N+6/7].                             (23)
```

Multiply by `H-2/7` and integrate. QED.

There is also a direct guard-complement form of the average margin:

```text
A(Q,h)=5/7-(1/7)integral_(E_h^c) N(8-N).              (24)
```

To see this, use
`sum_(i<j)X_iX_j=N(N-1)/2`, equation (7), and
`sum_i e_i=integral_(E_h)N-2/7`; then replace the guard integral using
`integral N=1`.

Equations (22)--(24) connect the relative-Hunter residual to the multiplicity
polynomials of THM-1122/1166. The obstruction is not an unlocated count of
relations; it is a signed correlation between guard position and the quadratic
danger-multiplicity profile.

## 5. Average-tree scope and hostile control

The average criterion is not equivalent to the maximum-tree gate. Put

```text
h=1,                  Q=(4,5,6,7,11,13,27).           (25)
```

This row is hereditarily primitive, though not divisor-complete. Exact
atomization gives

```text
Delta=181007/3783780,
A=(2/7)sum_(i<j)w_ij-Delta=-2467/26486460<0,           (26)
```

while the exact maximum spanning tree satisfies

```text
tau_h(Q)-Delta=86717/1891890>0.                        (27)
```

Thus the weighted graphic matroid can close a row that uniform Cayley
averaging misses. Equation (20) is a necessary condition for containment, not
a replacement for THM-2081 and not a characterization of survivors.

## 6. Transfer and assumption challenge

Candidate carriers included runners, global pairs, restricted events, Fourier
modes, relation triples, multiplicity levels, spanning trees, and proof
obligations. The Cayley average preserves the first moment of tree weight and
turns the genuine relation hyperedges into the scalar `D`. It destroys the
edge distribution that selects the maximum tree; (25)--(27) quantify that
loss.

A tournament orientation by larger edge weight is therefore not faithful.
The average uses each undirected edge with probability `2/7`; score histograms,
SCCs, directed cycles, and a tie Hamiltonian path cannot reconstruct either
`D` or the maximum weight. The faithful object remains THM-2081's weighted
graphic matroid, with `D` as a useful necessary sidecar.

The unique/coincident-channel split in the Gaussian Moment work has the same
warning. Counting channels resembles `R`; cancellation depends on their signed
aggregate, analogous to `D`. Unlike LRC, complex Gaussian channel sums have no
positive Hunter tree, so (20) does not transfer as a GMC proof. For the
Jacobian resonance dictionary, the similarity is likewise structural unless
one names a valuation or multiplicity map preserving the predicate.

## Exact referee

The companion exhausts the odd reduced charge spectrum, verifies every
constant in (6), (12)--(20c), and (22)--(24), replays the bounded terminal
bank, checks all modulo-seven charge floors, and checks the hostile
average/max-tree split (25)--(27) by exact rational atomization. It uses
explicit runtime checks; normal and `python -O` transcripts byte-match the
stored output and end in `PASS`. QED.

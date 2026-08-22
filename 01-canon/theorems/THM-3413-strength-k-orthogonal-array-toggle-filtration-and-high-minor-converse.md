---
id: THM-3413
title: "Strength-k orthogonal-array toggle filtration and high-minor converse"
status: PROVED + VERIFIED-EXACT + INDEPENDENTLY PROOF/REPLAY-AUDITED
source: root-2608-crouzeix-puzzle-2026-08-15
audit: combined independent Walsh-sign, Boolean-inverse, parity-converse, repetition, envelope, convolution-spectrum, regression, normal/-O/stored, hash, dependency, and scope audit clean
depends_on:
  - THM-3407-hadamard-core-multitoggle-response-plaquette-shells-and-trade-distance
related:
  - THM-3396-four-bit-pairwise-independent-fourier-cone
  - THM-3411-pairwise-independent-toggle-filter-and-sharp-high-minor-norm
verified_by:
  - 04-computation/hadamard_strength_k_toggle_filtration_thm3413.py
  - 05-knowledge/results/hadamard_strength_k_toggle_filtration_thm3413.out
script: 04-computation/hadamard_strength_k_toggle_filtration_thm3413.py
output: 05-knowledge/results/hadamard_strength_k_toggle_filtration_thm3413.out
script_sha256: dbdf97d32e4114abfd71393791ac8b956d4f1e46eef9a98ec14532107dd0a1a6
output_sha256: 55bbcec2736c72dbcadac16450fe326317971990693ab47e91676b3e51ec4172
semantic_sha256: 0db6229ea8922fb2f9c570103dd8aff7e30da424bcc196f29de19b54e2cab329
hash_basis: LF-normalized bytes
---

# THM-3413 -- strength-k orthogonal-array toggle filtration and high-minor converse

## 1. Statement

Let `H` be a normalized real Hadamard matrix of order `4m`, with sign core
`K` and binary maxdet core `B=(J-K)/2`.  Choose `t>=1` distinct core
positions

```text
e_a=(i_a,j_a),             delta_a=K_(i_a,j_a),       a in [t].       (1)
```

Rows or columns may repeat.  For `S subseteq [t]`, in inherited event order,
put

```text
M_S=(product_(a in S) delta_a)
    det K[(i_a)_(a in S),(j_a)_(a in S)],       M_empty=1.             (2)
```

Repeated selected rows or columns make the corresponding determinant zero;
no nonattacking hypothesis is imposed.  For `x in {+-1}^t`, toggle event `a`
exactly when `x_a=-1`, and write

```text
f(x)=det B_(toggle(x))/det B.                                      (3)
```

Fix `0<=k<=t`.  Let `nu` be any strength-`k` Rademacher law, meaning

```text
E_nu chi_T=0              for 1<=|T|<=k,
chi_T(x)=product_(a in T)x_a.                                      (4)
```

Thus a uniform `OA(N,t,2,k)` is one instance.  Let `u` be the uniform law on
the whole cube and define, for every nonempty `T subseteq [t]`,

```text
C_T=sum_(S superseteq T) (-1)^(|S|+|T|) M_S/(4m)^|S|.               (5)
```

Then the exact strength filtration is

```text
E_nu f-E_u f
 =sum_(|T|>k) (E_nu chi_T) C_T.                                    (6)
```

In particular, a strength-`k` law destroys every determinant interaction of
degree at most `k`; its complete surviving sidecar is the upper Boolean-zeta
transform `(5)` of the signed event minors of orders greater than `k`.

The converse is exact.  The following are equivalent:

1. `E_nu f=E_u f` for every strength-`k` Rademacher law `nu`;
2. the same equality holds for the finite parity-halfcube tests

   ```text
   nu_(T,sigma)(x)=2^(-(t-1)) 1_(chi_T(x)=sigma),
   |T|>k, sigma in {+-1};                                         (7)
   ```

3. `C_T=0` for every `|T|>k`;
4. `M_S=0` for every `|S|>k`.

Indeed the high minors are recovered explicitly from the high response by

```text
M_S=(4m)^|S| sum_(T superseteq S) C_T,             |S|>k.           (8)
```

Thus no lower-order Gram, plaquette, or minor data can replace the stated
high-minor packet in the universal strength-`k` problem.

## 2. Boolean-to-Walsh proof

Give the candidate toggles Boolean coordinates `z_a`.  THM-3407 proves the
complete event identity

```text
F(z)=det(B+sum_a z_a delta_a e_(i_a)e_(j_a)^T)/det B
    =sum_(S subseteq [t])(-1/(2m))^|S| M_S product_(a in S)z_a.     (9)
```

Set `z_a=(1-x_a)/2`.  For each `S`,

```text
product_(a in S)z_a
 =2^(-|S|)sum_(T subseteq S)(-1)^|T|chi_T(x).                      (10)
```

Therefore the Walsh coefficient of `f` at `T` is

```text
fhat(T)=(-1)^|T| sum_(S superseteq T)
                   (-1/(4m))^|S| M_S
       =C_T.                                                       (11)
```

The uniform law kills every nonconstant character, while `(4)` kills those
through degree `k`.  Pairing `(11)` with `nu-u` proves `(6)`.

For `|T|>k`, the law `(7)` has Fourier packet

```text
E_(nu_(T,sigma)) chi_U = sigma if U=T,
                         0     for nonempty U!=T.                  (12)
```

It is uniform on a `2^(t-1)`-run orthogonal array of strength `|T|-1`, hence
of strength `k`.  Equations `(6)` and `(12)` show that its response relative
to uniform is exactly `sigma C_T`.  This proves that the finite tests detect
every high coefficient.

Finally `(5)` is triangular on the reverse-inclusion poset, with nonzero
diagonal `(4m)^(-|T|)`.  Direct Boolean Mobius inversion gives `(8)`.  Hence
all high coefficients vanish exactly when all high signed minors vanish,
proving the four-way equivalence.  QED.

## 3. Integer orthogonal-array compiler

For an `OA(N,t,2,k)`, let

```text
A_T=sum_(rows x) chi_T(x).                                        (13)
```

Then `A_T=0` for `1<=|T|<=k`, and `(6)` becomes the exact integer average

```text
(1/N)sum_(rows x) f(x)-E_u f
  =(1/N)sum_(|T|>k) A_T C_T.                                     (14)
```

There is no matching or completion hidden in `(14)`: the OA rows choose
which toggle subsets are averaged, but they neither extend their columns to
the Hadamard core nor recover individual determinants from the average.

Two boundary specializations explain the earlier finite phenomena.

- For `t=3,k=2`, only the signed `3 by 3` event minor survives; the two
  parity halfcubes detect precisely THM-3407's oriented triple sidecar.
- For `t=4,k=2`, equations `(5)--(6)` give the four cubic plus one quartic
  packet of THM-3411.  THM-3396 supplies an additional sharp norm only in
  this four-bit case; no general high-moment-polytope norm is asserted here.

## 4. Convolution spectrum and the parity boundary

Let `nu^(star r)` be the law of the coordinatewise product of `r`
independent samples from `nu`, and put

```text
mu_T=E_nu chi_T.                                                   (15)
```

Characters multiply under this convolution, so

```text
E_(nu^(star r)) chi_T=mu_T^r.                                     (16)
```

In particular every convolution power remains strength `k`, and `(6)` gives
the exact response sequence

```text
Delta_r=E_(nu^(star r))f-E_u f
       =sum_(|T|>k) C_T mu_T^r.                                   (17)
```

Consequently its formal generating function is rational:

```text
sum_(r>=1) Delta_r z^r
 =sum_(|T|>k) C_T [mu_T z/(1-mu_T z)].                            (18)
```

Thus `Delta_r` is C-finite over the field generated by the moments.  It sees
only the grouped response

```text
D_lambda=sum_(T:mu_T=lambda) C_T,                                 (19)
```

not the individual `C_T` inside a repeated moment eigenvalue.  When the
active nonzero moments are distinct, their Vandermonde system recovers the
individual coefficients from finitely many consecutive responses; equality
of moment values is the exact collision boundary.

If

```text
rho=max_(C_T!=0)|mu_T|<1,
```

then

```text
|Delta_r|<=rho^r sum_(|T|>k)|C_T|.                                (20)
```

For a Rademacher character, `|mu_T|=1` exactly when `chi_T` is constant on
the support of `nu`.  Hence every nondecaying response is a sum of a fixed
part (`mu_T=1`) and a period-two part (`mu_T=-1`), plus exponentially
decaying modes.  The parity laws `(7)` attain this boundary sharply:

```text
Delta_r=sigma^r C_T.                                               (21)
```

This convolution statement concerns the response average, not matrix
multiplication or a Hadamard tensor product.

## 5. A universal envelope

Hadamard's determinant inequality gives

```text
|M_S|<=|S|^(|S|/2).                                               (22)
```

Since every character moment has absolute value at most one, `(5)--(6)`
imply the explicit, generally nonsharp estimate

```text
|E_nu f-E_u f|
 <=sum_(s=k+1)^t binom(t,s)
      [sum_(j=k+1)^s binom(s,j)] s^(s/2)/(4m)^s.                  (23)
```

The exact object is `(5)`, not this absolute-value envelope.  Cancellations
between nested minors can make `C_T` vanish even when its individual summands
do not, while the full converse requires all high `C_T` simultaneously and
then recovers all high `M_S` by `(8)`.

## 6. Connection and loss ledger

| field | content |
|---|---|
| source | `t` distinct sparse toggles in a normalized Hadamard core |
| target | their normalized determinant averaged by a strength-`k` sign law |
| map | THM-3407 Boolean event polynomial, followed by `z=(1-x)/2` and Walsh projection |
| preserved | every interaction above degree `k`, exact OA averages, parity-halfcube detection, repeated-row/column zeros, convolution eigenvalues |
| destroyed | all interactions through degree `k`, toggle localization, core equivalence, OA/Hadamard completion data, and separation inside equal moment eigenvalues |
| necessary sidecar | the signed minors `M_S` for `|S|>k`, equivalently their high zeta transform `C_T` |

The trigger is genuine strength-`k` whitening.  Biased or only partially
decorrelated laws retain lower Walsh layers, so `(6)` with its high-only
right side is then false.  This theorem does not construct a Hadamard matrix,
improve a determinant bound, or imply Grothendieck, Crouzeix, LRC(14), or
JC(2).

## 7. Verification status

Run

```text
python3 04-computation/hadamard_strength_k_toggle_filtration_thm3413.py
python3 -O 04-computation/hadamard_strength_k_toggle_filtration_thm3413.py
```

The standard-library companion checks the full transform and inverse on
abstract packets through `t=8`, all `697,004` character moments of the
parity-halfcube bank through `t=9`, and direct determinants for every one of
the `511` nonempty distinct-position sets in the Paley order-four core.  It
also checks five larger Paley packets, including repeated-row/column and
order-twelve controls, replays all four THM-3411 hostiles, and verifies a
four-eigenvalue convolution response through twelve powers and its exact
order-four recurrence.  Normal and optimized runs are byte-identical to the
frozen output.  The LF-normalized script/output hashes are respectively

```text
dbdf97d32e4114abfd71393791ac8b956d4f1e46eef9a98ec14532107dd0a1a6
55bbcec2736c72dbcadac16450fe326317971990693ab47e91676b3e51ec4172
```

and the semantic digest is
`0db6229ea8922fb2f9c570103dd8aff7e30da424bcc196f29de19b54e2cab329`.
The combined independent immutable-file proof/replay audit reconstructed the
Walsh signs, Boolean inverse, parity converse, repeated-index boundary,
envelope, convolution spectrum, and THM-3411 regression, and found no defect.
Normal, optimized, and stored outputs and every pinned hash agree.

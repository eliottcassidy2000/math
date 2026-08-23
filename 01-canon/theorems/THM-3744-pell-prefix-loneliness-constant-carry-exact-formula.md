---
id: THM-3744
title: "Pell-prefix loneliness constant-carry exact formula"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  For every initial
  prefix S_N={P_1,...,P_N} of the Pell numbers, its lonely-runner maximum is
  exactly A_N/(P_N+1), where A_N=(P_N-P_{N-1}+1)/2.  The maximizer is unique
  on [0,1/2], hence unique on the circle up to reflection.  Its runner
  residues have an exact constant-carry recurrence and a closed profile,
  symmetric on the two odd subsequences.  At odd square-triangular indices
  the numerator or reduced phase is a square/Pell factor.  This is a very safe
  structured family; it proves
  neither LRC(14) nor a planar-Jacobian statement.
source: root + cross_domain_skeptic / 2026-08-23
audit: >
  PASS.  A separate derivation supplied the interval-collapse upper proof.
  The exact companion checks 161,773 gates through N=200, independently
  exhausts every lower-envelope peak candidate through N=9, verifies both odd
  profile families, and passes normal/-O byte identity.  Single-speed and zero-speed
  hostiles destroy the claimed witness as predicted.
depends_on:
  - THM-3335-square-triangular-pell-markov-pythagorean-selector
related:
  - THM-778-centered-christoffel-endpoint-skew-product
  - THM-3570-universal-pell-conic-target-graph-factor-compiler
  - THM-3573-polynomial-target-graph-pell-parameter-descent-classification
  - THM-3575-coprime-pell-target-coordinate-euler-quotient-no-go
  - THM-3736-automorphic-cohn-complete-constant-sl2-polynomial-exposure-classification
  - THM-3740-automorphic-cohn-one-variable-right-shear-binomial-tower-classification
  - THM-3742-square-triangular-pell-mod13-central-sign-projective-cycle
script: 04-computation/pell_prefix_loneliness_constant_carry_thm3744.py
output: 05-knowledge/results/pell_prefix_loneliness_constant_carry_thm3744.out
script_sha256: 657aaf0a573a81588797e0c2d5984644901d81980ecc7a9f861e8f13e05c9860
output_sha256: 38cfd474d7b2c1c46c0b557969d710b2c53cfb423d11d307d83d736a9dd3934c
hash_basis: raw LF bytes
---

# THM-3744 -- every Pell prefix has an exact loneliness maximum

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

Because the legacy numeric namespace also contains `HYP-3744`, cite this
result as **THM-3744 plus its full theorem slug/file path**, never by the bare
number.

No literature-priority claim is made.  This is an elementary exact theorem
about a structured lonely-runner family, recovered from the square-triangular
Pell lane and proved here at all prefix lengths.

## 1. Inheritance, notation, and theorem

The closest proved mechanism is
[THM-3335](THM-3335-square-triangular-pell-markov-pythagorean-selector.md),
which identifies the square-triangular Pell ladder and the consecutive
coordinates `(a,s,x,b)`.  The corrected near miss is
[THM-3742](THM-3742-square-triangular-pell-mod13-central-sign-projective-cycle.md):
its mod-13 Pell clock has exactly the desired fourteen projective states, but
ordinary projectivization loses the central sign and its attractive offset
match lives on the wrong torsor.  The least-used sidecar here is not another
scalar residue set; it is the **integer carry at every recurrence step**.

Let

```text
P_0=0,  P_1=1,  P_(n+2)=2P_(n+1)+P_n                 (1)
```

be the Pell numbers, and put

```text
A_n=(P_n-P_(n-1)+1)/2,       U_n=A_n/(P_n+1).         (2)
```

The numerator in `(2)` is even.  For a finite positive speed set `S`, use

```text
M(S)=max_(t in R/Z) min_(v in S) ||v t||_(R/Z).       (3)
```

### Theorem 1.1 (all Pell prefixes)

For every `N>=1`, let `S_N={P_1,...,P_N}`.  Then

```text
M(S_N)=U_N=A_N/(P_N+1).                               (4)
```

The unique maximizer in `[0,1/2]` is `t=U_N`.  Thus the maximizers on the
circle are `U_N` and `1-U_N`, except that these coincide when `N=1`.

For `N=13`, the repository's LRC(14) convention gives thirteen nonzero
speeds, and `(4)` specializes to

```text
M({1,2,5,12,29,70,169,408,985,2378,5741,13860,33461})
  =99/338.                                             (5)
```

This is far above `1/14`.  It is a positive structured control, not a hard
LRC(14) candidate.

## 2. Two Pell inequalities

Write

```text
lambda=1+sqrt(2),              alpha=1-1/sqrt(2).      (6)
```

The Pell Binet formula gives the exact error identity

```text
A_n-alpha P_n=(1-(-1)^n lambda^(-n))/2.               (7)
```

Consequently `A_n=ceil(alpha P_n)`, and for every `n>=2`,

```text
(A_n-1)/(P_n-1) < alpha < A_n/(P_n+1)=U_n.            (8)
```

Indeed, the error in `(7)` lies strictly between `alpha` and `1-alpha` once
`n>=2`; rearranging gives `(8)`.

The candidate phases are nonincreasing.  A direct Cassini calculation gives

```text
A_n(P_(n+1)+1)-A_(n+1)(P_n+1)
  =(P_n-P_(n-1)+(-1)^(n+1))/2 >=0.                    (9)
```

Equality occurs only at `n=2`, so `U_2=U_3=1/3` is the only adjacent tie.

## 3. The phase U_N is safe

Fix `N` and write `delta=U_N`.  Monotonicity `(9)` and `(8)` imply, for every
`1<=i<=N`,

```text
A_i-1+delta <= P_i delta <= A_i-delta.                (10)
```

For the right inequality, `delta=U_N<=U_i`.  For the left inequality,
`i=1` is equality, while for `i>=2`,

```text
delta>alpha>(A_i-1)/(P_i-1).
```

Thus

```text
floor(P_i delta)=A_i-1,
delta <= {P_i delta} <= 1-delta.                      (11)
```

Every runner is at distance at least `delta`, so `(11)` proves the lower
bound `M(S_N)>=U_N`.

## 4. The safe interval collapses to one point

For `N=1`, the assertion is immediate.  Let `N>=2`, and suppose
`t in [0,1/2]` has every runner at distance at least `delta=U_N`.

Runner `P_1=1` gives `t>=delta`.  Runner `P_2=2` gives the initial upper
bound

```text
t <= (1-delta)/2=(A_2-delta)/P_2.                    (12)
```

We inductively prove

```text
t <= (A_i-delta)/P_i                 (2<=i<=N).       (13)
```

Define

```text
Delta_i=P_(i+1)A_i-P_i A_(i+1)
       =(P_(i+1)-P_i+(-1)^(i+1))/2.                  (14)
```

For `i>=2`, the Pell recurrence shows

```text
Delta_i/(P_(i+1)+P_i) < 1/4 < alpha < delta.          (15)
```

For the only nontrivial sign in the first inequality, `i` is odd and
`P_i-P_(i-1)>2`; the case `i=2` is direct.  Assume `(13)` at `i`.  From
`t>=delta`, `(8)`, and `(15)`,

```text
A_(i+1)-1+delta
  < P_(i+1)t
  < A_(i+1)+delta.                                    (16)
```

The upper inequality is exactly

```text
P_(i+1)(A_i-delta)/P_i < A_(i+1)+delta,
```

after moving terms and using `Delta_i<delta(P_(i+1)+P_i)`.  Since runner
`P_(i+1)` is safe, `(16)` cannot lie in the open `delta`-neighbourhood of
the integer `A_(i+1)`.  It must therefore satisfy

```text
P_(i+1)t <= A_(i+1)-delta,
```

which is `(13)` at the next index.  At `i=N`,

```text
(A_N-delta)/P_N=delta.                                (17)
```

Hence `t<=delta`, while runner one gave `t>=delta`.  Thus `t=delta`.
Reflection handles `[1/2,1]`, proving Theorem 1.1 and uniqueness.

## 5. Residue packet and constant carry

Put `D=P_N+1`.  Equation `(11)` lets us use the least nonnegative residue

```text
R_(i,N)=P_i A_N-(A_i-1)D.                             (18)
```

The standard d'Ocagne identity

```text
P_(i-1)P_N-P_i P_(N-1)=(-1)^i P_(N-i)
```

and `(2)` give the closed formula

```text
2R_(i,N)=D+P_(i-1)+(-1)^i P_(N-i).                   (19)
```

Therefore every runner distance is known without an optimization:

```text
||P_i U_N||
 =[D-|P_(i-1)+(-1)^i P_(N-i)|]/(2D).                 (20)
```

If `y_i={P_iU_N}`, the affine recurrence for `A_i`,

```text
A_(i+2)=2A_(i+1)+A_i-1,                              (21)
```

turns the Pell recurrence into

```text
y_(i+2)=2y_(i+1)+y_i-1.                              (22)
```

The carry is exactly `+1` at every internal step.  Unlike a geometric mean
of continued-fraction digits, `(18)--(22)` retain runner index, phase,
orientation, and temporal order.  This is why the Pell word is useful as a
model carrier even though the family is nowhere near extremal for LRC(14).

## 6. Where the squares and triangular numbers enter

For `k>=1`, define the consecutive Pell coordinates inherited from THM-3335:

```text
a_k=P_(2k-1),  s_k=P_(2k),
x_k=P_(2k)+P_(2k-1),  b_k=P_(2k+1).                  (23)
```

They obey

```text
x_k^2-2s_k^2=1,
T_((x_k-1)/2)=(s_k/2)^2,
a_k b_k=s_k^2+1.                                     (24)
```

Pell addition identities also give

```text
A_(4k-1)=2a_k^2,       P_(4k-1)+1=2a_k x_k,
A_(4k+1)=x_k^2,        P_(4k+1)+1=2x_k b_k.          (25)
```

Thus the two odd subsequences of `(4)` reduce to

```text
M(S_(4k-1))=a_k/x_k,
M(S_(4k+1))=x_k/(2b_k).                              (26)
```

The first rows are

| `k` | `(a_k,s_k,x_k,b_k)` | `M(S_(4k-1))` | `M(S_(4k+1))` |
|---:|:---|:---|:---|
| 1 | `(1,2,3,5)` | `1/3` | `3/10` |
| 2 | `(5,12,17,29)` | `5/17` | `17/58` |
| 3 | `(29,70,99,169)` | `29/99` | `99/338` |

The classical cannonball row is therefore not decorative:

```text
29*169=70^2+1,
T_49=35^2,
A_11=2*29^2,
A_13=99^2.                                           (27)
```

Both subsequences in `(26)` converge to `alpha=1-1/sqrt(2)`.

## 7. Symmetric profiles

Formula `(20)` and the Pell addition laws simplify further at the odd
square-triangular indices.  For `N=4k+1`, set `c=2k+1`.  Then

```text
||P_iU_N||=(P_c-P_|c-i|)/(2P_c),       1<=i<=N.       (28)
```

For `N=4k-1`, set `c=2k`,

```text
C_j=P_(j+1)+P_j,       d_i=max(|i-c|-1,0).
```

Then

```text
||P_iU_N||=(x_k-C_(d_i))/(2x_k),       1<=i<=N.       (29)
```

These follow by substituting the Pell sum/difference formulas into `(20)`.
The `4k+1` profile has one central peak; the `4k-1` profile has a three-runner
central plateau.  At `N=13`, the distance numerators over the common
denominator `338` are

```text
99,140,157,164,167,168,169,168,167,164,157,140,99.    (30)
```

This exact owner-labelled profile is much richer than the unordered square
or triangular residue images in THM-3742.

## 8. Khinchin, Jacobian, and hostile boundaries

The Pell recurrence comes from the eventually constant continued fraction

```text
alpha=[0;3,2,2,2,...].                                (31)
```

It is a quadratic, measure-zero orbit.  Khinchin's almost-everywhere theorem
does not explain `(4)`: after the initial digit, the geometric mean is
deterministically `2`, while the proof needs the ordered convergents and the
carry packet `(18)--(22)`.  This is a clean example of why “Khinchin content”
must be split into metric digit statistics, ordered Euclidean words, and
geometry-of-numbers flatness.

Two hostile changes locate the boundary:

- At the `N=13` phase `99/338`, replacing the central speed `P_7=169` by
  `170` drops that runner's distance to `35/169<99/338`.
- Adjoining `P_0=0` forces the maximum to zero.

Thus neither a nearby sequence nor a shifted indexing inherits the theorem.

The constant Pell matrices belong to the constant `SL_2` lane already
classified in
[THM-3736](THM-3736-automorphic-cohn-complete-constant-sl2-polynomial-exposure-classification.md).
They do not supply the interacting variable polynomial factors required by
[THM-3740](THM-3740-automorphic-cohn-one-variable-right-shear-binomial-tower-classification.md)
or construct a polynomial pair with constant Jacobian.  The apparent shared
matrix language therefore stops at a typed analogy: no planar-Jacobian
consequence follows.

## 9. Exact reproduction and scope

Run

```bash
python3 -B 04-computation/pell_prefix_loneliness_constant_carry_thm3744.py
python3 -B -O 04-computation/pell_prefix_loneliness_constant_carry_thm3744.py
```

The companion performs `161,773` exact gates.  It replays the symbolic upper
certificate through `N=200`, independently enumerates every lower-envelope
candidate through `N=9`, checks `(19)--(22)`, both profiles through `k=8`,
and the two hostiles.  Normal and optimized outputs are byte-identical.

The theorem supplies an infinite exactly solved LRC family, an owner-labelled
Pell carry model, and a lossless square-triangular interpretation on the two
odd subsequences.  It does not improve the `1/14` frontier, prove LRC(14), or
advance JC(2).

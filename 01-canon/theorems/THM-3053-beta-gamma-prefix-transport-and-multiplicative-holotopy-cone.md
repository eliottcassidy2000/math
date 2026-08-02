---
id: THM-3053
title: "Beta-Gamma prefix transport and multiplicative holotopy cone"
status: >
  PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT IMMUTABLE-FILE AUDIT
  REQUESTED.  A finite signed integer Gamma-shape inventory admits a
  cancellation-free factorization into Gamma moments and forward Beta ratios
  exactly when every prefix sum is nonnegative.  The canonical factorization
  is the oriented-path cut flow.  Applied to THM-3047, this enlarges the
  coordinatewise Gamma-flow cone to the exact inequalities d_0<=A and
  d_j<=I for j>=1.  Prefix positivity characterizes this multiplicative
  Beta-Gamma cone, not the full Stieltjes cone: a strict convex-mixture escape
  and sharp H2/H3 hostiles mark both sides.
source: kind-pasteur-2026-08-01-beta-gamma-prefix-transport
depends_on:
  - THM-3047-formal-corner-width-product-gamma-moment-and-strict-hankel-positivity
  - THM-3051-stieltjes-multiplier-gamma-flow-and-moving-lower-hankel-boundary
related:
  - THM-2828-lower-prefix-cone-factorial-moment-three-detection
  - THM-3021-the-hadamard-multiplier-question-is-false-and-sfc2-is-appell-squarefreeness
script: 04-computation/gmc_beta_gamma_prefix_transport_thm3053.py
output: 05-knowledge/results/gmc_beta_gamma_prefix_transport_thm3053.out
script_sha256: b61b007f0fa13a39e1c4d4847e53de91a65a060dee14defa44adbd0114db298b
output_sha256: bdc1fb1e495c53bfaeef42271e91eefff5d07289d3ef50a405ab38aa8cc6efd9
hash_basis: LF-normalized bytes
---

# THM-3053 -- prefix flow is the exact Beta-Gamma transport gate

**PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT AUDIT REQUESTED.**

THM-3051 gives an adjacent Gamma-shape preservation cone, but its own Beta
escape shows that coordinatewise nonnegative shape multiplicities are too
strict.  The missing coordinate is cumulative flow.  A negative exponent at
a high shape is lawful when it is the denominator of a Beta ratio fed by
earlier positive shape mass.

This produces an exact ordered transport theorem.  It also separates three
notions which should not be conflated:

```text
coordinatewise Gamma inventory
  proper subset of Beta-Gamma multiplicative prefix cone
  proper subset of the full Stieltjes cone.                             (1)
```

The last boundary is witnessed explicitly below: some prefix-negative ratios
are strict Stieltjes moments by convex mixture, while other first-prefix
violations fail at Hankel order two or three.

## 1. The prefix theorem

Fix `a>0`, `c>0`, an integer `N>=0`, and a finite signed integer inventory

```text
n=(n_0,...,n_N),
S_j=sum_(r=0)^j n_r.                                    (2)
```

Consider the positive sequence

```text
G_M=c^M product_(j=0)^N (a+j)_M^(n_j), M>=0.          (3)
```

The following are equivalent.

1. Every prefix satisfies `S_j>=0`.
2. The exponent vector `n` is a nonnegative integral sum of Gamma residues
   `e_j` and forward Beta edges `e_i-e_j`, `0<=i<j<=N`.
3. `(3)` has a cancellation-free representation as the moments of a scaled
   product of independent variables of the forms

```text
Gamma(a+j,1),
Beta(a+i,j-i), 0<=i<j<=N.                              (4)
```

There is a canonical adjacent representation:

```text
G_M=c^M (a+N)_M^(S_N)
       product_(j=0)^(N-1)
       [(a+j)_M/(a+j+1)_M]^(S_j).                     (5)
```

The quotient in `(5)` is the moment sequence of `Beta(a+j,1)`.  Thus `(5)`
constructs `(4)` directly whenever the prefixes are nonnegative.

Conversely, let `g_j` be the number of Gamma factors at level `j`, and let
`b_ij` count forward Beta edges.  Their exponent boundary is

```text
n_j=g_j+sum_(k>j)b_jk-sum_(i<j)b_ij.                  (6)
```

Summing `(6)` across the cut after `m` gives

```text
S_m=sum_(j<=m)g_j+sum_(i<=m<j)b_ij>=0.                (7)
```

This proves necessity and identifies the prefix as an oriented-path cut
invariant, not a coincidental telescoping trick.

## 2. Strict total positivity and its rank-one boundary

Under the equivalent conditions above, `(G_M)` is a Stieltjes moment
sequence.  It is strictly generalized-Hankel totally positive if and only if

```text
n !=0, equivalently some S_j>0.                        (8)
```

Indeed, any Gamma or Beta factor in `(4)` has a continuous positive density
on `(0,infinity)` or `(0,1)`, respectively.  A nonempty finite independent
product therefore has infinite positive support, and the
generalized-Vandermonde/Andreief argument applies.  If every `S_j=0`, then
`n=0` and `G_M=c^M`, represented by the point mass at `c`; all Hankel minors
of order at least two vanish.

This strictness statement concerns the moment sequence, not uniqueness of a
factorization.  Different edge flows can represent the same law.

## 3. Multiplicative holotopy by Beta subdivision

For `i<j<k`, independent Beta variables satisfy

```text
Beta(a+i,j-i) * Beta(a+j,k-j)
  has the Beta(a+i,k-i) law.                            (9)
```

Their moments telescope:

```text
(a+i)_M/(a+j)_M * (a+j)_M/(a+k)_M
 =(a+i)_M/(a+k)_M.                                    (10)
```

Both sides are supported on `[0,1]`, where moments determine the measure, so
`(10)` proves `(9)`.  Every long forward edge can therefore be subdivided
into adjacent edges.  The companion identity

```text
Gamma(a+j,1) * Beta(a+i,j-i)
  has the Gamma(a+i,1) law                              (11)
```

follows from the same Mellin cancellation and moves an unmatched Gamma
residue along the path.  Every feasible flow thereby expands to or contracts
from `(5)`.  We call this exact split/contract calculus the
**multiplicative holotopy**.  It is a measure-level equivalence generated by
`(9)` and `(11)`, not a claim that arbitrary moving-lower resultants are
topologically equivalent.

This is the multiplicative Mellin analogue of THM-2828's additive factorial
prefix cone.  Both use the oriented path incidence/cut map.  Addition of
adjacent factorial directions there is replaced by multiplication of Beta
ratios here; THM-2828's cubic detector is not transported.

## 4. Exact enlargement of THM-3047's Gamma-flow cone

Retain THM-3047's inventory

```text
e_0=A, e_1=B, e_j=0 (j>=2), I=A+B,                   (12)
```

and THM-3051's finitely supported integer transfers

```text
L_M=c^M product_(j>=0)T_j(M)^(d_j),
T_j(M)=(a+j+M)/(a+j),
n_j=e_j+d_(j-1)-d_j, d_(-1)=0.                       (13)
```

The prefix sums telescope much more strongly than the individual exponents:

```text
S_0=A-d_0,
S_m=I-d_m for m>=1.                                  (14)
```

Since `d_m=0` eventually, the terminal prefix is always `I>0`.  Hence

```text
F_M L_M has a strict Beta-Gamma product representation
  iff d_0<=A and d_m<=I for every m>=1.               (15)
```

Here “iff” refers exactly to the cancellation-free Gamma/forward-Beta
factorizations of Section 1.  It is not an iff for arbitrary Stieltjes
representability.

Condition `(15)` strictly enlarges THM-3051's coordinatewise Gamma cone
`n_j>=0`.  For example, at `k=2`, `d_0=-1` gives `n=(2,-1)`: the second
coordinate is negative, while the prefixes `(2,1)` are positive.  Formula
`(5)` is precisely the `Gamma(a)*Beta(a,1)` representation found in
THM-3051.

## 5. Why the prefix cone is not the full Stieltjes cone

The prefix condition is necessary for a **product** of the carriers `(4)`,
but convex addition creates new measures.  Take

```text
n=(1,-2,1), S=(1,-1,0).                               (16)
```

Although the middle prefix is negative,

```text
(a)_M(a+2)_M/(a+1)_M^2
 =a/(a+1)+[1/(a+1)] (a)_M/(a+1)_M.                   (17)
```

Thus `(17)` is represented by the probability measure

```text
[a/(a+1)] delta_1 + [1/(a+1)] Beta(a,1).             (18)
```

It has infinite support and is strictly generalized-Hankel totally positive.
So prefix nonnegativity is not necessary for Stieltjes representability.
The first operation outside the multiplicative cone is already a two-term
convex mixture.

There is nevertheless a useful necessary boundary for balanced ratios.  Put
`D=sum_j n_j`.  If `(3)` is Stieltjes, log-convexity makes

```text
q_M=G_(M+1)/G_M=c product_j(M+a+j)^(n_j)              (19)
```

nondecreasing.  Its asymptotic first forces `D>=0`.  If `D=0`, then
`q_M->c`, so `q_M<=c` and
`G_M^(1/M)->c`.  Every representing measure is therefore supported on
`[0,c]`; equivalently `G_M/c^M` must be a Hausdorff moment sequence and obey

```text
(-1)^r Delta^r(G_M/c^M)>=0 for all M,r>=0.            (20)
```

In the balanced case, let `P_r=sum_j n_j j^r` and take the first `r>=1` with
`P_r!=0`.  Expansion of
`log(q_M/c)` at infinity and `q_M<=c` give the necessary sign

```text
(-1)^r P_r>0.                                         (21)
```

For the mixture `(17)`, these inequalities hold strictly for `r>=1`.

## 6. Sharp failures outside the cone

Convex mixtures show that leaving the prefix cone need not be fatal.  But no
uniform theorem survives the first negative prefix.

For `n=(-1,1)`,

```text
G_M=c^M(a+1)_M/(a)_M=c^M(a+M)/a,
det[[G_0,G_1],[G_1,G_2]]=-c^2/a^2<0.                 (22)
```

The next hostile hides behind strict adjacent log-convexity.  At `a=c=1`
and `n=(-1,2)`,

```text
G_M=M!(M+1)^2,
G_(M-1)G_(M+1)>G_M^2 for every M>=1,                 (23)
```

but

```text
det[G_(i+j)]_(0<=i,j<=2)=-24.                        (24)
```

This is THM-3051's curvature hostile reinterpreted as the first-prefix cut
failure.  It proves that neither positivity of the terms nor every adjacent
Hankel inequality can replace the transport sidecar.

## 7. Scope and next operation

The theorem completely classifies **multiplicative** factorizations made from
the Gamma/Beta carriers `(4)`.  It does not classify all Stieltjes sequences
of the gamma-ratio form `(3)`, and it supplies no representation for the
non-Pochhammer low resultants in THM-3051's literal moving-lower hostiles.

The source-to-target contract is

```text
signed shape inventory n
  --oriented prefix/cut map--> nonnegative transport flow
  --Beta-Gamma Mellin factors--> strict Stieltjes width multiplier.       (25)
```

The map preserves the complete signed inventory and its cut feasibility.  A
mere moment sequence forgets the chosen product decomposition; convex mixing
in `(17)` leaves the transport category entirely.  The cheapest next test for
a moving low resultant is therefore to identify a dominant Beta-Gamma face
and bound the remaining faces, rather than demand an exact factorization of
the whole correction.

No physical-width, raw-chart, GMC, NC2, LRC, Jacobian, or tournament
consequence follows from `(15)` without such a map.

## 8. Exact companion

The dependency-free referee checks:

- all `19,607` signed inventories of lengths `1..5` with entries in
  `[-3,3]`, with `6,044` prefix-feasible cases, against an independent greedy
  long-edge transport;
- `312` materialized Beta-Gamma moment identities and `648` strict Hankel
  controls;
- `540` edge split/contract identities `(9)--(10)`;
- `9,604` exact THM-3047 transfer cells for `(13)--(15)`;
- `44` mixture identities, `120` strict Hausdorff differences, and both
  negative-prefix hostiles `(22)--(24)`.

Run

```text
python 04-computation/gmc_beta_gamma_prefix_transport_thm3053.py
python -O 04-computation/gmc_beta_gamma_prefix_transport_thm3053.py
```

Both modes equal the stored eight-line transcript after LF normalization.

**QED, pending independent immutable-file audit and status promotion.**

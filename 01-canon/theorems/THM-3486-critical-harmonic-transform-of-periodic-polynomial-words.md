---
id: THM-3486
title: "Critical harmonic transform of periodic-polynomial words"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.  For a nonzero complex
  period-p polynomial word of degree d,
  the critical sum through N is c H_N+C+O(1/N), where c is the mean leading
  lane coefficient, equivalently the top coefficient of the trivial Fourier
  colour.  After subtracting c n^d, the critical series is absolute exactly
  when the leading lane coefficients are all equal; otherwise it is
  conditional.  The constant is an explicit root-of-unity polylogarithm
  packet.  Periodic Boolean subsets are the degree-zero specialization.
source: codex-2026-08-16-critical-harmonic-transform
audit: >
  independently audited by death-star-2026-08-16: proof, trichotomy, and
  specializations sound after explicitly naming the Riemann zeta in (19);
  exact companion pins THM-3485 and THM-3455, checks four polynomial packets,
  the THM-3484 ternary word, exact Abel identities, both Berggren cap-seven
  periods, a K4/XOR Boolean-address control, normal/optimized replay, and
  security
depends_on:
  - THM-3485-periodic-polynomial-fourier-jordan-recurrence-classification
  - THM-3455-berggren-q-spine-cap-seven-atom-sieve-and-fibonacci-rank-spectrum
related:
  - THM-3484-ternary-weighted-determinant-minimal-recurrence-and-cubic-fourier-degree-drop
  - THM-3334-berggren-parabolic-spine-gaussian-collision-torsor
script: 04-computation/periodic_polynomial_critical_harmonic_transform_thm3486.py
output: 05-knowledge/results/periodic_polynomial_critical_harmonic_transform_thm3486.out
script_sha256: 0caa1fff76e4d18027f9e08efef3ac21cba0c4b6dd9964517aa93c6675f29202
output_sha256: f6b3b3afe749974fc45b1b70da6eeee28368e1193bd986c574a359cf931fed1d
semantic_sha256: d30150749708eee29dc4f1773fc1e757e6c8f0df0b7b707b509023e8e06db748
hash_basis: LF-normalized bytes
---

# THM-3486 -- only the trivial top colour carries logarithmic mass

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

This theorem connects two proved lanes without identifying their objects.
THM-3485 classifies the finite shift module of a periodic-polynomial word by
its Fourier/Jordan colours.  THM-3455 computes natural and harmonic
coefficients for two periodic Boolean subsets of recurrence indices.  The
missing general statement is that the same **top Fourier layer** completely
classifies critical harmonic divergence, conditional cancellation, and
absolute cancellation.

## 1. Setup

Fix `p>=1` and polynomials

```text
P_0(t),...,P_(p-1)(t) in C[t],                         (1)
```

not all zero.  Define, for `n>=1`,

```text
a_n=P_r(n) when n=r mod p.                             (2)
```

Let

```text
d=max_r deg P_r,
ell_r=[t^d]P_r,
c=(1/p) sum_(r=0)^(p-1) ell_r.                        (3)
```

Here `c` is the mean leading lane coefficient.  With a primitive `p`th root
`zeta`, use THM-3485's Fourier colours

```text
Q_j(t)=1/p sum_r zeta^(-jr)P_r(t),
a_n=sum_j zeta^(jn)Q_j(n).                            (4)
```

Then

```text
c=[t^d]Q_0.                                           (5)
```

Thus `c` is literally the top rung of the trivial shift eigenvalue, not an
average introduced after the recurrence calculation.

## 2. Theorem and sharp trichotomy

Write `H_N=sum_(n=1)^N 1/n`.  There is a finite constant `C(a)` such that

```text
boxed:
sum_(n=1)^N a_n/n^(d+1)
  =c H_N+C(a)+O(1/N).                                 (6)
```

Equivalently, the renormalized series

```text
sum_(n>=1) (a_n-c n^d)/n^(d+1)                       (7)
```

converges.  Its convergence type is exact:

```text
(7) converges absolutely
  iff ell_0=ell_1=...=ell_(p-1)=c,

(7) converges conditionally
  iff the leading coefficient vector is nonconstant. (8)
```

Consequently, because `d` is the actual maximum degree,

```text
c=0  iff the original critical series converges,
and then that convergence is conditional.             (9)
```

The last assertion excludes the zero word by hypothesis.  If the word is
zero, take `d=-infinity`, minimal recurrence one, and critical transform zero;
there is no analytic boundary to classify.

## 3. Abel-summation proof

Set

```text
beta_n=ell_(n mod p)-c.                                (10)
```

The period sum of `beta` is zero.  Hence its partial sums

```text
B_N=sum_(n=1)^N beta_n                                (11)
```

are bounded.  Lane by lane,

```text
a_n-c n^d=beta_n n^d+O(n^(d-1)),                      (12)
```

where the implied constant is uniform because there are only `p` lanes.
After division by `n^(d+1)`,

```text
(a_n-c n^d)/n^(d+1)=beta_n/n+O(1/n^2).                (13)
```

Finite Abel summation gives the exact identity

```text
sum_(n=M+1)^N beta_n/n
 = B_N/N-B_M/(M+1)
   +sum_(n=M+1)^(N-1) B_n(1/n-1/(n+1)).               (14)
```

Boundedness of `B_n` makes the tail in `(14)` `O(1/M)`.  The second term in
`(13)` is absolutely summable with the same tail order.  This proves `(6)`
and convergence of `(7)`.

If all `ell_r=c`, then `beta=0`, so `(13)` is `O(1/n^2)` and `(7)` is
absolute.  Conversely, if some `beta_r!=0`, then on that residue class

```text
|(a_n-c n^d)/n^(d+1)|=|beta_r|/n+O(1/n^2).            (15)
```

The reciprocal sum on one nonempty residue class diverges, so `(7)` is not
absolute.  This proves `(8)`.  If `c=0`, at least one `ell_r` is nonzero by
the definition of `d`, so `(9)` follows.  QED.

## 4. Exact Fourier/polylogarithm constant

Write `zeta_R(s)` for the Riemann zeta function, and write

```text
Q_j(t)=sum_(k=0)^d q_(j,k)t^k,                         (16)
```

allowing zero coefficients.  For an integer `s>=2`, let

```text
Li_s(lambda)=sum_(n>=1) lambda^n/n^s.                  (17)
```

For a nontrivial root of unity, use the Dirichlet/Abel boundary value

```text
Li_1(lambda)=sum_(n>=1)lambda^n/n=-Log(1-lambda),      (18)
```

where `Log` is the branch reached radially from the unit disk.  Substitution
of `(4)` into `(6)` gives the exact constant

```text
boxed:
C(a)=sum_(k=0)^(d-1) q_(0,k) zeta_R(d+1-k)
     +sum_(j=1)^(p-1) sum_(k=0)^d
          q_(j,k) Li_(d+1-k)(zeta^j).                 (19)
```

The omitted `j=0,k=d` term is exactly `c H_N`.  Formula `(19)` makes the
separation literal:

- the trivial top colour is the unique logarithmic channel;
- nontrivial top colours contribute finite `Li_1` boundary values; and
- every lower rung contributes an absolutely convergent polylogarithm.

For real lane coefficients, conjugate colours pair and make `(19)` real.
No claim that these constants are algebraic is made.

## 5. Exact connection to THM-3485's recurrence classifier

Finite Fourier inversion gives

```text
ell_0=...=ell_(p-1)
iff [t^d]Q_j=0 for every j!=0.                         (20)
```

Therefore the following three statements are equivalent:

1. every nontrivial Fourier colour loses its top Jordan rung;
2. all leading lane coefficients are equal; and
3. the renormalized critical series `(7)` converges absolutely.

If the leading vector is nonconstant, at least one nontrivial colour retains
degree `d`, and `(7)` is only conditionally convergent.  If its mean also
vanishes, the trivial colour loses the top rung and the **unrenormalized**
critical series converges conditionally.

The sharp hostile is the period-two word

```text
P_0(t)=t^d,            P_1(t)=-t^d.                    (21)
```

Then `a_n=(-1)^n n^d`, THM-3485 gives the sole primary factor
`(x+1)^(d+1)`, and the critical series is

```text
sum_(n>=1)(-1)^n/n=-log 2.                            (22)
```

Thus a full nontrivial Jordan block can carry no logarithmic mass while still
preventing absolute convergence.  Recurrence order and harmonic mass are
different projections of the same colour packet, not equivalent invariants.

## 6. Boolean subsets of the natural numbers

Let `A` be periodic modulo `p`, and put

```text
a_n=1_A(n),              delta=|A mod p|/p.            (23)
```

This is the degree-zero specialization.  Formula `(6)` becomes

```text
boxed:
sum_(n<=N,n in A) 1/n
  =delta H_N+C_A+O(1/N),                              (24)
```

and `(19)` becomes the finite cyclotomic-log packet

```text
C_A=-sum_(j=1)^(p-1) Q_j Log(1-zeta^j).               (25)
```

Hence every nonempty periodic subset has a divergent reciprocal subseries,
with both natural and harmonic coefficient `delta`.  If the subset is proper
and nonempty, the centered series

```text
sum_n (1_A(n)-delta)/n                                (26)
```

converges conditionally.  For the empty and full subsets it is identically
zero.

Every subset of `N` does select a subseries of the harmonic series, but
periodicity is load-bearing here.  An arbitrary subset need not have natural
or logarithmic density, and its reciprocal subseries may converge (powers of
two) or diverge.  Equations `(24)--(25)` are a finite-state theorem, not a
classification of all subsets.

## 7. Two proved specializations sharpen immediately

### THM-3484's ternary determinant word

Its three degree-seven lanes have the common leading coefficient `-16384`.
Therefore

```text
sum_(n=1)^N A_n/n^8
  =-16384 H_N+C_A+O(1/N),                             (27)
```

and the series after adding `16384 H_N` is absolutely convergent.  This is
the analytic shadow of the same two top-rung losses that lower the recurrence
order from 24 to 22.  It is not an LRC bispectrum or physical-current result.

### THM-3455's Berggren rank subsets

Let `E_t` indicate `rho_H(q_t)<=7` on the full parabolic spine.  THM-3455
proves period `1683` and density `76/187`; hence `(24)` sharpens its `O(1)`
statement to

```text
sum_(t<=T) E_t/t
  =(76/187)H_T+C_spine+O(1/T).                        (28)
```

For the Fibonacci-index sample `E'_n`, the proved period is `360` and the
density is `43/90`, so

```text
sum_(n<=N) E'_n/n
  =(43/90)H_N+C_Fib+O(1/N).                           (29)
```

Both constants have exact finite forms `(25)` at their displayed periods.
These are reciprocal sums over **indices**.  THM-3455's reciprocal sums over
the exponentially growing Fibonacci values and quadratic/exponential branch
labels converge for a different sparsity reason and are not instances of
`(24)` on those values.

## 8. K4, XOR, and tournament boundary

The eight transversals choosing one edge from each perfect matching of `K4`
split into four stars and four triangles by Boolean parity.  After declaring
any period-eight enumeration of those addresses, its star indicator has
harmonic coefficient `1/2` by `(24)`.  This is an address-word statement.

It does not orient the six edges, make the address order intrinsic, or turn a
static `K4` selector into a tournament process.  Likewise a periodic sequence
of tournaments gives periodic Boolean property indicators, but `(24)` sees
only their temporal frequency; it forgets vertices, arcs, ties, ancestry, and
every graph sidecar.

## 9. Exact companion and scope

Run

```bash
python -B 04-computation/periodic_polynomial_critical_harmonic_transform_thm3486.py
python -B -O 04-computation/periodic_polynomial_critical_harmonic_transform_thm3486.py
```

The companion checks exact top-layer packets, finite Abel identities, the
alternating-harmonic identity through 64 even horizons, THM-3484's common
leading layer, THM-3455's `1683/360` periods and `76/187,43/90` coefficients,
and the declared `4/4` K4 Boolean-address split.  These finite checks support
but do not replace the analytic proof.

This theorem proves no statement about arbitrary subsets, nonperiodic
automatic words, tournament classification, Berggren ancestry transport,
LRC currents, bispectra, Jacobian maps, or LRC(14).

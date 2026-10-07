# Spatial sidecars recover atoms at one fixed moment degree

2026-10-05. **PROVED:** a uniform, stable inverse and a finite signed
measurement certificate for the localized shell kernel at every integer
degree at least seven. **CONDITIONAL:** applying that certificate to the
actual Collatz law requires independently truthful, simultaneous measurement
intervals. **FINITE-EXACT:** the synthetic controls stated below. No new
actual Collatz atom positivity or ROOT coverage is claimed.

Script: `04-computation/experiments/collatz_specialization_sidecar_20261005.py`.
Output: `05-knowledge/results/collatz_specialization_sidecar_20261005.out`.

## 1. Inheritance and the coordinate that was discarded

The prior localized-moment packages fix a source index and increase moment
degree. In particular, `collatz_localized_resolvent_floor_20261005.md` and
`collatz_refinement_energy_dual_20261005.md` separate exact selector bounds
from the provenance of the measured moments. Their canonical hostile is a
positive smoothed value at a source whose actual atom is zero.

THM-4255,
`01-canon/theorems/THM-4255-specialization-kernel-and-transverse-hasse-jet-repair.md`,
proves that evaluation on a displayed formal arc has a kernel, and explains
how transverse coordinates can repair a restricted-source loss. Its current
correction is essential: that algebra alone does **not** refute the external
p-adic-zeta draft, whose alleged specialization map was not exhibited.
We use the principle only after displaying the actual map below. No formal
arc or zeta arithmetic is being asserted for the Collatz moment operator.

The board is **fixed source / shell reflection / spatial address / moment
degree / signed inverse / measurement error**. The new operation keeps degree
fixed and varies the measurement center. Neighboring addresses are measured
coordinates of one law, not replacement Collatz inputs.

For a probability law `(p_j)_(j>=0)`, define

```text
h_d = 4*2^d/(1+2^d)^2,             d>=0,
H_k(m) = sum_(j>=0) p_j h_|m-j|^k.                         (1)
```

Thus `h_0=1`, `h_1=8/9`, `h_2=16/25`. The intended application has
`p_j=lambda(6j+3)`, the inherited probability law on odd multiples of three.
The following proofs use only nonnegative masses summing to one.

The pointwise coordinate `h_|m-j|` identifies the two addresses `m-d,m+d`
when both exist. At `m=1`, indices 0 and 2 already collide. Let

```text
t=2^(m-j),   a=4t/(1+t)^2,   b=8t/(1+2t)^2.
```

Here `a=h_|m-j|` and `b=h_|m+1-j|`. Direct expansion gives

```text
t = (8/b-4/a-2)/3.                                        (2)
```

So two neighboring **pointwise** coordinates retain the lost orientation.
This is not an assertion that two averaged measurements determine a law.
The law-recovery theorem below uses the full spatial field and then bounds
finite approximations to its inverse.

## 2. A bounded inverse on the half-line

For bounded sequences on the nonnegative integers, define

```text
(T_k f)(m)=sum_(j>=0,j!=m) h_|m-j|^k f(j).
```

The half-line boundary is retained; negative indices are omitted rather than
wrapped or reflected. Equation (1) is the exact linear equation

```text
H_k=(I+T_k)p.                                             (3)
```

Since `h_d<=4*2^(-d)` for `d>=3`,

```text
||T_k||_infinity <= eta_k
 =2[(8/9)^k+(16/25)^k+1/(2^k-1)].                         (4)
```

The final term sums the geometric majorant for distances at least three.
Every summand decreases with `k`, and exact rational evaluation gives

```text
eta_7 = 3635701142046192706/3707501605224609375 <1,
eta_8 <=27/32,
eta_12 <1/2.                                              (5)
```

Consequently, for every integer `k>=7`,

```text
p = sum_(i>=0) (-T_k)^i H_k,
||(I+T_k)^(-1)||_infinity <=1/(1-eta_k).                   (6)
```

The series converges in operator norm, so it reconstructs every bounded
solution, not only probability laws. In particular it recovers the entire
probability law without a target-positivity premise. Supremum-norm error in
the measured field is amplified by at most `1/(1-eta_k)`, at most `32/5`
for degree eight and less than two for degree twelve.

Seven is the least integer degree for this uniform off-diagonal-norm
criterion. For `m>=2`, just the four neighbors at distances one and two give

```text
2[(8/9)^6+(16/25)^6]
 =145832200896512/129746337890625 >1.                      (7)
```

Smaller degrees have still larger contributions. This is **not** a minimum
information degree, nor an obstruction to inversion by another norm or
method at smaller degrees. The positive Fourier-symbol approach can address
that different question; it is unnecessary for the finite rational bounds
proved here.

## 3. A finite, exact signed readout

Let `T=T_k` and let `T_R` retain only distances `1,...,R`, where `R>=0`.
The omitted positive operator satisfies

```text
||T-T_R||_infinity <= gamma_(k,R)
 =2*4^k*2^(-k(R+1))/(1-2^(-k)).                           (8)
```

For an inverse depth `N>=0`, define the row and two exactly computable bills

```text
c = e_m^T sum_(i=0)^N (-T_R)^i,
C = sum_j |c_j|,
r = (T_R^(N+1) 1)(m).                                    (9)
```

The row is supported on `j>=0`, `|j-m|<=NR`, so at most `2NR+1` spatial
addresses are needed. Positive row propagation computes `r` without
knowing the law. The uniform bounds are

```text
C <= sum_(i=0)^N eta_k^i <=1/(1-eta_k),
0<=r<=eta_k^(N+1).                                       (10)
```

Write `E=T-T_R` and `S_N=sum_(i=0)^N(-T_R)^i`. Polynomial multiplication
in the single operator `T_R` gives the exact identity

```text
S_N H_k = p-(-T_R)^(N+1)p+S_N E p.                       (11)
```

There is no commutation assumption between `E` and `T_R` in (11).
Since `0<=p_j<=1`, the last term at `m` has absolute value at most
`beta=C gamma_(k,R)`. The power term is nonnegative before its displayed
sign and is at most `r`.

Suppose **one common law** gives truthful intervals
`l_j<=H_k(j)<=u_j` at every nonzero coefficient address. Form `A_-` by
using `l_j` for positive `c_j` and `u_j` for negative `c_j`; form `A_+`
by reversing those choices. Then

```text
N odd:   A_- - beta       <= p_m <= A_+ + beta + r,
N even:  A_- - beta - r   <= p_m <= A_+ + beta.            (12)
```

Intersect with `[0,1]`. Thus an explicitly positive left endpoint is a
genuine atom floor, conditional only on the stated simultaneous interval
premise. Merely passing the interval-format checks does not authenticate
that premise; mutually incompatible or fabricated oracle values do not
become valid measurements.

If every measurement interval has width at most `omega`, the resulting
atom interval has width at most

```text
C omega + 2 C gamma_(k,R) + eta_k^(N+1).                  (13)
```

This separates measurement error, finite spatial radius, and inverse-series
depth. It is uniform in the law and in the selected source.

## 4. A terminating precision compiler and retained floors

Given exact `0<epsilon<=1` and integer `k>=7`, the script chooses an odd
`N>=1` and `R>=0` such that

```text
eta_k^(N+1) <=epsilon/3,
2 gamma_(k,R)/(1-eta_k) <=epsilon/3,
omega=(1-eta_k)epsilon/3.                                 (14)
```

Successive exact rational multiplication and integer powers suffice;
no numerical logarithm is used. Both searches terminate because their
bases are strictly below one. Equations (10)--(14) guarantee an atom
interval of width at most `epsilon` from the requested finite measurement
packet. The compiler does not produce that packet.

For `epsilon=1/256`, its conservative parameters are:

| Fixed degree | Spatial cutoff R | Inverse depth N | Address radius NR |
|---|---:|---:|---:|
| 7 | 4 | 339 | 1356 |
| 8 | 3 | 39 | 117 |
| 12 | 3 | 9 | 27 |

This quantifies why degree seven, though sufficient, can be expensive.
Increasing degree is optional; it is a conditioning and address-count
tradeoff, not a logical requirement that degree diverge.

If truthful packets can be obtained to the prescribed accuracies, the
algorithm eventually certifies every actually positive `p_m`: any interval
containing it with width less than `p_m` has a positive lower endpoint.
It gives no proof that an arbitrary `p_m` is positive and no independent
method for producing the actual measurement intervals.

At each refinement retain the largest previous lower bound and smallest
previous upper bound for the **same index m**. Their intersection remains
valid and every established floor survives. The raw newly computed bound
need not improve by itself; changing the source index would also change
the proposition and is rejected by the retained-bound API.

## 5. Positive structure and two decisive hostiles

For the synthetic law

```text
p_0=1/16,             p_2=15/16,
```

the single scalar `H_8(0)` also comes from a zero-target law on indices
one and three. Indeed its mass at one can be chosen as

```text
w=(H_8(0)-h_3^8)/(h_1^8-h_3^8),       0<w<1.              (15)
```

This is an exact same-measurement, missing-atom hostile for the one-scalar
readout. A spatial packet separates the laws. For `(m,k,R,N)=(0,8,3,31)`,
the 94 exact neighboring measurements give a lower bound exceeding
`65437/1048576`, while the true atom is `1/16`. Measurement widths of
`1/100000` still give a lower bound exceeding `32713/524288`. These are
synthetic, independently specified laws; no Collatz target receipt enters.

Conversely, for `p=delta_0`,

```text
H_k(m)=h_m^k>0 for every finite m,
but p_m=0 for every m>0.                                  (16)
```

The inverse in (6) recovers those zero atoms exactly. The finite control
at `(m,k,R,N)=(3,12,3,15)` gives lower bound zero and upper bound less than
`1/1000`. Thus even an everywhere-positive measured field does not supply
positive source weights everywhere; the signed cancellation matters.

The transfer is now explicit:

| Source | Map | Preserved predicate | Lost information or premise |
|---|---|---|---|
| One point index j | `j -> h_|m-j|` | Shell distance from m | Left/right orientation; repaired pointwise by (2) |
| Probability law p | `p -> H_k`, all spatial centers | Entire law, through (6) | Measurement provenance is not produced by inversion |
| Full field | Finite signed row (9) | Interval for the original p_m | Spatial and series tails require the bills in (12) |
| Successive truthful packets | Intersect atom intervals | Every earlier source floor | Does not allow replacing the source or assuming a floor |

## 6. Exact finite universe and reproduction

The standard-library script uses `Fraction` throughout and explicit checks
that remain active under `-O`. It contains no ROOT encoder, orbit iterator,
completed weight bank, or actual moment oracle. Its universe comprises:

* contraction checks for degrees 7 through 16 and the exact degree-six
  hostile;
* precision plans for degrees 7, 8, 12 at tolerances `2^-1,2^-4,2^-8,2^-16`;
* the pointwise sidecar for `m=0..6`, `j=0..12`;
* ten singleton laws and four specified two-point laws, five source indices,
  three degrees and four radius/depth pairs: **840 exact atom brackets**;
* the positive, noisy, same-scalar ghost and everywhere-positive-field
  controls, plus retained-bound and malformed-input controls.

The synthetic positive readout has large exact rational coefficients. The
output prints rigorous 20-bit dyadic enclosures and numerator/denominator
bit sizes instead of many thousands of digits. This is an exact display
compression, not a floating-point substitute for the checks.

The first implementation exposed a cached-type alias: Boolean source keys
could reuse integer cache entries before domain validation. The cache now
retains argument types, and explicit Boolean, floating, negative and
malformed-packet controls cover the repaired boundary.

Reproduce from the repository root:

```text
python -B 04-computation/experiments/collatz_specialization_sidecar_20261005.py
python -B -O 04-computation/experiments/collatz_specialization_sidecar_20261005.py
```

Both modes perform **3,519 exact checks** and agree with the saved output.

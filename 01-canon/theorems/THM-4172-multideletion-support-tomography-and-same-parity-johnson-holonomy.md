---
id: THM-4172
title: "Multideletion support tomography and same-parity Johnson holonomy"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED. The complete deletion
  tower of a fixed exposure-capacity coordinate is the binomial transform
  of its tagged OCF outside-support histogram, and exact inversion recovers
  that histogram. More generally, deletion acts diagonally on the
  vertex-support grading of every edge polynomial. Since the Johnson
  numerator and denominator have support degrees three and four, every
  deletion depth obeys an exact normalized transport law. Even deletion
  depth preserves parity and makes the parent tilt a genuine positive
  barycenter of the corrected restrictions. FINITE-EXACT complete checks
  certify all 12,155 symmetric prime order-eleven classes and all 1,002
  prime attachments over the frozen asymmetric parent; a 100,000-row stream
  is a deterministic VERIFIED SAMPLE, not a census.
source: codex-frontier-synthesis-creative-20260826au
depends_on:
  - THM-002-ocf
  - THM-4128-johnson-slice-support-envelope-and-exposure-centrality-criterion
  - THM-4167-tournament-exposure-capacity-deletion-support-moment-and-parity-holonomy
related:
  - THM-3372-multiaffine-deletion-transform-variance-and-skew-join-current
  - THM-4133-strong-cyclic-substitution-johnson-centrality-counterexample
  - THM-4135-strong-tournament-centrality-complete-order-nine
  - THM-4168-prime-order-eleven-nontrivial-automorphism-johnson-centrality
  - THM-4169-prime-parent-one-vertex-augmentation-and-quartic-johnson-transfer
  - THM-4177-root-split-johnson-current-and-source-sink-boundary-census
script: 04-computation/tournament_multideletion_support_tomography_thm4172.py
output: 05-knowledge/results/tournament_multideletion_support_tomography_thm4172.out
independent_audit_script: 04-computation/tournament_multideletion_support_tomography_thm4172_independent_audit.cpp
independent_audit_output: 05-knowledge/results/tournament_multideletion_support_tomography_thm4172_independent_audit.out
finite_certificate_source: 04-computation/tournament_order11_two_deletion_certificate_thm4172.cpp
finite_certificate_outputs:
  - 05-knowledge/results/tournament_order11_two_deletion_certificate_thm4172_symmetric.out
  - 05-knowledge/results/tournament_order11_two_deletion_certificate_thm4172_cube.out
  - 05-knowledge/results/tournament_order11_two_deletion_certificate_thm4172_random100k.out
  - 05-knowledge/results/tournament_order11_two_deletion_certificate_thm4172_exhaustive.out
  - 05-knowledge/results/tournament_order11_two_deletion_certificate_thm4172_hostiles.out
script_sha256: 720aee79d95eaf1ab0049dcf3d7bdd099592be2786716f52e6c399994b65db4d
output_sha256: 8090ebd4969ffbb2ff822623fc0d98e498e3a46426fd71be50e90e097bcb1f80
independent_audit_script_sha256: 39ade4819410c0226aa9a27d40c4290eecffa1fc6255cf638739d194b03811aa
independent_audit_output_sha256: 9112310d8a42cda65e5a1dc7862184584fdbb9013074043ddd83e102b0f3704e
finite_certificate_source_sha256: c2113784df1dfabef739f50bd93eae6d2d9aec9f1fb9cdf96dd4d90ff975f350
finite_symmetric_output_sha256: 1d4bd17481c6babdba112fec5b6a6f45a1186621b7d967f3b867690304319125
finite_cube_output_sha256: b502f90f54b431e54304ad9a11254ac4e95463811fbc40bc75809c623c62fa4f
finite_random100k_output_sha256: f6d9ecea2b148009f169427c6e5cf72d20b23a71b3062e1319490bc3a058feb6
finite_exhaustive_output_sha256: a6a7a15aa0028c0faed3ea1f90f584bb68a1f9dae857d43f9cb8a78b38cf97f9
finite_hostiles_output_sha256: 9d4b14eaefc49da6fe960e84e3d030b037f9b943a794dfa0403d6e6bc5320f96
hash_basis: raw LF bytes
primary_audit: >
  PASS. Literal tagged odd-cycle packings check 410 capacity coordinates and
  1,202 deletion-layer identities through order four. Direct exact
  polynomial evaluation checks 15 all-depth rows and 632 restrictions
  through order nine, then evaluates the named order-eleven and order-twelve
  controls. Normal, optimized, and hash-seeded streams byte-match.
independent_audit: >
  ACCEPT. A clean-room C++17 referee derives capacities from exposed words
  through order five (10,650 coordinates and 42,162 layers), independently
  evaluates response capacities by ear-state Hamilton DPs for the named
  controls, and repeats the all-depth polynomial algebra. Clang O0/O3 and
  GCC O3 streams byte-match under warnings-as-errors builds.
finite_certificate_audit: >
  PASS. A separate exact C++ evaluator checks all 668,525 corrected
  restrictions over the complete THM-4168 symmetric-prime bank, all 55,110
  restrictions over the 1,002 prime children in the frozen THM-4169 cube,
  complete labelled orders six and seven, a deterministic 100,000-row
  order-eleven sample, and corrected-versus-actual hostile controls. Fresh
  O3 outputs byte-match all frozen streams; O0 also byte-matches the
  symmetry, cube, and hostile streams.
---

# THM-4172 -- multideletion support tomography and same-parity Johnson holonomy

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

THM-4167 found the first outside-support moment of the positive tagged OCF
atoms and the one-deletion parity law. The full deletion tower is stronger:
it recovers the entire outside-support histogram, while arbitrary polynomial
restriction diagonalizes by vertex-support size. Applied to the quadratic
Johnson packet, this gives an exact transport law at every depth. Two
deletions are especially useful because the parity factor disappears.

## 1. Capacity deletion tomography

Let `T` be a tournament on `V`, `|V|=n`, fix an edge `e={i,j}`, and retain
THM-4167's tagged positive atom family `P_e(T)`. Put `N=n-2` and

```text
W_(e,s)(T)
 =sum_(Gamma in P_e(T),
       |U(Gamma) intersect (V minus e)|=s) 2^|Gamma|,
                                                    0<=s<=N.       (1)
```

Thus `W_(e,s)>=0` and `c_e(T)=sum_s W_(e,s)`. For `0<=r<=N`, define the
labelled symmetric deletion layer

```text
A_(e,r)(T)=sum_(R subseteq V minus e, |R|=r)c_e(T-R).    (2)
```

> **Theorem 1 (binomial deletion transform).** For every `r`,
>
> ```text
> boxed: A_(e,r)=sum_(s=0)^N binom(N-s,r)W_(e,s).        (3)
> ```
>
> Equivalently,
>
> ```text
> A_e(t):=sum_(r=0)^N A_(e,r)t^r
>        =sum_(s=0)^N W_(e,s)(1+t)^(N-s).                (4)
> ```

Consequently the complete layer vector recovers the outside-support
histogram exactly:

```text
boxed:
W_(e,s)
 =sum_(r=N-s)^N (-1)^(r-N+s) binom(r,N-s)A_(e,r).       (5)
```

### Proof

By THM-4167's tagged deletion correspondence, an atom `Gamma` contributes to
`c_e(T-R)` exactly when `R` avoids its outside support. An atom with outside
support size `s` survives exactly `binom(N-s,r)` deletion sets of size `r`.
Double-counting `(R,Gamma)` proves `(3)`, and summing in `r` proves `(4)`.
After the substitution `x=1+t`, the coefficient of `x^(N-s)` in
`A_e(x-1)` is `W_(e,s)`, which is `(5)`. **QED.**

THM-4167's first moment is the `r=1` shadow. Indeed, with

```text
Omega_e=sum_s s W_(e,s),
```

equation `(3)` gives

```text
A_(e,1)=N c_e-Omega_e.                                  (6)
```

The alternating expressions in `(5)` are therefore nonnegative integers.
They recover support **sizes**, not the individual tagged atoms, their
supports, or overlaps between the atom families of different edges.

## 2. Deletion diagonalizes vertex-support degree

Let `P(z)` be any polynomial in symmetric edge variables `z_ab` on `V`.
The vertex support of a monomial is the set of endpoints incident to one of
its variables, regardless of exponent. Write

```text
P=sum_(s=0)^n P^[s]                                    (7)
```

for its decomposition by exact vertex-support size. Restriction to `V-R`
sets every variable incident to `R` equal to zero. Define

```text
Del_r P(z)=sum_(R subseteq V, |R|=r)P(z|_(V-R)).        (8)
```

> **Theorem 2 (support-degree diagonalization).** For every `r`,
>
> ```text
> boxed: Del_r P=sum_(s=0)^n binom(n-s,r)P^[s].         (9)
> ```

Each monomial of support `s` survives exactly the deletion sets disjoint
from its support, proving `(9)` term by term. Thus the deletion tower is the
binomial spectral transform of the support grading:

```text
sum_r (Del_r P)t^r=sum_s P^[s](1+t)^(n-s),              (10)
```

and coefficient extraction after `t=x-1` recovers every `P^[s]`. This is
the polynomial analogue of `(3)--(5)`; neither statement requires a
tournament quotient or an isomorphism-class census.

## 3. The Johnson packet at every deletion depth

Retain THM-4167's packet for an arbitrary symmetric edge tensor `z` on the
fixed oriented complete graph:

```text
d_i=sum_(j!=i)z_ij,
h_i=sum_(i->j)z_ij-sum_(j->i)z_ij,
C(z)=sum_i h_i d_i,
D(z)=sum_(e<f, e intersect f=empty)z_ez_f.              (11)
```

The cancellation in THM-4167 gives

```text
C(z)=2(sum_(e<f with common tail)z_ez_f
       -sum_(e<f with common head)z_ez_f).              (12)
```

Every surviving monomial of `C` has vertex support three, while every
monomial of `D` has support four. Theorem 2 therefore gives, for
`0<=r<=n-4`,

```text
boxed:
sum_(|R|=r) C(z|_(V-R))=binom(n-3,r)C(z),
sum_(|R|=r) D(z|_(V-R))=binom(n-4,r)D(z).               (13)
```

This is an all-depth theorem for arbitrary tensors; no capacity positivity
is used in `(12)--(13)`.

## 4. Multideletion parity holonomy

Let `m=n-r>=4`, assume `D(z)>0`, and assume every restricted denominator

```text
D_R:=D(z|_(V-R))>0,              |R|=r.                 (14)
```

Put

```text
kappa_k=2 if k is even,          kappa_k=4 if k is odd,
tau_k(z)=(k-3)C(z)/(kappa_k D(z)),
lambda_R=D_R/sum_(|S|=r)D_S.                            (15)
```

Then the `lambda_R` are positive and sum to one. Equations `(13)` and

```text
(m-3) binom(n-3,r)/binom(n-4,r)=n-3                   (16)
```

give the exact law

```text
boxed:
tau_n(z)=(kappa_m/kappa_n)
             sum_(|R|=r)lambda_R tau_m(z|_(V-R)).       (17)
```

The factor has the complete parity table

| deletion depth `r` | parent parity | factor `kappa_m/kappa_n` |
|:---|:---|:---:|
| even | either | `1` |
| odd | odd | `1/2` |
| odd | even | `2` |

Thus every **even** deletion depth is a genuine positive barycenter:

```text
boxed:
tau_n(z)=sum_(|R|=r)lambda_R tau_(n-r)(z|_(V-R)),
                                                     r even.       (18)
```

In particular the parent tilt lies in the convex hull of the corrected
restriction tilts. If all of them satisfy `|tau_(n-r)|<1`, then so does the
parent. This is a sufficient local certificate; parent centrality does not
imply that every child is central.

## 5. Exact order-eleven two-deletion certificate

For a tournament capacity tensor `z=c(T)`, all hypotheses in `(14)` hold.
Iterated capacity monotonicity from THM-4167 gives

```text
c(T)|_(V-R)>=c(T-R) coordinatewise,                     (19)
```

and THM-4128 gives positive disjoint-edge energy to the actual card tensor
when `|V-R|>=4`. Write

```text
b^(R)=c(T)|_(V-R).                                      (20)
```

These are corrected restrictions of the **parent** tensor, not the actual
capacities `c(T-R)`.

At `11 -> 9` with `|R|=2`, equations `(13)` and `(18)` become

```text
lambda_R=D(b^(R))/(21D(c(T))),
boxed: tau_11(c(T))=sum_(|R|=2)lambda_R tau_9(b^(R)).   (21)
```

Hence the exact strict centrality criterion is

```text
|sum_(|R|=2)lambda_R tau_9(b^(R))|<1.                  (22)
```

Requiring all 55 corrected pair restrictions to be central is stronger than
`(22)` but sufficient. This reframes the open asymmetric order-eleven bank
as a same-parity order-nine restriction problem. It does **not** inherit
THM-4135 automatically, because `b^(R)` need not equal the capacity tensor
of the actual order-nine card.

## 6. Aggregate quartic gates and exact evidence

At order eleven, `(13)` specializes further to

```text
sum_(|R|=2) C_R=28C,                 sum_(|R|=2) D_R=21D,

boxed:
sum_(|R|=2)(2D_R+3C_R)=42(D+2C),
sum_(|R|=2)(2D_R-3C_R)=42(D-2C).                    (22a)
```

Thus the parent strict rational gate is exactly positivity of both aggregate
signed child walls. Two convenient stronger sufficient loads are

```text
L_local=max_R 3|C_R|/(2D_R),
L_abs=3 sum_R |C_R|/(2 sum_R D_R).                    (22b)
```

Either `L_local<1` or the weaker `L_abs<1` implies `(22)`. The first demand
checks every corrected order-nine restriction; the second permits signed
cancellation only after paying the triangle-inequality tax.

The separate exact evaluator gives the following finite banks.

| universe | target rows | corrected restrictions | local/absolute failures | sharp local load | sharp absolute load |
|:---|---:|---:|---:|:---|:---|
| complete THM-4168 symmetric primes | `12,155` | `668,525` | `0 / 0` | `92494474/97222707` | `201109364/381606521` |
| prime children in the frozen THM-4169 cube | `1,002` | `55,110` | `0 / 0` | `1038483552/1131093929` | `571044751/1214691687` |

Both are **FINITE-EXACT** complete checks of the displayed universes. The two
uniform nonstrong cube attachments are the only local-certificate failures
in the full 1,024-pattern cube; they are not prime.

A deterministic SplitMix64 stream of `100,000` labelled order-eleven rows
contains `88,346` prime and `97,900` strong rows. No prime or strong row fails
the parent, local, or absolute gate. All `2,100` local failures in the whole
stream are nonstrong. This is a **VERIFIED SAMPLE**, not an isomorphism-class
census and not evidence by itself for the remaining universal claim.

Complete labelled small-order controls also have zero strong local failures:

```text
n=6: 22,320 strong rows,    334,800 restrictions, max load 70/87;
n=7: 1,677,488 strong rows, 35,227,248 restrictions, max load 570/1231.
                                                               (22c)
```

These computations make all 55 walls a plausible next certificate for the
asymmetric bank, but they do not prove that every strong order-eleven tensor
satisfies them.

## 7. Sharp controls

The prime order-eleven code from THM-4167,

```text
3169369058263173,
```

has

```text
tau_11=1055017002/11090656697.
```

All 55 corrected two-deletion restrictions are strictly central. The largest
absolute child tilt occurs at `R={7,10}` and equals

```text
1629373665/9484374388<1.                                (23)
```

Yet its **actual** card at vertex `10` is nonstrong and fails the order-ten
rational gate, as recorded in THM-4167. Thus actual cards cannot replace
the corrected restrictions in `(21)`.

The THM-4133 order-twelve hostile has

```text
tau_12=-53092739331/40435524866<-1.                     (24)
```

All 12 one-deletion corrected restrictions are individually central, so the
odd-depth factor `2` produces a false negative if it is forgotten. At depth
two, where parity is restored, exactly 65 of the 66 corrected restrictions
are noncentral. The unique central pair is

```text
R={9,11},    tau_10=-4461429/4143186457,                (25)
```

while the largest absolute tilt occurs at

```text
R={9,10},    tau_10=-770665315/458646486.               (26)
```

So same-parity deletion exposes the known failure almost everywhere. The
one-deletion parity factor and the corrected-versus-actual distinction are
both load-bearing.

There is also a corrected-versus-actual hostile on the THM-4133 row. For
`R={1,2}`, the actual order-ten card is strong and central:

```text
(C,D)=(-34428472,168571592),
7|C|/(2D)=4303559/6020414<1.                           (26a)
```

The corrected parent restriction is noncentral:

```text
(C_R,D_R)=(-51221650864,152533479120),
7|C_R|/(2D_R)=22409472253/19066684890>1.              (26b)
```

Its coordinatewise deletion defect is nonzero on all 45 edges, with total
`705,854` and maximum coordinate `35,962`. Thus even **strong central actual
cards** cannot substitute for corrected restrictions; the full edgewise
deletion defect is the missing sidecar.

## 8. Boundaries and replay

- The normalized law stops when fewer than four vertices remain because the
  disjoint-edge denominator vanishes.
- Denominator positivity is required for arbitrary tensors. The polynomial
  identities `(9)` and `(13)` do not require it.
- Capacity tomography reconstructs only each edge's marginal support-size
  histogram. Cross-edge atom correlations, which can affect `C` and `D`, are
  additional information.
- The strict wall is `|tau|<1`. Equality retains THM-4128's distinction
  between a central optimizer existing and every optimizer being central.
- Neither `(21)` nor the positive finite banks prove the full asymmetric
  order-eleven bank.

Run the primary exact audit with

```text
python3 -B 04-computation/tournament_multideletion_support_tomography_thm4172.py
python3 -B -O 04-computation/tournament_multideletion_support_tomography_thm4172.py
PYTHONHASHSEED=271828 python3 -B \
  04-computation/tournament_multideletion_support_tomography_thm4172.py
```

Compile the independent exposed-word/ear-response referee with

```text
clang++ -O3 -std=c++17 -Wall -Wextra -Werror \
  04-computation/tournament_multideletion_support_tomography_thm4172_independent_audit.cpp \
  -o /tmp/thm4172_independent_audit
/tmp/thm4172_independent_audit
```

Replay the finite banks with

```text
clang++ -O3 -std=c++17 -Wall -Wextra -Werror \
  04-computation/tournament_order11_two_deletion_certificate_thm4172.cpp \
  -o /tmp/thm4172_certificate

git show origin/main:05-knowledge/results/tournament_prime_nontrivial_automorphism_order11_thm4168.labels \
  | /tmp/thm4172_certificate --stdin
/tmp/thm4172_certificate --cube
/tmp/thm4172_certificate --random 11 100000
/tmp/thm4172_certificate --exhaust 6
/tmp/thm4172_certificate --exhaust 7
/tmp/thm4172_certificate --thm4133
/tmp/thm4172_certificate --root-hostile
```

Clang `-O0/-O3` and GCC `-O3` produce the same frozen stream. The finite
audits are hostile controls; the deletion formulas themselves are proved at
all orders. **QED.**

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
  barycenter of the corrected restrictions.
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
script: 04-computation/tournament_multideletion_support_tomography_thm4172.py
output: 05-knowledge/results/tournament_multideletion_support_tomography_thm4172.out
independent_audit_script: 04-computation/tournament_multideletion_support_tomography_thm4172_independent_audit.cpp
independent_audit_output: 05-knowledge/results/tournament_multideletion_support_tomography_thm4172_independent_audit.out
script_sha256: 720aee79d95eaf1ab0049dcf3d7bdd099592be2786716f52e6c399994b65db4d
output_sha256: 8090ebd4969ffbb2ff822623fc0d98e498e3a46426fd71be50e90e097bcb1f80
independent_audit_script_sha256: 39ade4819410c0226aa9a27d40c4290eecffa1fc6255cf638739d194b03811aa
independent_audit_output_sha256: 9112310d8a42cda65e5a1dc7862184584fdbb9013074043ddd83e102b0f3704e
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

## 6. Sharp controls

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

## 7. Boundaries and replay

- The normalized law stops when fewer than four vertices remain because the
  disjoint-edge denominator vanishes.
- Denominator positivity is required for arbitrary tensors. The polynomial
  identities `(9)` and `(13)` do not require it.
- Capacity tomography reconstructs only each edge's marginal support-size
  histogram. Cross-edge atom correlations, which can affect `C` and `D`, are
  additional information.
- The strict wall is `|tau|<1`. Equality retains THM-4128's distinction
  between a central optimizer existing and every optimizer being central.
- Neither `(21)` nor the positive named control proves the full asymmetric
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

Clang `-O0/-O3` and GCC `-O3` produce the same frozen stream. The finite
audits are hostile controls; the deletion formulas themselves are proved at
all orders. **QED.**

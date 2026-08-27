---
id: THM-4263
title: "Moving multigraph filtered jets and finite-factor density transport"
status: >
  PROVED + VERIFIED-EXACT. Moving multigraph evaluation has the graph ideal
  as kernel; total and rectangular Hasse jets recover precisely the expected
  powers of that ideal. Compatible truncated sections form a split filtered
  inverse system, while all-level short-jet injectivity is strictly stronger
  than limit injectivity. For sequences of finite factors, normalized fibre
  weights are uniformly integrable iff every uniform-target density-one
  compiler transports to uniform-source probability one. Migrated reset
  channels obey the exact conditional-hazard product, with no independence
  hypothesis. The Rule-30 consequence uses only THM-4204's fibre cap and does
  not solve a Rule-30 prize. MISTAKE-527 remains controlling: the chosen-section
  hostiles here do not refute an external universal-torsor identity.
source: root/cross-frontier-overnight/2026-08-27
depends_on:
  - THM-4204-rule30-debruijn-reset-and-dyadic-prefix-saturation
  - THM-4255-specialization-kernel-and-transverse-hasse-jet-repair
  - THM-4257-fixed-prime-exponent-orbit-density-one-compiler
mistake_firewall:
  - MISTAKE-527
script: 04-computation/graph_factor_density_transport_audit_thm4263.py
output: 05-knowledge/results/graph_factor_density_transport_audit_thm4263.out
script_sha256: 50272a096333874983000e18638a90dccd8747924896efa41e8c2eebff443020
output_sha256: a08aba5cc977bd7382649bad57b72fdf1e4bab7e2f60338aa1e3e9c10d204243
hash_basis: raw bytes
audit: >
  PASS. A dependency-free standard-library companion checks a moving
  two-coordinate graph, rectangular Hasse recovery, compatible and hostile
  filtered sections, all subsets of a regular finite linear factor, fresh
  and merged conditional hazards, a non-uniformly-integrable heavy fibre,
  periodic Rule-30 fibres through N=12, and the Cartier-channel hostile.
  Normal, optimized, and fixed-hash-seed streams are byte-identical.
---

# THM-4263 -- moving multigraph filtered jets and finite-factor density transport

**PROVED + VERIFIED-EXACT.**

## 1. Moving multigraph kernels

Let `A` be a nonzero unital commutative ring and put

```text
R=A[u_1,...,u_d][[f]],
g=(g_1,...,g_d) in A[[f]]^d,
ev_g(P)=P(f,g_1(f),...,g_d(f)),
I_g=(u_1-g_1,...,u_d-g_d)R.                            (1)
```

Then for every `N>=0`,

```text
ker(ev_g)=I_g,
ev_g^(-1)(f^N A[[f]])=I_g+f^N R.                      (2)
```

For a multi-index `alpha`, let `D_u^[alpha]` be the multivariate Hasse
derivative. The locally finite Hasse--Taylor formula is

```text
P=sum_alpha D_u^[alpha]P(f,g)(u-g)^alpha.              (3)
```

Consequently the total jet and rectangular jet kernels are

```text
ker((D_u^[alpha]P(f,g))_(|alpha|<=r))=I_g^(r+1),       (4)

ker((D_u^[alpha]P(f,g))_(0<=alpha_i<=r_i))
 =((u_1-g_1)^(r_1+1),...,(u_d-g_d)^(r_d+1)).           (5)
```

On the coordinate-degree box `deg_(u_i)P<=d_i`, the full rectangular jet
`0<=alpha_i<=d_i` is an isomorphism onto

```text
A[[f]]^(product_i(d_i+1)).                             (6)
```

### Proof

For a polynomial in the `u_i`, replace the variables by `g_i` one at a time.
The resulting telescoping difference is a sum of
`(u_i-g_i)q_i(u,g)`. Apply this coefficientwise in `f`. Regrouping is valid:
only the finitely many original coefficients of order at most `M` can affect
the coefficient of `f^M`, and each is a polynomial. This proves the first
identity in `(2)`; the second follows by subtracting the evaluated tail.

Formula `(3)` follows coefficientwise. Monomials of total normal degree at
least `r+1` are exactly `I_g^(r+1)`, proving `(4)`. A normal monomial lies
outside the rectangular box exactly when it is divisible by a generator in
`(5)`. On the bounded degree box, `(3)` terminates and its triangular inverse
is `b_alpha -> sum b_alpha(u-g)^alpha`, proving `(6)`.

Thus several torsor coordinates lose a graph **ideal**, not a cosmetically
chosen single element such as `u-f`.

## 2. Filtered sections and the short-jet criterion

Put `R_N=R/f^N R`. Truncated sections

```text
g^(N) in (A[[f]]/f^N)^d
```

come from one formal section iff

```text
g^(N+1) mod f^N=g^(N)                                  (7)
```

for every `N`. Under `(7)`, the graph quotients form a split inverse system:

```text
0 -> I_(g,N) -> R_N -> A[[f]]/f^N -> 0,
lim_N R_N/I_(g,N)=A[[f]],
lim_N I_(g,N)=I_g.                                     (8)
```

For an `A`-submodule `V<=R`, the finite short-jet map

```text
J_(V,N):V/(V intersect f^N R) -> A[[f]]/f^N            (9)
```

is injective iff

```text
V intersect (I_g+f^N R)=V intersect f^N R.            (10)
```

Indeed its kernel is the quotient of the left side of `(10)` by the right
side. All-level filtered order reflection requires `(10)` for every `N`.
This is strictly stronger than limit injectivity, which asks only

```text
V intersect I_g=0.                                    (11)
```

The boundary is sharp. Over `F_2`, under `ev(u)=0`, the span of `u+f` maps
injectively in the formal limit because its image is `f`, but its level-one
target jet is zero although `u+f` is not in `fR`. Also, choices such as
`g^(1)=0` and `g^(2)=1 mod f^2` violate `(7)` and define no filtered graph
quotient.

## 3. Exact finite-factor density law

Let `q:X->Y` map finite nonempty sets, let `mu_X,mu_Y` be uniform, and define

```text
w(y)=|q^(-1)(y)|,
rho(y)=|Y|w(y)/|X|,
nu=q_*mu_X=rho mu_Y.                                  (12)
```

Zeros in `w` are allowed. A source predicate is the pullback of a target
predicate iff it is a union of fibres. For a group homomorphism this is the
congruence condition `G+ker(q)=G`. For every `H subset Y`,

```text
mu_X(q^(-1)H)=nu(H)=|Y|^(-1)sum_(y in H)rho(y).        (13)
```

This equals uniform target density for every `H` iff the fibres are constant.
In particular, a finite group homomorphism restricted to its image has
regular fibres.

Now let `q_n:X_n->Y_n` be a sequence. The following are equivalent:

1. for every `B_n subset Y_n`,
   `mu_Yn(B_n)->0` implies `q_(n* )mu_Xn(B_n)->0`;
2. the normalized fibre weights are uniformly integrable:

```text
lim_(M->infinity) sup_n |Y_n|^(-1)
    sum_(rho_n(y)>M)rho_n(y)=0.                        (14)
```

For sufficiency, split `(13)` at `rho_n<=M`; the light part is at most
`M mu_Yn(B_n)` and the heavy tail is uniformly small. Conversely, failure of
`(14)` supplies a subsequence and heavy-tail sets
`B_(n_j)={rho_(n_j)>M_j}` with pushforward mass bounded below. Since every
`rho_n` has uniform mean one, Markov gives
`mu_Yn_j(B_(n_j))<=1/M_j->0`, a contradiction. A uniform fibre bound is a
convenient sufficient condition, not the exact boundary.

The sharp hostile gives one point of an `n`-point target a fibre of weight
`2^n` and every other point weight one. Its uniform bad density is `1/n`, but
its source mass is `2^n/(2^n+n-1)->1`.

## 4. Migrated channels and exact reset hazards

On `(Y_n,nu_n)`, let `E_(n,1),...,E_(n,m_n)` be migrated reset events and set

```text
A_(n,0)=Y_n,
A_(n,j)=A_(n,j-1)\E_(n,j),
p_(n,j)=nu_n(E_(n,j) | A_(n,j-1)).                    (15)
```

If the preceding survivor has zero mass, set `p=1`. By definition,

```text
nu_n(A_(n,m_n))=product_j(1-p_(n,j)).                 (16)
```

Therefore the compiler has density one iff

```text
sum_j -log(1-p_(n,j))->infinity,                      (17)
```

where `-log(0)=infinity`. No independence assumption occurs. Repeating one
event of mass `1/8` four times illustrates the quotient loss: the exact
hazards are `(1/8,0,0,0)` and the defect stays `7/8`, not `(7/8)^4`.

A source certificate is exactly the pullback compiler generated by these
events iff

```text
C_n=q_n^(-1)(union_j E_(n,j)).                         (18)
```

Together, `(14)`, `(17)`, and `(18)` give the precise firewall:

```text
fibre saturation / predicate congruence
 + correct pushforward weights (or uniform integrability for all compilers)
 + divergent conditional hazard.                     (19)
```

## 5. Cross-frontier consequences

### Graph specialization and Cartier channels

After finite coefficient and `f`-adic truncation over a finite ring, graph
evaluation is linear and has regular fibres on its image. Its kernel alone
does not bias a fibre-saturated density. The danger in THM-4255 is instead
that a source Cartier channel need not be constant on graph fibres: under
`u->f`, `u` and `f` have the same value while the explicit-`f` channel
distinguishes them. The first `u`-Hasse jet restores rank two on their span.

MISTAKE-527 is controlling. The abstract pivots `38=23+15` and `46=29+17`
remain chosen-section hostiles outside the respective short windows. They are
not counterexamples to the external Proposition 6.2/12.3 universal-torsor
arguments unless an actual chosen-section map is exhibited there.

### THM-4257 factorial cylinders

The exponent-class map to the fixed odd-prime orbit modulo `2^ell` is a
bijection, hence `rho=1`. THM-4257's complete disjoint all-one blocks have
exact fresh hazard `2^-L`; `(16)` recovers its defect bound

```text
(1-2^(-L))^floor((ell-c_P)/L).                         (20)
```

If a later quotient merges blocks, the conditional hazards replace the old
marginals.

### Rule 30

THM-4204 proves that every cyclic Rule-30 output fibre has size at most three.
Both finite spaces have size `2^N`, so `rho_N=w_N<=3` and `(14)` holds.
Consequently every uniform-output density-one language pulls back to a
uniform-input probability-one language. Exact finite densities can differ:
at `N=3` the output-fibre histogram for weights `0,1,2,3` is `(3,3,1,1)`,
and the all-one output has three predecessors.

This is a lawful density bridge only. It proves no named center-orbit theorem
and solves none of the three Rule-30 prizes.

## 6. Audit and scope

The standard-library companion verifies the finite models and all stated
hostiles. Reproduction is

```bash
python3 -B 04-computation/graph_factor_density_transport_audit_thm4263.py
python3 -B -O 04-computation/graph_factor_density_transport_audit_thm4263.py
PYTHONHASHSEED=4263 python3 -B \
  04-computation/graph_factor_density_transport_audit_thm4263.py
```

The theorem supplies exact transport criteria. It does not turn a
non-fibre-saturated source channel into a target predicate, replace weighted
pushforward measure by uniform counting, prove an external density manuscript,
or solve a Rule-30 prize. **QED.**

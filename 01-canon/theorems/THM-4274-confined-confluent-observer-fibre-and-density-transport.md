---
id: THM-4274
title: "Confined confluent observer fibre and density transport"
status: >
  PROVED ABSTRACT LEMMA + VERIFIED FINITE CONTROLS. On a nonempty finite
  candidate universe, a complete first-stage observer followed by a linear
  quotient has exactly the original certificate fibres iff the quotient
  kernel misses the finite observed pair-difference shell. Monic recurrence
  prefixes and rectangular normal-Hasse jets supply two complete observers;
  cleared Bezout scalars, resultant torsion, fixed-basis modular separation,
  total observer capacity, and normalized fibre weights are mandatory
  sidecars. Uniform integrability is exactly equivalent to transporting every
  target-negligible family. This theorem proves no Jacobian realization,
  p-adic density claim, Rule-30 prize, or LRC(14) case.
source: codex-creative-frontiers-20260827
depends_on: []
related:
  - THM-4255-specialization-kernel-and-transverse-hasse-jet-repair
  - THM-4258-w0-three-sample-attachment-recurrence-and-two-torsion-sidecar
  - THM-4263-moving-multigraph-filtered-jet-and-finite-factor-density-transport
  - THM-4264-w0-visible-incidence-two-edge-attachment-observer
script: 04-computation/confined_confluent_observer_transport_audit_thm4274.py
output: 05-knowledge/results/confined_confluent_observer_transport_audit_thm4274.out
script_sha256: 6a9307c753edc86cb862b07f69decd26c28c1cd9e798c31200244ccc9fd12bc9
output_sha256: 4ba8a9ecd01b5d84f3d0bddb555657b7e7319eb0040aad00ef7e3dabd14b519f
hash_basis: raw bytes
audit: >
  PASS. A dependency-free standard-library companion checks 340 field-gcd
  recurrence pairs and all their common words, direct equality of 324
  candidate fibre pairs, 348 sharp integer-box cells, all 512 subsets of a
  bounded-fibre model, capacity and balanced-overload controls, and explicit
  empty-universe, nonmonic, integral-torsion, bad-prime-no-degeneration,
  specialization, unbounded-order, merged-hazard, and heavy-fibre hostiles.
  Normal, optimized, and fixed-hash-seed streams are byte-identical.
---

# THM-4274 -- confined confluent observer fibre and density transport

**PROVED ABSTRACT LEMMA + VERIFIED FINITE CONTROLS.**

## 1. Exact two-stage fibre criterion

Let `R` be a unital commutative ring, let `C` be a nonempty finite set, let
`Z,A,B` be `R`-modules, and let

```text
F:C -> Z,             O:Z -> A,             lambda:A -> B.          (1)
```

The certificate map `F` need not be linear. Assume that `O` and `lambda`
are linear and that `O` is injective on the finite set `F(C)`. Put

```text
W=O(F(C)),
Delta W=W-W={w-w':w,w' in W},
q=lambda O F.                                                        (2)
```

Then the following are equivalent:

```text
(i)   F(c)=F(c') iff q(c)=q(c') for all c,c' in C;
(ii)  lambda restricted to W is injective;
(iii) ker(lambda) intersect Delta W={0}.                            (3)
```

If `0 in F(C)`, there is also the exact zero-test criterion

```text
q(c)=0 iff F(c)=0 for every c in C
    iff ker(lambda) intersect W={0}.                                (4)
```

The nonempty hypothesis in `(3)` is real: it ensures `0 in Delta W`. If
`C` is empty, the correct third condition would be the containment
`ker(lambda) intersect Delta W subset {0}`, not equality.

### Proof

Because `O` is injective on `F(C)`, equality of full certificates is
equivalent to equality of their observations. Thus `(i)` is equivalent to
the assertion that `lambda` creates no collision inside `W`, which is `(ii)`.
For `w,w' in W`,

```text
lambda(w)=lambda(w') iff w-w' in ker(lambda) intersect Delta W.     (5)
```

Since `0 in Delta W`, this proves `(ii) iff (iii)`. If `0 in F(C)`, then
`0=O(0) in W`. Condition `(4)` follows by applying the same argument to the
distinguished observed certificate zero. This proves the fibre and zero-test
claims. `square`

The theorem is intentionally about full fibres. A zero test alone does not
control multiplicities and therefore cannot by itself transport measure.

## 2. Two complete first-stage observers

### 2.1 Monic recurrence prefixes

Let `M` be an `R`-module, let `L>=1`, let `S` be the `R`-linear cyclic shift
on `M^L`, and let `D` be a monic polynomial of degree `0<=r<=L`. For `r>=1`,
write

```text
D(T)=T^r+d_(r-1)T^(r-1)+...+d_0.                                  (6)
```

On

```text
K_D={x in M^L:D(S)x=0},                                            (7)
```

the prefix observer

```text
pi_r(x)=(x_0,...,x_(r-1))                                         (8)
```

is injective. Indeed, the recurrence `(6)` determines successively
`x_r,...,x_(L-1)` from the first `r` coordinates; no division is used because
the leading coefficient is one. If `r=L`, the prefix is the whole word. If
`r=0`, monicity says `D=1`, so `K_D=0` and the empty observer is complete.

There is a typed Bezout gate. Let `X subset M^L` be the actual set of
candidate differences. Suppose every `x in X` is killed by `P(S)` and
`Q(S)`, and suppose

```text
A(T)P(T)+B(T)Q(T)=a D(T)                                           (9)
```

for `A,B,D in R[T]`, `a in R`, with `D` monic of degree `r`. If multiplication
by `a` is injective on the subset `D(S)X`, then every `x in X` lies in `K_D`.
Consequently `r` prefix samples distinguish the candidate certificates.

If additionally `D` divides both `P` and `Q`, and multiplication by `a` is
injective on

```text
D(S)(ker P(S) intersect ker Q(S)),                                 (10)
```

then one has the exact operator-kernel identity

```text
ker P(S) intersect ker Q(S)=ker D(S).                              (11)
```

For the forward inclusion, `(9)` gives `aD(S)x=0` and the scalar-injectivity
sidecar gives `D(S)x=0`. For the reverse inclusion, divisibility of `P,Q` by
`D` suffices. Over a field, take `a=1` and `D` to be the monic gcd. Over a
general ring, deleting the cleared scalar `a` is invalid.

### 2.2 Resultant torsion firewall

Let `R` be a domain with fraction field `K`, and let monic `P,Q in R[T]` be
coprime in `K[T]`. The Sylvester adjugate gives polynomials `A,B in R[T]`
with

```text
AP+BQ=Res(P,Q).                                                     (12)
```

Hence, for every `R`-module `M`,

```text
ker P(S) intersect ker Q(S)
    subset (M^L)[Res(P,Q)],                                        (13)
```

where the right side is the submodule killed by the resultant. In particular,
an injective resultant action forces the common kernel to vanish. After base
change to a residue field `k(p)`, only prime ideals `p` containing the
resultant can acquire a new common kernel.

That locus is permissive, not deterministic. For `P=T-1,Q=T+1`, the bad
prime two produces a nonzero common constant word over `F_2`. In contrast,
`P=T,Q=T+5` has resultant five, but modulo five both polynomials act as the
invertible cyclic shift, so the common kernel remains zero. Over a PID and a
fixed integral Sylvester presentation, its Smith factors give a finer torsion
atlas; no such Smith-normal-form claim is made for arbitrary domains.

### 2.3 Rectangular normal-Hasse observers

Let `A_0` be a commutative ring, let

```text
R_0=A_0[u_1,...,u_d][[f]],        g in A_0[[f]]^d,                  (14)
```

and fix coordinate bounds `e_1,...,e_d>=0`. On the module

```text
V_e={P in R_0:deg_(u_i)P<=e_i for every i},                         (15)
```

define the rectangular normal-Hasse observer

```text
J_(g,e)(P)=
  (D_u^[alpha]P evaluated at u=g)_(0<=alpha_i<=e_i).                (16)
```

Then

```text
J_(g,e):V_e -> A_0[[f]]^(product_i(e_i+1))                         (17)
```

is an isomorphism. The multivariate Hasse--Taylor identity is

```text
P=sum_(0<=alpha_i<=e_i)
    (D_u^[alpha]P)(f,g)(u-g)^alpha,                                (18)
```

and gives the inverse explicitly. This is the rectangular form of the normal
jet mechanism in THM-4263. It works over arbitrary characteristic because
Hasse derivatives use no factorial division.

A scalar specialization `u=g` is not a replacement for `(16)`: its kernel is
the graph ideal and can meet the candidate difference shell. Once `(16)` has
made the first stage complete, a later quotient is faithful on a confined
shell exactly when the kernel test `(3)` passes on its observed differences.

## 3. Quantitative finite confinement

Let a fixed free integral lattice `Lambda` have a chosen basis, and suppose
the observed shell `W` lies in the coordinate box `[-H,H]^s`, with integer
`H>=0`. Coordinatewise reduction

```text
Lambda -> Lambda/m Lambda                                          (19)
```

is injective on `W` whenever the positive integer modulus satisfies

```text
m>2H.                                                              (20)
```

The proof is the pair-difference test: every nonzero coordinate difference
has absolute value at most `2H<m` and cannot be divisible by `m`. For the
full box and `s>=1`, the threshold is sharp. If `m<=2H`, then
`-floor(m/2)` and `ceil(m/2)` both lie in `[-H,H]`, differ by `m`, and collide.

This statement requires the fixed integral basis and the coordinate quotient
`Lambda/mLambda`. No separation theorem using only the norm of an arbitrary
prime ideal is asserted. A fixed modulus cannot separate an unbounded moving
shell. The lawful order is

```text
characteristic-zero confinement and coordinate-height bound
 -> choose m>2H
 -> solve the finite reduced observer.                              (21)
```

### 3.1 Total observer capacity

For any finite total observer codomain `Omega`, faithfulness on `N>=1` distinct
certificates requires `N<=|Omega|`. In particular, if the **entire** codomain
is exactly `B^r`, with `|B|=q>=2`, then

```text
N<=q^r,              r>=ceil(log_q N),                             (22)
```

and some fibre always has size at least `ceil(N/q^r)`. Extra labels, torsion
flags, branches, or other sidecars enlarge `Omega` and must be counted;
`(22)` cannot be applied to selected coordinates while ignoring them.

The abstract recurrence modules in THM-4258 and THM-4264 enumerate respectively
all `4^3=64` ambient words and all `4^2=16` visible-gated words. Those are
abstract module counts, not geometric-realization claims.

Capacity overload forces a collision, but not density bias. If each of the
`q^r` values occurs with the same multiplicity, normalized fibre weight is
identically one even while raw multiplicities diverge.

## 4. Exact finite-factor density law

Return to the setup of Section 1. Let `Y=q(C)` and give `C,Y` their uniform
probability measures. Define

```text
w(y)=|q^(-1)(y)|,
rho(y)=|Y|w(y)/|C|.                                                 (23)
```

The function `rho` is nonnegative and has target mean one. For every
`H subset Y`,

```text
mu_C(q^(-1)H)
  =(1/|Y|) sum_(y in H) rho(y).                                    (24)
```

Consequently the sharp constant in the domination inequality over all
subsets is

```text
K_sharp=max_(y in Y)rho(y),
mu_C(q^(-1)H)<=K_sharp mu_Y(H).                                   (25)
```

Sharpness follows by taking a singleton at a maximizer. Under the fibre
criterion `(3)`, the numbers `w(y)` are precisely the multiplicities of the
full certificates. If those raw multiplicities are at most `K`, then

```text
K_sharp<=K|Y|/|C|<=K.                                              (26)
```

If `F` is injective, then `q:C->Y` is a bijection and `rho=1` exactly.

### 4.1 Uniform integrability is the exact sequential gate

For a sequence of finite systems indexed by `n`, write `rho_n` for `(23)` and
`mu_n` for uniform measure on `Y_n`. The following are equivalent:

```text
(a) lim_(M->infinity) sup_n integral_(rho_n>M) rho_n dmu_n=0;
(b) for every H_n subset Y_n with mu_n(H_n)->0,
    mu_(C_n)(q_n^(-1)H_n)->0.                                     (27)
```

For `(a) => (b)`, split the integral of `rho_n` over `H_n` at height `M`:
the bounded part is at most `M mu_n(H_n)` and the tail is uniformly small.
Conversely, if `(a)` fails, choose `M_k->infinity` and an unbounded subsequence
`n_k` whose tail integral is at least some fixed `epsilon>0`. The sets

```text
H_(n_k)={rho_(n_k)>M_k},                                           (28)
```

with all other `H_n` empty, satisfy `mu_(n_k)(H_(n_k))<=1/M_k` because
`rho_n` has mean one, but their source masses stay at least `epsilon`. This
contradicts `(b)` and proves the equivalence.

Thus completeness and modular separation preserve identities; they do not
manufacture uniform integrability. Likewise, if events on `Y` are merged or
reused, survivor mass is the product of their **conditional** hazards under
`q_*mu_C`; marginal event sizes alone do not imply independence.

## 5. Sharp failure controls

1. **Graph specialization:** in `A[u][[f]]`, `u` and `f` collide under
   `u->f`, while their first normal-Hasse jets `(f,1)` and `(f,0)` separate.
2. **Cleared scalar:** over `Z`, `P=T-1,Q=T+1` have rational gcd one but only
   a Bezout identity with cleared scalar two. A nonzero constant word over
   `Z/2` is killed by both.
3. **Nonmonic recurrence:** over `Z/4`, `D(T)=2T` kills `(0,0)` and `(0,2)`,
   which have the same one-coordinate prefix.
4. **Bad prime without degeneration:** `T,T+5` have resultant five, yet their
   cyclic-shift common kernel modulo five is zero.
5. **Empty universe:** for `C=empty`, `Delta W=empty`, not `{0}`.
6. **Moving shell:** every fixed modulus eventually collides on `[-H,H]`.
7. **Unbounded order:** for a proposed fixed prefix `r`, zero and the
   length-`r+1` word supported at coordinate `r` are killed by `S^(r+1)-1`
   and share their first `r` samples.
8. **Merged hazards:** repeating one event of mass `1/8` gives conditional
   hazards `(1/8,0,...,0)` and survivor `7/8`, not `(7/8)^m`.
9. **Heavy fibre:** one target of weight `2^n` among `n` targets has target
   mass `1/n` and source mass tending to one.
10. **Balanced overload:** equal repeated multiplicity causes collisions but
    leaves `rho=1`; identity recovery and measure distortion have distinct
    boundaries.

## 6. Reproduction and scope

Run

```text
python -B 04-computation/confined_confluent_observer_transport_audit_thm4274.py
python -B -O 04-computation/confined_confluent_observer_transport_audit_thm4274.py
```

The companion directly tests candidate-pair fibre equivalence, not merely a
histogram. It is finite corroboration; the proofs above establish the abstract
statements.

To apply the theorem, a lane must name:

```text
candidate shell; full certificate; complete first observer;
quotient kernel; observed pair-difference shell; all codomain sidecars;
certificate multiplicities; and conditional event filtration.             (29)
```

Nothing here constructs an off-fibre Jacobian Hom dictionary, proves that the
abstract recurrence words are geometrically realizable, validates an external
p-adic density specialization, transports a Rule-30 prize language, or closes
an LRC(14) edge. Those remain separate proof obligations.

---
id: THM-2069
title: "The k-deletion CRT wheel is the low-weight evaluation code and matroid cogirth"
status: >
  PROVED. Integer rows generating Z^r give, at every prime p, an injective
  rank-r evaluation code. A primitive covector fails some k-deletion gcd
  modulo p exactly when its codeword has weight at most k. Thus the initial
  weight distribution is the exact local wheel and its first failure radius
  is matroid cogirth. When every k-deletion complement has full rational
  rank, its determinantal indices confine all bad primes and CRT gives the
  exact vector, projective, and primitive-density products. Rank deficiency
  is exactly the no-finite-wheel branch. The theorem proves the Paley-e8
  deletion/design bridge and reframes, but does not solve, the [72,36,16]
  support-realization problem or LRC(14).
source: codex-2026-07-21-code-cogirth-wheel
script: 04-computation/k_deletion_code_cogirth_wheel_codex_20260721.py
result: 05-knowledge/results/k_deletion_code_cogirth_wheel_codex_20260721.out
script_sha256: f67736997dfba6973c077c1f1b952287ddf28281c29079f7341107e71019b1a1
result_sha256: 4532964632a4f7d772d5470567dc41df1935982bc236976dfbea7e3ebb89bd3f
hash_basis: normalized repository blobs (LF)
depends_on:
  - THM-2062
  - THM-481
related:
  - THM-211
  - THM-480
  - THM-487
  - HYP-2430
  - HYP-2764
---

# THM-2069 -- The deletion wheel is a low-weight code

Let `rho_1,...,rho_n` generate `Z^r` as an integer lattice. For an integer
covector `d in Hom(Z^r,Z)` and a set `D subset {1,...,n}`, put

```text
g_D(d)=gcd_(i notin D)|d(rho_i)|.                         (1)
```

The covector is **primitive** when its coordinate gcd is one. For a prime
`p`, let `bar d` denote reduction modulo `p` and define the evaluation code

```text
C_p={c_p(phi)=(phi(rho_1),...,phi(rho_n)):
                    phi in (F_p^r)^*} subset F_p^n.       (2)
```

Here `(F_p^r)^*` means the dual vector space, including zero. Write

```text
A_w(p)=#{c in C_p:wt(c)=w},
B_(p,k)=sum_(w=1)^k A_w(p),
G_(p,k)=p^r-1-B_(p,k).                                   (3)
```

The word `c_p(d)` below is nonzero: a primitive integer covector is nonzero
modulo every prime, and integral generation makes evaluation injective.

## 1. Exact local deletion theorem

For `0<=k<=n`, using `gcd(empty)=0` when `k=n`,

```text
there exists D, |D|=k, with p|g_D(d)
    iff wt(c_p(d))<=k.                                   (4)
```

Consequently the exact local counts are

```text
bad nonzero covectors       = B_(p,k),
bad projective directions  = B_(p,k)/(p-1),
good nonzero covectors      = G_(p,k),
good projective directions = G_(p,k)/(p-1).              (5)
```

### Proof

If `p|g_D(d)`, all word coordinates outside `D` vanish, so the support has
size at most `k`. Conversely, a support of size at most `k` can be enlarged
to a set `D` of size exactly `k`; every coordinate outside it then vanishes.
This proves (4).

Multiplication by `F_p^x` acts freely on nonzero covectors, preserves word
support, and has orbits of size `p-1`. This proves (5). QED.

The boundary cases are honest. At `k=0`, the full row family generates the
lattice, so every primitive specialization has gcd one. At `k=n`, every
direction is bad under the empty-gcd convention; this formal boundary has no
finite-index deletion theorem.

## 2. Cogirth, coloops, and exact cosimplicity

Let `M_p` be the rank-`r` column matroid over `F_p` represented by the matrix
whose columns are the `rho_i mod p`. Minimal nonzero supports in the row space
of a representing matrix are exactly the cocircuits of its column matroid.
Therefore

```text
d_min(C_p)=cogirth(M_p),                                 (6)
B_(p,k)=0 iff k<cogirth(M_p).                            (7)
```

Thus cogirth is literally the first deletion radius at which a local wheel
loses a direction. Two useful specializations are

```text
A_1(p)/(p-1)=#{coloops of M_p},                          (8)
B_(p,2)=0 iff M_p is cosimple.                           (9)
```

For (8), a singleton cocircuit is a coloop and has precisely `p-1` scalar
codewords. Equation (9) says that there is no cocircuit of size one or two,
equivalently that the dual matroid is simple. One must not interpret
`A_w/(p-1)` as the number of cocircuit supports for arbitrary `w`: a
nonminimal codeword support need not be a cocircuit.

At `k=1`, this recovers THM-2062's higher-rank coloop count. The weight-code
form is what survives for `k>=2`, where a bad deletion annihilator can have
dimension greater than one and different deletion kernels can overlap.

## 3. Finite determinantal wheel versus rank deficiency

Fix `0<=k<n`. For `|D|=k`, let

```text
L_D=sum_(i notin D) Z rho_i.                             (10)
```

Assume first that every `L_D` has rational rank `r`, and put

```text
I_D=[Z^r:L_D]
   =gcd of the absolute r by r minors of the retained rows,
Q_k=rad(lcm_(|D|=k) I_D).                               (11)
```

Then for every primitive `d`,

```text
g_D(d) divides I_D.                                     (12)
```

Indeed `d:Z^r->Z` is surjective. Reduction modulo `g_D(d)` vanishes on
`L_D`, so `Z^r/L_D` surjects onto `Z/g_D(d)Z`; its order is divisible by the
specialized gcd. It follows that

```text
p|Q_k iff B_(p,k)>0.                                    (13)
```

The forward direction says that some retained matrix loses rank modulo `p`,
so a nonzero covector annihilates it and has support inside its deleted set.
The reverse direction uses (4), or equivalently the same rank loss.

If some `L_D` has rational rank below `r`, every maximal minor is zero. Its
reduction has rank below `r` for **every** prime, so `B_(p,k)>0` for every
prime. This is the exact no-finite-wheel branch, not a large-index case.

## 4. Exact CRT and primitive-density products

Continue in the full-rank branch. A primitive covector is good for every
`k`-deletion precisely when

```text
wt(c_p(d))>k for every p|Q_k.                            (14)
```

Ordinary CRT gives exactly

```text
product_(p|Q_k) G_(p,k)                                 (15)
```

locally unimodular vector residues modulo `Q_k`. The corresponding product
of local projective choices is

```text
product_(p|Q_k) G_(p,k)/(p-1).                          (16)
```

Every factor is integral because scalar multiplication acts freely. A zero
factor means exact emptiness; unlike the `k=1` coloop wheel, positivity is not
automatic at larger deletion radius.

For `r>=2` and a bounded Jordan-measurable region `Omega subset R^r`, with
`d` and `-d` counted separately, Mobius inversion and equidistribution in
the fixed residue classes give

```text
#{d in T Omega : d primitive and every k-deletion is primitive}
 =vol(Omega)T^r/zeta(r)
    * product_(p|Q_k) G_(p,k)/(p^r-1) + o(T^r).          (17)
```

Equivalently the absolute Euler product is

```text
vol(Omega)T^r
 * product_(p not|Q_k)(1-p^(-r))
 * product_(p|Q_k) G_(p,k)/p^r.                         (18)
```

Formula (17) replaces the primitive local factor `(p^r-1)/p^r`; multiplying
`G/p^r` on top of `1/zeta(r)` would count that condition twice. The density
statement excludes `r=1`, where primitive covectors are only `+1,-1`.

## 5. A sharp one-bad-line family

For a prime `ell`, take the four rows in `Z^2`

```text
(1,0), (1,ell), (1,2ell), (0,1),       k=1.             (19)
```

Their deletion indices are `1,1,1,ell`. Modulo `ell`,

```text
C_ell={(a,a,a,b):a,b in F_ell},
W_C(z)=1+(ell-1)z+(ell-1)z^3+(ell-1)^2 z^4.             (20)
```

There are `ell-1` bad vectors but exactly one bad projective line. Hence

```text
G_(ell,1)=ell(ell-1),
relative primitive factor=ell/(ell+1).                  (21)
```

For a primitive integer parameter `(N,M)`, the first three deletion gcds are
one and the fourth is

```text
gcd(N,N+ell M,N+2ell M)=gcd(N,ell).                     (22)
```

Thus goodness is exactly `ell` not dividing `N`, and (17) gives density
`zeta(2)^(-1) ell/(ell+1)`. This family catches both the vector/projective
normalization and the replacement of the primitive Euler factor.

## 6. LRC leverage: robustness of a relation template

In a saturated multi-anchor LRC template, let the `rho_i` be the coefficient
rows and let `d` be the primitive anchor parameter. Then `d(rho_i)` are the
specialized speeds. THM-2062 computes hereditary primitivity under one
deletion. The present theorem computes every higher deletion layer without
inclusion-exclusion over deletion kernels:

```text
prime p survives k runner deletions
    iff the parameter codeword has weight > k.           (23)
```

The first fragile layer is the modular cogirth. This is a lossless atlas
filter for relation-plane and multi-anchor searches, and `k=2` identifies
exactly when the modular coefficient matroid is cosimple.

It does not finish LRC(14). A fixed-coordinate slice such as `N=constant`
meets projective bad sets nonuniformly and can be empty even when (17) is
positive; THM-2062 Section C performs that affine count separately. More
importantly, the wheel retains primitivity but forgets positivity, collision
walls, phase height, endpoint owners, and the persistent marked circuits left
by THM-2065. Those sidecars remain the structured LRC residual.

## 7. Paley `P_7`, `e8`, and the two Fano planes

Let `Q={1,2,4}` in `F_7`, and orient `i->j` when `j-i in Q`. For each vertex
`v`,

```text
N^+(v)=v+Q,                  N^-(v)=v-Q.                 (24)
```

Both triples are directed cycles. The two translation orbits are distinct:
if `Q+a=-Q+b`, summing their elements gives `3a=3b`, hence `a=b`, contrary
to `Q!=-Q`. Each of `Q` and `-Q` is a `(7,3,1)` difference set because its
six ordered nonzero differences are all of `F_7^x`. Their fourteen translates
therefore form two distinct block-disjoint cyclic Fano planes and together a
simple `2-(7,3,2)` design. They exhaust the directed triangles because

```text
C(7,3)-7 C(3,2)=35-21=14.                               (25)
```

Moreover `N^+(v)` and `N^-(v)` are disjoint and cover `F_7\{v}`, producing
the seven near-parallel classes. This proves the structural part of THM-211
and corrects its former explanation: a tournament gives one orientation to a
fixed cyclic vertex triple, not two orientations of one Fano line.

THM-481 identifies the binary loss-gauge code of `border(P_7)` with `e8`.
Its seven gauge rows and their complements have supports

```text
G_v={infinity} union N^-(v),
H_v={v}        union N^+(v).                             (26)
```

They are the fourteen weight-four words in

```text
W_e8(y)=1+14y^4+y^8.                                    (27)
```

The triangle bridge is a **marked deletion**: delete the fixed infinity
marker from `G_v`, but delete the varying row label `v` from `H_v`. Forgetting
that marker would identify objects that the construction keeps distinct.

Applying (3)--(7) at `p=2`, all fifteen nonzero/projective code directions
survive deletion radii `k<=3`. At `k=4`, fourteen fail and exactly the
all-ones direction survives; at `k=8`, none survives. Thus the Fano
near-resolution is the first cocircuit shell of the `e8` deletion wheel.

## 8. The `[72,36,16]` support-realization gate

If an extremal binary Type II `[72,36,16]` code exists, its prime-two wheel
is full through `k=15`. At `k=16`, exactly

```text
A_16=249849                                                (28)
```

of its `2^36-1=68719476735` nonzero projective directions become bad, leaving
`68719226886`. Their supports must be the blocks of the forced
`5-(72,16,78)` design; indeed

```text
78 C(72,5)/C(16,5)=249849.                               (29)
```

This sharpens HYP-2430: seek a rank-36 binary column matroid whose deletion
wheel is full through radius fifteen and whose first cocircuit shell realizes
the prescribed design. It does not construct that matroid. The scalar weight
count cannot recover block incidence, self-duality, or doubly-evenness.
THM-481/487's Paley `eQR(72)` route has minimum distance twelve, so its wheel
already fails at radius twelve; this is the same obstruction in code-wheel
language.

## 9. Finite-depth constant-term selector

There is a rigorous but deliberately finite by-product. Suppose the product
in (17) is positive, and let `Delta subset Z^r\{0}` be finite. Each exact
hyperplane

```text
H_delta={d:d(delta)=0}
```

contains only `O(T^(r-1))` lattice points in `T Omega`. Removing their finite
union leaves the leading asymptotic (17) unchanged. Hence a primitive good
covector can simultaneously avoid every member of `Delta`.

For a Laurent polynomial

```text
F(X)=sum_i c_i X^(u_i)
```

and a fixed depth `D`, take `Delta_D` to be the nonzero exponent sums
`sum_i k_i u_i` with `k_i>=0` and `1<=sum k_i<=D`. A good covector avoiding
`Delta_D` satisfies

```text
CT(F(z^d)^m)=CT(F^m),              1<=m<=D.              (30)
```

Indeed no nonzero exponent occurring through depth `D` can specialize to
zero. This is an exact finite-depth one-variable reduction, useful for
bounded constant-term audits and finite reconstruction atlases.

It is not a Gaussian Moment proof. The selector depends on `D`; no scalar
map is injective on all of `Z^r` for `r>1`; and Gaussian Wick weights depend
on the full vector occupation, not only its scalar charge. Physical angular
rotation is even less useful on a balanced channel: its total charge is
already zero. Arbitrary radial coefficients and global channel
noncancellation remain untouched.

## 10. Assumption challenge and tournament analysis

Three tempting quotients are invalid:

1. Rational spanning is insufficient. The rows `(2,0),(0,1)` span `Q^2`, but
   at `p=2` the nonzero covector `(1,0)` evaluates to the zero word. Integral
   generation is load-bearing.
2. A rational-rank-deficient deletion is not an index-zero member of a finite
   sieve; it is bad at every prime.
3. Code weights count parameter directions, not labelled support incidence.
   The latter is exactly the missing `[72,36,16]` design information.

For the general wheel, the faithful vertices are cocircuit supports or
deletion obligations, not runners; forcing an orientation destroys the
linear incidence and supplies no theorem. In the Paley application an honest
tournament is present. Its pairwise observable is the quadratic character of
`j-i`; the residue set `Q` fixes the gauge, and multiplication by a
nonresidue reverses it. There are no ties, so no artificial tie path is used.
The exact fingerprint is score histogram `{3:7}`, fourteen directed
three-cycles, one strongly connected component, and `189` Hamiltonian paths.

## Computational audit

The companion is integer and finite-field exact and uses explicit runtime
checks rather than Python assertions. It checks fifty integrally spanning
rank-two/rank-three templates over `p=2,3,5`, comprising `855` direct
deletion/weight equivalences, `855` projective orbit counts, `150`
cocircuit/cogirth comparisons, `404` index-prime equivalences, and `113`
direct CRT residues. The sharp family (19) is checked at four primes and on
`4096` primitive integer parameters. The same referee reconstructs the
Paley design, seven near-classes, `e8` enumerator and wheel, all `189`
Hamiltonian paths, and the exact length-72 arithmetic. Normal and `python -O`
runs byte-match the frozen output and end in `PASS`.

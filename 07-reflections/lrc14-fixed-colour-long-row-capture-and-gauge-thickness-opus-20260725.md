---
source: codex-2026-07-25-lrc-fixed-colour-thickness
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  Under an
  explicit post-transfer total-step hypothesis, THM-2323 closes every
  physical residue row whose order u is at least the product of the
  actual word/bare jump counts.  Thus complementarity is confined to a
  finite short-order regime.  A shifted u=2 interval hostile satisfies
  THM-2323's rational positive-comb hypotheses, has an XOR-aligned common
  endpoint, and has common support in every four gauges of every
  primitive 91- and 169-colour, yet its two prescribed diagonal samples
  remain complementary.  Post-summation cyclotomic idempotents cannot
  repair a zero product.  The surviving short-row sidecars are
  pre-summation private-address separation or genuinely two-coordinate
  target mixing together with collision-multiplicity separation;
  nonconstant transport on one active axis is ruled out by the
  strengthened THM-2344, while the audited pending THM-2353 handoff
  shows that two active axes alone can still be twist-constant.  This
  is not a canonical word stratum, excludes no scalar profile, and does
  not prove LRC(14).
depends_on:
  - THM-2304-deepest-boundary-cyclotomic-current-separation
  - THM-2323-primitive-fixed-colour-cross-correlation-and-same-gauge-word-alignment
  - THM-2344-correlation-inverse-rigidity-and-aligned-tooth-twist-hostile
related:
  - THM-2305-canonical-blocker-word-handoff-hypergraph
  - lrc14-word-current-parity-zero-divisor-opus-20260725
script: 04-computation/lrc14_fixed_colour_gauge_thickness_hostile_probe.py
output: 05-knowledge/results/lrc14_fixed_colour_gauge_thickness_hostile_probe.out
script_sha256: b336381b1db1703bc8eb699f31cdba480492975fb9392bb5a321cca88b6f7326
output_sha256: 5e1bd9d7cb518f1ae8b64c236a9b88d6efd7f89aa81849e10f502d5a92ab8d3f
hash_basis: LF-normalized working-tree bytes
---

# Fixed colour closes long rows, but not a thin diagonal

> **INHERITANCE UPDATE.**  Promoted THM-2349 (`4389031a2`) gives every
> repeated `(1,1,c)`, `5<=c<=19`, positive shallow-owner mass, a literal
> THM-2305 word after finite BV delay, and a marked primitive-`91`
> deepest edge.  Thus repeated-row marked incidence is no longer a
> separate debt: the THM-2331/2334/2343/2344 current boundary now reaches
> all `165` rows.  The theorem below still requires its explicit
> post-transfer total-pair hypotheses; THM-2349 supplies incidence, not
> the missing collision-separated phase.

## 1. The exact long-row transport lemma

Let `E_Q subset E` be a total word/bare pair after every bank belonging
to the chosen word has already been summed.  Fix integers

```text
b>lambda>=0,

N=13^(b-lambda)>1,

n_k=13^lambda(kappa+Nk),             0<=k<u,

gcd(kappa,N)=1.                                      (1)
```

Assume the legal post-transfer factors are

```text
H_Q=P^lambda 1_(E_Q),

H_E=P^lambda 1_E,                                   (2)
```

with the exact Fourier ancestry

```text
(H_Q)_hat(q)=(1_(E_Q))_hat(13^lambda q),

(H_E)_hat(q)=(1_E)_hat(13^lambda q).                (3)
```

Suppose further that:

```text
0<=H_Q<=H_E;

H_Q,H_E are rational-valued step functions with rational breakpoints;

integral H_Q>0;

support(H_E) subset D_a for some 13 not_divides a.  (4)
```

Write `J_Q,J_E` for their actual nonzero jump counts and put

```text
L=J_Q J_E.                                          (5)
```

Since `gcd(a,N)=1`, the arithmetic-comb form of THM-2323 applies at
modulus `N` and primitive colour `kappa mod N`.  If `kappa_0` is the
canonical representative and `kappa=kappa_0+NH_0`, then the physical
indices `0<=k<u` are exactly the consecutive gauge block
`H_0,...,H_0+u-1`.  The block form of THM-2323 says that every `L`
consecutive gauge indices contain a `k` for which

```text
(H_Q)_hat(kappa+Nk)
(H_E)_hat(kappa+Nk) !=0.                            (6)
```

If `u>=L`, one such index lies in the physical row (1).  Equation (3)
transports (6) back to the two original Fourier coefficients at `n_k`.
Because `n_k!=0`, distributional differentiation multiplies both by a
nonzero scalar.  Hence both **total** derivative currents are nonzero at
the same row coordinate.  In the cyclic row algebra this is

```text
u>=J_QJ_E  implies  F_jW_j!=0.                       (7)
```

This is a transport theorem, not a bankwise heuristic: (2) must be the
post-summation total pair and (3) must be the actual transfer ancestry
for the physical row.

For the exact THM-2323 Section 6 word/bare pair,

```text
J_E<=2S,                     J_Q<=6S,

L<=12S^2.                                           (8)
```

Thus `u>=12S^2` is a convenient inherited sufficient condition on every
stratum where (1)--(4) have been verified.  The sharper local cutoff is
the actual product (5).  Any complementary-row obstruction within this
scope is forced into

```text
u<J_QJ_E.                                           (9)
```

## 2. A hostile with XOR alignment and fixed-colour overlap

Let

```text
J=[-1/4420,1/4420],

I=[-1/4420,1/2210],

E=I disjoint-union (I+1/2),

Q=T^2J,                         T(x)=13x.            (10)
```

The left endpoints of `I,J` coincide, while the right endpoint of `J`
lies strictly inside `I`.  The wrong `T^2` preimages are separated by
`1/169` on the same component and `1/338` from the half-shifted one.
The interval widths give strict clearance, so

```text
E intersect T^(-2)Q=J.                              (11)
```

The common left endpoint is an even-clock XOR-aligned source/target
jump.  The right endpoint of `J` is target-only.

After one Perron transfer,

```text
f=P^1 1_J
  =(1/13)1_[-1/340,1/340],

g=P^1 1_E
  =(1/13)(
      1_[-1/340,1/170]
      +1_([-1/340,1/170]+1/2)).                    (12)
```

Therefore

```text
0<=f<=g,

support(g) subset D_2,

gcd(2,91)=gcd(2,169)=1.                             (13)
```

These are exactly the abstract positive rational-comb hypotheses of
THM-2323.  For every nonzero integer `n`, direct interval Fourier
calculation gives

```text
f_hat(n)=0  iff  170 divides n;

g_hat(n)=0  iff  n is odd or 340 divides n.         (14)
```

Thus `f_hat(n)g_hat(n)!=0` exactly when `n` is even and not divisible by
`170`.  In every four consecutive gauges of every primitive colour
modulo `91` or `169`, two sampled integers are even; they differ by twice
the odd modulus and therefore cannot both be divisible by `170`.
Consequently every such four-gauge block contains common support.  This
is stronger than the general `J_fJ_g=8` bound.

## 3. The prescribed `u=2` diagonal still misses

Use the literal arithmetic packet

```text
lambda=1,                 b=3,

A'=13,

c_*=2*13^5=742586,

A=A'+170c_*=126239633.                              (15)
```

After dividing the two physical rows by `13`, their two addresses are

```text
row 0:  1,        170;

row 1:  9710741,  9710910=57123*170.                (16)
```

Both row bases are congruent to `1 mod 169`.  The four primitive
`91`-colours are

```text
1, 40, 79, 27.                                      (17)
```

At address zero, both entries in (16) are odd, so `g_hat=0` and
`f_hat!=0`.  At address one, both are `170` times an odd integer, so
`f_hat=0` and `g_hat!=0`.  The first common coefficient in either
`169`-row occurs only at the later gauge `k=3`; the physical order is
only `u=2`.

For each of the two rows, interpolation in `C[X]/(X^2-1)` therefore
gives

```text
F_j proportional to 1-X,

W_j proportional to 1+X,

F_jW_j=0.                                           (18)
```

At the even samples, the aligned left jump and the target-only right
jump are separated by `1/170`.  Their phases are both `-1`, while their
jump signs are opposite.  They cancel **inside the total word current**
in the same exact-order/conductor cell.  This is why a common XOR
endpoint and abundant nearby fixed-colour overlap do not force the thin
diagonal (16) to hit.

The construction is a local word control.  It is not a full canonical
THM-2305 Boolean stratum, a global zero-safe cover, or an LRC
counterexample.

## 4. Why post-summation idempotents cannot repair it

Let

```text
R=K[C_u]
```

over a characteristic-zero splitting field.  It is semisimple and
commutative.  If the total row elements obey `FW=0`, then for every
central or exact-order idempotent `e`,

```text
(eF)(eW)=e^2FW=eFW=0.                               (19)
```

Cyclotomic refinement can diagnose which components are complementary;
it cannot create overlap after summation.  Separation must act before
the cancelling endpoints are merged.  The clean sufficient sidecar is:

```text
one XOR-aligned common endpoint has a relative-conductor/address
class whose total source and total word class sums are both nonzero
and which is private from every cancelling endpoint.                 (20)
```

This is the private top-digit route suggested by THM-2304.

## 5. The strengthened THM-2344 stopping rule

The current THM-2344 rules out a tempting weaker repair.  For an
arbitrary nonconstant transported word on one active axis and
`13 divides R`, its translated factors have the form

```text
L_t(a)=zeta^(at)L_0(a),

F_t(b)=zeta^(bt) I_hat(b).                           (21)
```

On the inverse-character line `b=a+m`, the target correlation is

```text
K_t=constant*zeta^(-mt),
```

so the full amplitude is constant in all twists.  Its genuine hostile
uses `I=g`, `W=d`, `R=13`, `a=1`, `b=2`: the word is nonconstant, the
full amplitude is positive, and all `169` twists remain constant.

Therefore generic word nonconstancy is not a phase sidecar.  After the
long-order reduction (7), a short row must supply one of:

```text
1. a private pre-summation relative-address class as in (20);

2. a thick block of actual physical gauges, rather than one diagonal
   sample per colour;

3. at least two independently translated base coordinates, or
   cross-axis word/base mixing, followed by a retained 91-unit and
   terminal-phase certificate.                                    (22)
```

Only after (22) creates a genuinely non-diagonal target operation can
THM-2344's detecting-twist alternative be invoked.

> **Collision-multiplicity stopping rule (audited incoming THM-2353
> handoff; not used as a proved dependency here).**  Two active axes,
> an active positive even word, and `l<=f` can still leave all `169`
> deep twists constant when two left atoms collide in the same target
> residue and the right atom sees only their aggregate.  The first
> mod-`169` Bockstein recovers a vertical collision tax, but not a
> nonvertical detecting phase.  Thus item 3 is only a source of
> candidate variation: a closure must additionally separate collision
> multiplicity, for example through a nonvertical ancestry/Bockstein
> coordinate or an acute phase-cone difference certificate.  Generic
> cross-axis nonconstancy is no more sufficient than generic
> one-axis word nonconstancy.

THM-2352's pushed `q`-adic prefix/collision-abscissa theorem gives the
orthogonal scale boundary.  Every finite interaction complex sees the
full simplex of indexed collision profiles, while every positive
Dirichlet abscissa and eventual termination remain dense Haar-null
phenomena inside each finite cylinder.  Thus no fixed finite prefix bank
can certify that the collision tax eventually disappears.  A private
class intended to work uniformly down the tower must retain an
infinite-depth end marker or radix clock; finite short-row enumeration
is legitimate only after the physical depth has itself been fixed.

## 6. Exact reproduction and scope

The standard-library companion verifies:

- the shifted interval and pullback clearances;
- the common and target-only endpoints;
- both exact zero sets in (14);
- every four-gauge block for all `72` primitive `91`-colours and all
  `156` primitive `169`-colours;
- the physical rows, colours, and complementary zero patterns;
- the endpoint phase cancellation; and
- the zero product in `C[C_2]`.

Normal Python, optimized Python, and the stored transcript agree
byte-for-byte.  The long-row implication is the mathematical block form
of THM-2323 plus the explicit ancestry (3), not a finite-computation
assumption.

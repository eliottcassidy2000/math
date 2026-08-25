---
id: THM-4107
title: "GCD-normalized exponent tournament holonomy and LRC blindness"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY VERIFIED-EXACT. Raw a^b versus b^a is a scalar preorder. Pairwise
  gcd normalization instead creates a scale-invariant, tie-free tournament:
  for a<b the smaller loses exactly on divisibility edges and ratio 2:3.
  Its directed triangles are the exact curvature of these incompatible edge
  gauges. Nevertheless AP13 and infinitely many uniformly loose thirteen-
  speed rows have identical labelled tournaments, so the carrier does not
  preserve LRC tightness, margin, phase, or arrival.
source: codex-lrc14-abc-exponent-reciprocity-20260825
depends_on: []
related:
  - THM-4057-stern-brocot-depth-pullback-and-rational-edge-tournament-gauge
  - THM-4088-order-tournament-arithmetic-type-blindness-and-lrc-margin-density
  - THM-4095-exact-arithmetic-field-transport-gapped-pair-margins-and-order-tournament-blindness
  - THM-4105-primitive-reciprocal-phase-descent-and-quantitative-arrival
  - THM-4108-abc-conditional-reciprocal-power-radicals-and-lrc-gauge-obstruction
script: 04-computation/reciprocal_power_abc_lrc_thm4107_4108.py
output: 05-knowledge/results/reciprocal_power_abc_lrc_thm4107_4108.out
script_sha256: 977490728ed77150daf580dbc294f7355bc502cf85c750748e8d4ed532d5691b
output_sha256: 6e91f10ea6db2b117937d439bd7d7a9b300ce3fa31a10b92e3811dee79ecbe32
independent_audit_script: 04-computation/gcd_normalized_power_tournament_thm4107_independent_audit.py
independent_audit_output: 05-knowledge/results/gcd_normalized_power_tournament_thm4107_independent_audit.out
independent_audit_script_sha256: 152651161920f7b0fc603f345c2d0092e67f29b16a487cd59f32dc39355641b1
independent_audit_output_sha256: ac300d909e83bc4bda38240b390b70da525798338d22b468cc624591e9593ffe
hash_basis: raw LF bytes
audit: >
  PASS. Independent classification checked 19,900 pairs and 82,160 triples,
  the rational equality catalogue, 19 AP13 terminal-prime twins, exact LRC
  maxima for AP13/V26 and pair controls, and normal/-O/frozen transcript
  identity.
---

# THM-4107 -- GCD-normalized exponent tournament holonomy and LRC blindness

**PROVED + VERIFIED-EXACT + INDEPENDENTLY VERIFIED-EXACT.** There are two distinct reciprocal-power objects.
Raw comparison has one global scalar potential and no curl. Reducing each
edge by its own gcd is compatible with LRC common dilation, but the gauges no
longer glue globally; their incompatibility creates an exact sparse
tournament holonomy.

## 1. Raw reciprocal powers are a scalar preorder

For positive reals put

```text
phi(x)=log(x)/x.                                          (1)
```

Then

```text
a^b>b^a iff phi(a)>phi(b).                               (2)
```

Hence raw comparison is a total preorder and contains no directed triangle.
The stronger separable model

```text
A_i^(B_j)>A_j^(B_i) iff log(A_i)/B_i>log(A_j)/B_j        (3)
```

is equally transitive. Ordered exponent syntax alone does not create a
tournament obstruction.

The multiplicative edge values `R(a,b)=a^b/b^a` satisfy

```text
R(a,b)^c R(b,c)^a R(c,a)^b=1,                           (4)
```

the weighted zero-holonomy identity.

### Rational equality ray

All non-diagonal positive rational solutions of `x^y=y^x`, with `x<y`, are

```text
x_n=((n+1)/n)^n,
y_n=((n+1)/n)^(n+1),                    n>=1.            (5)
```

Indeed, if `y/x=m/n` in lowest terms, then

```text
x=(m/n)^(n/(m-n)).                                       (6)
```

Since `gcd(n,m-n)=1`, rationality forces `m` and `n` to be perfect
`(m-n)`-th powers. If `d=m-n>=2`, writing `m=M^d,n=N^d` gives
`d=M^d-N^d>=2^d-1>d`, a contradiction. Thus `m=n+1`, which gives `(5)`.

The ratio nodes `(n+1)/n=[1;n]` form a Farey-neighbor Stern ray because

```text
(n+1)^2-n(n+2)=1.                                        (7)
```

Reciprocal reflection swaps the two endpoints. Clearing the additive gap in
`(5)` gives

```text
n(n+1)^n+(n+1)^n=(n+1)^(n+1).                           (8)
```

After division by the common gcd `(n+1)^n`, this is only the unit ABC triple

```text
n+1=n+1,                  with slots (n,1,n+1).          (9)
```

Thus the exponential equality height is common-factor amplification, not ABC
leverage. The only distinct positive-integer tie is `{2,4}`.

## 2. The primitive exponent tournament

For distinct positive integers `a,b`, put

```text
g=gcd(a,b),              p=a/g,              q=b/g,
a -> b iff p^q>q^p.                                        (10)
```

This relation is invariant under common dilation and is tie-free.

### Exact edge classification

For `a<b`,

```text
a -> b unless a|b or 3a=2b;
b -> a exactly in those two cases.                       (11)
```

To prove this, take coprime `p<q`. If `p=1`, then `1<q`; this is the
divisibility exception. If `(p,q)=(2,3)`, then `8<9`. For `p=2,q>=5`,
`2^q>q^2`, starting from `32>25` and continuing by induction. If `p>=3`,
`log(x)/x` is strictly decreasing on `[3,infinity)`, so `p^q>q^p`.
These cases are exhaustive and translate exactly to `(11)`.

The reversed smaller-to-larger edges therefore occupy the reduced rational
set

```text
{1/k:k>=2} union {2/3}.                                  (12)
```

The first lane is divisibility; `2/3` is the unique sporadic classic-power
exception. The orientation bit forgets which mechanism fired, so that tag is
a necessary sidecar in structural applications.

## 3. Exact triangle curvature

For `u<v`, let

```text
e(u,v)=1 iff u|v or 3u=2v.                               (13)
```

For every sorted triple `a<b<c`, its tournament is cyclic exactly when

```text
boxed: e(a,b)=e(b,c) != e(a,c).                          (14)
```

Start with the ascending transitive tournament and reverse the exceptional
edges. The two cyclic reversal words are precisely `(0,1,0)` and `(1,0,1)`
in the coordinate order `(ab,ac,bc)`, proving `(14)`.

The tournament on `[5]` is transitive. At order six the first cycles appear;
a minimal witness is

```text
2 -> 5 -> 6 -> 2.                                       (15)
```

Thus the curl does not come from transcendence or from the raw exponent
function. It is the holonomy of incompatible pairwise multiplicative gcd
charts placed on an additive total order.

Every primitive positive additive packet `A+B=C` with three distinct terms
gives a transitive triple under `(10)`: its terms are pairwise coprime, so the
pairwise gauges glue to one global potential order. The repeated-term boundary
`1+1=2` is not a three-vertex tournament. A cyclic exponent triple is
therefore not an ABC quality signal.

## 4. Higher power levels isolate the 2--3 exchange path

For an integer `k>=1`, define

```text
a ->_k b iff p^(q^k)>q^(p^k),
equivalently log(p)/p^k>log(q)/q^k.                      (16)
```

For every `k>=2`, the potential is strictly decreasing on integers at least
two, so the only reversed smaller-to-larger edges are divisibility edges.
Consequently

```text
T_1 triangle T_k={{2h,3h}:h>=1},             k>=2.      (17)
```

Write a vertex as `r 2^i 3^j` with `gcd(r,6)=1`. An edge in `(17)` transfers
one unit between the `2`-adic and `3`-adic exponents, preserving `r` and
`i+j`. Each component is the finite path

```text
r 3^K -- r 2*3^(K-1) -- ... -- r 2^K.                  (18)
```

This is a precise bridge to recurrent mod-six mechanisms in the repo. It is
only an exchange carrier until residues, owners, and a physical time are
attached.

## 5. Exact LRC(14) blindness fibre

Let

```text
A={1,2,...,13},
B_p={1,2,...,12,p},                 p>13 prime.          (19)
```

Under the identity map on the first twelve labels and `13 -> p`, the two
primitive exponent tournaments are identical. Internal edges agree. Both
terminal vertices lose to every `a in [2,12]` and beat `1`, because none of
the other exceptional conditions in `(11)` occurs.

Nevertheless,

```text
M(A)=1/14,                                                (20)
```

whereas at `t=1/13` every speed in `B_p` has a nonzero residue modulo `13`.
Therefore

```text
F_(B_p)(1/13)>=1/13,
M(B_p)>=1/13>1/14.                                       (21)
```

There are infinitely many such primes. One fixed **labelled** tournament
fibre therefore contains tight AP13 and infinitely many rows with a uniform
strict margin. The tournament preserves common-dilation classes and sparse
divisibility curvature, but destroys the modulus-13 residue, simultaneous
minimum, phase, arrival, and optimized LRC margin.

This hostile is stronger than a mere collision of `c3` or Hamiltonian-path
counts: even the complete labelled orientation agrees.

## 6. Generated next tasks

1. Attach two edge colours—divisibility and `2:3` exchange—to existing LRC
   maximizing trees. Test the full coloured tree against the mandatory fibre
   `(19)` before using it as a decoder.
2. On live mod-six ledgers, retain the path address `(r,v_2+v_3)` from `(18)`
   and test whether endpoints or exchange distance control an actual blocker.
3. Replace the scalar tournament by a modulus-indexed collision deck
   `p^q=+-q^p mod N`, retaining `N`, sign, owner, and zero cards. Compare this
   with the signed-pair deck of AP13 and the canonical V26 hostile.
4. Derive divisor-sum recurrences for the initial-segment cycle count directly
   from `(14)`; use the exact census only as a verifier, not as proof.
5. Treat the Stern equality ray `(5)` as an ABC normalization canary: any
   apparent exponential-height gain must survive primitive gcd removal.

## 7. Exact audit

Reproduce with

```text
python -B 04-computation/reciprocal_power_abc_lrc_thm4107_4108.py
python -B -O 04-computation/reciprocal_power_abc_lrc_thm4107_4108.py
```

The normal and optimized streams equal the frozen transcript. The exact audit
checks `7,140` sorted edges, `34,220` triples, `3,540` valuation packets,
`639` common-dilation rays, the AP13/V26 isomorphism control, and several
strict LRC hostiles. An independent path checks `19,900` edges, `82,160`
triples, `19` terminal-prime twins, and exact LRC maximum controls. LRC(14)
remains **OPEN**.

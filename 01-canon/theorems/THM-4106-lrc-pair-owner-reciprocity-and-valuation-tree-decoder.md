---
id: THM-4106
title: "LRC pair-owner reciprocity and valuation-tree decoder"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY VERIFIED-EXACT. For a primitive two-speed LRC observer, the exact
  owner imbalance is reciprocal to the product of the sum and difference
  clocks. Together with the optimized pair margin, it reconstructs every
  mixed-parity labelled primitive pair. For a primitive positive row that is
  not all odd, N-1 cross-v2 pair summaries on a spanning tree reconstruct the
  entire labelled row. The representation is lossless arithmetically but
  discards simultaneous phase and arrival, so LRC(14) remains open.
source: codex-lrc14-abc-exponent-reciprocity-20260825
depends_on:
  - THM-4095-exact-arithmetic-field-transport-gapped-pair-margins-and-order-tournament-blindness
related:
  - THM-668-pair-sum-ruler-witness-structure
  - THM-3742-square-triangular-pell-mod13-central-sign-projective-cycle
  - THM-4057-stern-brocot-depth-pullback-and-rational-edge-tournament-gauge
  - THM-4105-primitive-reciprocal-phase-descent-and-quantitative-arrival
script: 04-computation/lrc_pair_owner_reciprocity_thm4106.py
output: 05-knowledge/results/lrc_pair_owner_reciprocity_thm4106.out
script_sha256: 4266285916d0fbe725d89bb888f52212b63ca5dba9722de2d14fe4ae4bac2182
output_sha256: 471f97276a99243c86ff83d420882e435137c45d8590f6ab2540ce9287d899fb
independent_audit_script: 04-computation/lrc_pair_owner_reciprocity_thm4106_independent_audit.py
independent_audit_output: 05-knowledge/results/lrc_pair_owner_reciprocity_thm4106_independent_audit.out
independent_audit_script_sha256: 42035f1d1ce2032d374b8dda6dfc65361e3615975bf13955ec8b69f86b9f1693
independent_audit_output_sha256: 6e80fd339472037f3663bc3a76df351504b5803e367bb3276883d7abe43bcc1b
hash_basis: raw LF bytes
audit: >
  PASS. An independent cell-decomposition and common-grid audit checked all
  4,950 scaled pairs through 100, both directed decoders, 24,930 labelled
  five-row orders, and 3,600 Fourier correlations, with normal/-O/frozen
  transcript identity.
---

# THM-4106 -- LRC pair-owner reciprocity and valuation-tree decoder

**PROVED + VERIFIED-EXACT + INDEPENDENTLY VERIFIED-EXACT.** Addition and multiplication meet in an exact
two-speed identity. The sum and difference give the two equality clocks;
their product is a Pythagorean leg; its reciprocal is the owner imbalance.
Adding the optimized pair margin recovers the complete primitive ratio.

## 1. Exact pair-owner law

Let `0<m<n`, `gcd(m,n)=1`, and define

```text
d_k(t)=||kt||,
F_(m,n)(t)=min(d_m(t),d_n(t)),
rho(m,n)=lambda{t in R/Z : d_m(t)<d_n(t)},
M(m,n)=max_t F_(m,n)(t).                                (1)
```

Write `G_h={j/h mod 1:0<=j<h}`. Then the complete equality set is

```text
{d_m=d_n}=G_(n-m) union G_(n+m).                        (2)
```

Its cardinality is

```text
2n-1,                         m,n of opposite parity;
2n-2,                         m,n both odd.             (3)
```

The exact owner measure is

```text
rho(m,n)=1/2+1/[2(n^2-m^2)],  m,n of opposite parity;
rho(m,n)=1/2,                 m,n both odd.              (4)
```

Combining `(4)` with THM-4095's exact pair optimum gives the unified identity

```text
boxed: 1/2-M(m,n)=(n-m)(rho(m,n)-1/2).                  (5)
```

More explicitly,

```text
M(m,n)=1/2,                         m,n both odd;
M(m,n)=1/2-1/[2(m+n)],              opposite parity.    (6)
```

In mixed parity, put `q=m+n`. The complete maximizer set consists of the two
reflected points `+-k/q`, where

```text
m k=(q-1)/2 mod q.                                      (7)
```

In the odd--odd case the unique maximizer is `t=1/2`.

### Proof of the owner formula

Cosine is strictly decreasing as a function of circle distance on
`[0,1/2]`, so

```text
d_m(t)<d_n(t)
 iff cos(2 pi m t)>cos(2 pi n t)
 iff sin(pi(n+m)t) sin(pi(n-m)t)>0.                     (8)
```

Equality proves `(2)`. The intersection of the two clock groups has order
`gcd(n-m,n+m)`, which is one in primitive mixed parity and two in primitive
odd parity, proving `(3)`.

On `[0,1]`, let `epsilon_h(t)=sgn(sin(pi h t))`. Its square-wave Fourier series is

```text
epsilon_h(t)=(4/pi) sum_(j>=1, j odd) sin(pi j h t)/j.  (9)
```

If `s=n+m`, `r=n-m`, and `h=gcd(s,r)`, orthogonality gives

```text
integral_0^1 epsilon_s(t)epsilon_r(t) dt
 =h^2/(sr),  if s/h and r/h are both odd;
 =0,         otherwise.                                (10)
```

A term survives precisely when `js=ell r`. After division by `h`, odd
solutions exist exactly when both residual frequencies are odd; summing
`1/k^2` over odd `k` supplies the displayed constant. Since the integral is
`2rho-1`, `(4)` follows. Equations `(5)--(7)` follow from the exact optimum
and its equality conditions.

## 2. The sum/difference and Pythagorean decoder

In mixed parity put

```text
D=1/2-M(m,n),                 I=rho(m,n)-1/2.           (11)
```

Then

```text
q=m+n=1/(2D),
r=n-m=D/I,
m=(q-r)/2,
n=(q+r)/2.                                               (12)
```

Thus `(M,rho)` is a complete invariant of an ordered primitive mixed-parity
pair. Moreover,

```text
qr=n^2-m^2                                               (13)
```

is the first leg of the primitive Pythagorean triple

```text
(n^2-m^2, 2mn, n^2+m^2),                                (14)
```

and

```text
rho-1/2=1/[2 * (first Pythagorean leg)].                 (15)
```

The bridge is literal but sharply scoped: sum and difference create the two
phase clocks, multiplication creates one Pythagorean leg, and the owner
imbalance is its reciprocal. It transfers no Berggren ancestry or common
multi-runner time.

The parity boundary is exact. Every primitive odd pair has
`(M,rho)=(1/2,1/2)`. For example, `{1,3}` and `{1,5}` have the same two
summaries, while

```text
M({1,2,3})=1/4,                 M({1,2,5})=1/3.          (16)
```

Pair locations, not only pair scalars, matter after a third runner is added.

## 3. Common scale and the complete pair geometry

For distinct positive speeds `u<v`, write

```text
g=gcd(u,v),              (m,n)=(u/g,v/g).               (17)
```

Multiplication by `g` is onto on `R/Z`, so the pair maximum and owner measure
are those of `(m,n)`. The scalar summary intentionally loses `g`, but the
full geometry is:

- if the primitive core is mixed, the tie count is `2v-g` and the complete
  maximizer set has `2g` points, the `g` inverse images of each point in
  `(7)`;
- if the primitive core is odd--odd, the tie count is `2v-2g` and the
  complete maximizer set is
  `{(2j+1)/(2g):0<=j<g}`.

Hence a pair summary decodes a labelled primitive ratio exactly when
`v_2(u)!=v_2(v)`. To state this without first ordering the labels, define

```text
rho_(u,v)=lambda{t:||ut||<||vt||},
D=1/2-M(u,v),
I=rho_(u,v)-1/2.                                        (18)
```

For a cross-valuation pair, if `(p,q)=(u/g,v/g)` in the original label order,

```text
p+q=1/(2D),                    q-p=D/I.                 (19)
```

This recovers the labelled ratio `u/v=p/q`. Equal-valuation edges are exactly
the blind odd--odd shell.

## 4. A lossless valuation-tree representation of a row

Let

```text
S=(s_1,...,s_N)                                             (20)
```

be a labelled primitive row of distinct positive integer speeds. If all
speeds are odd, `t=1/2` gives distance `1/2` for every runner, so the row is
already safe for the standard Lonely Runner threshold, and indeed for every
threshold at most `1/2`.

Otherwise, form the graph whose edges join speeds with unequal `2`-adic
valuation. Primitivity supplies at least one odd speed, and the assumed even
speed supplies a second valuation class. The graph is a connected complete
multipartite graph. Choose any spanning tree. Its `N-1` directed summaries

```text
(M(s_i,s_j),rho_(s_i,s_j))                                (21)
```

recover every tree-edge ratio by `(19)`. Propagation recovers the labelled
projective row, and gcd-one normalization recovers `S` uniquely.

For the thirteen nonzero relative speeds in the LRC(14) normalization:

```text
all odd -> explicit antipodal witness;
mixed   -> twelve cross-v2 pair summaries encode the row exactly.          (22)
```

This is a representation theorem, not a loneliness inequality. Pair maxima
generally occur at incompatible times. The summaries discard the common-gcd
sheet and the locations of all lifted maximizers. The remaining anchor
problem is simultaneous signed-residue and arrival synchronization across the
twelve-edge tree.

Zero speeds, duplicate absolute speeds, and signed rows that become duplicate
after absolutization lie outside `(20)`; they must first be reduced to the
standard distinct positive hard row.

## 5. Generated next tasks

1. On each of the 17 live `11+2` types, select a cross-valuation spanning tree
   minimizing a synchronization cost such as `lcm(m_e+n_e)`, maximum clock,
   or `sum log rad(m_e n_e(m_e+n_e))`.
2. Retain, on every tree edge, its two signed maximizers and both equality
   clocks. Test whether a short tree cycle or one non-tree edge forces a
   common `1/14`-safe cell.
3. Enrich the existing support-two short-relation atlas with `(m+n,n-m)`, the
   Pythagorean leg, owner imbalance, and prime colours. A radical address is a
   colour, not a synchronization theorem.
4. Treat AP13 and the pair-summary twins in `(16)` as mandatory hostile
   controls for any compression of the phase locations.

## 6. Exact audit

Reproduce with

```text
python -B 04-computation/lrc_pair_owner_reciprocity_thm4106.py
python -B -O 04-computation/lrc_pair_owner_reciprocity_thm4106.py
```

The normal and optimized streams equal the frozen transcript. Independent
exact paths check all `1,965` primitive pairs `m<n<=80`, `1,770` scaled
pairs `u<v<=60`, `1,600` square-wave correlations, `136` exponent pairs,
and reconstruction of all `1,675` primitive mixed four-speed rows in
`[1,16]`, with `70` all-odd controls. An independent program checks `4,950`
scaled pairs, `24,930` labelled five-row orders, and `3,600` correlations.
LRC(14) remains **OPEN**.

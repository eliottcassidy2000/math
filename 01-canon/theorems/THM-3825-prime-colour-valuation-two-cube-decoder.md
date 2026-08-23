---
id: THM-3825
title: "Prime-colour valuations decode inert-scaled two-cube states"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  On the carrier
  where the primitive pair sum is inert and cube-free and the common scale
  is inert, quotient and remainder of every inert valuation modulo three
  recover the scale and primitive shell; the split valuations recover the
  Eisenstein cofactor.  One square test then recovers the coprime cube pair.
  Multiplication by an inert cube is the native scale law.  A 3-adic tag
  gives unordered and orientation-refined natural-number addresses for the
  5,855 support-two ratios.  These are fixed label-order gauges, not LRC(14)
  certificates; THM-3818's facet and residue packet remains RESERVED.
source: root + lrc_graver_address + cube_decoder_audit / incoming-signal extension session, 2026-08-23
audit: >
  INDEPENDENT HOSTILE AUDIT PASS (cube_decoder_audit, 2026-08-23).  The audit
  rederived the valuation decoder, found the orientation collision in the
  proposed 78-channel covector address, verified the repaired 156-channel
  tag, and froze ramified-prime, split-scale, exponent-three, relabelling,
  and singleton-failure boundaries.  The assertion-free companion has
  2,066,834 active gates; normal and optimized raw LF streams match its
  frozen transcript.
depends_on:
  - THM-463-two-cube-representations-are-a-divisor-property-on-the-split-axis
  - THM-3793-inert-prime-sum-all-scale-two-cube-singleton
related:
  - THM-3743-lonely-runner-polyhedron-khinchin-flatness-relation-reduction
  - THM-3818-scaled-inert-cubeclass-support-two-pair-packet
script: 04-computation/two_cube_prime_colour_decoder_thm3825.py
output: 05-knowledge/results/two_cube_prime_colour_decoder_thm3825.out
script_sha256: a3e2a1950a2c684e55cbb2c15ff3f8ea90d2b7dad042a589c077dc6ef27e971e
output_sha256: e002a15b000ad8bd00011da54ac5b1b5b4f18ffdca20f68af66f5e44970b9e8b
hash_basis: raw LF bytes
---

# THM-3825 -- inert valuations are scale digits plus shell digits

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  This is an
arithmetic decoder theorem and a typed natural-address theorem.  Its finite
LRC corollary does not assert loneliness, arrival, or a bad row.

## 1. Restricted carrier and exact decoder

Let `a<b` be coprime positive integers, put

```text
s=a+b,                       q=a^2-ab+b^2,             (1)
```

and assume

```text
every p|s satisfies p=2 mod 3 and v_p(s)<=2.           (2)
```

Let `g>=1` have only prime divisors congruent to `2 mod 3`, and define

```text
m=g^3(a^3+b^3)=g^3 s q.                               (3)
```

Then the prime factorization of `m` recovers `(g,a,b)` without divisor
enumeration.  For a positive integer `M` with `3` not dividing `M`, set

```text
G(M)=product_(p=2 mod 3) p^floor(v_p(M)/3),
S(M)=product_(p=2 mod 3) p^(v_p(M) mod 3),
Q(M)=product_(p=1 mod 3) p^v_p(M).                    (4)
```

On `(3)`, exactly

```text
(G(m),S(m),Q(m))=(g,s,q).                             (5)
```

Moreover `(4)` gives an exact image test for this restricted carrier.  Put

```text
D=(4Q(M)-S(M)^2)/3.                                  (6)
```

Then `M` has the form `(3)` if and only if `D` is the square `delta^2` of an
integer satisfying

```text
0<delta<S(M),             delta=S(M) mod 2.           (7)
```

In that case the unique decoded state is

```text
g=G(M),
a=(S(M)-delta)/2,          b=(S(M)+delta)/2.           (8)
```

The parity clause in `(7)` is automatic once `(6)` is an integer square, but
it is retained in the executable decoder as a hostile gate.

### Proof

THM-463 applied to the primitive pair `(a,b)` says that `q` is odd, every
prime divisor of `q` is `1 mod 3` except possibly `3`, and

```text
v_3(q)=1 iff 3|s.                                     (9)
```

Hypothesis `(2)` makes `3` divide neither `s` nor `q`; the inert scale `g`
is also 3-free.  Hence the inert and split factors of `(3)` are disjoint and

```text
v_p(m)=3v_p(g)+v_p(s)        for p=2 mod 3.           (10)
```

Because `0<=v_p(s)<=2`, quotient and remainder in `(10)` prove `(5)`.
Finally

```text
4q-s^2=3(b-a)^2,                                      (11)
```

which proves the forward square test.

Conversely, `(4)` gives `M=G^3SQ` by unique factorization.  Equations
`(6)--(8)` give positive distinct integers with sum `S`, and

```text
a^2-ab+b^2=(S^2+3delta^2)/4=Q.                       (12)
```

If a prime divided both `a` and `b`, it would divide `S` and its square
would divide `Q`; this is impossible because the prime supports of `S` and
`Q` are disjoint.  Thus `(a,b)` is primitive, and `(12)` proves `(3)`.
This proves the image test and uniqueness inside the restricted carrier.

THM-3793 supplies the stronger global consequence

```text
r_+(m)=1,                                             (13)
```

so no positive distinct two-cube representation outside the restricted
carrier competes with the decoded one.

## 2. Natural-number tags and their operation law

The bare value `m` is already a sparse injective natural address for
`(g,a,b)`.  To attach a fixed unordered coordinate pair `1<=i<j<=13`, put

```text
kappa(i,j)=binom(j-1,2)+(i-1) in {0,...,77},
A_un=3^kappa(i,j) m.                                  (14)
```

Since `3` does not divide `m`, `v_3(A_un)` recovers `kappa`, hence `{i,j}`;
stripping that power and applying `(4)--(8)` recovers `(g,a,b)`.  At
primitive scale this gives exactly

```text
78*5,855=456,690                                      (15)
```

unordered labelled states.

An unordered coordinate pair does **not** say which coordinate receives the
smaller speed.  Retain an orientation bit `epsilon`, where `epsilon=0`
assigns `(ga,gb)` to `(i,j)` and `epsilon=1` assigns `(gb,ga)`.  The repaired
tag is

```text
lambda=2kappa(i,j)+epsilon in {0,...,155},
A_or=3^lambda m.                                      (16)
```

This recovers the assigned pair and therefore the primitive support-two
covector, up to global sign:

```text
epsilon=0:  b e_i-a e_j,
epsilon=1:  a e_i-b e_j.                              (17)
```

There are `913,380` primitive oriented placements.  Both `(14)` and `(16)`
use the chosen coordinate order as a gauge; neither is invariant under
`S_13` relabelling.

The native scale operation is closed.  If `c` has only inert prime divisors,

```text
(g,a,b)->(cg,a,b)       sends m->c^3m,
A_un->c^3A_un,          A_or->c^3A_or.                (18)
```

Thus the arithmetic address carries inert dilation exactly.  Dense ordinal
successor and Stern--Brocot branching are different operations and are not
asserted to descend through `(14)` or `(16)`.

## 3. Sharp losses and the reserved larger packet

The first orientation hostile is already

```text
A_un=28 on labels {1,2}:
(1,3) gives covector (3,-1),
(3,1) gives covector (1,-3).                           (19)
```

The two covectors are not equal up to sign.  This is why `(16)` needs its
extra bit.  The fixed gauge also sends address `28` on `{1,2}` to `252` on
`{2,3}` under relabelling.

The arithmetic scope gates have small witnesses:

```text
ramified shell:       9=1^3+2^3,
ramified scale/tag:   756=3^3(1^3+3^3),
split scale:          9604=7^3(1^3+3^3),
exponent-three loss:  152=3^3+5^3, shell 8=2^3.       (20)
```

The first high-inert-exponent singleton failure is

```text
65728=12^3+40^3=31^3+33^3,       shell 64=2^6,       (21)
```

and the first failure with primitive exponent exactly three is

```text
515375=15^3+80^3=54^3+71^3,      shell 125=5^3.      (22)
```

THM-3818 remains **RESERVED / UNPROVED**.  Even the oriented arithmetic tag
does not reconstruct the ambient residue schedule, endpoint owner, root,
phase, arrival, other eleven speeds, or a loneliness predicate.  No exposed
facet or complete pair-sum covering theorem is smuggled into this result.

## 4. Exact verification

The assertion-free companion directly checks:

- all `94` admissible shells and `5,855` primitive pairs with `a+b<=356`;
- all inert scales through `100` and both `456,690`/`913,380` tag counts;
- the decoder and orientation round trips; and
- the hostile minima in `(19)--(22)`, including ramified- and split-scale
  singleton failures.

Run

```bash
python3 -B 04-computation/two_cube_prime_colour_decoder_thm3825.py
python3 -B -O 04-computation/two_cube_prime_colour_decoder_thm3825.py
```

Both raw LF streams equal
`05-knowledge/results/two_cube_prime_colour_decoder_thm3825.out`.

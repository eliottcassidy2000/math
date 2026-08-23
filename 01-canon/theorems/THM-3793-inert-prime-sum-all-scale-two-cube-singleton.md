---
id: THM-3793
title: "Inert pair sums with cube-free primitive quotient force two-cube singleton fibres"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  Let x<y be
  positive, g=gcd(x,y), d=x+y, and s=d/g.  If every prime divisor of d is
  congruent to 2 modulo 3 and every v_p(s) is at most two, then x^3+y^3 has
  exactly one positive distinct two-cube representation.  Arbitrary inert
  prime powers may occur in the common scale g.  An explicit two-prime
  subfamily gives H(Z^6)>=(A(Z)^2-B(Z))/5 and
  liminf H(X)/(log log X)^2>=1/20.  More generally, products of every fixed
  number j of distinct inert primes give liminf
  H(X)/(log log X)^j>=2/(5*2^j*j!); therefore H(X) dominates every fixed
  real power of log log X.  No support asymptotic follows.  The finite LRC
  address sidecar is injective but has no loneliness consequence.
source: root / cross_frontier_live_scout inert-prime all-scale lane, 2026-08-23
audit: >
  INDEPENDENTLY HOSTILE-AUDITED by root, 2026-08-23.  The audit rederived the
  valuation invoice with arbitrary nonprimitive scale, checked that the
  exponent cap belongs on the primitive sum, verified the divisibility and
  short-multiple step, reconstructed the unordered-prime mass correction,
  and retained split-prime and exponent-three hostiles.  MISTAKE-450's
  downstream terminology repair now reports 456,690 unoriented supports and
  913,380 oriented assignments separately.  The assertion-free companion
  passes normally and under optimization against its frozen transcript.  A
  second hostile audit independently checked the all-fixed-order product
  amplification, including primitive/nonprimitive shells, empty and
  undersized prime banks, ordered collision union bounds, the j=1 endpoint,
  exact liminf constants and exponent-three hostile.  It caught and repaired
  the provisional undefined 0*B*A^(-1) display at j=1.  The repaired primary
  and independent companions have 14,452 and 5,402 active gates; normal,
  optimized and frozen streams agree after LF normalization.
depends_on:
  - THM-463-two-cube-representations-are-a-divisor-property-on-the-split-axis
  - THM-3730-positive-distinct-two-cube-support-abscissa
related:
  - THM-3743-lonely-runner-polyhedron-khinchin-flatness-relation-reduction
  - THM-3818-scaled-inert-cubeclass-support-two-pair-packet
script: 04-computation/two_cube_inert_cubefree_singleton_thm3793.py
output: 05-knowledge/results/two_cube_inert_cubefree_singleton_thm3793.out
script_sha256: fecde0faf77e919a4d643d2ebede6822fb5b2749b481d2691e43a913a0a84fe9
output_sha256: e40beb9a11f934846dbb44938f443462b14d40fb1fd09af3b89ceeb97e15230c
amplification_script: 04-computation/two_cube_all_fixed_prime_product_amplification_thm3793.py
amplification_output: 05-knowledge/results/two_cube_all_fixed_prime_product_amplification_thm3793.out
amplification_script_sha256: 3b6b8f4118ad4b383328cab7d5adf551fd43dda26e7d65adb4b6870acd47bd2a
amplification_output_sha256: 413ebb4f472d98c691d536c9df2f025066efea4e24c0d245c323c172acb2de54
amplification_semantic_sha256: 083431ef666f4c7a8435da7015cf19671d17e62108f95d89455252283aeb4d41
amplification_independent_script: 04-computation/two_cube_all_fixed_prime_product_amplification_independent_audit_thm3793.py
amplification_independent_output: 05-knowledge/results/two_cube_all_fixed_prime_product_amplification_independent_audit_thm3793.out
amplification_independent_script_sha256: e32535c2fa1efb95c0ad5252387edfb0ac5f1d3a5073b042b3388a1971d40964
amplification_independent_output_sha256: bc0fafc4ae54c9af3edc167e54bd0934a3c8678ac8d812121a3607316c870bcc
amplification_independent_semantic_sha256: ae6aa4dd88588a96cf9f4e4a4ebbce62ac4c51c27d8ce1c8c2cb117563028082
hash_basis: raw LF bytes
---

# THM-3793 -- inert primitive pair sums give all-scale singleton fibres

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  The theorem
below strengthens the singleton subfamily in THM-3730.  It is an arithmetic
support theorem, not an LRC theorem or a support asymptotic.

## 1. Statement

Define

```text
r_+(m)=#{(x,y) in Z_(>0)^2 : x<y and x^3+y^3=m}.
```

### The singleton theorem

Let `x<y` be positive integers and set

```text
g=gcd(x,y),              d=x+y,
s=d/g,                   m=x^3+y^3.
```

Assume, for every prime `p|d`, that

```text
p=2 mod 3,                     v_p(s)<=2.             (1)
```

Thus the entire pair sum has only inert prime factors, while only its
primitive quotient `s` must be cube-free.  The common scale `g` may
contain arbitrary powers.  Then

```text
r_+(m)=1.                                              (2)
```

The quantifiers include nonprimitive pairs: no assumption on `gcd(x,y)` is
made.

### Quadratic critical-mass amplification

Retain THM-3730's deduplicated critical partial sum

```text
H(X)=sum_(m<=X,r_+(m)>0) m^(-2/3).
```

For real `Z>=11`, put

```text
P(Z)={p prime : 5<=p<=Z and p=2 mod 3},
A(Z)=sum_(p in P(Z))1/p,
B(Z)=sum_(p in P(Z))1/p^2.
```

Then

```text
H(Z^6) >= (A(Z)^2-B(Z))/5.                            (3)
```

With the exact constant already proved in THM-3730,

```text
C_P=14/e+1/2+(log 3+2/9)/2,
A(Z)>=1/2 log log Z-C_P,          B(Z)<1/4,           (4)
```

equation (3) gives, for

```text
X>=exp(6 exp(2C_P)),
```

the explicit inequality

```text
H(X) > (1/5)[(1/2(log log X-log 6)-C_P)^2-1/4].       (5)
```

In particular,

```text
liminf_(X->infinity) H(X)/(log log X)^2 >= 1/20.      (6)
```

This is a lower bound from an explicit singleton subfamily.  It asserts no
support counting asymptotic, pole, residue, or collision-tax asymptotic.

### All fixed orders

For each fixed integer `j>=1`, let `e_j(Z)` be the degree-`j` elementary
symmetric sum of the numbers `1/p`, for `p in P(Z)`.  Then

```text
H(Z^(3j)) >= (2/5)e_j(Z).                             (6a)
```

For `j=1` this reads `H(Z^3)>=(2/5)A(Z)`.  For `j>=2`,

```text
H(Z^(3j))
 >=2/(5j!)[A(Z)^j-C(j,2)B(Z)A(Z)^(j-2)].             (6b)
```

Consequently

```text
liminf_(X->infinity) H(X)/(log log X)^j
 >=2/(5*2^j*j!)                                      (6c)
```

for every fixed integer `j>=1`.  In particular, for every fixed real `R`,

```text
H(X)/(log log X)^R -> infinity.                       (6d)
```

The quantifier is pointwise in `R`; no choice `j=j(X)` or uniform-in-`j`
claim is made.

## 2. Proof of the singleton theorem

Write

```text
x=gX,        y=gY,        s=X+Y,
q=X^2-XY+Y^2.
```

Fix a prime `p|d`.  By (1), `p=2 mod 3`.  THM-463's primitive split lemma
says that `p` cannot divide `q`; for `p=2`, this is its separate assertion
that the primitive cofactor is odd.  Put

```text
alpha=v_p(g),                    e=v_p(s).
```

Then `v_p(d)=alpha+e`, and hypothesis (1) says `e<=2`.

Since

```text
m=g^3sq,                         d=gs,
```

we have the exact valuation identity

```text
v_p(m)=3alpha+e.                                      (7)
```

Now take any competing positive distinct representation

```text
m=u^3+v^3,
```

and write

```text
h=gcd(u,v),       u=hU,        v=hV,
beta=v_p(h),      d'=u+v.
```

Applying the same primitive split lemma to `(U,V)` gives

```text
v_p(m)=3beta+v_p(U+V),
v_p(d')=beta+v_p(U+V)=v_p(m)-2beta.                  (8)
```

Because `h^3|m`,

```text
3beta<=v_p(m)=3alpha+e.                              (9)
```

If `beta>=alpha+1`, then (9) would imply

```text
3alpha+3<=3alpha+e,
3<=e,
```

which is impossible because `e<=2`.  Therefore `beta<=alpha`.  Combining
(7)--(8),

```text
v_p(d')=3alpha+e-2beta
       =(alpha+e)+2(alpha-beta)>=v_p(d).             (10)
```

Equation (10) holds for every prime power dividing `d`; hence

```text
d|d'.                                                 (11)
```

For the competing positive distinct pair,

```text
4m=d'(d'^2+3(v-u)^2)>d'^3.                           (12)
```

For the original positive pair,

```text
m=d(x^2-xy+y^2)=d(d^2-3xy)<d^3.                     (13)
```

Thus

```text
0<d'<(4m)^(1/3)<4^(1/3)d<2d.                        (14)
```

The only positive multiple of `d` in (14) is `d'=d`.  Finally, `d` and `m`
determine

```text
uv=(d^3-m)/(3d),
```

so `{u,v}` is the unique root pair of the resulting quadratic.  It equals
`{x,y}`.  This proves (2).

## 3. Proof of the critical-mass bound

Choose an **unordered** pair of distinct primes `p<q` in `P(Z)` and set

```text
d=pq.
```

For any pair `x,d-x`, its gcd divides `d=pq`; hence its primitive
sum is one of `pq,p,q,1`.  Thus the integer `d` and every row pair
satisfy (1).  Since it is odd, there are exactly

```text
(d-1)/2
```

**unordered positive distinct** pairs `{x,d-x}`.  The singleton theorem makes
all of their cube sums distinct within the row.  It also makes rows with
different `d` disjoint: equality of two values would exhibit two different
positive representations.

Every value in the row satisfies

```text
m<d^3<=(pq)^3<=Z^6,
m^(-2/3)>d^(-2).
```

Since the smallest distinct pair is `(p,q)=(5,11)`, every such `d` is at
least `55`, and

```text
sum_(1<=x<d/2)(x^3+(d-x)^3)^(-2/3)
  >(d-1)/(2d^2)
  >=2/(5d).                                           (15)
```

Summing (15) over **unordered** prime pairs gives

```text
H(Z^6)
 >=(2/5)sum_(p<q in P(Z))1/(pq)
 =(2/5)(A(Z)^2-B(Z))/2
 =(A(Z)^2-B(Z))/5.                                   (16)
```

The middle factor `1/2` is the ordered-to-unordered correction; no prime-pair
or cube-pair multiplicity remains hidden.  This proves (3).

For `X>=exp(6 exp(2C_P))`, take `Z=X^(1/6)`.  Then

```text
log log Z=log log X-log 6>=2C_P,
```

so the lower bound for `A(Z)` in (4) is nonnegative and may be squared.
Substitution into (16), followed by `B(Z)<1/4`, proves (5).  Dividing by
`(log log X)^2` and taking the lower limit gives (6).

## 4. Proof of the all-fixed-order amplification

Choose an unordered set of `j` distinct primes from `P(Z)` and put

```text
d=product_(p in S)p.
```

For every `1<=x<d/2`, its common scale is `g=gcd(x,d)` and its primitive
pair sum is `d/g`, which is squarefree with only inert prime divisors.
The singleton theorem applies.  It makes every value in the row distinct
and also makes rows from different prime products disjoint.  Moreover

```text
x^3+(d-x)^3<d^3<=Z^(3j),
sum_(1<=x<d/2)(x^3+(d-x)^3)^(-2/3)
 >(d-1)/(2d^2)>=2/(5d).                               (17)
```

Summing `(17)` over unordered prime sets proves `(6a)`.  Expand `A^j` as
the mass of all ordered `j`-tuples.  The distinct tuples have mass `j!e_j`.
For each pair of positions, the collision mass is `BA^(j-2)`; the union
bound therefore gives, for `j>=2`,

```text
j!e_j>=A^j-C(j,2)BA^(j-2),                            (18)
```

which proves `(6b)`.  The separate `j=1` identity avoids any formal
`0*B*A^(-1)` at an empty prime bank.

Fix `j`, put `Z=X^(1/(3j))`, and use

```text
log log Z=log log X-log(3j),
A(Z)>=1/2 log log Z-C_P,          B(Z)<1/4.           (19)
```

The collision correction in `(18)` is
`A^j[1-C(j,2)B/A^2]`, whose bracket tends to one.  This proves `(6c)`.
Given a fixed real `R`, choose a fixed integer `j>R` before taking the
limit; `(6c)` then proves `(6d)`.

## 5. Exact boundary hostiles

The split-prime condition cannot be discarded:

```text
1729=9^3+10^3=1^3+12^3,
9+10=19,                       1+12=13,
13=19=1 mod 3.
```

The primitive-sum exponent cap in (1) is sharp already at exponent three:

```text
515375=54^3+71^3=15^3+80^3,
54+71=125=5^3.
```

A smaller-valued high-power control is

```text
65728=31^3+33^3=12^3+40^3,
31+33=64=2^6.
```

The theorem's condition is sufficient, not necessary.  It makes no converse
claim and does not classify representations when `3`, a split prime, or an
inert exponent at least three enters the primitive pair-sum.

## 6. Finite support-two LRC address sidecar

This paragraph is a typed finite corollary, not an LRC theorem.  THM-3743's
support-two branch has `19,314` unordered coprime coefficient ratios
`a<b`, `a+b<=356`.  Restrict to ratios whose sum satisfies (1).  The exact
census is

```text
94 admissible sums,
5,855 unordered coprime ratios,
456,690 unoriented supports after C(13,2)=78 choices,
913,380 oriented assignments after 13*12 choices.
```

On this subset, `(a,b)->a^3+b^3` is injective and the singleton divisor fibre
recovers `(a,b)`, hence its ordered Christoffel word.  It does **not** preserve
the other speeds, gcd scale, owner, phase, arrival, or the loneliness
predicate.  It excludes no LRC(14) row.

## 7. Exact verification

The assertion-free companion performs three complementary finite views:

1. it precomputes every positive coordinate fibre capable of representing a
   value from pair-sum at most `1000`, then checks all `57,829` values from
   all `243` rows satisfying the stronger total-sum exponent cap are
   singleton and mutually disjoint;
2. it separately checks the theorem's primitive-sum condition for every
   pair-sum at most `1000`: all `61,434` values are singleton and mutually
   disjoint, including `3,605` values excluded by the total-sum cap; and
3. it independently checks the valuation invoice (7) prime by prime, all
   three hostiles, the complete `l1<=356` LRC ratio census, and exact rational
   ordered/unordered identities through inert primes at most `5000`.

The all-fixed-order primary and independent companions additionally check
the symmetric collision inequality, empty and undersized banks, complete
small primitive and nonprimitive fibres, the height conversion, every
constant through `j=12`, and the exponent-three hostile.

Run

```text
python -B 04-computation/two_cube_inert_cubefree_singleton_thm3793.py
python -B -O 04-computation/two_cube_inert_cubefree_singleton_thm3793.py
python -B 04-computation/two_cube_all_fixed_prime_product_amplification_thm3793.py
python -B -O 04-computation/two_cube_all_fixed_prime_product_amplification_thm3793.py
python -B 04-computation/two_cube_all_fixed_prime_product_amplification_independent_audit_thm3793.py
python -B -O 04-computation/two_cube_all_fixed_prime_product_amplification_independent_audit_thm3793.py
```

and compare with

```text
05-knowledge/results/two_cube_inert_cubefree_singleton_thm3793.out
05-knowledge/results/two_cube_all_fixed_prime_product_amplification_thm3793.out
05-knowledge/results/two_cube_all_fixed_prime_product_amplification_independent_audit_thm3793.out
```

The normal and optimized streams line-normalize exactly to the frozen
transcript.  Raw SHA-256 of the stored files:

```text
fecde0faf77e919a4d643d2ebede6822fb5b2749b481d2691e43a913a0a84fe9
  two_cube_inert_cubefree_singleton_thm3793.py
e40beb9a11f934846dbb44938f443462b14d40fb1fd09af3b89ceeb97e15230c
  two_cube_inert_cubefree_singleton_thm3793.out
3b6b8f4118ad4b383328cab7d5adf551fd43dda26e7d65adb4b6870acd47bd2a
  two_cube_all_fixed_prime_product_amplification_thm3793.py
413ebb4f472d98c691d536c9df2f025066efea4e24c0d245c323c172acb2de54
  two_cube_all_fixed_prime_product_amplification_thm3793.out
e32535c2fa1efb95c0ad5252387edfb0ac5f1d3a5073b042b3388a1971d40964
  two_cube_all_fixed_prime_product_amplification_independent_audit_thm3793.py
bc0fafc4ae54c9af3edc167e54bd0934a3c8678ac8d812121a3607316c870bcc
  two_cube_all_fixed_prime_product_amplification_independent_audit_thm3793.out
```

The proof, not the finite range, carries the all-scale quantifier.

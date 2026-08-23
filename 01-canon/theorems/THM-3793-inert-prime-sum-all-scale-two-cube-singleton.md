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
  liminf H(X)/(log log X)^2>=1/20.  The finite LRC address sidecar is
  injective but has no loneliness consequence.
source: root / cross_frontier_live_scout inert-prime all-scale lane, 2026-08-23
audit: >
  INDEPENDENTLY HOSTILE-AUDITED by root, 2026-08-23.  The audit rederived the
  valuation invoice with arbitrary nonprimitive scale, checked that the
  exponent cap belongs on the primitive sum, verified the divisibility and
  short-multiple step, reconstructed the unordered-prime mass correction,
  and retained split-prime and exponent-three hostiles.  MISTAKE-450's
  downstream terminology repair now reports 456,690 unoriented supports and
  913,380 oriented assignments separately.  The assertion-free companion
  passes normally and under optimization against its frozen transcript.
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

## 4. Exact boundary hostiles

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

## 5. Finite support-two LRC address sidecar

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

## 6. Exact verification

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

Run

```text
python -B 04-computation/two_cube_inert_cubefree_singleton_thm3793.py
python -B -O .scratch/cross_frontier_live_scout_20260823/two_cube_inert_cubefree_singleton_probe.py
```

and compare with

```text
05-knowledge/results/two_cube_inert_cubefree_singleton_thm3793.out
```

The normal and optimized streams line-normalize exactly to the frozen
transcript.  Raw SHA-256 of the stored files:

```text
4b5e7052d0698a86e43aab4f6c99619aa6014abe2315456122ce02ba6db1adea
  two_cube_inert_cubefree_singleton_thm3793.py
3c4d8cd4d375ea823bc6649ce79730ff4e78b04179884034a1c82d13cd37459a
  two_cube_inert_cubefree_singleton_thm3793.out
```

The proof, not the finite range, carries the all-scale quantifier.

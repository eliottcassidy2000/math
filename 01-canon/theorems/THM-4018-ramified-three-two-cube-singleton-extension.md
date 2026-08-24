---
id: THM-4018
title: "Ramified-three extension of the two-cube singleton criterion"
status: >
  PROVED + CITED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  Let 0<x<y,
  g=gcd(x,y), d=x+y, and s=d/g.  If every prime divisor of d is either 3 or
  congruent to 2 modulo 3, with v_3(s)<=1 and v_p(s)<=2 at every inert prime,
  then x^3+y^3 has exactly one positive distinct two-cube representation.
  The cap at 3 is sharp: 4104=2^3+16^3=9^3+15^3 has v_3(18/gcd(2,16))=2.
  Pairing the squarefree-inert rows d and 3d sharpens THM-3793 to
  liminf H(X)sqrt(log log X)/sqrt(log X) at least
  (2/3)sqrt(2/3)E_32=0.4276598..., exactly 5/3 of its recorded coefficient.
  No support asymptotic or critical residue follows.
source: all-frontiers two-cube ramified-prime lane + audit_thm4018_ramified_three, 2026-08-24
audit: >
  PASS.  The primary verifier checks every admissible pair sum through 2,500
  against the complete height-forced competitor-sum universe, the primitive
  3-adic law, the first relaxed hostile, and an exact two-row Euler bank at
  815,204 gates.  An import-independent verifier instead enumerates every
  competitor by its coordinate cube bound, rederives all local valuation
  branches, checks 510,009 strict values with no hostile and 547,824 relaxed
  values with 172 hostile values, proves 4104 globally least, audits exact
  row disjointness and the all-real X/27 normalization, and passes 446,309
  gates.  Both normal/optimized pairs match their frozen LF outputs.
depends_on:
  - THM-463-two-cube-representations-are-a-divisor-property-on-the-split-axis
  - THM-3730-positive-distinct-two-cube-support-abscissa
  - THM-3793-inert-prime-sum-all-scale-two-cube-singleton
related:
  - THM-3825-prime-colour-valuation-two-cube-decoder
references:
  - "K. S. Williams, Mertens' Theorem for Arithmetic Progressions, Journal of Number Theory 6 (1974), 353--359, DOI 10.1016/0022-314X(74)90032-8."
script: 04-computation/two_cube_ramified_three_singleton_thm4018.py
output: 05-knowledge/results/two_cube_ramified_three_singleton_thm4018.out
independent_script: 04-computation/two_cube_ramified_three_singleton_independent_audit_thm4018.py
independent_output: 05-knowledge/results/two_cube_ramified_three_singleton_independent_audit_thm4018.out
script_sha256: 534cccdd382bc26212dd7229c16035adb804b2951c060492c64624c2bb6ce084
output_sha256: 35bbef26a27304a8b42049583de23fbfb7b900abe80377e5245e3a1f40db8c4c
independent_script_sha256: 455570b7dae6a23e07b3fa87d32476948915bc179e10bd5f0fb368ef5fee97b3
independent_output_sha256: 1b56b8ce29d0c21573ba15c9a44848bc60eb9a918499e2544fa9a37669dbd367
semantic_sha256: b59f0e8e093c8a98f77c778dfcf11d87667c8d917673ed623308f5eb43763119
independent_semantic_sha256: 9b6ae2ff591739574eb65ad52808832656d0fd07ba961fdbbe24ad9a06e8e005
hash_basis: raw LF bytes
---

# THM-4018 -- the ramified-three singleton extension

**PROVED + CITED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  This is a
sufficient all-scale singleton criterion.  It does not classify all two-cube
collisions involving three.

## 1. Statement

Define

```text
r_+(m)=#{(x,y) in Z_(>0)^2:x<y and x^3+y^3=m}.
```

Let `0<x<y` be integers and put

```text
g=gcd(x,y),             d=x+y,
s=d/g,                  m=x^3+y^3.                    (1)
```

Assume every prime divisor of `d` is either the ramified prime three or is
congruent to two modulo three, and impose the primitive-sum caps

```text
v_3(s)<=1,
v_p(s)<=2 for every p|d with p=2 mod 3.                (2)
```

Then

```text
r_+(m)=1.                                              (3)
```

Arbitrary powers of the allowed primes may occur in the common scale `g`.

The new cap at three is sharp:

```text
4104=2^3+16^3=9^3+15^3,
d=18, g=2, s=9 for (2,16).                             (4)
```

Thus the exponent cap cannot be raised uniformly from one to two.

## 2. Inert-prime invoice

Take a competing positive distinct representation

```text
m=u^3+v^3,
h=gcd(u,v), u=hU, v=hV, d'=u+v.                        (5)
```

Fix an inert prime `p|d), including `p=2), and write

```text
alpha=v_p(g), e=v_p(s),
beta=v_p(h),  e'=v_p(U+V).                             (6)
```

THM-463 says the primitive Eisenstein cofactor is prime to `p).  Hence

```text
v_p(m)=3alpha+e=3beta+e'.                              (7)
```

Since `h^3|m`, if `beta>=alpha+1` then `3<=e`, contrary to (2).
Therefore `beta<=alpha), and

```text
v_p(d')=beta+e'
       =3alpha+e-2beta
       >=alpha+e=v_p(d).                               (8)
```

## 3. The ramified-prime invoice

For coprime positive `X,Y`, THM-463 gives

```text
v_3(X^2-XY+Y^2)=
  0 if 3 does not divide X+Y,
  1 if 3 divides X+Y.                                  (9)
```

Put `alpha=v_3(g)`, `e=v_3(s)`, `beta=v_3(h)`, and
`e'=v_3(U+V)`.

If `e=0`, then `v_3(m)=3alpha`.  When `e'=0`, one has
`beta=alpha`.  When `e'>0`, equation (9) gives

```text
3beta+e'+1=3alpha,
e'=3(alpha-beta)-1.                                   (10)
```

Here `delta=alpha-beta>=1`, so

```text
v_3(d')-v_3(d)=beta+e'-alpha=2delta-1>=0.              (11)
```

If `e=1`, then `v_3(m)=3alpha+2).  A competitor with `e'=0` would
force `3beta=3alpha+2` and is impossible.  Thus `e'>0), and

```text
e'=3(alpha-beta)+1,
v_3(d')-v_3(d)=2(alpha-beta)>=0.                       (12)
```

Equations (8), (11), and (12) for every prime divisor of `d` imply

```text
d|d'.                                                  (13)
```

## 4. Short-multiple closure

For every positive distinct representation,

```text
4m=d'(d'^2+3(v-u)^2)>d'^3,                             (14)
```

whereas the original pair gives `m<d^3).  Consequently

```text
0<d'<(4m)^(1/3)<4^(1/3)d<2d.                          (15)
```

The only positive multiple of `d` in this interval is `d'=d`.  The common
pair sum and cube sum determine

```text
uv=(d^3-m)/(3d),                                       (16)
```

so `{u,v}` and `{x,y}` are the roots of the same quadratic.  This proves
(3).

## 5. Exact hostile boundary

The two independent finite oracles cover the same candidate rows in different
ways.  They find

```text
strict rows=1087,
strict values=510009,
strict external hostiles=0,
relaxed values=547824,
relaxed hostile values=172,
globally first relaxed hostile=4104.                   (17)
```

The second oracle proves global minimality of 4104 because every
representation of a smaller value has pair sum below 26, well inside the
candidate universe.  In the forbidden `e=2` branch, (4) realizes exactly the
escape absent from (12): the competing common-scale exponent can rise and
`v_3(d')<v_3(d)).

## 6. Euler-product critical-mass sharpening

Retain THM-3793's notation

```text
H(X)=sum_(m<=X,r_+(m)>0)m^(-2/3),
P(Z)={p prime:5<=p<=Z, p=2 mod 3},
e_j(Z)=sum_(|S|=j)1/product_(p in S)p.                 (18)
```

For each nonempty prime set `S`, put
`d=product_(p in S)p`.  Every unordered positive distinct pair in the rows
of pair sum `d` and `3d` satisfies (2).  The singleton theorem makes all
values globally disjoint, both between these two rows and across different
subsets and orders.

Both pair sums are odd.  The exact row counts and the lower bound
`m^(-2/3)>(pair sum)^(-2)` give the combined surrogate

```text
(d-1)/(2d^2)+(3d-1)/(18d^2)
 =2/(3d)-5/(9d^2).                                    (19)
```

Every value is below `27d^3`.  Therefore, for every real `Z>=11` and
integer `J>=1`,

```text
H(27Z^(3J))
 >=(2/3)sum_(1<=j<=J)e_j({1/p})
   -(5/9)sum_(1<=j<=J)e_j({1/p^2}).                   (20)
```

The negative term is uniformly bounded because

```text
sum_(j>=1)e_j({1/p^2})
 =product_(p in P(Z))(1+1/p^2)-1
 <=product_p(1+1/p^2)-1<infinity.                     (21)
```

THM-3793, using Williams' fixed-modulus Mertens theorem, proves

```text
E(Z)=product_(p in P(Z))(1+1/p)
    =E_32 sqrt(log Z)(1+o(1))                          (22)
```

and that its floor-sensitive first `J` Bernoulli-normalized layers contain
`E(Z)(1-o(1))).  For the all-real limit here, set

```text
Y=X/27,
L=log log Y,
J=floor(L/2-(1/2)log L+L^(2/3)),
Z=Y^(1/(3J)).                                          (23)
```

Then `27Z^(3J)=X`.  The factor 27 changes neither limiting logarithmic
ratio, while (21) becomes negligible after normalization.  Equations
(20)--(23) yield

```text
liminf_(X->infinity)
 H(X)sqrt(log log X)/sqrt(log X)
 >=(2/3)sqrt(2/3)E_32
 =0.4276598... .                                       (24)
```

This is exactly `5/3` times THM-3793's coefficient
`(2/5)sqrt(2/3)E_32`.  The gain is the sum of two transparent changes: the
exact inert row contributes asymptotic coefficient `1/2` rather than the
uniform `2/5` surrogate, and the ramified row contributes `1/6`.

Equation (24) is a lower bound for the deduplicated critical mass.  It is not
a support asymptotic, a residue, or a collision-tax law.

## 7. Reproduction

Run

```text
python -B 04-computation/two_cube_ramified_three_singleton_thm4018.py
python -B -O 04-computation/two_cube_ramified_three_singleton_thm4018.py
python -B 04-computation/two_cube_ramified_three_singleton_independent_audit_thm4018.py
python -B -O 04-computation/two_cube_ramified_three_singleton_independent_audit_thm4018.py
```

and compare with the frozen outputs named in the frontmatter.  The proof, not
the finite ranges, carries the all-scale quantifier. **QED.**

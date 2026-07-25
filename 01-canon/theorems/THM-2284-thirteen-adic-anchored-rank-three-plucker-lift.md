---
id: THM-2284
title: "Thirteen-adic anchored rank-three Plucker lift"
status: >
  PROVED + INDEPENDENTLY AUDITED. A general centered p-adic
  pivot-extension lemma turns a rationally new relation into a bounded row
  extending any fixed mod-p relation pivot. Applied to THM-2282 and
  THM-2283, every one of the 120 interior first-depth-one scalar profiles
  has three relations of heights at most 9841, 4921, and 7381 whose
  three-by-three coefficient minor is nonzero modulo thirteen, contains
  c_1, and is bounded by 54991358114. The fixed-section original minor is
  at most 109982716228. The terminal pivot is either blocker-complete or
  guard/unit-marked. It supplies a full mod-13^n address bank but no prescribed
  Fourier lift, ancestry return, profile exclusion, or proof of LRC(14).
source: codex-2026-07-25-thirteen-adic-rank-three-lift
depends_on:
  - THM-2203-fixed-dyadic-coordinate-section-and-covector-intersection
  - THM-2282-thirteen-adic-saturation-and-unit-anchored-minor
  - THM-2283-mixed-rank-two-safe-torus-floor-and-scalar-rank-three-harvest
related:
  - THM-2069-k-deletion-code-cogirth-crt-wheel
  - THM-2269-marked-expiration-root-spectrum-and-branch-state-no-go
  - THM-2278-two-shallow-proper-root-spectrum-and-gap-ancestry-activation
---

# THM-2284 -- thirteen-adic anchored rank-three Plucker lift

**PROVED + INDEPENDENTLY AUDITED.**

THM-2282 gives a bounded rank-two relation plane with a
`c_1`-anchored minor that is a unit modulo thirteen. THM-2283 independently
forces three-dimensional scalar relation rank at height `3540`. The two
statements do not immediately combine: the third short relation may reduce
modulo thirteen into the old plane.

The repair is a centered carry operation:

```text
new row lies in the old mod-p plane
  -> cancel its centered residue address
  -> divide the complete relation by p
  -> retain rational independence
  -> repeat.                                                   (1)
```

The operation works at every odd prime. Its height fixed point depends only
on half the sum of the old row heights, not on the prime.

## 1. Centered p-adic pivot extension

We first isolate the general lemma.

> **Pivot-extension lemma.** Let `p` be an odd prime and let
>
> ```text
> Lambda subset Z^n
> ```
>
> be a saturated sublattice. Suppose
>
> ```text
> rho_1,...,rho_d in Lambda
> ```
>
> have row rank `d` over `F_p`, with
>
> ```text
> ||rho_i||_infinity<=H_i.                              (2)
> ```
>
> If
>
> ```text
> s_0 in Lambda minus span_Q(rho_1,...,rho_d),
> ||s_0||_infinity<=S_0,                                (3)
> ```
>
> then there is
>
> ```text
> s_* in Lambda minus span_Q(rho_1,...,rho_d)           (4)
> ```
>
> such that the reductions of
>
> ```text
> rho_1,...,rho_d,s_*
> ```
>
> have row rank `d+1` over `F_p` and
>
> ```text
> ||s_*||_infinity
>   <=max(S_0,ceil((H_1+...+H_d)/2)).                   (5)
> ```

### Proof

Start with `s=s_0`. If its reduction extends the old row rank, stop.
Otherwise there are residues `a_i mod p` such that

```text
s=sum_i a_i rho_i mod p.                               (6)
```

Choose the centered representatives

```text
-(p-1)/2<=a_i<=(p-1)/2                                (7)
```

and put

```text
s'=(s-sum_i a_i rho_i)/p.                             (8)
```

The numerator is coordinatewise divisible by `p`. It belongs to `Lambda`,
and saturation of `Lambda` in `Z^n` therefore gives

```text
s' in Lambda.                                         (9)
```

It remains rationally new. Indeed,

```text
s' in span_Q(rho_1,...,rho_d)
  -> s=p s'+sum_i a_i rho_i
       in span_Q(rho_1,...,rho_d),                     (10)
```

contrary to the induction hypothesis.

Write

```text
G=Lambda/(Z rho_1+...+Z rho_d).                       (11)
```

The class of `s_0` has nonzero image in the free part of the finitely
generated abelian group `G`, because (3) says it is not torsion. Every
descent step gives

```text
[s]=p[s'] in G.                                       (12)
```

If the process never stopped, the nonzero free component of `[s_0]` would
belong to

```text
intersection_(m>=1) p^m Z^r={0},                      (13)
```

which is impossible. Notice that the fixed row lattice need not be
primitive. Prime-to-`p` torsion in `G` is harmless because the class in
(3) has a nonzero free component.

Finally, if `S=||s||_infinity`, equations (2), (7), and (8) give

```text
||s'||_infinity
 <=(S+((p-1)/2)sum_i H_i)/p.                         (14)
```

Every

```text
M>=max(S_0,ceil((sum_i H_i)/2))                       (15)
```

is invariant under (14), since

```text
M+((p-1)/2)sum_i H_i<=pM.                            (16)
```

Termination and (15) prove the lemma. QED.

This is the relation-to-carry operator behind the theorem. The cancellation
radius `(p-1)/2` and the division by `p` cancel exactly in the fixed point.

## 2. The scalar input

Use the scalar row

```text
w_*=(H,q_1,...,q_5,c_1,c_2,c_3)                     (17)
```

in one of the `120` interior first-depth-one profiles

```text
c_1=13u_1,       c_2=13^b u_2,       c_3=13^c u_3,

3<=b<=c-2,                    5<=c<=19.              (18)
```

Let

```text
Lambda_*={x in Z^9:x.w_*=0}.                         (19)
```

This kernel is saturated in `Z^9`.

THM-2282 supplies independent relations

```text
p,r in Lambda_*                                      (20)
```

with

```text
||p||_infinity<=9841,
||r||_infinity<=4921.                                (21)
```

There are labels `z` and `k` such that

```text
supp(p)={c_1,z},
k notin {c_1,z},
13 divides p_z,
13 does not divide p_(c_1)r_k,
|p_(c_1)|<=757,
p_k=0.                                               (22)
```

Here `z` is either `H` or one of the `q_i`, and its scalar value is a
thirteen-unit. Thus the two-by-two minor on `(c_1,k)` is

```text
Delta_2=p_(c_1)r_k!=0 mod 13.                        (23)
```

THM-2283 proves

```text
dim_Q span{x in Lambda_*:||x||_infinity<=3540}>=3.   (24)
```

If every relation in (24) lay in `span_Q(p,r)`, their span would have
dimension at most two. Hence there is

```text
s_0 in Lambda_* minus span_Q(p,r),
||s_0||_infinity<=3540.                              (25)
```

## 3. The third row and its exact height

Apply the pivot-extension lemma with

```text
p=13,       d=2,
rho_1=p,    rho_2=r.                                 (26)
```

The notation overload between the prime and the relation row is confined to
this display; below `p` again denotes the relation row.

Equations (21) and (25) give the invariant

```text
M=max(3540,ceil((9841+4921)/2))
 =7381.                                              (27)
```

Indeed the worst recurrence closes with equality:

```text
7381+6(9841+4921)=13*7381.                           (28)
```

We obtain a relation

```text
s_* in Lambda_* minus span_Q(p,r),
||s_*||_infinity<=7381,                              (29)
```

such that

```text
rank_(F_13)(p,r,s_*)=3.                              (30)
```

The proof uses the free component of

```text
Lambda_*/(Zp+Zr),                                    (31)
```

not the false assertion that this raw quotient must be torsion-free. The
unit minor (23) excludes thirteen-primary saturation debt, but it need not
exclude torsion at other primes.

## 4. Extending the anchored pivot

Over `F_13`, first eliminate the `c_1` entries of both `r` and `s_*` with
`p`. The reduced `r` still has the unit `k` entry `r_k`, because `p_k=0`.
Use this reduced row to eliminate the `k` entry of the reduced `s_*`, and
call the result `t`. These invertible row operations preserve all relevant
minors. Equations (23) and (30) give

```text
t_(c_1)=t_k=0,
t!=0.                                                (32)
```

The row `t` is still a relation modulo thirteen. It cannot be supported
only at `z`: otherwise

```text
0=t.w_*=t_z z mod 13,                                (33)
```

and the thirteen-unit value of `z` would force `t_z=0`, contradicting
(32). Therefore some

```text
l notin {c_1,k,z}                                    (34)
```

satisfies `t_l!=0`. The three-by-three minor on `(c_1,k,l)` is consequently
a thirteen-unit:

```text
Delta_3
 =det [[p_(c_1),p_k,p_l],
       [r_(c_1),r_k,r_l],
       [(s_*)_(c_1),(s_*)_k,(s_*)_l]]
 !=0 mod 13.                                         (35)
```

Because `p` is supported on `{c_1,z}`, equations (22) and (34) give

```text
p_k=p_l=0,

Delta_3=p_(c_1)(r_k(s_*)_l-r_l(s_*)_k).             (36)
```

The sharp inherited bound is

```text
0<|Delta_3|
 <=2*757*4921*7381
 =54991358114.                                       (37)
```

Both `k` and `l` lie in

```text
B=[9] minus {c_1,z}.                                 (38)
```

There are only two other blocker labels. Hence the pivot has the exact
dichotomy

```text
{k,l}={c_2,c_3}
  -> the pivot is blocker-complete on (c_1,c_2,c_3);

otherwise
  -> at least one of k,l is H or a q_j,
     so the pivot is guard/unit-marked.              (39)
```

This is a label dichotomy, not an ancestry return.

## 5. Smith and mod-13^n address sidecars

Let

```text
U=Zp+Zr+Zs_* subset Z^9.                             (40)
```

For a rank-three integer row lattice, the gcd of the three-by-three minors
is the index of `U` in its rational saturation. Equation (35) therefore
shows that this index is prime to thirteen. Equivalently:

```text
U has no thirteen-primary Smith obstruction.         (41)
```

For every `n>=1`, the pivot in (35) remains invertible modulo `13^n`.
Consequently the coefficient map

```text
(Z/13^n Z)^3 -> (Z/13^n Z)^9,
(a,b,c) |-> ap+br+cs_*                               (42)
```

is injective. It gives exactly `13^(3n)` distinct relation addresses.

Choose centered representatives

```text
|a|,|b|,|c|<=(13^n-1)/2.                             (43)
```

Using (21) and (29), every address in (42) has an integer representative
of height at most

```text
B_n
 =((9841+4921+7381)(13^n-1))/2
 =22143(13^n-1)/2.                                   (44)
```

Thus the result is not merely a mod-thirteen rank statement: it supplies a
lossless finite address bank at every thirteen-adic depth. It still does
not identify an integer Fourier lift within any residue class.

## 6. Fixed-section original-row lift

THM-2203 lifts scalar relation coefficients by

```text
(x_H,x_rest) |->(2x_H,x_rest)                        (45)
```

on the fixed nine-coordinate original-row section. This map is injective
and preserves relation rank. Among the three pivot columns in (35), at most
one is `H`, so the determinant is multiplied by either one or two. Hence

```text
Delta_3^original!=0 mod 13,

0<|Delta_3^original|
 <=2*54991358114
 =109982716228.                                      (46)
```

The lifted rows have heights at most

```text
2*9841=19682,
2*4921=9842,
2*7381=14762.                                        (47)
```

## 7. Exact-frequency stopping boundaries

The outside labels in (39) cannot be prescribed. For example, on the toy
scalar row

```text
(1,1,2,3,4,5,13,169,2197),
```

the relations

```text
13e_H-e_(c_1),
e_(q_1)-e_H,
e_(c_2)-169e_H                                      (48)
```

have a unit pivot on `(c_1,q_1,c_2)`, and the third direction lands on the
deep blocker `c_2`. This is an algebraic pivot control only; its depth-two
`c_2` lies outside the live profile range (18), and no cover claim is made.

More importantly, even a named nonzero root character does not select an
integer Fourier lift. Let `F` be a centered interval of length `1/14`
inside the unit danger arc `D_1`. Since

```text
1/14<1/13,
```

the map `T(x)=13x mod 1` is injective on `F`. Every point of `T(F)` has
exactly one active root sheet, so every nonzero `C_13` root character has
amplitude one. Nevertheless, for `f=1_F`,

```text
f_hat(14)
 =phase*sin(pi*14/14)/(14pi)
 =0,                                                 (49)
```

even though

```text
14=1 mod 13.                                         (50)
```

A length-`1/q` interval gives the same obstruction at a prescribed lift
`q` whenever `q>13` and `13` does not divide `q`; these hypotheses retain
one-sheet injection and a nonzero root character. Thus a unit Plucker
coordinate, the full address bank (42), or local
root-character activation cannot by themselves force the exact coefficient
needed by an owner crossing. An integer-lift/ancestry sidecar or a signed
whole-residue-class argument is genuinely necessary.

## 8. Frontier effect

The proved chain is now

```text
120 interior scalar profiles
  -> rank two with a c_1-anchored F_13 pivot
  -> scalar rank at least three at height 3540
  -> centered thirteen-adic pivot extension
  -> a bounded c_1-anchored rank-three flag over F_13
  -> blocker-complete or guard/unit-marked pivot.     (51)
```

The connection and loss ledger is

```text
source:
  THM-2282's anchored rank-two packet and THM-2283's independent
  bounded rank-three relation space;

target:
  a bounded anchored rank-three Plucker coordinate over F_13 and
  a full mod-13^n relation-address bank;

map:
  cancel a centered residue address and divide the complete relation
  until the third row extends the old pivot;

preserved:
  integrality, the scalar relation equation, rational independence,
  the c_1 anchor, a uniform height invariant, and all 13-adic addresses;

destroyed:
  the original third row, descent length, terminal labels k,l,
  integer Fourier lift, root history, and owner ancestry;

cheapest hostile probes:
  prime-to-thirteen torsion in the raw quotient, the deep-blocker toy
  pivot (48), and the one-sheet zero Fourier lift (49);

needed sidecar:
  couple a terminal unit column or a signed whole-residue aggregate to
  the exact root-history frequency used by the owner transition.        (52)
```

No scalar profile is excluded. The other `45` first-depth-one profiles lie
outside THM-2282's interior scope. LRC(14) remains open. QED.

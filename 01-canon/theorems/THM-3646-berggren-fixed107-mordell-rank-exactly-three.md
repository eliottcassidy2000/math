---
id: THM-3646
title: "Berggren fixed-107 Mordell rank exactly three"
status: >
  PROVED + CITED + FINITE-EXACT; PENDING INDEPENDENT HOSTILE AUDIT.
  The Mordell curve y^2=x^3+1225041 has rational rank exactly three.
  The lower bound uses an explicit third point and three finite reductions;
  the upper bound applies Bandini's explicit 3-isogeny Selmer-support theorem
  with exact class numbers 137 and 1.  No integral-point classification is
  claimed.
source: kps-s189 / 3-isogeny descent frontier, 2026-08-21
depends_on:
  - THM-3620-berggren-fixed-107-mordell-rank-two-and-local-collision
  - THM-3643-berggren-fixed107-real-quadratic-class-number-one
related:
  - THM-3640-berggren-positive-cube-slope-atlas-through-401
script: 04-computation/berggren_fixed107_mordell_rank_three_thm3646.py
output: 05-knowledge/results/berggren_fixed107_mordell_rank_three_thm3646.out
script_sha256: e88875bf71b655bbf37eb02202c9d4cbe7817c55637ce273e474cda22cffeb2f
output_sha256: 9fde201367be2f00b806ba99f8031c162b4a57f2b783b84e3f91ad43238e66eb
semantic_sha256: 30cf28dee9ac8842062cddb843db2c43e7dd66f2ba3d23d7bd70143aa2e74ae9
hash_basis: raw LF bytes for files; canonical JSON for semantic ledger
---

# THM-3646 -- fixed-`107` Mordell rank exactly three

**PROVED + CITED + FINITE-EXACT; PENDING INDEPENDENT HOSTILE AUDIT.**  Put

```text
B=107^3-2=1225041=3*408347,
E:  y^2=x^3+B.                                           (1)
```

Then

```text
rank E(Q)=3.                                             (2)
```

This corrects the tempting continuation of THM-3620: its two-point subgroup
is `3`-saturated, but it is not the full Mordell--Weil group.

## 1. The rational `3`-isogeny

Let

```text
E': Y^2=X^3-27B=X^3-33076107.                            (3)
```

The standard rational `3`-isogeny and its dual are

```text
phi:E -> E',
phi(x,y)=((x^3+4B)/x^2, y(x^3-8B)/x^3),

psi:E' -> E,
psi(X,Y)=((X^3-108B)/(9X^2),
          Y(X^3+216B)/(27X^3)).                         (4)
```

Direct substitution gives `psi phi=[3]` and `phi psi=[3]`.  Their geometric
kernels have nonzero points with ordinates `+-sqrt(B)` and
`+-sqrt(-27B)`, respectively; neither kernel has a nonzero rational point.

The imported descent result is Theorem 3.5 of Andrea Bandini,
[*3-Selmer groups for curves y^2=x^3+a*](https://dml.cz/bitstream/handle/10338.dmlcz/128267/CzechMathJ_58-2008-2_9.pdf),
*Czechoslovak Mathematical Journal* 58 (2008), 429--445.  In the notation used
here, it embeds each isogeny Selmer group in a cubic Kummer group `H(S')` in
the quadratic field `Q(sqrt(-3a))`, with the explicit support set

```text
p!=3 belongs to S'(Q)
  iff v_p(4a) is 2 or 4 and -3a is a square in Q_p,

3 belongs to S'(Q)
  iff [v_3(a)=1 and a/3=2 mod 3]
   or [v_3(a)=5 and a/243=2 mod 3].                     (5)
```

The Selmer classes additionally have cubic rational norm.  Equations `(4)`
and `(5)` are the only imported mathematical ingredients beyond standard
elliptic-curve and quadratic-form facts.

## 2. An explicit third direction

The rational point

```text
U=(50913/16,-11482065/64) in E'(Q)                       (6)
```

gives

```text
R=psi(U)-P
 =(8279053120/216766729,
   -3611785597108493/3191456551067) in E(Q),             (7)
```

where `P=(232,3703)` is the first THM-3620 point.  The denominators in `(7)`
are exactly `14723^2` and `14723^3`.

Together with

```text
Q=(4960,349321),                                        (8)
```

the three points are independent modulo `3E(Q)`.  Indeed, exact finite-group
enumeration at the good primes `5,11,13` gives, for `a,b,c in F_3`,

```text
aP+bQ+cR in 3E(F_5)   iff 2a+b+c=0,
aP+bQ+cR in 3E(F_11)  iff  a+b+c=0,
aP+bQ+cR in 3E(F_13)  iff 2a+2b+c=0.                   (9)
```

The coefficient matrix in `(9)` has determinant `2 mod 3`.  Reduction of a
rational triple-divisible point is triple-divisible at every good prime, so
`(9)` forces `a=b=c=0`.  THM-3620 proves `E(Q)_tors=0`; hence

```text
rank E(Q)>=3.                                           (10)
```

## 3. The forward Selmer bound

For `phi:E->E'`, Bandini's quadratic field is

```text
K_-=Q(sqrt(-3B))=Q(sqrt(-408347)).                      (11)
```

The fundamental discriminant is `-408347`.  Exact enumeration of primitive
reduced positive binary quadratic forms

```text
(A,C,D),  C^2-4AD=-408347,
|C|<=A<=D,
C>=0 on the reduction boundary                         (12)
```

gives exactly `137` forms.  Therefore

```text
Cl(K_-) has order 137, so Cl(K_-)[3]=0.                 (13)
```

Applying `(5)` to `a=B` leaves only the rational prime `3`: the prime `2` is
excluded because `-3B=5 mod 8`, the prime `408347` occurs to exponent one,
and `v_3(B)=1` with `B/3=2 mod 3`.  Since `3` splits in `K_-`, the set `S'`
has two prime ideals.

The exact `S'`-unit sequence and `(13)` give

```text
dim_F3 H(S')=2.                                         (14)
```

There is no unit-rank contribution because `K_-` is imaginary quadratic and
is not `Q(sqrt(-3))`.  Rational norm sends the two split-prime valuation
coordinates to their sum; it has rank one modulo cubes.  The cubic-norm
condition therefore yields

```text
dim_F3 Sel^phi(E/Q)<=1.                                 (15)
```

## 4. The dual Selmer bound

Apply `(5)` to `a'=-27B=-33076107`; its isogenous partner is `(1)`.  The
quadratic field is

```text
K_+=Q(sqrt(-3a'))=Q(sqrt(B)).                           (16)
```

THM-3643 proves its ordinary class number is `1`.  This time `(5)` leaves
only the rational prime `2`: `v_2(4a')=2` and
`-3a'=81B=1 mod 8` is a `2`-adic square, while `v_3(a')=4` and
`v_408347(a')=1` contribute no support.  The prime `2` splits because
`B=1 mod 8`, so again `S'` has two prime ideals.

The real field has unit rank one.  Thus

```text
dim_F3 H(S')=1+2=3.                                     (17)
```

The fundamental unit has rational norm `1` by THM-3643, while norm again sums
the two split-prime valuations.  The norm map has rank one, and hence

```text
dim_F3 Sel^psi(E'/Q)<=2.                                (18)
```

The use of the **ordinary** class number in `(17)` is correct: the narrow to
ordinary quotient has order `2` and is invisible to `3`-torsion.

## 5. Rank upper bound and conclusion

Because both rational isogeny kernels are trivial, `(4)` gives the exact
Mordell--Weil quotient sequence

```text
0 -> E'(Q)/phi E(Q) -> E(Q)/3E(Q)
  -> E(Q)/psi E'(Q) -> 0.                               (19)
```

Each outer quotient embeds in its corresponding Selmer group.  Using
torsion-freeness and `(15),(18)`,

```text
rank E(Q)=dim_F3 E(Q)/3E(Q)<=1+2=3.                     (20)
```

Together, `(10)` and `(20)` prove `(2)`.

## 6. Reproduction and strict boundary

Run

```text
python3 04-computation/berggren_fixed107_mordell_rank_three_thm3646.py
python3 -O 04-computation/berggren_fixed107_mordell_rank_three_thm3646.py
```

The companion verifies the isogeny identities on `P,Q,U`, every point and
denominator in `(6)--(8)`, the complete finite groups and all `27` coefficient
triples at each prime in `(9)`, all `137` reduced forms in `(12)`, the two
support calculations, the THM-3643 hashes, and an assertion-free AST.  Normal
and optimized modes reproduce the stored transcript byte for byte.

Rank three does not classify `E(Q)` as a lattice, determine its index over a
displayed basis, or classify integral points.  In particular it neither proves
that `(1)` has another positive integral Berggren point nor supplies a new
admissible slope or two-cube ray.

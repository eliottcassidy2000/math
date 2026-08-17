---
id: THM-3540
title: "Fixed Keller depth-two newest-prime residue saturation and four-packet descent"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.  At the depth-two newest divisor H, the
  predecessor cubic factors over kappa(H) as its unique birational ancestry
  section times an irreducible quadratic.  The residue decomposition image
  on the three predecessor blocks is therefore the full marked-point
  stabilizer S2.  It has the full point-and-unordered-pair orbits required by
  THM-3539, so THM-3538's six raw depth-two factors descend to exactly four
  valuation packets.  This proves neither the full decomposition group nor
  saturation at depth at least three.
source: codex/turnpike-atlas/2026-08-16
audit: >
  An independent derivation accepted the field identification, factor
  division, arbitrary-leading-coefficient discriminant identity, DVR
  nonsquare gate, residue orbit action, packet multiplicities, and replay.
  MISTAKE-421 records its repair of one unnecessary etale-to-coordinate-
  derivative shortcut; the axis value d=4 is the lawful nonvanishing gate.
depends_on:
  - THM-2473-sporadic-keller-branch-tower-depressed-trisection-anatomy
  - THM-2570-jelonek-cusp-cylinder-normalization-and-conductor
  - THM-2576-composite-jelonek-image-divisor-and-two-component-nonproperness-law
  - THM-3530-fixed-keller-all-level-image-prime-and-component-tower
  - THM-3538-fixed-keller-newest-prime-prescribed-coordinate-index-criterion
  - THM-3539-fixed-keller-newest-prime-decomposition-centralizer-and-lca-packet-floor
related:
  - THM-3533-fixed-keller-newest-prime-reduced-different-and-index-square
  - THM-3535-fixed-keller-full-wreath-and-all-level-linear-primitivity
script: 04-computation/keller_depth_two_newest_prime_residue_saturation_20260816.py
output: 05-knowledge/results/keller_depth_two_newest_prime_residue_saturation_20260816.out
script_sha256: a92e5237f4ceed6b086e0f29c89536c494bbd60588630e4331b3b5443fe2d968
output_sha256: 6477e75c328111f4cbac497b82c1f369bbb845bfe7b66a34dcc90773d5e5f0a8
semantic_sha256: 81e36ca0659eb75dc7a0f660f4c8a40ae70bf0d4fba56c08f4ef17241349488e
hash_basis: LF-normalized bytes
---

# THM-3540 -- depth two has the full residue orbit gate

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

Retain the fixed sporadic Keller map `F`, the old Jelonek prime `L=P_0`,
and its image prime `H=P_1` up to the harmless rational normalization in
THM-3530.  At the generic point of `V(H)`, this theorem determines the
residue action on the three predecessor blocks.  The result is exactly the
first positive instance of THM-3539's decomposition-saturation gate.

## 1. The birational ancestry section gives one rational root

THM-2570 normalizes `V(L)` by `A^2_(tau,lambda)`:

```text
X=lambda^2(3-tau lambda)/27,
Y=lambda(4-tau lambda)/3,
Z=tau.                                                   (1)
```

The normalization is finite birational, so

```text
kappa(V(L))=Q(tau,lambda).                              (2)
```

THM-2576 and THM-3530 prove that

```text
F:V(L) -> V(H)                                          (3)
```

has generic degree one.  Therefore `(2)` is also `kappa(H)`, and the generic
point

```text
q_*=(X,Y,Z) in V(L)                                     (4)
```

is a `kappa(H)`-rational point in the predecessor fibre over the generic
point `eta_H`.

Write the target coordinates of `eta_H` as `(a,b,c)` and abbreviate

```text
L_o=L(a,b,c),       T_o=4-3bc,
S_o=27ac^2-9bc+8.                                      (5)
```

The primes `L` and `H` are distinct, so `L_o` is nonzero in `kappa(H)`.
THM-2473's inverse `x`-core is

```text
E(W)=L_o W^3+T_o W-2c.                                 (6)
```

The point `(4)` makes `W=X` a root.  Exact division gives

```text
E(W)=(W-X) Q(W),
Q(W)=L_o W^2+L_o XW+(L_o X^2+T_o).                    (7)
```

Thus the predecessor factorization has one rational ancestry root and one
residual quadratic.  The question is whether that quadratic splits.

## 2. The residual quadratic has square class `[-L_o]`

Put

```text
d=E'(X)=3L_o X^2+T_o,
delta=Disc(Q)=(L_o X)^2-4L_o(L_o X^2+T_o).             (8)
```

The discriminant of a linear factor times a quadratic gives the exact
identity

```text
Disc(E)=delta d^2.                                     (9)
```

THM-2473's cusp identity supplies the independent expression

```text
Disc(E)=-4 S_o^2 L_o.                                  (10)
```

The axis calculation `(13)` below gives `d|_(lambda=0)=4`, so `d` is a
nonzero rational function.  Combining `(9)` and `(10)` yields

```text
delta=-L_o (2S_o/d)^2,
[delta]=[-L_o] in kappa(H)^*/kappa(H)^(2).             (11)
```

This is a square-class statement, not yet a nonsquare statement.  The next
axis specialization supplies the missing valuation.

## 3. The normalization axis proves nonsquareness

On the prime divisor `lambda=0` of the normalization plane, formulas `(1)`
and the exact map give

```text
q_*=(0,0,tau),             F(q_*)=(tau,0,0).            (12)
```

Consequently

```text
L_o=16tau,       S_o=8,       d=4,
E(W)=4W(4tau W^2+1),
Q(W)=4(4tau W^2+1),
delta=-256tau.                                         (13)
```

Suppose `-L_o` were a square in `Q(tau,lambda)`.  Its `lambda`-valuation is
zero because its residue in the `lambda`-adic DVR is the nonzero function

```text
-16tau in Q(tau).                                      (14)
```

A square root would therefore be a unit in that DVR, and reduction modulo
`lambda` would make `(14)` a square in `Q(tau)`.  This is impossible: its
`tau`-valuation is one.  Hence

```text
-L_o is not a square in kappa(H).                      (15)
```

Equations `(11)` and `(15)` make `delta` nonsquare.  In characteristic zero,
a nondegenerate quadratic is reducible exactly when its discriminant is a
square.  Therefore `Q` is irreducible over `kappa(H)`.

The proof uses the axis only as a divisor-valued nonsquare detector.  It does
not replace the generic point of `H` by the special target `(tau,0,0)`.

## 4. The residue action is the full marked-point stabilizer

The predecessor cubic over `kappa(H)` now has factorization type

```text
1+2.                                                    (16)
```

Its linear factor is the marked ancestry block `b_*`; its irreducible
quadratic forces the residue Galois action to swap the other two blocks.
Therefore the image `bar D` of the newest-prime decomposition group on the
three predecessor blocks is

```text
bar D=Stab_(S_3)(b_*) isomorphic to S_2.               (17)
```

At depth two, this is exactly THM-3539's maximal marked-leaf action `H_1`.
Its point and unordered-pair orbits have sizes

```text
blocks:             1,2;
unordered pairs:    2,1.                               (18)
```

Thus the point-and-pair saturation gate `G_(2-orb)` is proved for `n=2`.
Notice the typing: `(17)` determines the **image on predecessor blocks**.
It does not determine the kernel of the residue action on last-step roots and
therefore does not prove that the full decomposition group equals the full
centralizer from THM-3539.

## 5. The exact four-packet descent

For any THM-3538 integral observation `theta=y,z,u=1/x`, write `A_*` for the
internal factor on the marked block, `A_o` for either off-ancestry block,
`R_(*,o)` for either marked/off resultant, and `R_(o,o')` for the resultant
between the two off blocks.  THM-3539's orbit formula specializes to

```text
i_(theta,2)
 =v(A_*)+2v(A_o)+2v(R_(*,o))+v(R_(o,o')).              (19)
```

There are exactly

```text
2 internal packets + 2 pair packets = 4=2^2            (20)
```

instead of the raw list of three internal factors and three resultants.
THM-3538 separately proves that all these factors are units for `y_2,z_2,u_2`.
The new result is the lawful residue descent of that list, not another route
to its unitness.

## 6. Equality and failure boundary

THM-3540 proves no statement at depth `n>=3`.  In particular, it does not
prove:

1. full decomposition-centralizer equality even at depth two;
2. the marked-LCA point/pair saturation gate at depth three or higher;
3. unitness of any new THM-3538 representative at `n>=5`;
4. an all-level prescribed-coordinate index formula beyond the conditional
   packet identity;
5. transport to the old-`L` prime, where THM-3537's predecessor is ramified;
   or
6. an arbitrary-Keller classification, `JC(2)`, a current, or LRC(14).

The next residue target is now precise: at `n=3`, factor the degree-nine
predecessor algebra over `kappa(J)` far enough to recover its orbits on the
nine blocks and their unordered pairs.  Full field identification is more
than THM-3539 needs; point-and-pair saturation remains the economical gate.

## 7. Exact companion

Reproduce with

```text
python -B 04-computation/keller_depth_two_newest_prime_residue_saturation_20260816.py
python -B -O 04-computation/keller_depth_two_newest_prime_residue_saturation_20260816.py
```

The companion verifies the normalization, map, root division, both
discriminant identities, square ratio, axis valuation, and all point/pair
orbits symbolically.  It contains no executable `assert`; ordinary and
optimized transcripts match the stored output exactly.

**QED.**

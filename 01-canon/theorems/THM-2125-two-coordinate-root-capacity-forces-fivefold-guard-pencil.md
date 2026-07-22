---
id: THM-2125
title: "Two-coordinate root capacity forces a fivefold guard pencil modulo thirteen"
status: >
  PROVED. In a rank-two guard/eight-terminal torus cover, assume the guard is
  nonzero modulo 13 and at least one terminal is a 13-blocker. If between one
  and four nonblocker terminals are projectively parallel to the guard modulo
  13, an almost-everywhere cover is impossible. Divide the true blockers and
  make them safe by settled LRC. On the 13^2 root fiber, guard danger forbids
  at most four guard-coordinate columns and each parallel nonblocker at most
  two; a surviving column exists. Each remaining transverse terminal forbids
  at most two points in that column, and there are at most six. A 169-sheet
  measurable pigeonhole gives positive-measure escape. With THM-2123, every
  non-guard-blocker rank-eight cover needs at least five aligned nonblockers.
source: codex-2026-07-22-LRC-two-coordinate-root-capacity
depends_on:
  - THM-2114
  - THM-2122
  - THM-2123
related:
  - THM-2124
  - THM-2130
external: Settled Lonely Runner Conjecture for at most eight integer speeds (LRC(9) in the total-runner convention).
---

# THM-2125 -- two-coordinate root capacity

Let `Gamma` be a saturated rank-two character lattice,

```text
K=Hom(Gamma,R/Z),                                      (1)
```

and take a nonzero guard `g` and eight nonzero terminal characters `c_i`.
Assume `g mod 13` is nonzero. Split the terminal labels into

```text
B={i:c_i=0 mod 13},                         b=|B|>=1;
R={i:c_i in F_13^* g mod 13},               r=|R|;
T={1,...,8}\(B union R),                    t=|T|.       (2)
```

Thus `B` consists of true terminal 13-blockers, `R` of nonzero projective
guard-parallel reductions, and `T` of transverse reductions. Suppose

```text
1<=r<=4.                                               (3)
```

We prove that the terminal danger bands cannot cover

```text
C={X:||g.X||>1/7}                                     (4)
```

even almost everywhere.

## 1. An open phase where every true blocker is safe

For `i in B`, write `c_i=13u_i`. THM-2122 Section 1 applies verbatim: choose
a generic integral cocharacter separating the nonzero sign classes of the
`u_i`, and apply settled LRC for at most `b<=8` distinct speeds. It produces
`Y_0 in K` with

```text
||u_i.Y_0||>=1/9>1/14                    for all i in B. (5)
```

Hence there is a nonempty open neighborhood `U` of `Y_0` on which all these
inequalities remain strict.

## 2. Coordinates on all 169 thirteenth roots

Multiplication by thirteen is a finite covering `[13]:K->K`. Choose a local
inverse sheet `sigma:U->K`. Since `g mod 13` is a nonzero linear functional on
the two-dimensional space `K[13]`, choose a basis `z,v of K[13]` with

```text
g.z=1/13,                         g.v=0.               (6)
```

All 169 local root sheets above `U` are

```text
X_jk(Y)=sigma(Y)+jz+kv,             (j,k) in F_13^2.  (7)
```

For a blocker `i in B`,

```text
c_i.X_jk(Y)=13u_i.X_jk(Y)=u_i.Y,                       (8)
```

so (5) makes it safe on every sheet.

## 3. First choose a safe guard column

Fix `Y in U`. As `j` varies, the guard values in (7) form a translate of the
thirteen-grid, independent of `k`. Its closed radius-`1/7` unsafe arc has
length `2/7<4/13`, so at most four values of `j` fail strict guard safety.

For `i in R`, choose `lambda_i in F_13^*` with

```text
c_i=lambda_i g mod 13.                                 (9)
```

Equation (6) gives `c_i.v=0` and `c_i.z=lambda_i/13`. Therefore `c_i.X_jk`
is independent of `k`, and as `j` varies it is another complete translate of
the thirteen-grid. Its radius-`1/14` danger arc has length `1/7<2/13`, so it
forbids at most two columns.

The guard and all `r` parallel nonblockers together forbid at most

```text
4+2r<=12<13                                            (10)
```

columns. Choose a column `j` on which the guard and every label in `R` are
safe.

## 4. Then choose a safe point inside the column

For `i in T`, transversality says `c_i.v!=0`. Holding the chosen `j` fixed and
varying `k`, its thirteen values form a complete translated grid, so its
danger band hits at most two points.

Because `b>=1` and `r>=1`,

```text
t=8-b-r<=6.                                             (11)
```

Thus all transverse terminals forbid at most

```text
2t<=12<13                                               (12)
```

points of the surviving column. Choose `k` outside their union. Equations
(5), (8), (10), and (12) show that `X_jk(Y)` is guard-safe and safe for every
terminal.

## 5. The 169-sheet measurable pigeonhole

The pair `(j,k)` may depend on `Y`. For every root index define

```text
A_jk={Y in U:X_jk(Y) is guard-safe and terminal-safe}. (13)
```

The 169 measurable sets `A_jk` cover the positive-measure open set `U`.
Hence one has positive measure. Its image under the corresponding local root
sheet also has positive measure and lies in `C` outside all open terminal
danger bands. Removing the finitely many equality hypersurfaces leaves a
positive-measure strict escape set.

This contradicts an almost-everywhere cover of (4), proving the theorem. QED.

## 6. Fivefold projective invoice

Now assume the rank-two span hypotheses of the live rank-eight cover. If the
guard is nonzero modulo thirteen, THM-2114 forces `b>=1`. THM-2123 excludes
`r=0`, while the theorem above excludes `1<=r<=4`. Therefore every surviving
cover with a non-13-blocker guard must satisfy

```text
r>=5: at least five nonblocker terminal reductions lie in F_13^* g.          (14)
```

Since a true terminal blocker is also present, at most two terminal labels
remain outside this guard pencil. This is a located multiplicity statement,
not merely the determinant existence of THM-2123.

THM-2124 treats the complementary branch where the **guard** is the
13-blocker. THM-2130 extends the two-coordinate capacity count to prime eleven
and to the applicable higher-rank residual counts. Neither extension is used
in the proof above.

## 7. Assumption challenge and Tournament Analysis

The challenged assumption was that all 169 roots should be treated as an
unstructured finite-plane cover. The projective relation to the guard makes
the root fiber triangular: first choose the guard coordinate `j`, then the
kernel coordinate `k`. This two-stage capacity is stronger than a flat union
bound, which would total all forbidden root points and lose the constant-band
structure.

Candidate tournament vertices were terminals, the 169 roots, thirteen
columns, projective directions, and coverage obligations. The faithful
vertices are the thirteen guard columns with inner thirteen-point fibers.
Orienting columns by forbidden-label count and using cyclic `j` order as the
tie Hamiltonian path preserves search priority, but score histograms, cycles,
SCCs, edge flips, and path counts erase which labels are constant on columns.
The certificate must retain the two-level incidence matrix and the local-sheet
positive-measure sidecar.

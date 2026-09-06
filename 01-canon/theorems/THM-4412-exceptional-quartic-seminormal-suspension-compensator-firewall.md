---
id: THM-4412
title: "Exceptional-quartic seminormal suspension and compensator firewall"
status: >
  PROVED + VERIFIED-EXACT, relative to the exact THM-4381/4404 exceptional
  curve. If r is its unique seminormal generator and D=z+r, then the embedded
  suspension B_susp=S[D] in K[x,z] satisfies B_susp intersect K[x]=S and does
  not contain z. Every constant-z section has image S[r]=S^sn: it preserves
  all conductor point fibres, changes nothing at the 86 nodes, and adds
  exactly the missing first-jet direction at the plane triple. The target
  section D=0 instead requires the non-descended source graph z=-r. This is
  an algebraic one-coordinate model, not stable cancellation, a planar
  embedding, chart/seam entry, a Keller pair, JC(2), or DC(2).
source: root / JC2 and arXiv continuation session, 2026-09-05
audit: >
  PASS relative to the independently audited dependencies. THM-4381 supplies
  S^sn=S+K*r, r notin S, the radical conductor and the one-triple-plus-86-node
  local classification; THM-4404 supplies the nonzero retained period of r'.
  The companion checks the translated-coordinate sections, the all-degree
  leading-coefficient pattern behind the intersection lemma, the exact
  three-row jet determinant, defect bookkeeping, and the contrast with
  Long's independently audited 1+2 value idempotent. Normal, optimized, and
  fixed-hash-seed runs byte-match the frozen output.
depends_on:
  - THM-4381-exceptional-quartic-seminormalization-and-conductor-fibre-classification
  - THM-4404-exceptional-quartic-descended-two-form-seminormal-cokernel
related:
  - THM-4408-rank-two-poisson-hamiltonian-primitive-nondescent-compensator-firewall
  - THM-4411-first-order-collision-transgression-seminormal-tradeoff
  - THM-4397-rank-two-poisson-counterexample-symplectic-gauge-equivalence
script: 04-computation/jc2_exceptional_quartic_seminormal_suspension_formal_audit_s616.py
output: 05-knowledge/results/jc2_exceptional_quartic_seminormal_suspension_formal_audit_s616.out
script_sha256: 82b774b1e32701b63f109228bf92587bfae8789257342022d1bdf2d8d6e73dbb
output_sha256: c7310cf5b048f21deba398b9f2927177f58ac9b708661f848a68af1de1cd0bef
hash_basis: raw LF bytes
---

# THM-4412 -- one extra coordinate realizes the seminormal jet but cannot descend

**PROVED + VERIFIED-EXACT, RELATIVE TO THM-4381/4404. THIS IS AN EMBEDDED
ALGEBRAIC SUSPENSION OF ONE RESTRICTED CURVE. IT IS NOT A SYMPLECTIC MAP,
STABLE CANCELLATION, AN AFFINE-PLANE EMBEDDING, CHART OR SEAM ENTRY, A KELLER
PAIR, `JC(2)`, OR `DC(2)`.**

## 1. The exceptional curve and its missing class

Work over the characteristic-zero exceptional quartic field `K`. Retain the
normalization and restriction algebra of THM-4381:

```text
S=K[B,C,E] subset N=K[x],
normalization(S)=N.                                    (1)
```

With `h_172` from the conductor calculation, put

```text
r=x(x^2-1)h_172.                                       (2)
```

THM-4381 proves

```text
r notin S,
S^sn=S+K r=S[r],
dim_K(S^sn/S)=1,
r^2 in S.                                              (3)
```

It also classifies the reduced conductor fibres as 86 ordinary nodes and one
ordinary plane triple with normalization points `x=(-1,0,1)`. The polynomial
`r` vanishes simply at every conductor branch. At each node it already lies
in the local conductor and hence in the local ring. At the triple its first
derivative is the unique seminormal jet missing from `S`.

## 2. The translated suspension and exact intersection

Let `z` be an indeterminate over `N` and define

```text
D=z+r,
B_susp=S[D] subset N[z].                               (4)
```

Since `D` is monic of degree one in `z`, it is transcendental over `N` and
`B_susp` is abstractly a polynomial extension `S[D]`. Nevertheless its
position inside `N[z]` remembers the non-descended translate. Precisely,

```text
boxed: B_susp intersect N=S,
       z in B_susp iff r in S.                         (5)
```

To prove `(5)`, write an element of `B_susp` as

```text
a_0+a_1D+...+a_mD^m,                  a_i in S.        (6)
```

If `m>0`, its coefficient of `z^m` after substituting `D=z+r` is `a_m`, so
the element cannot lie in `N`. This gives `B_susp intersect N=S`. Since
`r=D-z`, membership of `z` implies `r in B_susp intersect N=S`; conversely
`z=D-r` belongs to `B_susp` if `r in S`. By `(3)`, the actual conclusion is

```text
boxed: z notin B_susp.                                 (7)
```

This is an embedded-ring statement. Abstractly `B_susp` is just `S` with one
polynomial variable; `(7)` is not a cancellation obstruction.

## 3. Source sections realize the seminormalization

Evaluation on the source section `z=0` sends

```text
D |-> r,
```

so `(3)` gives

```text
boxed: image(B_susp -> N, z=0)=S[r]=S^sn.              (8)
```

The same holds on every constant source section `z=c`, `c in K`, because

```text
D |-> c+r,                    S[c+r]=S[r].             (9)
```

Thus the extra output coordinate realizes exactly the seminormalization, not
the full normalization. The length ledger stays

```text
dim_K(N/S)=89,
dim_K(S^sn/S)=1,
dim_K(N/S^sn)=88.                                     (10)
```

At the 86 nodes, adjoining `r` makes no local change because the nodes are
already seminormal and `r` is in their local conductor. All of the change in
`(8)` is concentrated at the retained plane triple.

## 4. The point collision remains, but its missing jet is restored

At `x=(-1,0,1)`, THM-4381 gives

```text
(B,C,E)=(0,0,-3),                  r=(0,0,0).          (11)
```

Hence on the constant section `z=c`, all three branches still have

```text
(B,C,E,D)=(0,0,-3,c).                                (12)
```

The point collision is not split. The change occurs one conductor layer
later, in first derivatives. The old tangent rows are

```text
C'=(3,3,3),                    E'=(-9,4,9),            (13)
```

and their unique relation is

```text
ell=(5,-18,13).                                        (14)
```

Writing `r'=(a,b,q)` on the ordered branches, direct expansion gives

```text
det [ C' ; E' ; r' ]=3(5a-18b+13q).                   (15)
```

THM-4404 proves that the period on the right of `(15)` is nonzero. Therefore
the three branch tangents in the augmented `(C,E,D)` coordinates are
independent. The added row is exactly the one-dimensional jet missing from
the planar tangent image. This is the local differential manifestation of
the global equality `S[r]=S^sn` in `(8)`.

The distinction is sharp:

```text
point fibre before and after adjoining D:     still three branches,
first-jet span before adjoining D:             dimension 2,
first-jet span after adjoining D:              dimension 3.             (16)
```

## 5. The target section moves the sidecar into the source graph

Imposing the descended target hyperplane `D=0` in `(4)` forces

```text
z=-r.                                                     (17)
```

The restricted output algebra is

```text
B_susp/(D) ~= S.                                        (18)
```

The source graph `(17)` is polynomial in the normalization coordinate `x`,
but it is not descended through the core algebra because `r notin S`.
Equations `(7)--(8)` show that the auxiliary coordinate has not been
eliminated over the core: choosing a constant source graph exposes the
seminormal coordinate in the target, while choosing a constant target graph
hides it in the source.

This is the exact jet-level analogue of THM-4408's compensator-intersection
firewall.

## 6. Value layer versus jet layer: the Long comparison

On Long's reduced three-point core fibre, THM-4408 computes

```text
Hbar=(-1/48,-1097/192,-1097/192),                     (19)
```

whose normalized residue is the idempotent `(0,1,1)`. Thus a constant source
slice for Long's independent coordinate makes `D=D_0+H` take two values and
splits the point fibre as `1+2`. That primitive carries **branch-value data**.

The exceptional class `r` instead has values `(0,0,0)` in `(11)` and a
nonzero missing derivative period in `(15)`. A constant `z` slice preserves
the threefold point collision but turns the planar two-dimensional tangent
span into a three-dimensional one. It carries **branch-jet data**.

These are adjacent but distinct layers of branch information:

```text
Long H:       visible modulo branch maximal ideals; value idempotent;
quartic r:    zero modulo branch maximal ideals; visible modulo their squares.
                                                                    (20)
```

Here the second line refers to the retained triple's first infinitesimal
neighbourhood. Both constructions require a non-descended translate in their
fixed embedded cores. The difference explains why the Long slice separates a
branch while the
seminormal suspension repairs the infinitesimal embedding without separating
the underlying points.

## 7. Consequence and boundary

The positive result is a precise stabilized model: one auxiliary coordinate
is sufficient to realize the unique seminormal jet, and `(5)` identifies the
exact obstruction to eliminating that coordinate by a core-descended graph.
Together with THM-4411, it shows that the missing seminormal branch jet can
be carried after stabilization, but not while simultaneously retaining a
planar descended two-form and a collision-preserving normal motion in the
fixed compiler.

Nothing here supplies a polynomial endomorphism of `A^2`, much less a
noninjective Keller map. The target on a constant source section has the
extra coordinate `D`, and the target section `D=0` uses the non-descended
source graph `(17)`. No planar chart or seam entry follows; `JC(2)` and
`DC(2)` remain open.

## 8. Reproduction

Replay the formal exact audit from the repository root:

```text
python3 -B 04-computation/jc2_exceptional_quartic_seminormal_suspension_formal_audit_s616.py
python3 -B -O 04-computation/jc2_exceptional_quartic_seminormal_suspension_formal_audit_s616.py
PYTHONHASHSEED=1729 python3 -B 04-computation/jc2_exceptional_quartic_seminormal_suspension_formal_audit_s616.py
```

All modes byte-match the frozen LF output and perform `25` dynamic checks.
The new proof is formal relative to the exact dependencies; the companion
does not rerun their degree-178 conductor computation.

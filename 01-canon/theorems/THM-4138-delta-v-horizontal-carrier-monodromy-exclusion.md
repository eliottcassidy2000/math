---
id: THM-4138
title: "Delta-V horizontal-carrier orbit-merger exclusion"
status: >
  PROVED RELATIVE TO THM-4103/4130/4134 + VERIFIED-EXACT. The degree-16
  and degree-15 finite horizontal BC carriers left by THM-4134 are
  impossible. Together with THM-4134's exclusion of degrees 20 and 19,
  this empties only the theta-only exact-M=8 Delta_V=0 wall. The other
  collision walls, other reduced cells, JC(2), and DC(2) remain OPEN.
source: codex-planar-jacobian-cycle6-63-20260825
depends_on:
  - THM-4103-jc23-theta-boundary-ramification-and-degree-response
  - THM-4130-theta-only-extremal-seam-critical-monodromy-obstruction
  - THM-4134-delta-v-collision-wall-strict-transform-and-high-branch-exclusion
related:
  - THM-3996-etale-node-address-balance-cycle-and-nonproperness-dichotomy
  - THM-4120-jc23-extremal-target-degree-twenty-one-response
  - THM-4122-planar-keller-asymptotic-width-and-resonant-shear-contraction
script: 04-computation/jc23_delta_v_horizontal_carrier_monodromy_thm4138.py
output: 05-knowledge/results/jc23_delta_v_horizontal_carrier_monodromy_thm4138.out
script_sha256: fc6216331dd19181706a3ec51f1b9f29985c5550390e4d60ed08cb94e218da97
output_sha256: 00e4466e097cc9108eb066103697d738bffd6355a3b3a75cdc847f1bff99b190
semantic_sha256: e1593d26e2a4b6931c1550868dbc3c85c9cf9430813bc5ba28c095dd72b28b66
hash_basis: raw LF bytes
audit: >
  PASS. The exact checker freezes both wall ledgers, the sharp permutation
  index budget, distinct-value and collided-value carrier controls, two
  minimal failure hostiles, and the identity between THM-4134's displayed
  polynomial section and 3P on the quadratic base change. Normal and
  optimized runs byte-match the frozen output.
---

# THM-4138 -- horizontal carriers cannot merge enough sheet orbits

**PROVED RELATIVE TO THM-4103/4130/4134 + VERIFIED-EXACT; JC(2) OPEN.**
Retain the theta-only exact-`M=8` collision wall

```text
Delta_V=0,                    Phi!=0.                    (1)
```

THM-4134 leaves precisely two low-degree alternatives:

```text
generic wall:       n=16, genus=8, #Crit_aff=19,
secondary wall:     n=15, genus=7, #Crit_aff=18.         (2)
```

In each alternative the two index-two `BC` punctures have finite target
image. This theorem proves that neither alternative is realizable by a
planar Keller pair.

## 1. Inheritance and the missing invariant

The closest proved mechanism is THM-4130's injection from affine preimages
of the two target nodes into fixed sheets of their vanishing monodromies.
The canonical hostile is THM-4134's quadratic horizontal carrier: it removes
the full-boundary hypothesis, can pass through a target node, and inserts
one or two additional meridians into the punctured torus. Thus the smooth
commutator proof cannot simply be copied.

The least-used sidecar is not the Mordell--Weil label of the horizontal
point, but the **total permutation index** of its meridians. The two `BC`
points have ramification indices `2,2`, hence total finite-carrier defect

```text
(2-1)+(2-1)=2.                                          (3)
```

This remains `2` whether their target images are distinct or collide. The
proof compares that exact orbit-merger budget with the fixed sheets forced
by `(2)`.

## 2. The generic covers and their finite branch set

Fix a smooth generic pencil value `q_*`. Let

```text
bar(varphi): C -> E                                      (4)
```

be the corresponding finite map of smooth projective curves, and let `O`
be the target origin. In either row of `(2)`, let `B` be the set of finite
images of the two geometric `BC` punctures. Thus `|B|` is one or two.

The packets from THM-4134 split as

```text
generic:    over O  (7,5,3,1),     over B  (2,2),
secondary:  over O  (7,3,2,2,1),   over B  (2,2).       (5)
```

Their ramification defects are respectively

```text
generic:    12+2=14=2*8-2,
secondary:  10+2=12=2*7-2.                              (6)
```

Hence Riemann--Hurwitz is saturated by the displayed boundary points. In
particular, if `Z_b` is a small-meridian permutation for `b in B`, then

```text
sum_(b in B) ind(Z_b)=2,                                (7)
ind(sigma):=n-#Cycles(sigma).
```

If the two images are distinct, the two `Z_b` are transpositions. If they
coincide, the one local branch permutation is a product of two disjoint
transpositions. Distinct ramification points give disjoint local cycles, so
both cases have the same index ledger `(7)`.

After removing `O union B` and its complete inverse image, `(4)` is an
unramified cover. Its source is a connected compact Riemann surface with
finitely many points deleted, hence remains connected. Therefore its
monodromy action on the `n` sheets is transitive.                         (8)

## 3. Punctured fixed-sheet transport, including a node collision

Let `o_0,o_1` be the two affine nodes of the target pencil and put

```text
r_i=#F^(-1)(o_i) in A2.                                  (9)
```

The gradient identity and reduced critical calculation inherited in
THM-4130/4134 give

```text
r_0+r_1=#Crit_aff=
  19=n+3  on the generic wall,
  18=n+3  on the secondary wall.                        (10)
```

We need the fixed-sheet injection in the punctured, rather than the
full-boundary, cover.

> **Punctured fixed-sheet lemma.** There are loops `delta_0,delta_1` in
> `E-(O union B)` whose images after filling `B` are THM-4130's two standard
> vanishing generators of `pi_1(E-O)`. If
>
> ```text
> X=mon(delta_0),                 Y=mon(delta_1),         (11)
> ```
>
> then
>
> ```text
> #Fix(X)>=r_0,                  #Fix(Y)>=r_1.            (12)
> ```

### Proof

For every affine `z in F^(-1)(o_i)`, the Keller identity makes `F` a local
biholomorphism. In target holomorphic Morse coordinates one has

```text
E-q_i=uv,                                                (13)
```

and pulling `u,v` back at `z` gives the same equation, with no ramified
base replacement. Choose pairwise disjoint inverse neighborhoods of the
finitely many points in `(9)`, all mapping biholomorphically to one small
target neighborhood.

For a nearby nonzero `tau=q-q_i`, the local target fibre in `(13)` is a
Milnor annulus. Its intersection with the moving finite branch set is
finite. A core circle

```text
|u|=c |tau|^(1/2),             v=tau/u                  (14)
```

can therefore be chosen to avoid every such point. It has one closed
degree-one lift in each inverse neighborhood. These lifts are distinct.

Resolve the compactified graph of the pencil map, and remove from the
`q`-line its finite additional Hurwitz discriminant. Starting at a small
`tau!=0`, choose a generic base path to the common reference value in this
complement. Proper smooth transport and isotopy extension carry the
punctured core `(14)` together with all its degree-one lifts. The induced
bijection of sheet labels **conjugates** the local monodromy permutation;
it does not postmultiply it by a branch permutation. Hence all `r_i` closed
lifts remain distinct fixed sheets at the reference fibre. This proves
`(12)`.

This argument also covers the case in which a horizontal branch value tends
to `o_i`. Only the central fibre `tau=0` contains the collision. On every
nearby annulus the branch set is finite and a parallel core `(14)` avoids
it. The resulting reference loop may differ from another punctured lift of
the same filled vanishing generator by based `B` meridians. We retain the
actually transported loop and its conjugate monodromy; we never multiply
its permutation by those meridians. Thus a branch braid does not alter the
fixed-sheet injection.

Finally, filling the punctures `B` induces a surjection

```text
pi_1(E-(O union B)) -> pi_1(E-O),                       (15)
```

whose kernel is normally generated by the `B` meridians. THM-4130 proves
that the filled images of `delta_0,delta_1` are the two geometric handle
generators of `pi_1(E-O)`. Choose disjoint based arcs to the points of `B`.
Cutting the punctured torus along the two transported handle loops and those
arcs gives a disk. Equivalently, the standard punctured-surface presentation
gives

```text
pi_1(E-(O union B))=<delta_0,delta_1, meridians of B>.  (16)
```

This is exactly the generating set used below. QED.

The local conclusion uses only etaleness and the unramified base identity
in `(13)`. It neither assumes that the node is outside the nonproperness
set nor that the horizontal point is a constant section.

## 4. The orbit-merger inequality

For a permutation `sigma`, write `supp(sigma)` for its nonfixed letters.
Equations `(10),(12)` imply

```text
|supp(X)|+|supp(Y)|
  =2n-#Fix(X)-#Fix(Y)
  <=2n-(n+3)=n-3.                                      (17)
```

If a nonidentity permutation has support size `s`, then

```text
ind(sigma)<=s-1;                                        (18)
```

each nontrivial cycle saves one more than the crude support bound. Thus,
if at least one of `X,Y` is nonidentity,

```text
ind(X)+ind(Y)<=n-4.                                     (19)
```

If both are the identity, the left side is zero and the conclusion below
is even stronger.

For any permutations `g_1,...,g_k` on `n` letters, form a graph by adding a
spanning tree on each nontrivial cycle of every `g_j`. It has
`sum_j ind(g_j)` edges, and its connected components are the orbits of the
generated group. Consequently a transitive generated action requires

```text
sum_j ind(g_j)>=n-1.                                    (20)
```

Apply `(20)` to the complete generating set `(16)`. Equations `(7),(19)`
give

```text
ind(X)+ind(Y)+sum_(b in B) ind(Z_b)
                         <=(n-4)+2=n-2<n-1.             (21)
```

If `X=Y=1`, the total is only `2`. Thus the action has at least two orbits
in every case, contradicting transitivity `(8)`. This excludes both rows of
`(2)`.

## 5. Sharp controls and the `3P` collision hostile

The inequality is sharp. On `n` letters let `X` be one `(n-3)`-cycle,
let `Y=1`, and use two transpositions to attach two of the three fixed
letters successively. The generated action has exactly two orbits and total
index `n-2`, attaining `(21)`.

Both numerical margins are load-bearing:

1. with only `#Fix(X)+#Fix(Y)=n+2`, an `(n-2)`-cycle and two attaching
   transpositions generate a transitive action with total index `n-1`;
2. with the actual `n+3` fixed-sheet sum but carrier defect `3`, an
   `(n-3)`-cycle and three attaching transpositions do the same.

The horizontal branch can genuinely meet a target node. On the quadratic
base change

```text
q=a^3/2+rho^2,
E_q: V^2=U^3-(3a^2/4)U+q-a^3/4,                        (22)
```

put `P=(a/2,rho)`. The group law gives

```text
2P=(-a,-rho),
3P=(a/2+16rho^2/(9a^2),
    -rho-64rho^3/(27a^3)).                              (23)
```

Thus THM-4134's displayed polynomial section is exactly `3P`. At `rho=0`
it meets the node `(a/2,0)`. Writing `xi=U-a/2` gives the exact local form

```text
q-a^3/2=V^2-(3a/2)xi^2-xi^3,                           (24)
```

while on `3P`, `xi=O(rho^2)` and `V=-rho+O(rho^3)`. This
is the explicit collision hostile handled by the parallel-core argument
`(14)`. The identity `(23)` is not used to classify every possible BC
image, and no Mordell--Weil completeness claim is needed.

## 6. Consequence and boundary

THM-4134 already excludes the full-boundary degrees `20` and `19`.
Equation `(21)` now excludes its remaining horizontal degrees `16` and
`15`. Therefore the theta-only exact-`M=8` `Delta_V=0` wall is empty.

This proof does not cross the `Delta_D=0` or `delta+theta=0` walls, treat
maximum residual weight at least nine, leave the inherited reduced `(2,3)`
cell, or prove `JC(2)` or `DC(2)`. The `3P` calculation is a hostile control,
not a bridge from elliptic divisibility dynamics to arbitrary Keller maps.
**QED.**

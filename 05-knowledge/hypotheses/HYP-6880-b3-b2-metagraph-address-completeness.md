---
id: HYP-6880
title: Joined B3-face and folded-B2 sidecars form a recursively useful metagraph address
status: OPEN GENERAL/TRANSPORT CLAIM; finite exact Omega+B2 injectivity proved through n=8 by THM-801/809
source: codex-2026-07-15-S12
depends_on: [THM-553, THM-796, THM-801]
related: [THM-805, HYP-2685, HYP-3234, HYP-6825, HYP-6870, HYP-6885]
---

# HYP-6880 — joined B3/B2 metagraph address

The full three-subtriangle chart preserves overlap ownership and the missing
gap-contraction face.  The mirror-folded two-coordinate chart preserves the
reflection-orbit phase that face-node projection forgets.  The hypothesis is
that their join is a useful recursive address for complement lines and their
node-pair fibres beyond the current finite atlas.

## Finite result now proved

THM-801 completed the exact audit:

```text
n       lines       Xi cells    Omega cells    Omega+B2 cells
4           4              4              4                 4
5          32             32             32                32
6         512            509            510               512
7      16,384         16,031         16,308            16,384
8   1,048,576               -              -         1,048,576
```

Thus `Omega+B2` is injective on literal lines through `n=8`.  At `n=7`, the
gap face removes 277 of `Xi`'s 353 collision excess and `B2` removes the
remaining 76.  Inside coloured node-pair fibres, `B2`, `B3`, and their join
leave respectively 172, 1,368, and 16 excess collisions; the two charts
preserve genuinely different information.

The same computation found exact compact companions:

- the gap-face support row resolves five of endpoint deletion's eight support
  twin pairs at `n=7`;
- the boundary-curvature polynomial `K_u(z)` separates 238/272 nodes, or
  249/272 with total `C3`;
- the common-core node pair improves `Xi_7` to 16,110 cells but is weaker than
  the gap face;
- the joined address remains a static codec, not a proved recursive Markov
  state.

The finite question is exact: disintegrate literal complement lines first by
their `Omega_n` upper/three-face node tensor, then by raw mirror crossing-layer
counts, retaining loop multiplicity only at the endpoint-incidence marginal.
Record every remaining collision rather than replacing a multigraph fibre by
simple adjacency.

THM-809 proves the `n=8` row by the strictly weaker lower-only key

```text
Lambda=((A,A'),(B,B'),(C,C'),UABC,S2),
```

where all six nodes are `n=7` face nodes.  No upper `n=8` node label is used.
Its collision-excess ladder is

```text
B3 lower nodes + colour                 418
+ tau=3                                 252
+ tau=4                                 148
+ tau=5                                  74
+ tau=6                                  52
+ tau=7                                   0
+ fixed tau=8                             0.
```

Every collision before the last step is a double collision.  Thus the new
size-three crossing layer `tau=7` is the first complete separator, while the
fixed layer is redundant after it.  The literal B3 lower-line triple is also
injective on all `2^20` lines, independently checking the Cech descent.

The all-size claim remains deliberately open.  Even an injective finite codec
need not be a continuation-complete state: THM-796 already proves failure of
strong node lumpability.  A recursive theorem would have to transport the
sidecars and their action on the rooted Hamiltonian-path fibre, not merely
name a static signature.

## Exact completion and minimization target

For a literal core `c`, write THM-796's reconstruction variables as
`(p_L,p_H,a)` and define

```text
I_(u,c)(p_L,p_H,a)=1[reconstruct(c,p_L,p_H,a) in F_n(u)].
```

The full anchored Boolean Möbius transform of this function is invertible and
therefore continuation-safe if its substitution action is retained.  The
research question is whether `Omega+B2` identifies a small closed sector of
these coefficients.  Static equality must be tested separately from equality
under future deletion, complement, reflection, and lift operations.

The first undecided finite size is now `n=9`.  The exact next question is not
merely whether the codec collides, but whether the decisive `tau=n-1` layer
continues to remove every lower-node collision and whether its positional
moment, rather than raw counts, first becomes necessary.

## Continued-fraction connection

THM-778 supplies the exact analogy and the guardrail.  THM-808 now supplies an
actual transported coefficient: on a degree-one prime-sheet fibre, a centered
Christoffel owner-count block acts on the duplicate root by
`d'=d-sum c_a w_a^(-1)`.  The same block and same least-path mask can still
have two different target masks, so the root/owner lift is genuinely needed.
A continued-fraction
digit is useful only with the induced substitution on its labelled token
fibre; a metagraph address is useful only with the induced action on its
path/core stalk.  HYP-6880 therefore asks for a transported address, not a
larger tuple of static scalars.  MPA-38 proposes the first direct test: apply
centered Christoffel substitutions to the low/high leg variables and their
Möbius coefficients.

THM-812 now performs that direct test for the first nontrivial centered word
`d=(1,2)`.  Its coordinate-copy action `X_5->X_6` is complement/reflection
equivariant.  Bare nodes are not functorial, but every one of the 20 projected
coloured edge cells has a unique, distinct target cell.  On the coefficient
side, the induced general law is subset-image pushforward; degree at most
three loses one fixed-core target-node pair and three explicit quartic
coefficients repair it.  The remaining open transport claim is therefore not
whether a CF digit can act at all, but whether a finite coefficient sector
containing those quartics is closed under arbitrary Euclidean words and the
`Omega+B2` recursion.

Assumption challenge: the useful third face is not induced deletion of a
tournament vertex.  Its vertices are gap-contracted interval roots.  It
preserves a valid lower tournament tiling but destroys literal induced-
subtournament ancestry; that loss must be stated whenever the face is pushed
to merged nodes.

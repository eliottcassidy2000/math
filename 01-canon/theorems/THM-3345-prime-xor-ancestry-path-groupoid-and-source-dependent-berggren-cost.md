---
id: THM-3345
title: "Prime-XOR ancestry path groupoid and source-dependent Berggren cost"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  In the first
  spine-selected affine V4 fibre, the two edges of every prime-XOR matching
  have different Berggren ancestry costs and reduced path labels, so neither
  the prime nor its folded weight defines a source-independent ancestry word.
  The exact endpoint paths instead form a flat source-dependent groupoid;
  ambient-tree holonomy and graph H^1 vanish.  Rooted automorphisms and fixed
  descendant words are excluded, while arbitrary source-reading finite-state
  transducers, LRC transport, and JC flux remain open.
source: codex-2026-08-12-prime-xor-ancestry-groupoid
audit: >
  An independent hostile audit rederived all six longest-common-prefix
  lengths, path distances, and depth jumps; checked the constant-label,
  rooted-automorphism, flat-square, ambient-holonomy, frozen-basepoint, and
  finite-state scope boundaries; and reproduced the normal, optimized, and
  stored transcripts after LF normalization.  Both recorded hashes match.
depends_on:
  - THM-3334-berggren-parabolic-spine-gaussian-collision-torsor
  - THM-3336-primitive-gaussian-multiplication-content-curved-farey-triangulation
related:
  - THM-2606-affine-v4-parity-channels-partial-cubes-and-feuerbach-origin
  - THM-2622-affine-torsor-holonomy-fixed-section-spectrum-and-v4-c13-dictionary
  - THM-2753-six-edge-parity-erasure-and-three-matching-resolvent-restoration
  - THM-3339-fibonacci-three-ray-berggren-transplant-and-moving-owner-obstruction
script: 04-computation/prime_xor_ancestry_groupoid_thm3345.py
output: 05-knowledge/results/prime_xor_ancestry_groupoid_thm3345.out
script_sha256: 98e19c42a2c893e3d44b47b5dff0d1f5c6202e17a9156726289b4d22132e8a11
output_sha256: 3ba3d22f42df234752fa037086a8410b17c0c929b9b00ff6b82576b0b7e9b71b
hash_basis: working-tree bytes (LF)
---

# THM-3345 -- prime-XOR ancestry paths are source-dependent

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

No literature-priority or global-novelty claim is made.  The payload is a
typed comparison of the external Gaussian collision torsor with paths in the
ambient Berggren tree.

## 1. Inheritance and the two carriers

[THM-3334](THM-3334-berggren-parabolic-spine-gaussian-collision-torsor.md)
proves that the first spine-selected affine `V4` collision fibre occurs at

```text
c=1105=5*13*17.
```

Choose the factor-allocation chart

```text
X={000,100,010,001} isomorphic to F_2^3/<111>,            (1)
```

where the three singleton directions are labelled by `5,13,17`.  The four
parents and their unique root-to-child Berggren addresses are

| bit | parent | word | depth |
|---|---|---|---:|
| `000` | `(943,576,1105)` | `DAUD` | 4 |
| `100` | `(817,744,1105)` | `UDUA` | 4 |
| `010` | `(1073,264,1105)` | `DADDD` | 5 |
| `001` | `(47,1104,1105)` | `U^22` | 22 |

The external `K4` joins equal-hypotenuse parents.  It is not a subgraph of
the Berggren tree: every nonempty Berggren child word strictly increases the
hypotenuse.  An external edge must therefore be transported by the unique
undirected path that rises to a common ancestor and descends again.

[THM-3336](THM-3336-primitive-gaussian-multiplication-content-curved-farey-triangulation.md)
gives the intrinsic folded weights of the three XOR directions:

```text
K_1105(5) ={5,221},
K_1105(13)={13,85},
K_1105(17)={17,65}.                                      (2)
```

These are symmetric edge colours.  They do not yet specify an ancestry path.

## 2. Exact ancestry paths and costs

Let `w_x` be the word of `x`, and let `r(x,y)` be the length of its longest
common prefix with `w_y`.  Write `bar(v)` for the reverse of a word with every
letter inverted in the free group `F(U,A,D)`.  The unique tree path is

```text
P(x,y)=bar(w_x[r:]) w_y[r:],
L(x,y)=|P(x,y)|=|w_x|+|w_y|-2r(x,y).                     (3)
```

Applying (3) to both edges of each matching gives the complete table:

| prime | folded weight | edge | `r` | path cost `L` | absolute depth jump |
|---:|---|---|---:|---:|---:|
| 5 | `{5,221}` | `000--100` | 0 | 8 | 0 |
| 5 | `{5,221}` | `010--001` | 0 | 27 | 17 |
| 13 | `{13,85}` | `000--010` | 2 | 5 | 1 |
| 13 | `{13,85}` | `100--001` | 1 | 24 | 18 |
| 17 | `{17,65}` | `000--001` | 0 | 26 | 18 |
| 17 | `{17,65}` | `100--010` | 0 | 9 | 1 |

The first spine collision `c=85`, with words `U^5` and `AD`, is a positive
control: its unique ancestry distance is `7`.

Three consequences are immediate.

1. A prime-XOR direction, even together with its folded weight, does not
   determine ancestry distance: the two costs in every matching differ.
2. It does not determine a fixed free-group path label, even up to reversal,
   because equal or inverse reduced labels have equal length.
3. No rooted Berggren-tree automorphism realizes an entire prime matching.
   Such an automorphism preserves depth, while every matching above contains
   an edge with unequal endpoint depths.

Still less can a nonempty positive Berggren word act on a fibre edge: it would
send one endpoint to a strict descendant with larger hypotenuse.  These are
obstructions to a constant path label, descendant-word action, and rooted-tree
automorphism.  They do **not** exclude an arbitrary finite-state transducer
reading the source address; one finite fibre cannot decide that global
language question.

## 3. The exact lift is a flat source-dependent path groupoid

For a direction `p in {5,13,17}` define the oriented arrow

```text
a_p(x)=P(x,x+e_p).                                       (4)
```

Reversing the external edge inverts the path, so

```text
a_p(x) a_p(x+e_p)=1.                                     (5)
```

For distinct directions `p,q`, uniqueness of paths in a tree gives the flat
square identity

```text
a_p(x)a_q(x+e_p)
 =a_q(x)a_p(x+e_q)
 =P(x,x+e_p+e_q).                                        (6)
```

The exact companion checks all `12` directed involutions and all `24` ordered
affine squares.  It also checks all `60` directed closed walks whose two,
three, or four nonterminal vertices are distinct in the external `K4`; after
substitution of the source-dependent paths, every one freely reduces to the
identity.

This is structural, not a finite accident.  The path groupoid of a tree is
thin: there is exactly one reduced path between any two vertices.  Therefore
every closed ancestry path is trivial and

```text
H^1_graph(Berggren tree; A)=0                            (7)
```

for every constant abelian coefficient group `A`.  The external `K4` itself
has first Betti number three, but its exact ancestry-path realization sends
all three graph cycles to null paths in the ambient tree.

## 4. Why freezing one arrow per colour creates false holonomy

If one forgets the source and freezes only the three arrows based at `000`,
the quotient relation `e_p+e_q=e_r` is no longer respected.  The three reduced
defect lengths are

```text
(p,q;r)=(5,13;17): 39,
        (5,17;13): 39,
        (13,17;5): 37.                                   (8)
```

These nonzero words are hostile controls for a colour-only representation of
`V4` in `F(U,A,D)`.  They are not intrinsic curvature or cohomology classes:
the products in (8) identify arrows with different sources before composing
them.  Restoring the source-dependent second arrow gives (6) and zero
holonomy exactly.

Thus the lawful structure is

```text
external XOR torsor --source-dependent flat functor--> Berggren path groupoid,
```

not a translation-invariant action by three fixed ancestry words.

## 5. Tournament, cohomology, and transfer boundary

The ancestry path orients only after an ordered source and target are chosen;
reversing them gives the inverse path.  It therefore adds no intrinsic
tournament orientation to the symmetric XOR matching.

Equation (7) also separates two cohomologies that cardinality can blur.
THM-3336's `tau` and `lambda_p` are genuine `H^1` classes on the primitive
Gaussian associate multiplication group.  They are not Berggren-tree classes.
Any proposed word-current or D5 flux class must retain a nontrivial quotient,
current, boundary, or other sidecar; ancestry in the tree alone cannot carry
it.

Nothing here proves an LRC current or owner, a JC flux, a finite-state
prime-toggle transducer, or a global bridge between the external `K4` and the
internal four-box of THM-3339.  The theorem identifies the exact missing
coordinate: source-dependent ancestry path, not another scalar edge colour.

## 6. Exact controls and reproduction

The companion independently rebuilds the `177` Berggren vertices through
hypotenuse `1105`, verifies the four addresses, all six path rows, the `c=85`
positive control, the three constant-label and depth hostiles, equations
(5)--(6), all closed-walk reductions, and the three frozen-source defects.

Reproduce with

```bash
python3 04-computation/prime_xor_ancestry_groupoid_thm3345.py
python3 -O 04-computation/prime_xor_ancestry_groupoid_thm3345.py
```

Both modes must match
`05-knowledge/results/prime_xor_ancestry_groupoid_thm3345.out` after LF
normalization.

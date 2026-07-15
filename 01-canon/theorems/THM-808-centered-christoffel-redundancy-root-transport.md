---
id: THM-808
title: Centered-Christoffel blocks act affinely on the prime-sheet redundancy root, while the metagraph mask quotient is not continuation-complete
status: PROVED (general prime-sheet root cocycle) + FINITE-EXACT (five block actions and mask no-go on the THM-778 ten-wall movie)
source: codex-2026-07-15-S13
depends_on: [THM-773, THM-778]
related: [THM-781, THM-801, HYP-6880]
verification:
  - 04-computation/continued_fraction_redundancy_root_transport_codex_S13.py
  - 05-knowledge/results/continued_fraction_redundancy_root_transport_codex_S13.out
  - 05-knowledge/results/continued_fraction_redundancy_root_transport_codex_S13.json
---

# THM-808 — the redundancy root is the transported CF stalk

Let `p` be an odd prime and let `p+1` labelled owner tokens occupy
`k_a in F_p`.  Suppose a chamber covers every residue and has exactly one
duplicated residue `d`.  If a centered-Christoffel block contains `c_a`
events of owner `a`, whose speed is nonzero modulo `p`, then after the block,
provided the terminal chamber again has one duplicate, its duplicate is

```text
d'=d-sum_a c_a w_a^(-1)  (mod p).                        (1)
```

Thus the centered mechanical word has a genuine affine action on the
degree-one redundancy root.  This action need not descend to the least-path
tiling mask or its metagraph isomorphism-class node.

## Proof of the cocycle

In a covered `p+1`-token chamber the token multiset is `F_p` plus one extra
copy of `d`.  Since

```text
sum_(x in F_p) x=p(p-1)/2=0  (mod p),
```

the duplicated residue is intrinsically

```text
d=sum_a k_a.                                             (2)
```

THM-778's prime-sheet skew product says that one event of owner `a` sends
`k_a` to `k_a-w_a^(-1)`.  Summing that update over all events in the block
and applying (2) at its endpoints proves (1).  Notice that no assumption is
made about coverage in the intervening chambers.

The same proof works for any odd modulus for which all owner steps are units
and the sum of the complete residue system is zero.  Primality is retained in
the statement because THM-773/778's present sheet application is over
`F_p`.

## The five exact actions in the ten-wall movie

For

```text
W=(108,169,143,213,206,197,30,162),       p=7,
```

THM-778 reduces the ten covered walls to five owner-count blocks and their
reflection.  Applying (1) gives:

| ordinary walls | individual events | occurrences | root translation |
|---:|---:|---|---:|
| 57 | 58 | `0,8` | `+4` |
| 301 | 306 | `1,7` | `+3` |
| 3 | 3 | `2,6` | `+3` |
| 24 | 25 | `3,5` | `+4` |
| 329 | 336 | `4` | `+3` |

The verifier computes each translation twice: from the adjacent duplicated
roots and from `-sum c_a w_a^(-1)`.  It also reconstructs every adjacent
chamber and checks (2).  All failures are zero.

This is the first exact centered-continued-fraction action on a nontrivial
metagraph stalk in the repository.  The partial quotients and centered parity
cocycle generate the event block; its owner-count vector is the abelianized
substitution, and (1) is the induced affine representation on `F_p`.

## A sharp quotient obstruction

The three-wall block has owner-count vector

```text
(0,1,1,1,0,0,0,0).                                    (3)
```

It occurs twice.  Both occurrences start from the same least-path mask
`31115`, but their target masks are respectively `14635` and `615`.  Hence

```text
(source mask, centered-CF block) -> target mask          (4)
```

is not a function even on this ten-wall skeleton.  The two source duplicate
roots are `4` and `6`; formula (1) sends them to `0` and `2`.  Adding the
owner-labelled redundancy root repairs this witnessed collision.

This is stronger than merely observing that the isomorphism node is constant.
It exhibits equal literal mask and equal CF block data with unequal futures.
The missing information lives above the mask in the owner-labelled token
fibre.

The repair is proved for the root observation, not claimed as a complete mask
automaton: `(mask,root,block)` happens to determine the nine recorded target
masks, but those cells are mostly singletons.  A general continuation theorem
must retain owner assignment, simultaneous blocks, inverse steps, global
carry, and the metric core-safe component as required by THM-778.

## Relation to the metagraph Möbius stalk

THM-801's corewise Boolean Möbius stalk and the present root cocycle are typed
differently:

```text
marked tournament path/core stalk     owner-labelled prime-sheet stalk
        |                                        |
        v                                        v
   Omega+B2 line address                 redundancy root d
```

A centered-Christoffel substitution acts canonically on the right-hand stalk
by (1).  It can act on the left only after a specified lift relates owners and
token assignments to a marked Hamiltonian-path presentation.  Appending a CF
digit directly to an isomorphism node is therefore not a well-defined action.
The correct combined object is a fibre product over the literal path orbit,
not a tuple of unrelated scalar invariants.

## Tournament Analysis and preservation boundary

The computation treats candidate preservation carriers—block length, owner
count, mask, root, and their joins—as tournament vertices.  The pairwise
observable is how many of the nine literal between-wall gaps each carrier
separates; retention and retention per description bit are the two gauges.

For the LRC application the richer vertices are owner events, Euclidean
blocks, token obligations, and path-orbit lifts, not runners or arcs.  The
root quotient preserves degree-one prime-sheet coverage and its affine block
transport.  It destroys which owner occupies each sheet, the least-path mask,
intermediate cover failures, metric wall positions, and the core-safe
component.  Those losses are why (1) is a rigorous transport coordinate but
not by itself an LRC(14) proof.

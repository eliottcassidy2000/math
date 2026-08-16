# The missing bridge is a primitive harmonic K4 face, not a matching tournament

**Historical research reflection, 2026-08-16.** The proof source is
[THM-3509](../01-canon/theorems/THM-3509-reduced-fraction-harmonic-k4-face-and-fibonacci-unit-cassini-ray.md),
not this note.

## Inheritance pass and no-novelty audit

The assigned archive was read theorem-by-theorem before deriving a new
statement. The exact inheritance ledger is:

| source | inherited mechanism | boundary left open for this session |
|---|---|---|
| THM-3333 | every reduced `m/n` gives Euclid's Gaussian square and a Farey light-cone point | no primitive harmonic `K4` face converse or matching-orientation location |
| THM-3339 | `W=(x,y,x+y,x+2y)`, the three Berggren spinor maps, edge products, Fibonacci unit Cassini, six edge-order states | treats the Fibonacci ray and its owner packet; does not classify all fractions by a primitive edge face |
| THM-3357 | branch Walsh/Horn collapse and the `A/C` leg-swap near miss | no global two-tree parity involution on all reduced fractions |
| THM-3379 | the Fibonacci `T4` mod-three projection | a projection, not the four/six carrier or its lost antipodal sign |
| THM-3382 | dual index/address harmonic bifurcation | no fraction decoder and no ancestry recovery |
| THM-3454 | the U-spine `t/(t+1)`, Lorentz/Farey metric, and root tie | a one-parameter boundary, not all of `(0,1)` or the golden interior ray |
| THM-3455 | sampled rank spectrum and cap obstruction | ranks do not reconstruct the arithmetic source |
| THM-3487 | two 24-state bundles and their `H1` repair classes | finite periodic connection data, not a full fraction-tree equivalence |
| THM-3497 | exact ancestry-language densities and `A,C` versus `B` boundary | density forgets the individual word and does not supply the harmonic edge face |

Two sidecars were also load-bearing. THM-2596 distinguishes the binary
Stern--Brocot tree from the ternary Berggren cross-section. THM-2753 proves
that the six-edge `S4` action is faithful while the three-matching quotient
has kernel `V4`. The closest corrected near misses were MISTAKE-394 (a tied
root is not a tournament) and MISTAKE-352 (a modular shadow is not an
action).

The genuinely new content is therefore narrow but structural:

```text
all reduced fractions
  <-> primitive positive recurrence windows
  <-> primitive positive harmonic K4 faces,
```

with an exact gcd decoder, a linear Pythagorean-current formula, an exact
two-tree parity/leg quotient, and a proof that Cassini sign is erased by the
three-matching quotient.

## Provisional tetrahedral atlas audit

The recent tetrahedral representation atlas was also audited read-only on
`origin/codex/session-turnpike-atlas-20260815-codex2`: the base atlas, the
Haar/XOR square-cycle refinement, and the missing/double directed-channel
refinement. At the final pre-promotion fetch used for this note, its tip was
not an ancestor of `origin/main`, so it is not a proved dependency here.

Its useful provisional decomposition is compatible with THM-3509:

```text
4 K4 vertices;
6 edge coordinates = vertices of the octahedron L(K4);
3 perfect matchings = antipodal two-edge quotient;
cut/face/H1 sidecars retain data the matching quotient loses.
```

The atlas does not subsume the theorem. It is a representation-level
decomposition of general tetrahedral data. THM-3509 identifies a particular
multiplicative rank-one edge locus, singles out the face `(01,02,12)`, proves
the Diophantine converse on that face, and recovers the four recurrence
vertices by gcds. Conversely, THM-3509 does not prove the atlas's general
linear reconstruction or directed-channel claims.

## Incoming THM-3506 face-orbit sidecar

The user supplied a precise incoming connection after the first theorem
draft. THM-3506's formal state update

```text
(e,m) |-> (7e-2m,3e-2m)
```

becomes, under `(x,y)=(e-m,m)`,

```text
(x,y) |-> (4x+4y,3x+y).
```

This is a real connection. Every coprime odd THM-3506 pair becomes a
primitive even-`x` seed, its harmonic face is

```text
(m(e-m),e(e-m),em),
```

and its raw current is `(e^2-m^2,2em,e^2+m^2)`, of content two. Halving gives
exactly THM-3506's displayed primitive triples. The verified/exposed pairs
`(7,3),(43,15),(271,99)` give seeds `(4,3),(28,15),(172,99)` and primitive
triples `(20,21,29),(812,645,1037),(31820,26829,41621)`.

The proposed fixed-word obstruction required repair. In the present seed
chart, the branch determinants are

```text
det(A),det(B),det(C)=(1,-1,1),
```

so branch words are not all determinant `+1`; they are unimodular with sign
set by the number of `B` letters. The THM-3506 matrix has determinant `-8`.
If it were a rational projective scalar times a branch word, the scalar
square would be `-8` or `8`, neither a rational square. The true obstruction
is therefore the determinant square class `[-2]`, not orientation reversal.
The first verified jump is even sharper: ancestry words `B -> A^6B` show it
is not a descendant move.

Only the fixed verified/exposed pairs are used unconditionally. Iterating
the matrix as a full face-packet orbit remains conditional on THM-3506's two
renewal faces. The sidecar embeds that conditional orbit in the even tree's
vertex set; it does not turn the matrix into a Berggren action.

## Portfolio and concept board

The session portfolio was:

```text
Anchor:   classify the all-fraction Euclid/Farey/Berggren bridge exactly;
Niche:    determine the honest 4/6/3 tetrahedral carriers and quotient loss;
Wildcard: look for an Egyptian/harmonic face with an integral converse.
```

The wildcard became the theorem. The live concept board stabilized at seven
objects:

| object | representation/invariant | operation | lost coordinate | decisive test |
|---|---|---|---|---|
| reduced fraction `m/n` | seed `(x,y)=(n-m,m)` | mediant | none | recover `m/n` from edge gcds |
| four-term window | `(x,y,x+y,x+2y)` | componentwise addition | Farey flank labels | Stern--Brocot depth 12 |
| harmonic face | `(xy,x(x+y),y(x+y))` | reciprocal splitting | common scale | primitive converse through `v,z<=600` |
| Euclid triple | linear current `(u+v,2z,v+2z-u)` | normalize content | leg order and parity lane | `1/2` versus `1/3` |
| octahedral six-edge carrier | `L(K4)` plus numeric weights | edge comparison | rational vertex square class | all six weights equal to two |
| golden ray | `kappa=e03-e12` | Fibonacci mediant recurrence | sign under matching quotient | `k=3` versus `k=4` |
| THM-3506 face state | `(e,m)` or even seed `(e-m,m)` | determinant-`-8` update | renewal faces / branch action | projective determinant square class |

Each pull changed the other lanes. The harmonic identity

```text
vz=u(v+z)
```

initially looked like a pleasant Egyptian-fraction shadow. Its primitive
converse made it useful: writing `v=da,z=db` forces
`(u,v,z)=k(ab,a(a+b),b(a+b))`, and primitivity forces `k=1`. The three pairwise
gcds then recover `(x,y,x+y)` exactly. That decoder turns a visual tetrahedral
face into an arithmetic equivalence.

## Why the parity sidecar is exactly two trees

The first spinor coordinate `x=n-m` is invariant modulo two under all three
Berggren maps. Exact inverse descent partitions the positive primitive seeds
into roots `(1,1)` and `(2,1)`. The first tree gives raw primitive triples;
the second gives raw content two.

The involution

```text
(x,y) -> (2y,x)       for x odd,
(x,y) -> (y,x/2)     for x even
```

is the reduced form of `r -> (1-r)/(1+r)`. It swaps the two roots, conjugates
`A<->C` and fixes `B`, and swaps the two Pythagorean legs after normalization.
Thus all fractions map two-to-one onto primitive triples with unordered legs.
This is not a numerical coincidence at `(3,4,5)`; it is the exact quotient.

The ancestry loss has two levels. The unordered triple retains only the
orbit of two branch words under root swap and `A<->C`. If a word is further
abelianized to branch counts, even order is lost: chronological `AB` gives
`3/8`, while `BA` gives `5/8`.

## Why four and six are honest, but three loses Cassini

The carrier question has a clean answer.

- The four recurrence entries are the labelled `K4` vertices.
- The six pair products are vertices of `L(K4)`, the octahedron. Its
  share-an-endpoint relation is intrinsic even when numeric values tie.
- The three perfect matchings are the three antipodal fibres of that
  octahedron.

A comparison tournament is an additional structure, not the carrier itself.
Away from `1/2`, all four vertex values and all six edge weights are distinct,
so numerical `<` supplies transitive `T4` and `T6` relations. At `1/2`, ties
remain ties. This directly answers THM-3510's warning that six edge-objects
need fifteen comparisons before they can be called a tournament.

On the golden ray, the only changing edge comparison is

```text
e03-e12=(-1)^k.
```

But `{03,12}` is one perfect matching. Passing from six edges to three
matchings identifies the two terms whose difference carries the sign. The
matching quotient locates the relevant fibre but cannot orient it. The
recurrence labels provide the needed two-sheet gauge. An isolated `03/12`
swap fixing the other four edges is not induced by `S4`.

## Hostile controls and stopping boundaries

The hostile ledger is deliberately orthogonal:

1. **Parity:** `1/2` gives `(3,4,5)`; `1/3` gives `(8,6,10)`, which normalizes
   to the leg swap of the same primitive triple.
2. **Ancestry:** `AB` and `BA` have equal letter counts but give `3/8` and
   `5/8`.
3. **Orientation:** windows `(1,2,3,5)` and `(2,3,5,8)` have the same
   transitive `T4` order but Cassini signs `-1,+1`.
4. **Primitivity:** `(2,4,4)` is harmonic but scaled and has no primitive
   integral vertex decoder.
5. **Toric lift:** six edge weights all equal to two satisfy every opposite-
   product equation but force a rational vertex square to equal two.
6. **Operation:** mediants add vertex windows; edge products require the two
   mixed polarization terms.
7. **Conditional face orbit:** the update has determinant `-8`, whereas every
   Berggren word has determinant `+/-1`; neither determinant ratio is a
   rational square. The all-level orbit still requires THM-3506 renewal.

These controls also state the stopping boundary. A matching packet without
its oriented two-edge fibre cannot recover Cassini. Opposite products without
a square class cannot recover rational vertices. An unordered Pythagorean
triple cannot recover which parity tree supplied it.

## Exact verification and scope

The standard-library companion checks 17,543 reduced fractions, 19,682
generated Berggren nodes, 4,095 Stern--Brocot nodes, 2,392 Farey-neighbor
edges, 255 primitive harmonic triples in the stated box, Fibonacci indices
through 120, THM-3506's fixed/exposed pairs, and the hostile controls above.
Normal and optimized execution agree. Script/output/semantic LF SHA-256 are

```text
c12dc71df13e0c627740aaebf971c07342ce3302d674b276cc1600cf69266dd9
3f31f6f80989b360227c792d525e013ecff9b06a8df09d252d31af801d8fdb9a
7ca802f1f8706658b62fb330504544b9c46d473db6481429cb011b078695a46e
```

No META-PATTERNS promotion is proposed: the successful moves are already
covered by the repository's primitive-converse, quotient-loss, and hostile-
sidecar cards. Most importantly, this theorem gives no physical time, owner,
current, or exclusion map for LRC(14), and no polynomial map or flux map for
the Jacobian conjecture. The shared words *tree*, *current*, *Cassini*, and
*matching* are grammar only until such a map is written.

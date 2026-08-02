---
status: VERIFIED-EXACT RESEARCH REFLECTION; NO NEW PHYSICAL OR CANONICAL THEOREM
source: codex-2026-08-01-modular-v4-c3-point-direction-decoder
script: ../04-computation/modular_v4_c3_two_prime_point_direction_decoder_referee.py
output: ../05-knowledge/results/modular_v4_c3_two_prime_point_direction_decoder_referee.out
script_sha256: 402a6b562126abd3a28d026dc59b68c62e611b0dc556392d240c1d8e11f0844c
output_sha256: 6fffe9a9fdeb9d3be9af9e44ceed394b020bea2f21c397e6920cbe80d537cf06
hash_basis: LF-normalized bytes
---

# V4 points and C3 directions: an exact two-prime decoder

This is an exact scratch corollary/decoder, not canon.  Its closest proved
inheritance is THM-2595's affine `V4 semidirect S3` action and THM-3045's
six-edge integral binary/ternary clutch.  The new datum here is the literal
tensor decoder on four affine points times three matching directions; it
does not supersede either theorem.

Let `V=F_2^2`.  Its four affine points are the four prospective quartic
sections, while its three nonzero directions are exactly the three perfect
matchings of `K_4`: direction `d` pairs `x` with `x+d`.  Thus

```text
AGL(V)=V4 semidirect GL_2(F_2)=V4 semidirect S3=S4
```

acts equivariantly on the four points and three matching channels.  The
marked matrices from THM-2595 give the quotient `C2*C3 -> GL_2(F_2)`.
The tensor carrier below is

```text
V times (V\{0}) = twelve directed-edge flags,
```

and maps two-to-one to the six undirected edges used by THM-3045.  It must
not be identified with that six-edge lattice.

The integral Walsh decoder on the four points has Smith form

```text
diag(1,2,2,4), det=16.
```

Hence it is invertible after adjoining `1/2`, while modulo `2` it becomes the
all-ones rank-one matrix.  Point mass and point augmentation cannot be split
by this canonical Fourier direct sum in characteristic `2` (noncanonical
vector-space complements of course still exist).  These point-Walsh
coordinates also require the chosen affine origin.

Over the Eisenstein integers, the `C3` Fourier decoder on the three
directions has determinant an associate of `(1-zeta_3)^3`, with norm `27`.
It is invertible after adjoining `1/3`; modulo `(1-zeta_3)` it too becomes
all ones and has rank one.  Direction mass and direction augmentation cannot
be split by the canonical constant/augmentation direct sum in characteristic
`3`; this is not a denial of noncanonical complements.

This direction transform is canonical only after marking the order-three
generator, equivalently choosing a cyclic orientation of the matching
triple.  The full `S3` does not act linearly on that marked Eisenstein basis:
each reflection reverses the cycle and acts semilinearly by
`zeta_3 -> zeta_3^(-1)`.  The integral exceptional-prime statement is
`S3`-invariant; a labelled Fourier coordinate is not.

The Kronecker point-direction decoder therefore has determinant norm

```text
2^24 3^12.
```

Indeed `det(A tensor B)=det(A)^3 det(B)^4`; taking Eisenstein norms gives
`16^6*27^4=2^24*3^12`.  After inverting `6`, both the four-point/parity grammar and the
three-direction/matching grammar split simultaneously.  At the two bad
primes they fail for different, complementary reasons.  This is a precise
sense in which `2` and `3` are the two faces of the same affine `S4` object.

For a K4 hafnian amplitude, the three matching monomials live in the
direction module.  In characteristic `3`, `(1,1,1)` lies both in the trivial
line and in the augmentation kernel, so a scalar matching sum cannot choose
a complementary direction jet.  In characteristic `2`, the four-section
Walsh transform cannot choose a point-origin jet.  THM-3058's valuation-face
augmentation is compatible with this picture but does not itself provide
either missing decoder.

The prime-two point-Walsh collapse is not THM-2595's global
`H^1(C2*C3,V4)` compatibility bit, and it is not THM-3045's binary
`F2[three matchings]` opposition clutch.  Those are distinct modules and
distinct quotient losses which merely occur on the same affine object.

The perfect-matching source in
`05-knowledge/reference/CORE-PAPERS-COMPOSITIONAL-RELATIONS.md` makes one
external analogy literal: a coloring amplitude is a sum of products over its
perfect-matching fibre, and `K4` has exactly these three matching products.
The napkin-ring/Cavalieri analogy says only that a fibrewise invariant can
forget an ambient radius; partial-cube and graceful-labelling pictures say
only that binary cuts or edge labels can organize a carrier.  None supplies
the common amplitude-preserving map demanded below.

Scope is deliberately abstract.  THM-2596's boundary remains: the
isomorphisms `PSL_2(F_2)=S3` and `S4/V4=S3` do not canonically identify a
physical LRC, Keller, or quartic-origin carrier.  A transfer still needs a
map on one common object which preserves the required amplitude or owner.

## Exact verification and audit

Run

```text
python 04-computation/modular_v4_c3_two_prime_point_direction_decoder_referee.py
python -O 04-computation/modular_v4_c3_two_prime_point_direction_decoder_referee.py
```

Both modes LF-byte-match the stored transcript.  An independent reconstruction
recomputed the Walsh Smith form, Eisenstein determinant, Kronecker norm, all
24 affine point actions, all six induced direction permutations, all 72
matching-equivariance checks, and the two-to-one 12-flag/six-edge map.  Its
verdict is that this is a useful exact reflection/corollary: THM-2606 owns the
affine object and THM-3045 owns the integral six-edge clutch; the new invariant
here is the explicit oriented-flag tensor discriminant.

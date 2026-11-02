# Verification, scope, and implementation boundary

## Commands

```
python3 -B local_universality.py > normal.json
python3 -B -O local_universality.py > optimized.json
cmp normal.json optimized.json
python3 -B operators.py > probe.json
sha256sum -c SHA256SUMS
```

All code uses only the standard library, with no repository imports and no
assertions removed by `-O`. The normal and optimized main receipts match
byte-for-byte.

## Evidence layers

1. **Actual graph foundation:** independent subset/partition and literal
   permutation/cut engines check the bilinear cyclic polynomial identity
   on16 actual factor pairs at11 evaluation points, and check all four
   anchor H/Q2 pairs.
2. **Actual polynomial semantics:** 1,950 checks evaluate each gate and
   positive inverse through full polynomial functions modulo M=2..40,
   independently of its claimed small matrix. Reset sidecars are checked.
3. **Complete positive recipes:** 226 cases, including mixed moduli35,65,91,
   120,144, replay the whole constructed expression at the polynomial-function
   level. Every inverse is expanded into a nonnegative repetition of its
   genuine primitive operation, rather than applied as a formal inverse
   matrix.
4. **Compact controller arithmetic:** 24,688 simultaneous CRT landings cover
   every compatible target for M<=256 and selected larger primary/mixed
   moduli. These large checks do not literally expand the enormous graph or
   inverse words. The all-modulus theorem supplies that finite expansion.
5. **Hostile controls:** E kills the nonzero vector(2,3) modulo13; R/E fix
   (1,1) modulo5. The construction avoids the first kernel by preconditioning
   and avoids the second orbit by forcing an actual H5 seed factor.

The finite unit-group base has all24 units modulo72. The all-exponent lifting
argument is in the proof, not extrapolated from numerical tests. The width
constraint is compositional: every primitive graph operation is J or
C3[T,U,1], with all auxiliary factors built in the same grammar.

## Independent review

Two separate reviewers accepted the complete all-order proof, including the
unit lift, simultaneous CRT control, actual positive inverses, exceptional-
prime guards, final strong extraction and unbounded fibres. Both pin the
final proof hash `b7da98667ef7547dd8814d522363cde9c7209f39283bfe4bce7a41cb434d4bf2`.

- `independent/positive_gates/WIDTH_UNIVERSALITY_AUDIT.md` accompanies an
  independent implementation with121 literal factor pairs,4 literal seeds,
  3,562 simultaneous landings,132 full positive recipes,33,705 primitive
  applications and2,566 inverse requests. It imports no producer code.
  Its normal/optimized receipt hash is
  `523f754eeb2ed25a23918d7c064fec66c12eeedf78f531bc3a649768c0ea7987`.
  Both scripts were rerun after relocation here and reproduced that receipt.
- `independent/AUDIT.md` gives a second symbolic review, explicit positive
  products generating the mod72 unit base, and a separate scalar checker
  with24,832 landings and11,394 safe R/E steps. Its normal/optimized receipt
  hash is `910b9ec30ff8fd527b98e4dd9ee874c364545725cfdcc392f373edec4b7f83e9`.
  A fresh replay also matched. This checker explicitly does not claim actual
  graph replay; the first independent implementation supplies the positive-
  polynomial layer.

## Research lineage

The first unary P/C/D orbit probe found every residue for M<=128 but exhausted
its20,000-state limit at some larger powers. No negative conclusion was drawn.
The universal proof instead arose from the bilinear polynomial law, actual
closed two-coordinate gates, and explicit exceptional-prime control.

A first attempt to infer full mod5 transitivity from nonzero vectors was
rejected: diagonal vectors are fixed. Forcing H5 with positive exponent
repairs the source carrier. E's determinant-13 singularity is likewise kept
visible and avoided, not divided away. Intermediate outputs are not strong
because of neutral padding; only the final padding is removed after every
continuation has ended, preserving the exact H and revealing a strong root.

No claim is made about minimal witness order, practical expanded graph size,
all-context profile universality, exact odd-spectrum completeness, or global
boundedness of exact spectral substitution width.

All work remains in the assigned Phase2 directory, below1MB. Phase1 was not
changed. No git mutation, fetch, commit, or push was performed by this child.

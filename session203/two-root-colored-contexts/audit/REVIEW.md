# Independent audit of the two specified colored repair contexts

6 November 2026. **PASS: exact all-order identities and equality classification, restricted to the two displayed physical core pairs.** No claim about all two-root contexts, free target coverage, or an eventual odd-spectrum threshold follows.

## Frozen input and physical operation

The reviewed proof has SHA-256 9af8664c047d96410dd2395fa3ffb230185da4a0ae3da06e40525835fe201fa1. The producer verifier has SHA-256 1260bfd79686bd0c9c2e813d43b190c9331d10c951e6c59f2ec9775f1af9c96a. The producer's exact tensor receipt, ordered tensor receipt, and input manifest are retained. Both independent scripts use explicit failure checks that remain active under Python optimization.

For A, the parent core is (2,24,3,37,13,23), with body masks (15,63), and the child core is (2,8,67,101,15,87,19), with masks (79,127). Exactly old arc 1->4 is reversed, and vertex 6 is born. It loses to every body vertex. Old vertices 4 and 5 are the two color-sensitive vertices.

For B, the parent core is (32,5,1,39,47,70,31), with masks (0,5), and the child core is (0,5,129,39,47,199,159,27), with masks (0,5). Exactly old arc 0->5 is reversed, and vertex 7 is born. It beats every body vertex. Old vertices 0 and 2 are color-sensitive.

In both families every old body edge and every retained core/body incidence is unchanged. The count identities use B's full eight-vertex intermediate child. Old vertex 0 is not a universal sink when a color-0 body vertex is present, so the earlier monochromatic deletion cannot be silently applied. No strong-output assertion is needed or made.

## Why the arbitrary-order identity is certified

The variables x,y,z,w are functions counting nonempty body Hamilton paths with endpoint colors 00,01,10,11. Products are subset convolution. The multiplicative unit is supported only on the empty body. Thus H=1+x+y+z+w, including H(empty)=1. A repeated variable still denotes distinct labeled path slots, not an unlabeled partition or an ordinary scalar power.

Every Hamilton path in a context has a unique core order and a unique set of occupied gaps. Each occupied gap contains one nonempty directed path in the body. Its first and last colors determine its permissible core incidences. Conversely a permissible core order and a labeled spanning path cover reconstruct one Hamilton path. This bijection works for arbitrary internal body orientations, not just a transitive spine.

The independent tensor checker enumerates permitted ordered endpoint-type choices in all 6!+7!+7!+8!=51,120 core orders, before commutative aggregation. It reproduces every parent, child and excess coefficient, checks the separately supplied ordered C++ excess tensors, and verifies the two displayed polynomial factorizations exactly. This finite enumeration exhausts the fixed cores and therefore certifies the identity for every body order. It is not an interpolation from selected body counts.

## Cut positivity and equality

Adjoin an apex q which beats color 0 and loses to color 1. The exact defect is

    W = H(T+q)-H(T) = y-z+xy+xw+y^2+yw.

The empty body gives zero on both sides. For nonempty bodies, splitting at q gives the claimed endpoint-cover formula. If Q_ab counts ordered spanning two-path covers whose first path ends at a and whose second starts at b, then the path-flow identity gives

    W(T) = sum over a->b of (Q_ab+Q_ba) epsilon_a(1-epsilon_b),

where epsilon is 1 on color 0. Every capacity is strictly positive: a Hamilton path after removing a,b can be split immediately before the first vertex beaten by either endpoint. Attach that endpoint before the suffix and the other after the prefix. The no-such-vertex and empty-remainder cases also produce a required cover. Hence W is nonnegative, and is zero exactly when no edge points from color 0 to color 1. This applies to every induced colored subgraph, which is essential for convolution positivity.

Set F=y+xy+y^2 and J=y+y^2+yw. The independently checked identities are

    D_A = (4+x)W + (x+z)F + y(4+3x+4y+3z+4w),
    D_B = (21+42y+20w)W + 20(z+w)J + 22yH,

where D_i=H(C_i)-2H(P_i)-H. All other factors have nonnegative coefficients. Thus D_A>=4W and D_B>=21W. Any color-0 to color-1 edge makes the excess strictly positive. If all cross edges point from color 1 to color 0, both W and y vanish on every induced subgraph, so every summand vanishes. The claimed iff equality boundary is therefore sound, including empty bodies and empty color classes.

## Independent literal controls

The body checker uses nonmemoized directed-path enumeration, explicit labeled-bin assignments for convolution monomials, and literal permutation/cut enumeration for Q_ab. This differs from the producer's subset dynamic programs and recursive convolution evaluator.

- All 1,099 labeled colored tournaments through order 4 pass the cut identity and equality test. All 410 unordered endpoint pairs in the underlying uncolored universe have positive capacities.
- Canonical representatives are independently generated by testing every relabeling. The uncolored class counts at orders 0 through 4 are 1,1,1,2,4. All 87 colorings of these representatives are checked in both physical contexts.
- Twelve seeded order-5 inputs, two explicit nontransitive order-6 inputs, and the unsorted two-vertex spine are also checked.
- In total, 204 literal parent/child pairs agree with the tensor identity and exact equality test. There are 152 strict cases. The literal path counter visits 15,477,646 directed prefixes.
- The unsorted spine with rows (2,0) and colors (0,1) gives A: parent 133, child 279, excess 12; B: parent 115, child 295, excess 64.
- For body C3->C3, H(body)=9. Colors (1,1,1,0,0,0) give A: 1413->2835 and B: 765->1539, both with zero excess. Reversing the colors gives cut defect 90 and strictly positive excesses 864 and 15,192. Thus equality does not require either color block to be transitive.

Both independent normal/optimized receipts agree byte-for-byte. The producer verifier is also replayed normally and under optimization, with exact agreement with its frozen receipt.

## Scope and use

On a transitive spine, equality holds exactly for words 1*0*. In these surviving cases A's color-1 block is a global initial module and B's color-0 block is a global terminal module, explaining the ordered-product mechanism. Other words do not give the intended unit correction through these two fixed operations.

This result does not exclude other two-root cores, old body-edge surgery, target-dependent contexts, or a construction that controls the positive excess. It establishes no missing scalar count and no all-order threshold. The positive excess may be useful as an additional coupled coordinate, but this audit does not provide a coverage or iteration theorem for it.

## Reproduction

Run in this directory:

    python3 check_tensors.py > TENSORS_NORMAL.json
    python3 -O check_tensors.py > TENSORS_OPTIMIZED.json
    cmp TENSORS_NORMAL.json TENSORS_OPTIMIZED.json
    python3 check_bodies.py > BODIES_NORMAL.json
    python3 -O check_bodies.py > BODIES_OPTIMIZED.json
    cmp BODIES_NORMAL.json BODIES_OPTIMIZED.json
    python3 PRODUCER_VERIFIER.py > PRODUCER_REPLAY_NORMAL.json
    python3 -O PRODUCER_VERIFIER.py > PRODUCER_REPLAY_OPTIMIZED.json
    cmp PRODUCER_REPLAY_NORMAL.json PRODUCER_REPLAY_OPTIMIZED.json
    cmp PRODUCER_REPLAY_NORMAL.json PRODUCER_TENSORS.json

The producer also writes VERIFICATION.json; it is a redundant byproduct. No external tournament bank or completeness assumption is used. Publication and any status promotion remain separate from this audit.

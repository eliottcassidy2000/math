# Higher dyadic energy forms: proof and independent controls

October 23, 2026. **PROVED:** the finite-matrix filtration theorem below.
**FINITE-EXACT:** the stated independent small-matrix controls, replay of
the producing consumer, and an independent reconstruction consuming its
saved literal lifts. No all-order bound on the largest
critical-group exponent is inferred from the smaller graphs.

## An elementary proof of every layer

Let A be a nonsingular symmetric integer matrix of size m. For j >= 1 put

    L_j = {x in Z^m : A x is divisible coordinatewise by 2^j},
    V_j = image of L_j in (Z/2)^m.

For u,v in V_j, choose any lifts x,y in L_j and set

    q_j(u,v) = x^T A y / 2^j modulo 2.

This is a well-defined symmetric bilinear form. Integrality follows from
Ay being divisible by 2^j. Changing x to x+2z changes the quotient by
`2 z^T (Ay/2^j)`, an even integer. Changing y is the same calculation using
symmetry. Addition of binary vectors may change a chosen lift by an even
vector, so the same calculation proves bilinearity.

**Theorem.** `rad(q_j) = V_(j+1)`. Moreover `dim V_j` is the number of Smith
invariant factors of A divisible by 2^j, and `rank q_j` is the number whose
2-adic valuation is exactly j.

Proof of the radical assertion: one inclusion is immediate. If u has a
lift x in L_(j+1), then symmetry makes every numerator x^T A y divisible
by 2^(j+1), so u is in the radical.

Conversely take x in L_j and write Ax=2^j b. A lift x+2z belongs to
L_(j+1) exactly when

    A z = -2^(j-1) b modulo 2^j.                       (1)

Over the finite group `(Z/2^j)^m` with its perfect dot-product character
pairing, symmetry gives

    image(A) = annihilator(kernel(A)).                 (2)

Inclusion follows from symmetry; equality follows by counting, since both
subgroups have cardinality `2^(jm)/|kernel(A)|`. Therefore (1) is solvable
exactly when `b^T y = 0 modulo 2` for every y satisfying Ay=0 modulo 2^j.
Those reductions are precisely V_j. By symmetry the condition is exactly
`q_j(u,v)=0` for every v in V_j. This proves the converse.

For the dimension statement use an integer Smith decomposition. Its
unimodular right change is invertible modulo 2. In a scalar coordinate with
Smith factor d, an odd residue admits a lift annihilated modulo 2^j exactly
when 2^j divides d. Each such coordinate contributes one binary dimension.
Subtract consecutive dimensions and apply the radical identity. QED.

For a connected graph use a reduced Laplacian. Its nonsingularity follows
from connectivity. For a disconnected graph apply the result separately
to each component after deleting one root in that component. Literal loops
contribute zero to the Laplacian and are not inferred from these forms.

In the notation of Session36, q_1 is the half-energy form on the binary
blue-harmonic space modulo constants. Its radical is V_2. Choosing lifts
with Ax divisible by 4 gives q_2, the quarter-energy form, whose radical
counts factors divisible by 8. The same rule continues at every layer.

## Independent small-matrix controls

`dyadic_controls.py` enumerates **all** potential lifts modulo 2^j for
j=1,...,4 for seven explicit matrices. Cases include diagonal matrices
with Smith factors `(2,4,8)` and `(1,2,16)`, their nontrivial unimodular
congruences, and the matrices

    [[0,2^e],[2^e,0]],  e=1,2,3.

These last controls have identically zero diagonal energy but a nonzero
bilinear form at their active layer, preventing an accidental replacement
of the full form by its diagonal. Every available lift choice is compared
against the chosen lifts, the radical is compared with the next explicitly
enumerated lift space, and the dimensions are compared with the known Smith
factors. Normal and optimized Python each pass **489 explicit checks**.

```text
python -B 05-knowledge/results/metagraph_allsizes_20261023/dyadic_controls.py
python -B -O 05-knowledge/results/metagraph_allsizes_20261023/dyadic_controls.py
```

## Actual order-eight consumer audit

The code in `metagraph_order8_electrical_20261023/check.py` constructs the
half-form from the complete blue/paired-black graph. Blue-isolated unit
potentials and a kernel basis of the active reduced blue block together
give a basis of V_1. The root row is redundant by the zero row-sum law.

Every subsequent radical basis is lifted as a literal integer potential.
All corrections are even, so the previous binary class is retained. The
consumer explicitly checks divisibility of the full Laplacian image at
every accepted lift and computes the next pairing by an integer dot product.
Thus even a failed intermediate search for a correction cannot silently
produce a certificate: final divisibility is checked, and all constructed
corrections preserve the source class. The independent consumer below also
checks the exact parity span of every saved basis, without relying on how
the producer found its corrections.

The audited optimized run passes **439933 checks** and reproduces receipt
LF SHA256 `417e4cc3af2688df39bbb71701789a6967677ae82ea8d0f0badb8ff50d739b27`.
It gives dim V_1=3363, rank q_1=3358 and rank q_2=4. One mode remains
through levels 3,4,5,6,7; q_8 has rank 1 and no radical. Therefore the exact
2-primary part of this finite graph's critical group is

    (Z/2)^3358 direct-sum (Z/4)^4 direct-sum Z/256.

Its contribution to the tree-count valuation is 3358+2*4+8=3374.
The large final factor shows why the maximum exponent 4 observed through
order seven cannot be promoted to an all-order theorem.

This review independently checks the matrix mechanism and replays the
finite graph consumer. Exhaustive tournament-to-graph export coverage is
inherited from the separately authenticated order-eight export packet.

## Independent saved-lift consumer

`dyadic_saved_consumer.py` imports no code from the producing packet. It
reads the complete pinned edge export and constructs the sparse integer
reduced Laplacian directly. Highest-bit Gaussian elimination constructs
the full mod-two kernel; this differs from the producer's active-blue-block
construction and lowest-bit pivots. Direct integer column sums and binary
dot products recompute every entry of the entire 3363-dimensional half-form.
Its independently computed radical has dimension five.

For every saved level from two through eight, the consumer checks that:

- each potential is a literal integer vector with the root coordinate zero;
- its full Laplacian image is divisible by the asserted power of two;
- its mod-two classes are independent and span exactly the previous form's
  radical, expressed in the underlying vertex coordinates;
- every saved form entry is the actual integer energy quotient modulo two;
- independent elimination gives the saved rank and next radical dimension.

No higher correction is searched for or reconstructed. In particular, this
is a consumer of the literal certificates, not a replay of the producer's
correction ladder. Saved lifts need not be bounded least-residue vectors.
A positive control adds four to one coordinate of a level-two lift and
preserves the certificate. Hostile controls reject a float alias, an image
divisibility mutation, an omitted independent mode, and a changed form bit.

Normal and optimized consumer runs pass **445217 explicit checks** including
comparison with the independently saved receipt (445216 certificate-body
checks). They give the following independently reproduced filtration:

| j | dim V_j | rank q_j | dim V_(j+1) |
|---:|---:|---:|---:|
| 1 | 3363 | 3358 | 5 |
| 2 | 5 | 4 | 1 |
| 3 | 1 | 0 | 1 |
| 4 | 1 | 0 | 1 |
| 5 | 1 | 0 | 1 |
| 6 | 1 | 0 | 1 |
| 7 | 1 | 0 | 1 |
| 8 | 1 | 1 | 0 |

**PASS.** Together with the proved filtration theorem, this independently
establishes for the authenticated finite order-eight graph
`(Z/2)^3358 direct-sum (Z/4)^4 direct-sum Z/256` and tree-count two-adic
valuation **3374**. It does not assert that every larger graph has this
exponent, or that electrical data alone reconstruct all colored ports.

```text
python -B 05-knowledge/results/metagraph_allsizes_20261023/dyadic_saved_consumer.py
python -B -O 05-knowledge/results/metagraph_allsizes_20261023/dyadic_saved_consumer.py
```

LF-normalized SHA256 pins at audit completion:

| File | SHA256 |
|---|---|
| `metagraph_order8_20261023/merged_edges.tsv` (decompressed) | `6687671346fbf83ac73bb2f952ffee7705602036f1d083b59b524ea35f4543b5` |
| `metagraph_order8_20261023/merged_nodes.tsv` | `61e4e25a0e83bbde3cdb756bd4a71febb4d0f2b8c5237143184e782738d02835` |
| producer `check.py` | `587f40058a9a9296b93fa0159f610491edf97852492dacfc5e65b02f01eebf08` |
| producer `receipt.json` | `417e4cc3af2688df39bbb71701789a6967677ae82ea8d0f0badb8ff50d739b27` |
| `dyadic_controls.py` | `875654e16956d79d4666ae1d73341414da8e1a6673320975626ec849cc63e6c8` |
| `dyadic_saved_consumer.py` | `4668a9e27a965fbc2901848da28f793b558138bf86cc7cd5309bd504599fcf03` |
| `dyadic_saved_receipt.json` | `7b8c0ed9c67453029fa23fe19916b801bf0bf057cddeef529f3c0d66ba07b524` |

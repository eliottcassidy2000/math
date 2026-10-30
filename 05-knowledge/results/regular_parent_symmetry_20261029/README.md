# Tournament symmetry across one-vertex growth

**Status: PROVED all-size deductions with FINITE-EXACT controls.
This packet extends credited prior recurrences; it makes no literature-priority
claim and does not solve a major open conjecture.**

The strongest result is a regular-parent law: **adding one vertex to a regular
tournament cannot create an automorphism.** Every automorphism fixes the added
vertex and restricts to the stabilizer of its attachment subset.

This turns rooted orbit counts into exact unmarked extension counts:

    number of child classes = number of subset orbits - number of vertex orbits.

The subtraction accounts precisely for cloning an old vertex. Children of
nonisomorphic regular parents never overlap. This gives an exact component
of A000568 enumeration, with group information retained, and removes child
isomorphism comparisons once the parent action and its subset orbits are known.
It is not yet a benchmarked speedup of the repository's fastest full generator.

The canonical augmentation rule simply rejects attachment masks
N_T^+(w) union {w}. If L_n is the labelled regular-parent count, the
labelled child count is exactly (n+1)(2^n-n)L_n.

| Regular parent | Parent group order | Rooted subset orbits | Unmarked children |
|---|---:|---:|---:|
| C3 | 3 | 4 | 3 |
| R5 | 5 | 8 | 7 |
| R7 | 7 | 20 | 19 |
| other regular seven-tournament | 3 | 48 | 45 |
| Paley seven | 21 | 12 | 11 |
| C3[C3] | 81 | 24 | 23 |

Summing the three seven-vertex rows gives **75 eight-vertex classes with a
regular card**. The ternary tower has the all-size stabilizer polynomial

    E_k(z)=P_k(z)-z^((3^k-1)/2-k),

where the rooted polynomial P_k was already proved in the repo. Child counts
at parent sizes 1,3,9,27,81 are 1,3,23,4623,32955830623. There are exactly two
rigid unmarked children at every ternary depth k>=2.

Two further results explain why symmetry size alone is insufficient:

- **Random named ports have a critical scale.** Triangle towers require about
  (log_2 n)/2 ports; carousel-five and Paley-seven towers require about
  (log_2 n)/4. At the critical scale their residual groups converge to products
  of a Poisson number of small cyclic groups. The controlling quantity is
  the minimum number of color equalities forced by a nonidentity base
  automorphism, not simply the order of the base group.
- **Symmetry can disappear while the extremal obstruction survives.**
  Explicit rigid Paley-eleven-free tournaments U_k have
  |U_k|=7^k+1 and largest transitive subtournament 3^k+1. The recursive
  attachment has weights (7^k-3^k)/2 and (7^k+3^k)/2, from a two-state
  matrix with eigenvalues 7 and 3. This is a guarded Erdős–Hajnal test
  family: removing automorphisms does not remove the difficult extremal
  geometry.

The 81/27 hostile pair retains identical 67 labelled cyclic triples but has
1620 versus 1728 binary patterns with trivial **core** stabilizer. These
twelve-vertex parents are not regular, so these last figures do not assert
that the corresponding unmarked thirteen-vertex children are rigid.
The exact orbit action still matters after the leading probability law agrees.

## Proofs and boundaries

[PROOFS.md](PROOFS.md) contains R1--R6, including proofs, equality cases,
and the nonregular 3 -> 1 -> 5 group-order counterexample to unrestricted
prime inheritance. [RECOVERY.md](RECOVERY.md) records prior sources,
the concept board, information losses, and the literature. These deductions
are kept in this result packet rather than assigned unreserved canon IDs.

Named ports, literal orientation, and the relevant deletion guard must be
preserved. The order-168 switching symmetry of the source-plus-Paley-seven
class is a different action from the literal tournament's group of order 21.
No IID model is claimed for arithmetic ports. The Poisson law controls bounded
observables; it supplies no missing weighted LRC estimate or all-source
Collatz return certificate.

## Reproduce

From the repository root, with standard-library Python only:

    python -X utf8 -B 05-knowledge/results/regular_parent_symmetry_20261029/verify.py
    python -X utf8 -B -O 05-knowledge/results/regular_parent_symmetry_20261029/verify.py --output <temporary-output.json>

Compare the latter JSON byte-for-byte with [certificate.json](certificate.json).
The 35,661 explicit truth gates remain active under optimization. The universe
includes all 32,834 labelled tournaments at even orders 2,4,6; all regular
parent classes through order seven; every attachment of those parents and the
nine-vertex ternary tower; seed colorings for Q=2,3,4; exact probability
recurrences for Q=2,4,8,16,32 through depth four; both 4096-cut hostile controls;
all module tests for P11 and its eleven cards, and for all P7 cards; native
automorphism tests of the explicit rigid examples at orders 6,8,26,50; and
the rigid ten-vertex counterexample where adding a port creates R5.

[manifest.json](manifest.json) pins inherited Git blobs and packet bytes.
[HANDOFF.md](HANDOFF.md) gives the next decisive questions.

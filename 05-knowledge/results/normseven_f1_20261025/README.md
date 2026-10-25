# Norm-seven arithmetic fibres, integral wrap torsion, and the F1 boundary

2026-10-25. **PROVED** statements are in the linked proof files;
finite exhaustive controls are labelled separately. These results preserve
the marked base point, branch labels, characters and arithmetic twists.
They do not settle integer-square hourglasses, Collatz, plane colouring,
planar Jacobian, or all-order tournament connectivity.

1. [Arithmetic fibre reconstruction and F1 coordinates](KUMMER_AND_F1.md).
   The two complementary norm-seven quotient covers recover every quadratic
   twist over the original field, with an explicit inverse. The universal
   root-algebra statement holds over every coefficient ring, even when its
   characteristic divides the root order. The inherited six-character
   monomial matrix has determinant three; changing one character gives
   determinant one. More generally every invertible binary character matrix
   has a unimodular Laurent lift preserving its individual characters.
   This is an actual torus coordinate change over pointed monoids, with
   explicit limits at the branch boundary.
2. [Integral wrap theorem](INTEGRAL_WRAP.md).
   The cyclic incidence kernel is `R^d + R[D]` over every commutative ring.
   Its integer cokernel has just one possible nontrivial Smith factor D.
   At prime length p, D is exactly p when the inherited extra mode survives.
   For composite lengths prime to three, a boundary-balance gcd gives the
   additional possibilities: explicit supports at lengths8 and16 yield
   D=2 and4, while the full length49 support yields D=49.
3. [Twist average and positivity obstruction](TWIST_AVERAGE.md).
   The 64-twist average of branch-free point counts is exactly the base
   count `q^2-6q+12` for characteristic>3, and `q^2-6q+10` in characteristic3.
   The first has formal value seven at q=1 but a negative coefficient in
   powers of q-1, excluding a finite split-torus decomposition. The positive
   multiplicative carrier must retain the additive branch map separately.
4. [Ports under the arithmetic basis change](PORTS_AND_CHARACTERS.md).
   The same D-character repair sends old cubic parity to a new pair
   character. It preserves the complete root fibre but changes which
   observations are called pair marginals, and sends one cube edge to a
   distance-two pair. Retaining the old port observables resolves the issue.

## Inheritance, concept board, and connections

The anchor is arithmetic descent of the norm-seven fibre product; the niche
is integral rather than field-valued incidence; the wildcard is the formal
seven in twist averages. Source pins and hashes are in [SOURCES.json](SOURCES.json).

Closest proved mechanisms: the two norm-seven complementary pairs, the
all-prime anchored-path kernel, and the Session59 negative-torus-coefficient
test. Canonical hostiles: a split geometric cover with a nonsplit arithmetic
twist; cyclic constants killed by an integer sum; and a negative coefficient
after a seemingly attractive specialization. The corrected near miss is
treating the determinant of one chosen character lift as an invariant of its
quadratic cover. The least-used relevant sidecars were base units and the
full integral coefficient module.

| Live concept | Exact map and preserved predicate | Lost data and necessary sidecar |
|---|---|---|
| Norm-seven complementary covers | Root monomials, K1-K4; same-source cover points and all twists | An unmarked product forgets the common base and branch labels |
| F1 multiplicative carrier | Unimodular lattice map; torus and root-algebra coordinates | Addition in the branch map and behavior where a root vanishes |
| Wrapped incidence | W4; every physical partition equation | Field rank alone loses Smith exponents; retain D |
| Signed tournament cancellation | Incoming exterior-invariant/Burnside polynomial | Root multiplicity is not the rank of legal fixed toggles |
| Twist average | Sum over constant square classes at each base point | Mean count does not select a soluble fixed twist |

Incoming commit `9807ae1678` was fetched, read, and merged during the session.
Its [clone theorem](../session62_agent_theme_bridges_20261025/clone/README.md)
gives actual tournaments with toggle rank1 and cancellation order3^(k-1).
In its nine-gap example, division by (1+t)^3 leaves negative coefficients.
That and the negative q-1 coefficient above are exact positivity obstructions
to two proposed decomposition interpretations. They do not identify the
tournament and arithmetic geometries. Similarly, the Smith invariant of a
fixed incidence map cannot be removed by an integral basis change, whereas
the determinant of an optional lift can change when radicand corrections
change the ambient map. The permitted operation is part of each claim.

Incoming `1f804be357` was subsequently read and merged. Its
[projective-eight map](../hadwiger_six_seven_eight_20261025/projective_eight/RESULT.md)
and [complete boundary theorem](../hadwiger_six_seven_eight_20261025/cube_boundary/RESULT.md)
give a precise new comparison: the very shear used by our easy character
repair moves the missing cubic observation into a pair observation, but
changes the ports and cube metric. The all-rank fixed-character Laurent
repair avoids that relabeling while retaining the source-unit corrections.

The scoped [audit and hypothesis ledger](AUDIT.md) records boundaries,
rejected guesses, and the next useful obligations. None is labelled an
independent agent audit: the independent checks here are computational paths.

## Reproduction

From this directory, with Python3 and SymPy installed:

```bash
python3 check.py > normal.json
python3 -O check.py > optimized.json
python3 audit_smith.py > smith_normal.json
python3 -O audit_smith.py > smith_optimized.json
python3 check_twists.py > twists_normal.json
python3 -O check_twists.py > twists_optimized.json
python3 check_ports.py > ports_normal.json
python3 -O check_ports.py > ports_optimized.json
python3 verify.py
```

The first checker uses only the standard library. It exhausts all120
complementary labelled Fano pairs and 94,080 Boolean character matrices.
The Smith audit constructs physical matrices directly and covers2,052
declared supports, including all nonempty supports at lengths5,7,8,10.
The twist checker directly counts roots over six prime fields and F9,F25;
the port checker exhausts the two parity tables on the eight states.
The saved receipts state exact universes and counts; normal and optimized
outputs are required to agree byte for byte. Finite data support the stated
proofs and do not substitute for their unbounded quantifiers.

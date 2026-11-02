# Verification

Run from any working directory:

    PYTHONDONTWRITEBYTECODE=1 python3 /path/to/check.py
    PYTHONDONTWRITEBYTECODE=1 python3 -O /path/to/check.py
    PYTHONDONTWRITEBYTECODE=1 python3 /path/to/face_bases.py
    PYTHONDONTWRITEBYTECODE=1 python3 -O /path/to/face_bases.py

`check.py` uses direct adjacency subset-DP for Hamilton paths and ordered
covers, independently of the compiler being tested. It checks:

- the four actual lattice witnesses and the two transitive axis families
  through parent order10;
- every cut on transitive parents through order10 against the exact unit-height
  zero-interior face law;
- all 1,098 labelled local tournaments of orders2..5 with the fresh mark last, including 75 zero-cost
  sealing recipes, for positivity, integrality, height parity, native output
  coefficients and rank-one sealing;
- all 74 labelled marked sources with parent order1..3 against all 74 local
  tournaments of order2..4 with fresh mark last: 5,476 independent literal compiler checks;
- all 156 words of length0..3 in a five-operation keep/re-mark palette on four
  actual sources: 624 literal state and original-port collapse checks;
- the explicit F3 neutral-padding hostile and its successful F5 alternative;
- two independent all-length capped closures at173 for the inherited 32-matrix
  bank: 6,234 positive states versus 2,592 source-cost states, including85
  overflow tags, with identical support and missing odds7,21.

`face_bases.py` independently computes SCCs and applies the proved deletion
normal form to every labelled parent of orders2..5 and every singleton cut
with g=0: 1,598 literal reductions. It also exhausts every strong block through
orders3 and5, sufficient at those fixed prime heights, and records actual
base recipes. The minimal base lists are h1:{1}, h3:{3}, h5:{5,6,7,9,16}.
These finite bases become all-order classifications only through the proved
core reduction, not through a chosen graph-order cutoff.

Hostiles reject a negative-column map which passes the four small actual
calibration sources; a nonnegative integer parity-preserving map which fails
actual face membership; treating the actual H11 source as if its L coordinate
could be erased; and treating H-neutral sink padding as a harmless source
normalization. Assertions are not used as proof checks; failures raise explicit
exceptions, so `-O` retains every check.

The copied `migration.py` is unchanged from S134 (SHA-256
`91b0e21515f8f590f9093f26be650f6fc8a6e89276251920e75dda6c801f41d2`).
Its universal compiler theorem is an inherited dependency, not inferred from
these finite controls. No external graph atlas, network access, random seed,
or third-party package is required.

Normal/optimized receipts and fresh-directory replay receipts are retained.
The public manifest excludes exploratory face-matrix searches and their large
output; those are not proof inputs.

## Bounded-interior strengthening

Run `bounded_interior.py`, also with `-O`. It checks all 32,767 nonempty cuts
on every labelled parent of orders1..5; 5,084 literal adjacency reconstructions
after SCC suffix reordering; 300 unselected-block commutation controls; the
complete h1,L1 base universe through its proved order bound4; and all cuts on
transitive parents through order10 (81 h1,L1 instances). The maximal L in the
complete parent<=5 universe is45. All complete native states, not only H,
agree after normalization and restored tails.

The exact h1,L1 language is `{1,3} union {3*2^k-1,5*2^k-1:k>=1}`. Its two
canonical actual native bases are `(1,1,3,1)` and `(1,3,5,1)`. This all-order
claim uses the normal-form theorem and a complete finite base universe, not
extrapolation from the order10 control. Normal/-O and relocated normal/-O
outputs are byte-identical.

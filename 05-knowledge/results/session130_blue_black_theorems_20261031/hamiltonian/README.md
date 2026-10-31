# Session 130 Hamiltonian/parity bridge

The new refinement is

`H(T) = 1 + 2q + 4i (mod 8)`,

where `q` counts reflection-fixed odd directed cycles and `i` counts
intersecting nonfixed mirror-cycle pairs. See [PROOFS.md](PROOFS.md) for the
all-order proof, central insertion/cut bridge, automorphism-normalized black
degree, exact order-five graph table, and failure boundary.

The strongest hostile says that even the entire length-decorated odd-cycle
conflict graph with reflection action does not recover the blue fibre.
Two rigid five-vertex tournaments have `H=9` and identical such data, yet
blue counts one and three. The path/reflection sidecar is essential.

Run from this directory (Python standard library only):

```sh
python3 check_hamiltonian_bridges.py --max-blue-n 8 > normal.json
python3 -O check_hamiltonian_bridges.py --max-blue-n 8 > optimized.json
cmp normal.json optimized.json
python3 verify_hostile.py > hostile.normal.json
python3 -O verify_hostile.py > hostile.optimized.json
cmp hostile.normal.json hostile.optimized.json
sha256sum -c SHA256SUMS
```

All assertions are explicit exceptions and survive `-O`. The main checker
tests all 4,694 blue presentations at orders 3 through 8 and all 64 order-five
tilings. The independent hostile checker directly enumerates vertex words
and the full labelled tournament universes through order four.

`SOURCE_PINS.json` gives exact repository blob and SHA-256 pins; `sources/`
contains the inherited inputs. No current-checkout or network dependency is
needed to replay the finite checks. The package needs no source checkout to execute.

This packet is research staging, not a publication claim. The full retained
black-line twin theorem and centralizer formula are inherited and credited;
there is no claim of a new proof of Redei parity or OCF, no global forbidden
Hamiltonian count claim, and no arbitrary-order connectivity claim.

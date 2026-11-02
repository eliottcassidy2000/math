# Exact binary content of a reversal-closed path response

`PROOF.md` proves that the coefficient gcd equals the gcd of the binomial
moments, and gives the reversal/free-orientation lower bound
`2^(N-Umax+1)`. The exact contents are `2^276` at n25 and `2^296` at n26;
only the latter attains the reversal bound. A specified guard leaf gives
a concrete source/sink division of the reversal orbits.

This is a separate bounded refinement of the square-path packet, not a
new census or an external-priority claim. It also distinguishes maximum
response from maximum union support, and records the exact falling mean
and positive-response probability when passing from n25 to n26.

## Reproduction

Copy the selected packet including `sources/`, then run from any directory:

```sh
python3 /path/to/square_paths_content/verify.py > /tmp/content-normal.json
python3 -O /path/to/square_paths_content/verify.py > /tmp/content-optimized.json
cmp /tmp/content-normal.json /tmp/content-optimized.json
cmp /tmp/content-normal.json /path/to/square_paths_content/normal.json
```

Python standard library only. Checks remain active under `-O`. The checker
exhausts all 74 undirected guards at orders 2 through 4, explicitly verifies
902 nonempty compatible subsets and 112 leaf-half cases, and compares
with literal full tournaments in this small universe. Negative controls
reject n1 and a non-reversal-closed bank. Existing n25/n26 union profiles
are recompiled into moments, coefficients, gcds and exact probabilities;
no larger census is performed.

The two source snapshots are authenticated in `SOURCE_PINS.json` and come
from the terminology-corrected square-path core. They are local portable
inputs, not filesystem or network dependencies outside this packet.
`PORTABILITY.json` records isolated normal/-O replay. `PUBLIC_FILES.json`
and `SHA256SUMS` define and authenticate the release payload.

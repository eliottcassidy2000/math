# Local universality of the strong one-apex H-spectrum

**PROVED structural theorem; FINITE-EXACT implementation checks.**

For every modulus, strongly connected one-apex tournaments realize every
parity-compatible Hamiltonian-path-count residue, simultaneously across CRT
factors, and each residue contains arbitrarily large realized values. Thus
one-apex, strong, and unrestricted spectra all have the same profinite closure:

    1 + 2 Zhat.

This does **not** establish exact integer completeness. The one-apex family
misses 13 globally despite its local universality; the unrestricted spectrum
misses the inherited 7 and 21. No congruence-only finite observer on H can
explain these gaps. An upper bound on H or other noncongruential/graph sidecars
can still distinguish them.

## Files and reproduction

* `PROOF.md`: actual-word source, odd modular commutator, synchronized CRT,
  positive inverse compilation, polynomial numeric word bound, and limitations.
* `check.py`: compact actual-word constructor, separately implemented literal
  adjacency/subset-DP consumer, residue universe, and hostile controls.
* `normal.json`, `optimized.json`: matching finite receipts.

Run from any directory, without packages or network:

```sh
python3 -B /path/to/local_closure/check.py
python3 -B -O /path/to/local_closure/check.py
```

The receipts record 12,880 c-residue targets for every modulus 1..160,
70 selected larger mixed-modulus targets, 1,759 literal expanded words,
76 literal constructor adjacency counts, every binary word through length 9
(1,023 graphs), and 43,679 odd-modular affine states. The graph check is
algorithmically separate from the affine constructor; these two producer
runs alone are not described as an independent-agent audit.

The graph constructor is also callable:

```python
from check import construct, action, length, expand, adjacency
word, metadata = construct(M, z)   # guarantees c=(H-1)/2 == z modulo M
state = action(word, (0, 0), M)
assert sum(state) % M == z % M
# Expand only if its explicit length is manageable:
if length(word) <= 1000:
    bits = expand(word)
    outgoing_masks = adjacency(bits)
```

The compressed word uses only `A=0`, `B=10`, sequence, and nonnegative
repetition. `A` and `B` are actual appended apex-neighborhood bits, not
formal state updates. For larger output, this compact representation and
the documented adjacency rule specify the graph without allocating its full
adjacency matrix. The largest word in the finite check has 22,846,022 bits;
no claim is made that every large graph was expanded or subset-counted.

For H==h modulo m, solve 2z==h-1 modulo m and call `construct(m,z)`.
This equation has a solution exactly for parity-compatible residues.
The theorem uses no fixed-order graph enumeration and no all-target positive
inverse assumption. See the proof for the exact inherited source paths.

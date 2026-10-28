# Scope, validity checks, and reproduction

Source checkout: b38ec1a1eaf7, 2026-10-28, in an isolated session worktree.
No theorem IDs are reserved. The local S/P identifiers always refer to
THEORY.md. Recovered literature is credited; no general canonization or
self-converse representation discovery is claimed.

## Proof gates

S2 checks both necessary and sufficient blue conditions, empty rows/columns,
prefixes of size<=2, and the dependency of compressed updates on retained
entries only. S3 retains **ordered** end ports and uses native path complement,
not global converse. Fresh reset is class-constant; named/conditioned path
selection is outside the theorem. S4 supplies literal maps, a direct
geometric explanation, and an infinite strong family. Its minimal-order
claim depends on a replayed exhaustive result through6, not a search failing
to find smaller examples.

P1 keeps the common labeling and fixes the extra root. P2 is a class bijection,
not an equivalence preserving full automorphism groups. P3 distinguishes
independent pairs from compatible extensions and scalar counts from generated
classes. All claims of faster runtime are avoided: current runs overlap other
research computations and are not a matched benchmark. The reduction in stored
root records is exact; parent dictionaries and lookup storage remain additional
costs. Lazy recovery also costs time when requested.

## Finite universes and independent controls

* `check_theory.py`: every7-presentation grouped by exact endpoint-rooted types;
  both explicit witness isomorphisms;21845 contexts through seven insertions;
  all mixed-end/nonblue seeds through6 with contexts through three insertions;
  strong family orders7..30; full81/27 automorphism and attachment controls.
* `boundary_profile.py`: all presentations at1..7 and all three one-step
  updates; replay of inherited horizon7 partitions and all five horizon8
  separators at6. No strong-only or necklace filter is used.
* `pair_square.py`: all2^(n²) bit squares for0<=n<=3, all S_n permutations,
  explicit inverses, and independent complete ordinary self-converse banks
  through7. The channel prefix through100 uses inherited exact arithmetic.
* `lazy_root_switch.py`: all parent attachment orbits through9. Full-child
  canonicalization and known counts are audit-only. Every stored map is checked
  against literal arcs. Every lazy row through8 and a deterministic1026-row
  sample at9 is recomputed. The independent saved bank supplies all191536
  reference9-classes. No order10 enumeration is claimed.
* `count_extension.py`: exact integral degree-bucket and terminal-contraction
  kernels agree on0..350, including the saved300 prefix. Their common proof
  ancestry is disclosed. Nonempty strong coefficients are independently
  checked by ordered-component convolution; empty strong coefficient is0.
* `cross_domain.py`: exact symbolic determinant -2 and three rational images
  from existing THM-1300; blue syndrome/zero-test distinction after J_L.

Discarded probe: the two7-presentations have equal empty-context syndromes,
both5216, so they do not witness different syndromes until a left singleton
is added. The repaired values1820815/1982543 are distinct and nonzero. This
does not affect the closed-carrier proof or the all-size collision.

## Commands

From the repository root, with Python3 and SymPy:

```text
python -B 05-knowledge/results/singleton_pair_reconstruction_20261028/verify.py
python -B -O 05-knowledge/results/singleton_pair_reconstruction_20261028/verify.py
python -B 05-knowledge/results/singleton_pair_reconstruction_20261028/lazy_root_switch.py --max-order 9
python -B 05-knowledge/results/singleton_pair_reconstruction_20261028/count_extension.py --bound 350
python -B 05-knowledge/results/singleton_pair_reconstruction_20261028/count_extension.py --bound 350 --reference
```

The default verifier replays the small exhaustive mathematics, validates all
saved350 coefficients and their full equality, and checks source/receipt hashes.
It does not silently rerun the larger enumeration or timing jobs. Those are
the explicit final three commands. Fresh timings change the raw receipt hash;
compare semantic fields, not wall-clock times, when rerunning large jobs.
Hashes normalize LF for text and retain raw bytes for gzip. No assertion can
be disabled with `python -O`; checks raise explicit exceptions.

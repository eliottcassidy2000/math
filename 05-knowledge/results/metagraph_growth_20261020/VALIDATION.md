# Validation and provenance

20 October 2026. **PASS within the stated universes.** The all-size
statements rely on the proofs in growth.md and their named dependencies,
not on extrapolation from these finite computations.

## Exact replays

Run from the repository root:

```text
python -X utf8 -B 05-knowledge/results/metagraph_growth_20261020/audit.py
```

The auditor reruns enumeration.py and growth_check.py normally and with
Python optimization, compares each pair byte-for-byte after explicit
CRLF-to-LF normalization, and compares each with its retained literal
stdout. It does not regenerate or overwrite the lane outputs.

| Lane | Declared universe and independent controls | Result |
|---|---|---|
| enumeration | All classes and path tilings through n=7; all vertex-attachment patterns through parent n=6; full-permutation canonicalization through n=6; labelled-orbit mass, HP/Aut fibre, converse and line-incidence checks | PASS |
| growth_check | All 1099 labelled tournaments through n=5 versus all-subset maximal-module check; all path cubes through n=6; 152 inflation controls at k=2,3; 76 full-permutation class controls; 16 run-multinomial versus independent HP-DP checks; 33856 two-end colour controls; literal clone-allocation hostile | PASS |
| root audit | Independent all-labelled/full-permutation class counts through n=4; minimal two-rib triangle; isomorphic gap graphs with different suspension classes | PASS |

No load-bearing check uses a Python `assert`, so optimization does not
disable validation. The retained [audit.out](audit.out) reports:

- enumeration stdout LF SHA256:
  `a8df122d2295e8a009775c21a79dc9022614147a3236fb465f6400c6068aeb82`
- growth_check stdout LF SHA256:
  `a07ff82f53ae84ef804895b2f18018241550afb3e44480dcb22d66add3238860`

The packet manifest records normalized hashes of the final source and
notes as well. It excludes itself and generated Python cache files.

## Independent proof review

The inheritance reviewer and enumeration reviewer separately accepted:
the unique uniform-root criterion, canonical chain-gcd divisor
stratification, automorphism and converse preservation, and the exact
complement-defect support/cardinality. The root also checked these
arguments directly against THM-4463 and the coordinate definitions.

One pre-publication wording error was repaired: clone labels need to
increase inside each contiguous block visit, not globally across
separated visits. The repaired multinomial weight is
(k!)^n divided by the product of run-length factorials. The added ABBCCA
witness explicitly demonstrates the missing allocation coordinate.

The order-three growth hostile is limited to an unmarked source
metagraph; separate H and Aut distinguish that pair. The order-seven
hostile refutes the specified ten local coordinates, not arbitrary
information from the complete marked graph. The divisor law isolates
primitive counts but does not determine them from smaller orders.

## Inheritance and concurrent work

All-size additive pullback, phase torsor and colour transfer come from
THM-796. Canonical transitive modules come from THM-4463. The fixed-path
fibre and class-BFS methods come from THM-781/1430. Blowup construction
itself already occurs in THM-3235. These are credited as prior work.

The primary external count comparison is
[OEIS A000568](https://oeis.org/A000568), read 20 October 2026. Its Davis
formula is not used by the independent enumeration algorithm.

The session began from clean origin/main commit `013f47465485`, after
fetching and rebasing a newly managed isolated worktree. The shared
checkout had unrelated changes and was not modified. The targeted
processor peek could not run because the worktree had no `.machine-id`;
no registration was invented. Contemporary THM-4729/4730/4731 entries
were reserved empty stubs and were not dependencies.

The repaired historical explainer is outside the packet manifest.
Its minimal mathematical correction is reproduced by audit.py; exact
lineage is recorded in 01-canon/MISTAKES.md.

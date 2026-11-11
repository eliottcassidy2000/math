# Tournament aspects as sequences: metagraph, extensions, and seed geometry

**FINITE-EXACT:** 82 precisely defined sequences through order eight, using all
7,412 ordinary classes of orders 1–8. Selected sequences extend through all
191,536 ordinary classes at order nine. **PROVED IN NOTE:** a refined
11-presentation generating function and the sharp order of the preceding 33
unit-offset roots, conditional only on the explicitly inherited complete census
and strong-fibre certificates. General reconstruction and Chan–Pak remain OPEN.
Same-agent audit using independent algorithms; no external novelty claim.

The main deliverable is a reusable registry, not a list of unexplained numbers:
[sequence definitions and terms](sequences.tsv), [machine-readable registry](registry.json),
[per-class profiles](profiles.tsv), [full exact records](records.jsonl.gz), and
[query tool](query.py). Every sequence declares its universe and convention.
Blank order-nine cells mean not computed here; they never mean zero.

## What 33, 11, and 7 actually measure

The inherited path-presentation theorem is
`tc(T)=H(T)/|Aut(T)|`; see THM-781,
[`THM-781-hamiltonian-path-inverse-metagraph-fibre.md`](../../../01-canon/theorems/THM-781-hamiltonian-path-inverse-metagraph-fibre.md).
At order six, exactly five classes have H=33. Three are rigid, with tc=33;
two have automorphism order three, with tc=11. This is an actual quotient
mechanism behind 33=3*11. The line table also contains a black connection
between an H33 class and the H11 class at order six, with literal multiplicity
two. [Number roles, witnesses, and all-order formulas](NUMBERS.md) retain these
different meanings and recover the earlier H33 multiplicity work.

Seven is impossible as H or as a complete path-orbit count, but possible as
an endpoint-restricted path count, automorphism order, or odd seed coefficient.
Those exclusions and boundary examples are inherited from THM-4812 and
THM-4987/4988, not rediscovered by this bounded census.

## Aspects worth keeping

Each of twelve class signatures generates its number of values, collision-pair
count, and largest fibre. Other sequences track totals, maxima, and failures:

- named and rooted-orbit extension distributions;
- deletion H decks and ordered endpoint roles;
- endpoint-score flux, checked against the actual metagraph current;
- unit selections, their symmetry orbits, and the coefficients they produce;
- incompatible pairs of child targets inside a fixed parent-H fibre;
- line multiplicity, loops, equal-H lines, and failures of H/tc descent;
- how much the complete neighbour-H distribution distinguishes merged nodes.

Selected rows (orders 3 through 9):

| Aspect | 3 | 4 | 5 | 6 | 7 | 8 | 9 |
|---|---:|---:|---:|---:|---:|---:|---:|
| Ordinary classes |2|4|12|56|456|6880|191536|
| Converse-merged classes |2|3|10|34|272|3528|97144|
| Distinct H values |2|3|7|19|77|320|1520|
| Distinct one-vertex extension supports |2|3|10|34|272|3528|97144|
| Classes admitting positive unit selections |2|4|9|26|110|1053|17423|
| Unit-selection orbits |4|9|22|59|194|1300|18849|
| Distinct odd unit coefficients |2|3|7|24|89|636|5529|
| H-fibres with incompatible target pairs |0|0|3|10|62|309|—|

**Useful finite separator.** The unordered set of one-vertex extension H values
distinguishes every converse-merged class through order nine. Multiplicity is
unnecessary in this range. Equality is tested on complete integer vectors, not
hash digests. This is a bounded result, not an all-order reconstruction theorem.
The older order-12 81/27 decoration hostile and the order-13 uniform-incidence
hostile are also separated, as recorded in [hostiles.json](hostiles.json).

**Source identity is essential.** At order five, the three H3 sources collectively
offer targets 17 and 29, but no one source offers both. There are already 13 such
incompatible pairs over all H-fibres at that order, and 32,226,456 at order eight.
These counts refer to parents of exactly the stated order, not an all-height
nonexistence theorem. [Exact source/attachment witnesses](source_obstructions.json)
make the first failed gluing implication reviewable.

## A sharp continuation of the previous seed work

The two independent complete order-nine seed scans agree: among odd c<=4097,
exactly the preceding packet's 33 coefficients have no unit-offset base.
Adding a source outside the selected set preserves H, W, and U, so absence at
base order nine implies absence at every smaller base order. The preceding
[33 explicit bases of order ten](../odd_c_seed_symmetry_20261110/README.md)
therefore have **minimum root order eleven in the THM-5126 selected-set tail
construction**. This is not a lower bound for arbitrary tournaments realizing
one point of the ray. [Proof and exact boundary](PROOFS.md).

The complete scan additionally saves one explicit witness for each of its
5,529 odd coefficients in [seed_witnesses9.tsv](seed_witnesses9.tsv). This set is
not an initial interval; universal odd-c construction remains OPEN.

## Efficiency and validation

One pair of subset path tables supplies all endpoint-cover coordinates. All
2^n attachments are evaluated from one quadratic response in O(n^2 2^n)
arithmetic operations per parent, instead of rebuilding each child path table.
The order-nine screen processes 97,144 converse representatives and evaluates
49,737,728 attachments, representing all 98,066,432 named parent attachments.
The measured screen time is in [screen9-seconds.txt](screen9-seconds.txt); it is
a local timing, not a cross-machine performance claim. The seed audit separately
processes all 191,536 ordinary sources, including both converse directions.

Reproduce with `python3 -B reproduce.py` (Python standard library and g++).
Compilation and intermediate files use a temporary directory. The independent
audits use degree-permutation canonicalization of all 33,868 frames through
seven, literal child/tail counts, a second all-order-nine seed algorithm, and
relabeling/converse controls. Input files are pinned by SHA256. The order-eight
and order-nine class universes are inherited from their audited complete banks;
they are not claimed as newly generated here. No Lean proof or independent-agent
review is claimed.

`SHA256SUMS` fingerprints every payload file, normalizing text CRLF to LF and
leaving gzip bytes unchanged. The reproduction command refreshes it after the
new local timing is recorded. Source-bank text fingerprints use the same
newline convention; the compressed order-nine source fingerprint is literal.

Next questions are recorded in [the board](BOARD.md). Incoming seven-coordinate
compiler work is integrated there as a different, bounded-context response
observer; it does not establish all-order completeness of the extension-support
signature measured here.

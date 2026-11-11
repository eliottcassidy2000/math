# Reproduction and evidence scope

The mathematical claims are in the proof notes. Receipts distinguish
all-size proofs from bounded experiments and finite bases made complete by
proved cutoffs. No inherited theorem or outside source was broadly audited.
Python uses only the standard library; the two shell engines use C++17.

## Quick independent interface controls

From the repository root:

    python3 -B 05-knowledge/results/session260_quasiperfect_chanpak_20261111/review/third_wave_check.py
    python3 -B 05-knowledge/results/session260_quasiperfect_chanpak_20261111/review/bipartite_defect_tiny.py
    python3 -B 05-knowledge/results/session260_quasiperfect_chanpak_20261111/arithmetic/norm_check.py
    python3 -B 05-knowledge/results/session260_quasiperfect_chanpak_20261111/coverage/marked_module.py

The close-out replay runs these four under normal and optimized Python and
compares their JSON output with the associated receipts. See REPLAY.json.
This is a scoped interface replay, not a new full census.

## Author engines

All paths below are relative to this packet. Run each Python file with
`python3 -B <path>` from the repository root or its absolute path. Most
engines also rewrite their adjacent receipt deterministically.

- arithmetic/check.py → arithmetic/receipt.json: recovered source cases,
  cyclic products, explicit successor graphs, and fan controls.
- arithmetic/norm_check.py → arithmetic/norm_receipt.json:121 substitution
  pairs, all1099 labelled tournaments through order5 for the justified
  norm-switch base, and independent51/99 permutation controls.
- bridges/check.py → bridges/receipt.json: physical congruence controls and
  alternating-root tests. Its finite m range is not the universal proof.
- bridges/crossing_check.py → bridges/crossing_receipt.json:196 literal
  crossing-edge checks; all-size proof is in CROSSING_FLIP.md.
- local_check.py → local_receipt.json:154 single-edge controls.
- constructors/module_ladders.py and marked_position.py: exact finite
  skeleton coefficients plus physical module controls for all-length laws.
- constructors/relocation_atlas.py:1374 graph presentations in the declared
  two-core box,627 seed presentations,212 odd slopes, and literal heads.
- coverage/marked_module.py:108 marked-module controls. The later two-root
  coarsening/gap birth is deduced from the accepted formulas, not a census.
- bridges/domination_check.py:511 chain words and592 arbitrary bipartite
  scouts. The former OPEN question is superseded by BIPARTITE_DEFECT.md.
- bridges/bipartite_defect_check.py:5405 exhaustive roots through order5,
  plus600 order6 roots; every cut fibre is checked.
- constructors/two_edge_repair.py: all apex cuts for a=2,...,8 and the
  finite symbolic rank-template coefficient identity.
- constructors/bipartite_boundary.py: complete rooted-tree state compiler
  through order10, component products below225 using the spectrum receipt,
  and literal leaf/twin/hub controls.
- bridges/defect_quantitative_check.py:202013 exhaustive rooted cases
  through order6,1098 newborn cuts on75 small bases,6800 seeded roots at
  orders7–10. The quantitative inequalities remain CONJECTURAL despite
  zero recorded counterexamples. `--max-order 5` shortens the exhaustive
  portion; the command then has a different finite scope.

## C++ shell and spectrum engines

From the repository root, with P denoting this packet path:

    P=05-knowledge/results/session260_quasiperfect_chanpak_20261111
    g++ -O3 -std=c++17 "$P/constructors/two_edge_repair.cpp" -o /tmp/two_edge_repair260
    /tmp/two_edge_repair260 2 6 2
    /tmp/two_edge_repair260 2 6 3
    g++ -O3 -std=c++17 "$P/constructors/bipartite_spectrum.cpp" -o /tmp/bipartite_spectrum260
    /tmp/bipartite_spectrum260 9

The first shell is complete for two flipped edges at a=2,...,6. The second
is the explicitly selected three-edge union, not all triples. Fourth-edge
extensions are exactly the near-miss continuations described in the note.
The bipartite engine enumerates all matrices on every declared shore split
through order9, retaining duplicates. The rooted-tree cutoff, not the
finite enumeration alone, makes the below225 result valid at all orders.

## Reviews and publication

Focused proof/interface reports are in review/REVIEW.md, THIRD_WAVE.md,
BIPARTITE_DEFECT.md, BIPARTITE_COVERAGE.md, and FINAL_INTERFACES.md. These
are internal mathematical reviews, not external refereeing or priority claims.

SOURCE_PINS.json identifies inherited and concurrent source commits.
SHA256SUMS pins every other packet artifact except itself. Use
`sha256sum -c SHA256SUMS` from this packet directory. A manifest is an
integrity record, not a substitute for a proof.

Published checkpoints before this final synthesis:

- 6f495220a636cf3fa8413db82a095f0e763d3a9a: arithmetic and Carmichael bridge.
- 4315a74fe0b2a4305741429958804fafb58a130d: seed grammar and all-edge rigidity.
- 46b06e8a864c048cc62bc640447b9dfc23b79ba6: domination theorem and norm repairs.

Each was fetched back and its complete tree or explicit new blobs matched
the local publication state. Concurrent main commits were preserved; no
force push was used. Final publication adds the remaining proved interfaces,
reviews, synthesis, and navigation. Its commit is the version containing
this file, rather than an impossible self-referential stored commit hash.

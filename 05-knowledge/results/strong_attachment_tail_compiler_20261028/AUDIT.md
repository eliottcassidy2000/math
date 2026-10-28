# Audit, sources and reproduction

2026-10-28. Same-agent mathematical review with separate computational paths.
No independent human/agent certification or formal proof claim.

## Commands

Run from the repository root, using Python3.13 or another standard-library
Python with arbitrary precision integers:

```text
python -B 05-knowledge/results/strong_attachment_tail_compiler_20261028/verify.py
python -B -O 05-knowledge/results/strong_attachment_tail_compiler_20261028/verify.py
python -B 05-knowledge/results/strong_attachment_tail_compiler_20261028/verify.py --replay-large
python -B 05-knowledge/results/strong_attachment_tail_compiler_20261028/counting.py --bound 300
python -B 05-knowledge/results/strong_attachment_tail_compiler_20261028/root_switch_skip.py --max-order 9 --compare-small --output root-switch-replay.json
python -B 05-knowledge/results/strong_attachment_tail_compiler_20261028/benchmark.py --output timing-replay.json
python -B 05-knowledge/results/strong_attachment_tail_compiler_20261028/benchmark_roots.py root-timing-replay.json
```

The default verifier rebuilds the entire native through-seven experiment,
checks the residual order-eight witnesses, reruns small scalar kernels and
checks the complete large receipts. `--replay-large` also reruns the300 kernel.
Constructive order-nine regeneration is the separate root-switch command
or `verify.py --replay-roots`.
Output paths in the last two commands are caller-chosen new receipts; ordinary
verification does not overwrite saved evidence. Mathematical guards use
RuntimeError, not removable assertions. The imported June control has its own
assertion; our explicit equality checks remain active under `-O`.

## Universes and controls

| computation | exhaustive universe / observations | positive and hostile controls | separate path |
|---|---|---|---|
| native classes | all X_n for1<=n<=7, legal right-relocation generators | counts1,1,2,4,12,56,456; blue flag retained | all vertex permutations through5 |
| finite carrier | ordinary class and blue; tau, fresh reset, left/right singleton, max7 | full through5; five unresolved size-six pairs | order-eight exact isomorphism search separates every residual pair |
| group frames | six specified strong tournaments; all group elements; all candidate bases; all distinct ordered2/3-ports | free frames and root stars; P7 joint/pair collision;81/27 transporter | actual adjacency-preserving maps, explicit stabilizer checks |
| markers | both specified P7 triples, all six marker orders | all H pairs differ | independent Hamilton DFS for the4281/4315 pair versus subset DP |
| fields | all shifts on a chosen normal basis in F64/F256 | reversed-pair obstructions; actual order-three S9 action | inherited odd-cycle proof and direct pair-orbit traversal |
| scalar count | every n0..300, exact signed/integer weights | inherited whole-vector baseline to300; prefix extension250 | saved Fraction DP through200, direct odd partitions through40, OEIS through76 |
| pair-free core | every n0..50 | signed excess and nonnegative final counts | inherited Fraction receipt through50; independent polynomial substitution at every horizon0..20 |
| strong count | every n0..300, A051337 empty-object convention1 | ordered strong-component convolution | independent OEIS table through50 |
| timings | three alternating whole-kernel samples at100/200/250 and core50 | every vector compared, two-stage slower control | no historical cross-machine speed ratio |
| root-switch skip | all parent attachment orbits through the receipt's max order | disjoint complete row supports; post-generation count | same representatives as inherited generator through8 and full saved bank at9 |

Matched constructive timings at8 use the same producer with skipping on/off,
including parent setup, and omit the separate Hamilton/prime audits in both
modes. Three alternating samples give7.4959s/1.4112s, or5.31x. The47.72s
order-nine run is unpaired and is not used for a comparative timing claim.

The closure and frame universes are different. A base is an ordered tuple in
one fixed atom; it is not a representation of every Hamilton marking.
The marker operation is not asserted to be modular substitution. Strongness
is not assumed closed under complement: the closure computation keeps all
tournaments, including nonstrong ones.

`sources.json` pins imported code, mathematical dependencies and external
reference tables. `manifest.json` seals this packet except its own manifest.
Text SHA256 inputs normalize CRLF to LF, so Git's Windows line-ending policy
does not invalidate the seals; compressed binary banks use their raw bytes.
`validation.json` records exact checks, not elapsed time; `benchmark.json`
contains machine-sensitive times. The probe receipts are exploratory runs,
not additional independent benchmarks.

The [concurrent integration](INCOMING.md) adds a third exact ordinary kernel
through300 and independent selected-root class bijections through7. Its
single-run32.53s is not a paired speed ratio. The incoming algorithms keep
their own credit and distinguish representative selection from full rows.

## Literature checked

- [OEIS A000568](https://oeis.org/A000568): sequence convention, Davis formula
  and the linked Keith Briggs table through76. The table is inherited and
  pinned, not regenerated from our answers.
- [OEIS A051337](https://oeis.org/A051337): the prepended constant1 and
  2-1/A(x) convention. Andrew Howroyd's [b-file](https://oeis.org/A051337/b051337.txt)
  through50 was retrieved2026-10-28 using the web tool, with line breaks
  normalized for storage; direct urllib access returned403.
- Royle, Praeger, Glasby, Freedman and Devillers,
  [Tournaments and Even Graphs are Equinumerous](https://arxiv.org/abs/2204.01947):
  the graph property is automorphism-sign-even, not even vertex degrees.
  No operation-compatible bijection is inferred from equinumerosity.

The repository [CORE-PAPERS](../../reference/CORE-PAPERS.md), current guardrails
and research protocol were read before literature claims. The historical
source routes and correction precedence are listed in [HISTORY](HISTORY.md).

## Stopping boundary

The largest scalar run has1,683,206 live prefix states at its peak. The
terminal optimization leaves that bottleneck intact. A further size increase
alone would not test a new mechanism; the next useful experiment is adaptive
middle-stage contraction or a matched GMP comparison. Constructive generation
at10 would require processing tens of millions of pointed states, and has not
been run. All-order minimality for the literal singleton menu remains OPEN.

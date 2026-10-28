# Reproduction and validity audit

**Review:** same-agent mathematical review with alternate computational paths.
No independent-agent, external referee, or new Lean-build claim.

## Commands

From the repository root, with Python 3 and a C++17 compiler `g++`:

```sh
python3 -B 05-knowledge/results/consecutive_flags_order10_20261028/verify.py
python3 -B -O 05-knowledge/results/consecutive_flags_order10_20261028/verify.py
python3 -B 05-knowledge/results/consecutive_flags_order10_20261028/verify.py --full
```

The light verifier checks source pins, compressed and raw checksums, every
accepted address, the saved exact census agreements, every finite closed
classifier, and the independent literal-permutation/Hamilton-reset audit.
`--full` compiles and reruns the generator, full-bank alternate canonicalizer,
pointed-state producer, and even/odd double checker. It compares every exported
bank byte and every deterministic receipt field, ignoring only timings.
Full builds and raw data use a temporary directory and are removed afterwards.
The light verifier rewrites only deterministic `structure_audit.json`, checked
against its saved value. Both the complete rebuild and light replay passed.

`publish.py` is the separate republishing command: it builds the census and
audit, updates timing receipts and recompresses the bank. Use `verify.py` to
validate unchanged saved artifacts. The original timing run used GCC 15.2,
`-O3 -std=c++17 -march=native`, on Windows with Python 3.13. The portable full
replay omits `-march=native`; timing is not an assertion. Its generation took
13.51 seconds, with identical deterministic results and bank bytes.

## Exact universes and controls

| Experiment | Universe | Alternate path or decisive control |
|---|---|---|
| Ordinary bank | all parent-mask orbits, orders 1..10 | full-child alternate canonicalization of all 9,733,056 outputs; unique codes and labelled mass 2^45 |
| Closed rooted path quotient | all N*2^((N-1)(N-2)/2) states for N=2..7 | literal all-permutation pointed classes and actual uniform Hamilton resets through N=5 |
| Ordered deletion forest | all 6,504 path/root states for N=2..6 | native partner restriction versus separately constructed forest complement, including every root-mask edge |
| Self-converse doubles | every input for n=1..4, both parities | source-permutation orbits versus target canonical classes and full centralizers; 5,393,161 active checks |
| Signed B/D scores | all 599 inputs for n=1..3, both parities | direct outgoing degrees versus integer doubled root sums |
| 81/27 pair | two explicit twelve-vertex rows | same 67 labelled cyclic triples; direct full groups 81/27; Hamilton counts 462,915 each; named identity path true/false |

The last control imports a pinned existing direct isomorphism/backtracking
implementation, distinct from the new C++ canonicalizer. Production C++ also
checks C3, regular T5 and P7 groups 3,5,21 and an isomorphic 81/27 control made
by reversing one four-vertex module. The two implementations need not use the
same vertex labels. `structure_audit.json` records 43,962 active checks.

Explicit `need` exceptions replace Python assertions, so optimized Python
retains the gates. `source_pins.json` hashes twelve dependency files after
CRLF-to-LF normalization; it does not freeze unrelated repository history.
The bank manifest uses raw-byte SHA256 for plain and compressed artifacts.

## Mathematical gates

- **Types:** ordinary class, pointed class, fixed-path/root state, aligned pair,
  and marked/unmarked self-converse double are distinct. Matrix conjugation is
  not the induced tournament action.
- **Iff directions:** equal exact deleted-parent/mask keys extend to a root
  isomorphism. The invariant cell is only a preliminary filter.
- **Symmetries:** full parent actions and all attaining transports survive;
  rooted stabilizer times root-orbit size recovers the child group. A double's
  marked group is its anti-involution centralizer, not automatically its full group.
- **Quotient loss:** the minimum is for strongly lumpable random reset plus
  native complement with the declared observation. It is not a minimum for
  arbitrary attachments, named guards, or weak trace equivalence.
- **Actual consequence:** uniqueness plus exact orbit mass proves full bank
  completeness. Matching known counts or cyclic triples alone would not.
  Scalar counts, representative banks, and full metagraph outputs are separate.
- **Literature:** Eplett and Coxeter roots are credited classical inputs. The
  ETH document is a source-reported claim; conditional deductions and
  unconditional output/information bounds remain separate.

Order eleven is not attempted. A full order-ten native metagraph matrix is not
exported. Additional exact transports and transition rows must be computed when
requested. The weak observer's equality with the pointed-pair partition for
every N>=5 remains open beyond N=7. Session89's independent-pair/stream formulas
are integrated, but no order-ten independent-pair bank or arbitrary-stream
runtime is claimed. No major existence conjecture is closed by the signed
incidence connection or the conditional SAT discussion.

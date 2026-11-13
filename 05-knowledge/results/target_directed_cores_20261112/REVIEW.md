# Review, exact universes, and reproduction

2026-11-12/13. Proofs are in notes, not Lean formalizations. No numbered
canon namespace was reserved. The arithmetic and grammar inherited from
older work are credited in TARGET_INVERSE and INCOMING.

## Independent paths and scope

* Root and the ports reviewer audited the exact binary inverse, intrinsic
  P2 bound, weak-source reduction, and every pruning inequality. The
  history lane independently supplied the growth argument. A separate C++
  program includes neither producer source nor its generated bank: it uses
  depth-first core addresses and direct enumeration of every bounded outer
  word. It matches all core layers at both excluded targets. Unpruned
  Python grammar/word-automaton controls at881 and1997 additionally test
  pruning and binary conventions. See COMPLETE_INVERSE_AUDIT in ports.
* The43 positive405-ray sources retain full chronological addresses and
  literal outgoing rows. The six-state compiler supplies the full counts;
  three small sources also have literal subset-DP counts. Large adjacency
  certificates are not described as large subset censuses.
* The ports reviewer independently checks all seven moments of the
  two-entry defect-two core, including its joint correction. All16 short
  words and their NR children have independent literal four-state counts.
  Six zero-prefix extensions preserve old arcs and both ports. The exact
  all-word formula proves the27-vertex source and infinite family.
* The selective bridge is independently rebuilt from literal graphs,
  its six-state source/target arithmetic, and the retained map. It has75
  changed nonroot pairs, including8 old predecessor/module pairs. Every
  old vertex, root comparison and internal TT17 arc is retained.
* The arbitrary-module881 split is proved by reduced-word coefficients
  and the maximum three-run boundary, with small literal module controls.
  Its physical induced9-vertex overlap is independently audited. Bounded
  paired/triple-scale and pure-ternary menus retain explicit stopping scopes.
* The history reviewer checks the all-module462 proof, port change,
  resultant gate, and exponent5039 exclusions. Root reads both proofs.
  The general ordinal-square factory and explicit all-depth H35 family are
  independently reviewed by both other lanes. Ordinary/optimized controls
  agree for the root selector; they do not replace the all-parameter proof.
* The periodic-extension theorem and positive rewrite debt are independently
  checked by the ports reviewer; root reads both. Their finite checker
  covers1024 legal words,285 rooted sources and1140 literal children.
  Negative-Collatz conjugacy is independently derived by root and tested
  over49,152 arithmetic inputs plus literal23/45 controls. Its universal
  graph-lift obstruction uses inherited all-order7/21 exclusions.
* The ports reviewer checks SINGLETON_CONTINUATION's two physical extension
  formulas, complete ten-case divisor gate and positive mixed-block drift.
  The exclusion is only for fixed boundary words around one repeated block;
  it leaves coupled variable blocks and changing cores open.

## Rebuild the target producer and positive certificates

Run from the repository root. Generated core banks are intentionally
ignored; full deterministic source is retained.

```sh
python3 -B 05-knowledge/results/target_directed_cores_20261112/build_core_bank.py --order 12
python3 -B - <<'PY'
from pathlib import Path
from subprocess import run
from tempfile import TemporaryDirectory
p = Path('05-knowledge/results/target_directed_cores_20261112').resolve()
with TemporaryDirectory(prefix='.reproduce_', dir=p) as folder:
    inverse, complete = Path(folder)/'inverse', Path(folder)/'complete'
    for source, exe in [('target_inverse.cpp', inverse), ('complete_ready_inverse.cpp', complete)]:
        run(['c++', '-O3', '-std=c++17', str(p/source), '-o', str(exe)], check=True)
    run([str(inverse), str(p/'core_bank.txt'), '44', str(p/'target_inverse_receipt.json')], check=True)
    for target, output in [(53084231, 'complete_53084231.json'),
                           (159252697, 'constructors/ready_after_defect_two.json'),
                           (1997, 'complete_1997.json'),
                           (159252695, 'complete_159252695.json'),
                           (3887, 'complete_3887.json'), (4967, 'complete_4967.json'),
                           (8045, 'complete_8045.json')]:
        run([str(complete), str(target), str(p/output)], check=True)
PY
python3 -B 05-knowledge/results/target_directed_cores_20261112/certify_sources.py
```

The two negative targets use the complete intrinsic universe, not just
the small bank. Target53084231 retains5,640,026 cores, target159252697
retains16,075,330. The small bank has10,946 cores through12 and supplies
the finite positive head with arbitrary outer word length. The three
extra positive targets3887/4967/8045 belong to the constructor lane's
pure-ternary near misses, not to a new complete spectrum interval.

## Independent inverse and physical controls

```sh
python3 -B 05-knowledge/results/target_directed_cores_20261112/ports/audit_complete.py
python3 -B 05-knowledge/results/target_directed_cores_20261112/ports/audit_defect_two.py
python3 -B 05-knowledge/results/target_directed_cores_20261112/ports/audit_gap_bridge.py
python3 -B 05-knowledge/results/target_directed_cores_20261112/ports/module_lift.py
python3 -B 05-knowledge/results/target_directed_cores_20261112/ports/port_escape.py
python3 -B 05-knowledge/results/target_directed_cores_20261112/ports/intertwiner_probe.py
python3 -B 05-knowledge/results/target_directed_cores_20261112/constructors/selective_split.py
python3 -B 05-knowledge/results/target_directed_cores_20261112/constructors/defect_two_inverse.py
python3 -B 05-knowledge/results/target_directed_cores_20261112/constructors/singleton_continuation.py
python3 -B 05-knowledge/results/target_directed_cores_20261112/growing_family.py
python3 -B 05-knowledge/results/target_directed_cores_20261112/history/check.py
python3 -B 05-knowledge/results/target_directed_cores_20261112/history/negative_strata.py
```

The independent inverse automatically removes its temporary executable.
Its direct outer-address totals are74,845,805 and233,412,223, respectively.
The small881 positive and1997 hostile controls are included. The module
lift checks literal K1,TT2,TT3,C3 inputs and reduced run coefficients;
larger-module claims use the universal run proof. The root growing-family
selector checks624 parameter pairs and all31 family targets≤10^9, with
explicit parameter cutoffs. Its negative examples are only family
nonmembership, never Hamilton nonattainment.

Normal/optimized controls were run for the independent history and
constructor lanes and root growing-family selector. Heavy exhaustive
computations were not repeated merely for packaging after independent
agreement. Every saved claim retains either a proof or an explicit finite
universe. Packet links, whitespace and content hashes, plus the repository
documentation surface, are checked at close-out.

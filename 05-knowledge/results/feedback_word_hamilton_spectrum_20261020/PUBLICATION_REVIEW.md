# Publication review and reproduction

20 October 2026. The three frozen mathematical pieces were independently reviewed and accepted: the one-apex/bounded-feedback automaton, the fixed-run density theorem, and the exact intercyclic rays with their dyadic density bound. The [main audit](audit/AUDIT.md) and [intercyclic audit](audit/INTERCYCLIC_AUDIT.md) record the proof checks, primary-source qualifications, corrections, and hostile examples.

All twenty approved producer files are preserved exactly, together with eighteen audit/report/replay files. Superseded source snapshots are retained outside this publication, not presented as current bibliography or proof. The separate binary-product density candidate is not included.

## Portable checks

From this directory, run each command normally and again with `python3 -O`:

```sh
python3 verify.py
python3 verify_runs.py
python3 verify_intercyclic.py
python3 probe_affine.py
python3 probe_blocks.py
python3 audit/independent_checks.py
python3 audit/independent_intercyclic.py
```

All fourteen copied-layout runs passed and reproduced their archived output bytes. The producer's bounded matrix checks require NumPy (replayed with 2.3.5). The separate clean-room checks use only Python's standard library and do not import producer mathematics code. Import dependencies are preserved: verify_intercyclic uses verify; verify_runs and probe_blocks use probe_affine; both independent audit programs are colocated.

The three producer checksum files retain the original `session16-hamilton-spectrum/` prefix. For example, from this directory:

```sh
sed 's@session16-hamilton-spectrum/@@' SHA256SUMS | sha256sum -c -
sed 's@session16-hamilton-spectrum/@@' RUN_SHA256SUMS | sha256sum -c -
sed 's@session16-hamilton-spectrum/@@' INTERCYCLIC_SHA256SUMS | sha256sum -c -
```

Audit reports contain historical `source_snapshot/` replay commands. Those document the original audit layout; use the portable commands above in the published layout. The complete publication inventory also seals probe_affine.py, which the historical RUN_SHA256SUMS did not list.

## Scope retained

Raw Hamilton-path counts are not isomorphism-orbit counts. The bounded-feedback automaton is an exact all-length counting device, not an all-odd-value surjectivity theorem. Fixed-run polylogarithmic sparsity requires fixed feedback and fixed run complexity. The intercyclic exact rays and density deficit are separate statements, with the deficit claimed at dyadic cutoffs. H=29 is realized outside that restricted class, not a global spectrum omission.

The current complete odd tournament spectrum remains open. No unreviewed continuation or exploratory file is included in this packet.

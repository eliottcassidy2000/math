---
id: THM-4149
title: "Superseded minimum-boundary audit for the first-window universal odd-tail LRC(14) transfer"
status: >
  SUPERSEDED by THM-4148-first-window-width-universal-odd-tail-lrc14-transfer
  after concurrent integration; see MISTAKE-510. No distinct live theorem is
  asserted here. The independently audited alternative m=1,2 clock partition,
  the genuine (1,13) first-window trap, and its isolated 3/14 repair remain
  reproducible provenance. LRC(14) remains OPEN.
source: codex-lrc14-planar-jc-breakthrough-20260825
depends_on: []
related:
  - THM-4148-first-window-width-universal-odd-tail-lrc14-transfer
  - THM-4136-fixed-body-universal-odd-tail-lrc14-completion
script: 04-computation/lrc14_minimum_boundary_universal_odd_tail_transfer_thm4149.py
output: 05-knowledge/results/lrc14_minimum_boundary_universal_odd_tail_transfer_thm4149.out
independent_audit_script: 04-computation/lrc14_minimum_boundary_universal_odd_tail_transfer_thm4149_independent_audit.cpp
independent_audit_output: 05-knowledge/results/lrc14_minimum_boundary_universal_odd_tail_transfer_thm4149_independent_audit.out
script_sha256: cf146ae0b9a5fbbbdb6bf8597550c7923848cea1ae13ccba453b11c1b855e1c7
output_sha256: 84bca834c0ccbbb656c7d06ebc5f4fa30d4a9fb682d46464518dff400a0589d8
semantic_sha256: 90de52f9075568cead74de79315daf0494ad7f1f49e5ae7303eff0eae7881f6d
semantic_fnv64: 6784b7e0b01a759d
independent_audit_script_sha256: 283477210905e2ae4729ad978a9e67c7496581a413ac8a8fb6dc5b616a2df35b
independent_audit_output_sha256: b8457993e1ae0802aa6b9d65d36cdb6b0a501dfb6f9e9f8ad084ab3a1ea5f036
independent_semantic_fnv64: 6784b7e0b01a759d
hash_basis: raw LF bytes
---

# THM-4149 -- superseded minimum-boundary audit

**SUPERSEDED by
THM-4148-first-window-width-universal-odd-tail-lrc14-transfer; see
MISTAKE-510.**

Concurrent work incorporated the unrestricted `m=1,2` boundary directly
into THM-4148 before this local promotion was integrated. THM-4148 is
therefore the sole live theorem for

```text
W=13/(14 max H)-1/(14 min H)>=2/189
  =>  2H union {a,b} is 1/14-safe for all distinct positive odd a,b.
```

This file is retained only to route a genuinely different exact audit. Its
`m=2` residual ratios split into clock classes of sizes `58+8+2`. Its
`m=1` ratios split as `56+11+1`; the last class is `(p,q)=(1,13)`. For that
ratio, the first window lies strictly inside the open cross-comb component
`(6/91,8/91)`, so the endpoint method really fails. The isolated body-safe
phase

```text
y=3/14,                         x=y/2=3/28
```

repairs the row, with tail gaps `3/28` and `11/28`. The Python and
independent C++ artifacts freeze this alternate certificate and reproduce
the final `60,301,653,510` census.

No proved theorem should depend on THM-4149. Cite THM-4148 for the live
statement and this file only for the alternate audit/provenance.

# Independent audit of the eight later ready targets

**PASS — full source and primary read, normal/optimized/frozen replay, and
independent reconstruction of every actual graph and observer.** This audits
[TARGETS.md](TARGETS.md) as a finite positive witness table. It does not audit
an unprovided named transport or an order-minimality claim.

The eight counts are exactly

```
27*1429 + 9*d1 + 3*d2 + d3,  with d1,d2,d3 in {2,4}:
38609,38611,38615,38617,38627,38629,38633,38635.
```

I checked the complete source, its literal constructors, and its use of the
already audited ready-target and joint-response implementations. The source
reconstructs each matrix from its actual core, entry set and binary tail;
it does not read an asserted signature as an existence certificate. Every
auxiliary has its designated root as sole outgoing neighbour. Direct deleted-
graph counts establish L and K; the general sole-exit proof then establishes
the exact joint response G=F, with scalar gap2 in every row.

A separate verifier, `surgery_targets_audit.py`, imports no producer code.
It independently checks each literal tournament, strongness, the claimed
order, root and auxiliary, and sole-exit adjacency. It rebuilds every funnel,
indegree-one birth and converse directly from the definitions. It also checks
that all eight actual cores are strong, have orders7–10, and have total
starting response2 over the specified entries.

Its fresh backward endpoint recurrence recomputes H and every root-ending
response. The recurrence on the converse gives S; separate induced deletion
graphs give L and K. All nine reported quantities agree in every row:
H,F,S,L,K,G,root defect,sigma,tau. In particular sigma2 and tau0 are actual
joint input conditions, not deductions from scalar H alone. The resulting
witness orders are15,14,15,14,15,14,14,14, respectively.

The scope statements are correct. These are witness orders, with no proof of
global minimality. Canonical roots and auxiliaries have not been transported
from any particular earlier physical name set. A smaller target cannot
replace a larger child while retaining its old names. Prefix padding can
repair order in the weak category, but this table supplies no strong named
continuation and no return-class proof. The earlier finite bank's missing
responses are repaired by larger actual cores; the global spectrum and
unbounded closure do not follow.

## Reproduction and pins

From the session worktree:

```
python3 -B 04-computation/nonzero_preparation_20261103/surgery_targets.py --output /tmp/nonzero_targets_audit.json
python3 -B -O 04-computation/nonzero_preparation_20261103/surgery_targets.py --output /tmp/nonzero_targets_audit_optimized.json
python3 -B 04-computation/nonzero_preparation_20261103/surgery_targets_audit.py > /tmp/nonzero_targets_independent.json
```

Both producer modes match frozen `targets.json` and `targets.out` byte for
byte:56 always-active gates,8 witnesses. The independent verifier passes73
checks and reproduces [TARGETS_AUDIT.json](TARGETS_AUDIT.json). I checked every
length and hash in TARGETS_MANIFEST.json. No producer file was modified.

|Artifact|Bytes|SHA256|
|---|---:|---|
|surgery_targets.py|3560|efaf4f130d675d61e5b742758e6dc814c84808c2e41a3d8b35fb73bf56101f5a|
|targets.json|6393|6c2693e85bc1560ab4132b3eec5a30da8bb773fe6a689e39e85d78011b635a78|
|targets.out|137|f533d8f2cdc8cd472958fd6d7f52484b8d50a8acd07998239f8840aac56bccd6|
|TARGETS.md, audited prepromotion version|2550|c0b8f41e6408a15dc5530aa9d2502197ce5388fe2314a70fa133d3bf123f236b|
|surgery_targets_audit.py|3160|1e2fe78627b0fd096c9cde6dfc9dd9fb255da12ccd35f5a32cc3e811e5d0aad0|
|TARGETS_AUDIT.json|2522|3ca688a1f26cba533caed61e7bd7dacfb3d5057fcc0093df9def550b64002fcf|

No correction remains. Promotion and integration belong to the parent.

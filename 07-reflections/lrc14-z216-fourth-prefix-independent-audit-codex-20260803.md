# Independent audit of THM-3320

**Status: VERIFIED-EXACT independent audit; no new theorem claim.**

This artifact records an independent reconstruction of
`THM-3320-projected-k3-z216-fourth-ruler-prefix-and-affine-multicover-closure`.
The canonical theorem appeared concurrently after this calculation had already
reconstructed and screened the target prefix.  The audit was then pinned to
the incoming theorem, source, output, and semantic hashes and converted into a
field-by-field discrepancy check.

## Verdict

The discrepancy set is empty.  The independent calculation agrees with
THM-3320 on all ten explicit comparison classes:

1. the reconstructed complete-family queue and strict cost boundary;
2. selected indices, costs, bodies, component counts, and packet hash;
3. family and rowwise state/crude/status/residual counts;
4. the integral zero-constant modular support taxonomy;
5. the full affine common-table support taxonomy;
6. the three affine demand/cost/gap calculations;
7. the row-64 residual control;
8. the genuine row-138 support-three control;
9. the ledger, wall, family, and cap consequences; and
10. the projected-only logical scope.

The exact shared numerical packet is

```text
queue: gcd24/L76440 row 141, then gcd24/L30576 rows 133,219,359
screen: 230 = 172 crude + 58 status + 0 residual
rowwise: 141:10=8+2; 133:10=0+10; 219:189=147+42; 359:21=17+4
integral modular support counts: 55/1/2 at support 1/2/3
full affine-dual support counts: 55/3 at support 1/2
ledger: 373157 -> 373153; wall: 353 -> 349; families: 31 -> 29
```

The independently reproduced comparable hashes are

```text
selected packet: b8c63cd3f0b78cd83048158873dfdf8078f26ebd3ac6bab6c643b7292f7ae27d
selected screen: ea3f2619819dba921e3d874d1592d1e07b9792c14bb5d42c5b29bb7f991387ce
atlas rows:      53db9e1d3df2cf2b0398847682d909da81705e43a53ae2553d102fd152337649
```

## Independent checks

The audit does not import THM-3320's affine-majorant search.  It instead:

- rebuilds all `480` atlas rows through two inherited parsers;
- reconstructs all three prior cost prefixes from the original closure sets;
- reruns both the wrapper screen and the direct exact status engine;
- independently checks all `58` inherited Farkas vectors over the rationals;
- recompiles every status packet with THM-3308's finite integral template;
- verifies the three repaired affine inequalities directly on all sixteen
  Boolean patterns;
- constructs exact tables for all `18` singleton systems and all `32`
  empty-tail marginal systems;
- rejects one deliberately negative Farkas-alpha mutation on each selected
  row;
- replays row 64's eight residuals and exact feasible common table; and
- proves row 138 remains rank three by exact feasibility of all six singleton
  and fifteen pair systems.

The two half-affine circuits have exact data

```text
T=(1,5), a=1/2, w=(1/2,1/2,1/2,3/2), demand/cost/gap=3458/3276/182.
```

The third circuit has

```text
T=(2,5), a=0, w=(1,1,2,0), demand/cost/gap=3058/2548/510.
```

## Representation-level differences, not discrepancies

Some private audit hashes intentionally differ from THM-3320's hashes because
the persisted record schemas differ.  In particular, this audit's status
record omits the canonical affine-search metadata, and its row-64 packet omits
the full list of 57 one-layer kills.  Counts, addresses, certificates, exact
tables, and consequences agree wherever the records are semantically
comparable.  No mathematical discrepancy was found.

## Orbit connection for the next family

The two apparent integral support-three packets are exactly related by an
`S_4` permutation of the four status coordinates.  Canonicalizing their
`(q,marginals,capacities,histogram)` packets under coordinate permutations
reduces the three multi-layer addresses to two orbits, with multiplicities
`2` and `1`.  This makes the next natural operation precise: screen the
nineteen-row `gcd72/L7056` family with a cache keyed by coordinate-orbit
canonical capacity signatures, while retaining row labels as a sidecar.

The same pair also independently exposes the template/full-dual separation:
its zero-constant integral binary support is three, but its rational affine
tail support is exactly two.  Each singleton is feasible, whereas `(1,5)` is
infeasible.

## Reproduction and hashes

From the repository root:

```powershell
python 04-computation/lrc14_j7_k3_z216_fourth_prefix_independent_audit_20260803.py
python -O 04-computation/lrc14_j7_k3_z216_fourth_prefix_independent_audit_20260803.py
```

Both modes produce the same 23-line LF-normalized transcript.

```text
source_sha256:   bd4dac8af105dfd81c544fade471c3e296b31bbeb6bf0587824d97d46e9bef36
output_sha256:   db2a54804501d5bfc0fe30d8fc9fd6ad361668fa94a7a5dcf85985f845a0e209
semantic_sha256: 1ecc946ac1eaf9731b16b377ada123e3bb38437736a5897b6644fb415c587a35
normal/-O:       byte-identical
```

Canonical pins used by the audit:

```text
THM-3320 theorem:  be64218cace2d1b67b9be38e1f6e4405d192406c1e0957888f284cb21573a0b1
THM-3320 source:   b515c70174d58ad859a08c29949cdee36a4e04122451a66237732570ab5ee213
THM-3320 output:   d1845611f27d427a1d38afe349ed07bc964590fd39a3e88cd45e6ea34a86bc38
THM-3320 semantic: c201276eb84f71806a7eb683e42d723b930f7cd0f79862e8d6b1e0ad07d37dd9
```

Everything here remains inside the necessary projected `k=3,z1=216` atlas.
It yields no physical entry, arbitrary-`k` result, rung, or LRC(14) conclusion.

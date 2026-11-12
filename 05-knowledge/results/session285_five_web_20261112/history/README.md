# History lane: retrieval is broad; mathematical bridges are typed

Read [MECHANISM_ATLAS.md](MECHANISM_ATLAS.md) for the synthesis and proposed
experiments. `typed_connections.json` is a small mathematical connection graph;
`cooccurrence.json` is a different, purely lexical graph. Do not confuse them.

## Reproduce

From the repository root:

```sh
python 05-knowledge/results/session285_five_web_20261112/history/catalog.py --frozen-scope 05-knowledge/results/session285_five_web_20261112/history/searched_manifest.json.gz
python 05-knowledge/results/session285_five_web_20261112/history/history_catalog.py
python 05-knowledge/results/session285_five_web_20261112/history/check_mechanisms.py
python -O 05-knowledge/results/session285_five_web_20261112/history/check_mechanisms.py
```

The first command scans a pinned tree and the frozen path/blob selection.
It requires those Git objects locally; missing local blobs remain explicitly
omitted. Omitting `--frozen-scope` conducts a fresh availability-dependent scan
and may search more paths if additional objects have since been fetched.
It does not fetch every missing object or scan every historical file version.
The exact original metrics are in `summary.json`; all original exclusions in
`unavailable.json`; historical-subject scope in `commit_scope.json`.

`inventory.json.gz` records paths, Git blob IDs, query counts and initial matching
line excerpts. `expansion.json.gz` records nonseed paths touching at least two
chosen lateral topic vocabularies. These are leads, not claims of complete
reading. Separate deep consultation is described in the atlas. The path count
is not a count of unique sessions, mathematical ideas or proofs.

`lateral_sample.json` gives the deterministic random-like excursion rule,
selected paths, content hashes and section line numbers. It selects among
locally materialized reflections, not all repo objects. The prose records
which results remain historical hypotheses or finite evidence.

`check_mechanisms.py` uses explicit exceptions (not removable assertions).
The normal and optimized JSON outputs are identical in this run. These small
controls check the displayed elementary lattice mechanism and distinctions
between counting laws; they do not reverify all recovered research results.

Large catalogs are deterministic gzip (mtime zero); `catalog_hashes.json` records compressed and decompressed hashes. The script emits raw JSON when rerun.

Original scope: **34,864 of34,887 candidate paths** read; **7,509 seed-hit paths** (6,877 distinct blobs), **10,069 expansion-only paths**. All23 omissions are listed explicitly. Stress-query hits at the pinned baseline: **0**. Reachable shallow history: **4,541 commit subjects**,157 matching. These are retrieval metrics, not a full-session reading claim.

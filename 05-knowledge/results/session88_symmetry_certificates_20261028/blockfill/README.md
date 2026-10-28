# Session88: exact root-switch block filling

Start with [BLOCK_FILLING.md](BLOCK_FILLING.md).

- **Proved:** one complete **named** positive row identifies and fills exactly its
  equivalence fiber; aligned maps transport all full automorphism/root data.
- **Finite-exact through order eight:** same complete Q, raw bank, parent maps,
  every pointed root's complete transporter set and full automorphism group.
- **Order eight:** 434,048 → 55,040 root classifications; local matched full-output
  benchmark 12.255 s → 3.323 s. Not an asymptotic claim.
- **Second experiment:** certified subgroup propagation reduces queries further
  to 47,391, but did not produce a runtime win in its single illustrative pass.
- **Hostile:** anonymous `[7]` does not identify an object or joint ports. The
  two cyclic Z7 examples have different named keys, groups and pair orbits.
- **Factored API:** exact Qv in linear arithmetic work over pointed states,
  retaining the dictionary and transport sidecar.

## Files

| File | Role |
|---|---|
| `baseline_sources.py` | Hash-verified loader for unchanged Session87 source; no vendoring |
| `check_blockfill.py` | Matched full-output comparison, independent audits and optional repeated benchmark |
| `check_adaptive.py` | Optional certified subgroup propagation |
| `check_named_dictionary.py` | Direct S7 audit of anonymous-row hostile |
| `check_transport_replay.py` | Direct adjacency checks of every transporter and every sparse Q² row |
| `check_factored_operator.py` | Immutable row factors, exact Qv and Q²v |
| `verify.py` | Aggregate replay against five pinned semantic result hashes |

No raw receipt or baseline snapshot is needed. Each checker emits JSON locally;
redirect stdout to create a receipt if desired. `verify.py` ignores timing fields
and checks deterministic full-output hashes, audit results, operation counts,
scope metadata and hostile witnesses against compact expected hashes.

From this folder in the repository:

```sh
python3 -B verify.py
python3 -B -O verify.py
```

For an external staging folder, explicitly point to the source repository:

```sh
python3 -B -O verify.py --baseline-root /path/to/repository
```

The loader rejects missing or changed baseline bytes. Scope is ordinary
(non-converse-merged) tournaments, all masks, no hidden family filter.
Full labelled census is through six; independent full-child/group and complete
root-switch audits are through eight. No higher computational coverage claimed.

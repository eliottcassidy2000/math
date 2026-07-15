#!/usr/bin/env python3
"""Export the exact raw uint16 n=8 merged-node atlas used by THM-818.

The source is the certified 6,880-class exhaustive classifier.  Grid
reflection realizes converse, so an ordinary class and its reflected class
are sorted into a pair and the 3,528 distinct pairs are ranked.  The binary
format is deliberately minimal: 2^21 little-endian uint16 values in mask
order, with no header.  The C++ consumer rechecks all H8 and R8 census totals.
"""

from __future__ import annotations

import argparse
import contextlib
import hashlib
import io
import json
import runpy
from pathlib import Path

import numpy as np


LEGACY = Path("04-computation/metagraph_flow_n8_check_opus_S305.py")
EXPECTED_SHA256 = "30debad3387a4ea0ef51108ea132115efda2ac2fcdfcc2c5c1d4d23155095835"


def export(legacy: Path, output: Path, metadata: Path | None) -> dict:
    with contextlib.redirect_stdout(io.StringIO()):
        ns = runpy.run_path(str(legacy))
    bits = ns["bits_all"]
    cls = ns["cls_of"]
    refl = ns["refl"]
    n_masks, n_bits = int(ns["N"]), int(ns["m"])
    reflected = np.zeros(n_masks, dtype=np.int64)
    for source in range(n_bits):
        reflected |= ((bits >> source) & 1) << int(refl[source])
    pairs = np.stack((np.minimum(cls, cls[reflected]), np.maximum(cls, cls[reflected])), axis=1)
    unique, rank = np.unique(pairs, axis=0, return_inverse=True)
    atlas = rank.astype("<u2", copy=False)
    assert n_masks == 1 << 21 and n_bits == 21
    assert len(set(map(int, cls))) == 6880 and len(unique) == 3528
    assert int(atlas.min()) == 0 and int(atlas.max()) == 3527
    output.parent.mkdir(parents=True, exist_ok=True)
    atlas.tofile(output)
    digest = hashlib.sha256(output.read_bytes()).hexdigest()
    assert output.stat().st_size == 4_194_304
    assert digest == EXPECTED_SHA256
    result = {
        "schema_version": 1,
        "source_classifier": str(legacy),
        "format": "raw little-endian uint16 merged-node ranks in mask order",
        "masks": n_masks,
        "ordinary_classes": 6880,
        "merged_nodes": len(unique),
        "bytes": output.stat().st_size,
        "sha256": digest,
        "output": str(output),
    }
    if metadata:
        metadata.write_text(json.dumps(result, indent=2) + "\n")
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--legacy", type=Path, default=LEGACY)
    parser.add_argument("--output", type=Path, default=Path("/tmp/n8_merged_node_rank_u16.bin"))
    parser.add_argument("--metadata", type=Path)
    args = parser.parse_args()
    result = export(args.legacy, args.output, args.metadata)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()

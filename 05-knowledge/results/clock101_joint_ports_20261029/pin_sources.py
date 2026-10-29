"""Pin local historical sources linked by this packet, without changing them."""
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]


def read_bytes(path):
    # The repository contains canonical filenames exceeding legacy MAX_PATH.
    text = str(path.resolve())
    if os.name == "nt" and not text.startswith("\\\\?\\"):
        text = "\\\\?\\" + text
    return Path(text).read_bytes()


def main():
    paths = set()
    for name in ("README.md", "RECOVERY.md"):
        for link in re.findall(r"\]\(([^)]+)\)", (HERE / name).read_text(encoding="utf-8")):
            if link.startswith("../"):
                paths.add((HERE / link.split("#")[0]).resolve())
    records = []
    for path in sorted(paths):
        relative = path.relative_to(ROOT).as_posix()
        data = read_bytes(path)
        records.append({"path": relative, "sha256": hashlib.sha256(data).hexdigest(),
                        "bytes": len(data), "inspection": "targeted passages; not a full independent audit"})
    revision = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    output = {"source_revision": revision, "source_count": len(records), "sources": records}
    (HERE / "source_hashes.json").write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"source_count": len(records), "source_revision": revision}))


if __name__ == "__main__":
    main()

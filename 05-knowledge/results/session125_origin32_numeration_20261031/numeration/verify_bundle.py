"""Verify this exact portable payload from any working directory."""

import hashlib
import json
from pathlib import Path, PurePosixPath


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def main():
    root = Path(__file__).resolve().parent
    paths = (root / "WHITELIST.txt").read_text().splitlines()
    require(paths == sorted(set(paths)), "whitelist must be sorted and unique")
    for name in paths:
        path = PurePosixPath(name)
        require(name and not path.is_absolute() and ".." not in path.parts,
                "unsafe whitelist path: " + name)
    actual = sorted(str(p.relative_to(root)) for p in root.rglob("*") if p.is_file())
    require(actual == paths, "payload files do not exactly match whitelist")
    digest_lines = (root / "SHA256SUMS").read_text().splitlines()
    hashes = {}
    for line in digest_lines:
        digest, name = line.split("  ", 1)
        require(name not in hashes, "duplicate checksum: " + name)
        hashes[name] = digest
    require(set(hashes) == set(paths) - {"SHA256SUMS"}, "checksum coverage mismatch")
    for name, digest in hashes.items():
        require(hashlib.sha256((root / name).read_bytes()).hexdigest() == digest,
                "checksum mismatch: " + name)
    manifest = json.loads((root / "PUBLICATION_MANIFEST.json").read_text())
    expected = set(paths) - {"SHA256SUMS", "WHITELIST.txt", "PUBLICATION_MANIFEST.json"}
    require(set(manifest["files"]) == expected, "manifest coverage mismatch")
    for name, metadata in manifest["files"].items():
        data = (root / name).read_bytes()
        require(len(data) == metadata["bytes"], "size mismatch: " + name)
        require(hashlib.sha256(data).hexdigest() == metadata["sha256"],
                "manifest digest mismatch: " + name)
    for name in paths:
        require(not (root / name).is_symlink(), "symlink in payload: " + name)
    print(json.dumps({"status": "PASS", "files": len(paths),
                      "checksum_entries": len(hashes), "strict_whitelist": True},
                     sort_keys=True))


if __name__ == "__main__":
    main()

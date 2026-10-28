"""Pinned, reproducible search coverage; occurrence is not endorsement or review."""
from pathlib import Path
import collections
import json
import re
import subprocess

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
PATTERNS = {
    "66": r"(?<![\w])66(?:th)?(?![\w])|sixty[ -]?six",
    "5": r"(?<![\w])5(?![\w])|\bfive\b",
    "15": r"(?<![\w])15(?![\w])|\bfifteen\b",
    "platonic_fermat": r"platonic|fermat.prime|five.solids|icosahed|dodecahed",
    "moonshine": r"supersingular|baby.monster|moonshine|\bmonster\b",
    "divisor_balance": r"F\s*=\s*U\s*\+\s*S|F\s*=\s*S\s*\+\s*U|divisor.balance|p\^?2\s*q\s*r",
    "graph_triple": r"petersen|K[_ ]?\{?3[, ]+3|kuratowski",
    "graceful_square": r"graceful|square.sum|square.sigma|aliquot",
    "berggren": r"berggren",
    "seven_twentyone": r"7\s*[,/–-]\s*21|seven.twenty.?one|7.and.21",
    "idoneal": r"idoneal|one.class.per.genus|one class in each genus",
}

def git(*args):
    return subprocess.check_output(["git", *args], cwd=ROOT)

def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--revision", required=True)
    args = ap.parse_args()
    rev = git("rev-parse", args.revision).decode().strip()
    # Read committed prose, including narrow historical matches, not generated
    # numeric tables, copied agent inboxes, binary artifacts, or this packet.
    entries = []
    for row in git("ls-tree", "-rz", rev).split(b"\0"):
        if not row:
            continue
        meta, path = row.split(b"\t", 1)
        mode, kind, oid = meta.decode().split()
        name = path.decode()
        if (kind == "blob" and name.endswith((".md", ".tex", ".org"))
                and not name.startswith(("agents/", ".", str(HERE.relative_to(ROOT))))):
            entries.append((name, oid))
    unique = list(dict.fromkeys(oid for _, oid in entries))
    proc = subprocess.Popen(["git", "cat-file", "--batch"], cwd=ROOT,
                            stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    compiled = {k: re.compile(v, re.I) for k, v in PATTERNS.items()}
    hits = {}
    for oid in unique:
        proc.stdin.write((oid + "\n").encode())
        proc.stdin.flush()
        head = proc.stdout.readline().split()
        size = int(head[2])
        body = proc.stdout.read(size).decode("utf-8", errors="replace")
        proc.stdout.read(1)
        d = {}
        for lineno, line in enumerate(body.splitlines(), 1):
            for key, pattern in compiled.items():
                if pattern.search(line):
                    d.setdefault(key, []).append(lineno)
        hits[oid] = d
    proc.stdin.close()
    if proc.wait() != 0:
        raise RuntimeError("cat-file failed")
    records = [{"path": name, "blob": oid, "hits": hits[oid]}
               for name, oid in entries if hits[oid]]
    data = {"revision": rev, "scope": "tracked md/tex/org, excluding agents, dot paths and this packet",
            "patterns": PATTERNS, "files_scanned": len(entries),
            "unique_blobs_scanned": len(unique), "matched_files": len(records),
            "files_by_theme": {k: sum(k in r["hits"] for r in records) for k in PATTERNS},
            "records": records}
    (HERE / "occurrences.json").write_text(json.dumps(data, separators=(",", ":")) + "\n")
    print(json.dumps({k: v for k, v in data.items() if k not in ("records", "patterns")}, indent=2))

if __name__ == "__main__":
    main()

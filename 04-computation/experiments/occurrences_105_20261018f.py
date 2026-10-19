"""Frozen-Git census of literal 105, with separate PDF/archive extraction.

Standard library plus pypdf. Reads Git blobs, never a moving working copy.
Every literal hit is retained as a one-based (line, column) address. Numeric
token hits forbid adjacent decimal digits; semantic relevance is NOT inferred.
The compressed manifest is deterministic (gzip mtime=0). No network access.
"""
from __future__ import annotations

import argparse
import collections
import gzip
import hashlib
import io
import json
from pathlib import Path
import re
import subprocess
import sys
import tarfile

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_REV = "1352f3b3bb4d7b95248fc1373d57711c222d69dd"
STEM = ROOT / "05-knowledge/results/occurrences_105_20261018f"
TOKEN = re.compile(r"(?<![0-9])105(?![0-9])")
ALIASES = re.compile(
    r"one hundred(?: and)?[ -]+five|\b7!!(?!\!)|"
    r"(?<![0-9])3\s*[*×·]\s*5\s*[*×·]\s*7(?![0-9])",
    re.IGNORECASE,
)


def git(*args: str) -> bytes:
    return subprocess.check_output(["git", *args], cwd=ROOT)


def candidates(revision: str, pattern: str, extended: bool = False) -> set[str]:
    cmd = ["git", "grep", "-I", "-l", "-z", "-i",
           "-E" if extended else "-F", "-e", pattern, revision, "--"]
    p = subprocess.run(cmd, cwd=ROOT, capture_output=True)
    if p.returncode not in (0, 1):
        raise RuntimeError(p.stderr.decode(errors="replace"))
    prefix = revision.encode() + b":"
    return {
        item[len(prefix):].decode("utf-8", "surrogateescape")
        for item in p.stdout.split(b"\0") if item
    }


class Blobs:
    def __init__(self):
        self.p = subprocess.Popen(
            ["git", "cat-file", "--batch"], cwd=ROOT,
            stdin=subprocess.PIPE, stdout=subprocess.PIPE,
        )

    def read(self, sha: str) -> bytes:
        self.p.stdin.write((sha + "\n").encode())
        self.p.stdin.flush()
        header = self.p.stdout.readline().split()
        if len(header) != 3 or header[1] != b"blob":
            raise RuntimeError(header)
        data = self.p.stdout.read(int(header[2]))
        assert self.p.stdout.read(1) == b"\n"
        return data

    def close(self):
        self.p.stdin.close()
        self.p.wait()


def scan(text: str) -> dict:
    records, examples, aliases = [], [], []
    literal_count = token_count = 0
    # Match Git/editor LF line addresses; splitlines would also split U+2028,
    # form feeds and vertical tabs inside a source line.
    for n, line in enumerate(text.split("\n"), 1):
        all_cols = [m.start() + 1 for m in re.finditer("105", line)]
        token_cols = [m.start() + 1 for m in TOKEN.finditer(line)]
        if all_cols:
            records.append([n, all_cols, token_cols])
            literal_count += len(all_cols)
            token_count += len(token_cols)
            # A short sample for EACH matching line, not only the first hit.
            c = (token_cols or all_cols)[0] - 1
            examples.append([n, line[max(0, c-100):c+180]])
        alias_cols = [[m.start()+1, m.group()] for m in ALIASES.finditer(line)]
        if alias_cols:
            aliases.append([n, alias_cols])
    return dict(literal_count=literal_count, token_count=token_count,
                lines=records, examples=examples, aliases=aliases)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--revision", default=DEFAULT_REV)
    ap.add_argument("--output-stem", type=Path, default=STEM)
    ap.add_argument("--find", help="Regex over saved matching-line excerpts")
    ap.add_argument("--path", default="", help="Regex over source paths")
    ap.add_argument("--limit", type=int, default=40)
    ap.add_argument("--skip", type=int, default=0)
    ap.add_argument("--list-files", action="store_true")
    ap.add_argument("--verify-index", action="store_true")
    ap.add_argument("--write-concordance", action="store_true")
    args = ap.parse_args()
    if args.find is not None or args.list_files or args.verify_index or args.write_concordance:
        sys.stdout.reconfigure(encoding="utf-8")
        raw = gzip.decompress(args.output_stem.with_suffix(".jsonl.gz").read_bytes())
        rows = [json.loads(line) for line in raw.splitlines()]
        if args.write_concordance:
            summary = json.loads(args.output_stem.with_suffix(".json").read_text())
            prose = [r for r in rows if r["kind"] == "git-text"
                     and r["path"].endswith(".md") and r["token_count"]]
            lines = [
                "# Repository 105 prose concordance", "",
                "Generated locator index, **not a proof or semantic classification**.",
                "Includes metadata, decimals, references and historical copies.",
                f"Frozen source: {summary['revision']}.",
                f"All {len(prose)} Markdown files with digit-delimited105 are listed.",
                "The compressed manifest additionally includes every other matching text",
                "file and all matching extracted PDF/gzip units, with every literal",
                "line and column address. Source links below resolve in this checkout;",
                "use the manifest's blob hashes for the immutable source.", "",
                "[Mathematical synthesis](105_research_atlas_20261018f.md) ·",
                "[Full manifest](occurrences_105_20261018f.jsonl.gz) ·",
                "[Scope and counts](occurrences_105_20261018f.json)", "",
                "In each row: numeric-token count; source line numbers. Long line lists",
                "are abbreviated after24 entries; the manifest never abbreviates them.", "",
            ]
            for root in sorted({r["path"].split("/")[0] for r in prose}):
                lines.extend([f"## {root}", ""])
                for r in prose:
                    if r["path"].split("/")[0] != root:
                        continue
                    ns = [str(n) for n,_,cols in r["lines"] if cols]
                    short = ",".join(ns[:24]) + ("; more in manifest" if len(ns)>24 else "")
                    lines.append(f"- [{r['path']}](../../{r['path']}): "
                                 f"{r['token_count']} hits; lines {short}.")
                lines.append("")
            dest = args.output_stem.parent/"105_prose_concordance_20261018f.md"
            dest.write_text("\n".join(lines).rstrip()+"\n",encoding="utf-8")
            print(f"Wrote {len(prose)} source-file locators to {dest.name}")
            return
        if args.verify_index:
            summary = json.loads(args.output_stem.with_suffix(".json").read_text())
            if hashlib.sha256(raw).hexdigest() != summary["manifest_sha256_uncompressed"]:
                raise RuntimeError("manifest digest mismatch")
            revision = summary["revision"]
            # An independent Git/PCRE pass checks EVERY file's numeric hit count.
            p = subprocess.run(
                ["git", "grep", "-I", "-o", "-z", "-P", "-e",
                 r"(?<![0-9])105(?![0-9])", revision, "--"],
                cwd=ROOT, capture_output=True, check=True,
            )
            observed = collections.Counter()
            for line in p.stdout.splitlines():
                key, value = line.rsplit(b"\0", 1)
                if value != b"105":
                    raise RuntimeError("unexpected grep format")
                observed[key[len(revision)+1:].decode("utf-8","surrogateescape")] += 1
            expected = {r["path"]:r["token_count"] for r in rows
                        if r["kind"] == "git-text" and r["token_count"]}
            if observed != expected:
                mismatch = [(p,expected.get(p,0),observed.get(p,0))
                            for p in sorted(set(expected)|set(observed))
                            if expected.get(p,0) != observed.get(p,0)]
                raise RuntimeError(f"independent Git/PCRE mismatch: {mismatch[:20]}")
            for row in rows:
                if sum(len(cs) for _,cs,_ in row["lines"]) != row["literal_count"]:
                    raise RuntimeError("literal address count mismatch")
                if sum(len(cs) for _,_,cs in row["lines"]) != row["token_count"]:
                    raise RuntimeError("numeric address count mismatch")
            print(f"Verified digest, every address count, and {sum(observed.values())}"
                  f" numeric hits in {len(observed)} files by independent Git/PCRE.")
            return
        selected = (r for r in rows if re.search(args.path,r["path"],re.I))
        found = []
        for row in selected:
            address = row["path"] + ("#"+row["unit"] if "unit" in row else "")
            if args.list_files:
                found.append(f"{address}\tliteral={row['literal_count']}\tnumeric={row['token_count']}")
            else:
                for n, excerpt in row["examples"]:
                    if re.search(args.find,excerpt,re.I):
                        found.append(f"{address}:{n}: {excerpt}")
        print("\n".join(found[args.skip:args.skip+args.limit]))
        print(f"Displayed {args.skip}..{min(args.skip+args.limit,len(found))} of {len(found)} rows."
              " Excerpts are search aids; use the recorded Git blob for full context.")
        return
    revision = git("rev-parse", args.revision).decode().strip()
    tree = {}
    for item in git("ls-tree", "-r", "-z", "-l", revision).split(b"\0"):
        if not item:
            continue
        metadata, path = item.split(b"\t", 1)
        mode, kind, sha, size = metadata.split()
        if kind == b"blob":
            tree[path.decode("utf-8", "surrogateescape")] = {
                "blob": sha.decode(), "bytes": int(size),
            }
    matches = candidates(revision, "105")
    alias_candidates = candidates(
        revision,
        r"one hundred( and)?[ -]+five|7!!|3[[:space:]]*[*×·][[:space:]]*5[[:space:]]*[*×·][[:space:]]*7",
        True,
    )
    derived = {p for p in tree if p.endswith((".pdf", ".gz", ".zip"))}
    blobs, output, extraction = Blobs(), [], []
    root_totals = collections.defaultdict(lambda: [0, 0, 0])
    for i, path in enumerate(sorted(matches | alias_candidates | derived), 1):
        data = blobs.read(tree[path]["blob"])
        if path in matches or path in alias_candidates:
            row = scan(data.decode("utf-8", "replace"))
            row.update(path=path, kind="git-text", **tree[path])
            output.append(row)
            counts = root_totals[path.split("/")[0]]
            counts[0] += 1
            counts[1] += row["literal_count"]
            counts[2] += row["token_count"]
        if path in derived:
            info = dict(path=path, **tree[path], extracted_units=0, error=None)
            try:
                units = []
                if path.endswith(".pdf"):
                    from pypdf import PdfReader
                    reader = PdfReader(io.BytesIO(data))
                    units = [(f"page:{j}", page.extract_text() or "")
                             for j, page in enumerate(reader.pages, 1)]
                elif path.endswith(".tar.gz"):
                    with tarfile.open(fileobj=io.BytesIO(data), mode="r:gz") as tf:
                        units = [(m.name, tf.extractfile(m).read().decode("utf-8", "replace"))
                                 for m in tf.getmembers() if m.isfile()]
                elif path.endswith(".gz"):
                    units = [("gunzip", gzip.decompress(data).decode("utf-8", "replace"))]
                else:
                    raise ValueError("Unsupported archive format: explicit exclusion")
                info["extracted_units"] = len(units)
                for unit, text in units:
                    row = scan(text)
                    if row["literal_count"] or row["aliases"]:
                        row.update(path=path, unit=unit, kind="extracted", **tree[path])
                        output.append(row)
            except Exception as exc:
                info["error"] = f"{type(exc).__name__}: {exc}"
            extraction.append(info)
        if i % 500 == 0:
            print(f"Scanned {i} candidate blobs", flush=True)
    blobs.close()
    raw = b"".join(
        (json.dumps(row, ensure_ascii=True, separators=(",", ":")) + "\n").encode()
        for row in output
    )
    target = args.output_stem.with_suffix(".jsonl.gz")
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_bytes(gzip.compress(raw, mtime=0))
    direct = [r for r in output if r["kind"] == "git-text"]
    extra = [r for r in output if r["kind"] == "extracted"]
    summary = {
        "scope": "All tracked Git text blobs at the exact revision; all tracked PDF and gzip payloads separately extracted",
        "revision": revision, "tracked_blobs": len(tree),
        "literal_text_files": len(matches),
        "numeric_token_text_files": sum(r["token_count"] > 0 for r in direct),
        "literal_text_occurrences": sum(r["literal_count"] for r in direct),
        "numeric_token_text_occurrences": sum(r["token_count"] for r in direct),
        "alias_text_files": sum(bool(r["aliases"]) for r in direct),
        "matching_extracted_units": len(extra),
        "extracted_literal_occurrences": sum(r["literal_count"] for r in extra),
        "extracted_numeric_token_occurrences": sum(r["token_count"] for r in extra),
        "manifest_sha256_uncompressed": hashlib.sha256(raw).hexdigest(),
        "manifest_compressed_bytes": target.stat().st_size,
        "line_address_format": "[one_based_line, all_literal_columns, numeric_token_columns]",
        "roots": dict(sorted(root_totals.items())),
        "paths_containing_105": sorted(p for p in tree if "105" in p),
        "extraction": extraction,
        "limits": [
            "Numeric token means no adjacent decimal digit, not mathematical relevance.",
            "Metadata, decimals, code, outputs, references and copied historical text are retained.",
            "PDF text extraction is not OCR; graphical-only numerals may be missed.",
            "Images, executable payloads, untracked/ignored files, Downloads attachments and historical-only deleted Git blobs are outside the text census.",
            "Aliases are a deliberately stated finite search family, not all algebraic expressions evaluating to105.",
        ],
    }
    args.output_stem.with_suffix(".json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="utf-8"
    )
    print(json.dumps({k: v for k, v in summary.items()
                      if k not in ("roots", "paths_containing_105", "extraction")}, indent=2))


if __name__ == "__main__":
    main()

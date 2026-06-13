#!/usr/bin/env python3
"""
Math Research Inbox Processor
==============================
Run this to pre-process files dropped into the inbox/ directory.
Generates PROCESSING-REPORT.md for Claude to read and integrate.

Usage:
    python3 inbox/processor.py

Output:
    inbox/PROCESSING-REPORT.md   — Claude reads this to integrate content
    inbox/processed/YYYY-MM-DD/  — Processed files are moved here
    inbox/.hash-db.json          — Hash database of all known files
"""

import os
import sys
import hashlib
import json
import difflib
import shutil
from datetime import datetime
from pathlib import Path

# ─── Paths ────────────────────────────────────────────────────────────────────
SCRIPT_DIR   = Path(__file__).resolve().parent
REPO_ROOT    = SCRIPT_DIR.parent
INBOX_DIR    = SCRIPT_DIR
PROCESSED_DIR = INBOX_DIR / "processed"
HASH_DB_PATH  = INBOX_DIR / ".hash-db.json"
REPORT_PATH   = INBOX_DIR / "PROCESSING-REPORT.md"

# ─── Config ───────────────────────────────────────────────────────────────────
SIMILARITY_THRESHOLD = 0.85     # above this → near-duplicate, show diff
PREVIEW_LINES        = 40       # lines of new-file content to show in report
MAX_DIFF_LINES       = 80       # max diff lines to include in report

# Files that are part of the system itself, never treated as contributions
SYSTEM_FILES = {
    "processor.py", "README.md", "PROCESSING-REPORT.md",
    ".hash-db.json", ".gitkeep", ".DS_Store"
}

# Extensions that can be read as text
TEXT_EXTENSIONS = {
    ".md", ".txt", ".tex", ".py", ".js", ".ts", ".html",
    ".json", ".csv", ".rst", ".ipynb"
}

# Directories to skip when indexing the repo
SKIP_DIRS = {"inbox", ".git", "__pycache__", "processed"}


# ─── Utilities ────────────────────────────────────────────────────────────────

def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while chunk := f.read(65536):
            h.update(chunk)
    return h.hexdigest()


def read_text(path: Path) -> list[str] | None:
    if path.suffix.lower() not in TEXT_EXTENSIONS:
        return None
    try:
        with open(path, encoding="utf-8", errors="replace") as f:
            return f.readlines()
    except Exception:
        return None


def index_repo() -> dict[str, str]:
    """Walk the repo (excluding inbox) and build {hash: relative_path}."""
    db: dict[str, str] = {}
    for path in REPO_ROOT.rglob("*"):
        if path.is_file():
            # Skip system / inbox / git dirs
            parts = set(path.relative_to(REPO_ROOT).parts)
            if parts & SKIP_DIRS:
                continue
            if path.name in SYSTEM_FILES:
                continue
            try:
                file_hash = sha256(path)
                db[file_hash] = str(path.relative_to(REPO_ROOT))
            except Exception:
                pass
    return db


def find_most_similar(
    new_lines: list[str], repo_index: dict[str, str]
) -> tuple[str, float, list[str]] | None:
    """Return (rel_path, similarity_ratio, existing_lines) for best match."""
    best_ratio = 0.0
    best_match = None
    for _, rel_path in repo_index.items():
        full_path = REPO_ROOT / rel_path
        if not full_path.exists():
            continue
        existing_lines = read_text(full_path)
        if not existing_lines:
            continue
        ratio = difflib.SequenceMatcher(None, new_lines, existing_lines).ratio()
        if ratio > best_ratio:
            best_ratio = ratio
            best_match = (rel_path, ratio, existing_lines)
    return best_match if best_match and best_match[1] > 0 else None


def archive(path: Path, disposition: str) -> Path:
    """Move processed file to inbox/processed/YYYY-MM-DD/<disposition>/."""
    dest_dir = PROCESSED_DIR / datetime.now().strftime("%Y-%m-%d") / disposition
    dest_dir.mkdir(parents=True, exist_ok=True)
    dest = dest_dir / path.name
    # Avoid overwriting if same filename already exists
    if dest.exists():
        stem = path.stem
        suffix = path.suffix
        ts = datetime.now().strftime("%H%M%S")
        dest = dest_dir / f"{stem}-{ts}{suffix}"
    shutil.move(str(path), str(dest))
    return dest


# ─── Main ─────────────────────────────────────────────────────────────────────

def process_inbox():
    PROCESSED_DIR.mkdir(parents=True, exist_ok=True)

    # Collect inbox files (exclude system files)
    inbox_files = sorted([
        f for f in INBOX_DIR.iterdir()
        if f.is_file() and f.name not in SYSTEM_FILES
    ])

    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    lines = [
        f"# Inbox Processing Report\n\n",
        f"**Generated:** {timestamp}  \n",
        f"**Files found:** {len(inbox_files)}  \n\n",
        "---\n\n",
    ]

    if not inbox_files:
        lines.append("No files found in inbox. Nothing to process.\n")
        _write_report(lines)
        print("No files to process.")
        return

    # Build repo index
    print("Indexing repository for duplicate detection...")
    repo_index = index_repo()
    print(f"  Indexed {len(repo_index)} files in repo.")

    # Load persistent hash db (hashes of previously processed files)
    persist_db: dict[str, str] = {}
    if HASH_DB_PATH.exists():
        try:
            with open(HASH_DB_PATH) as f:
                persist_db = json.load(f)
        except Exception:
            pass

    # Merge: everything in repo + previously processed
    all_known_hashes = {**repo_index, **persist_db}

    new_persist = dict(persist_db)  # will be updated with processed files

    for path in inbox_files:
        lines.append(f"## `{path.name}`\n\n")
        lines.append(f"*Size: {path.stat().st_size:,} bytes | Extension: {path.suffix or 'none'}*\n\n")

        try:
            file_hash = sha256(path)
        except Exception as e:
            lines.append(f"⚠️ **ERROR reading file:** {e}\n\n")
            continue

        # ── Exact duplicate ──
        if file_hash in all_known_hashes:
            existing_path = all_known_hashes[file_hash]
            lines.append(f"**STATUS: EXACT DUPLICATE ✓**  \n")
            lines.append(f"Identical to: `{existing_path}`  \n")
            lines.append("No action needed. File archived.\n\n")
            dest = archive(path, "duplicate")
            new_persist[file_hash] = str(dest.relative_to(INBOX_DIR.parent))
            lines.append("---\n\n")
            continue

        # ── Text analysis for near-duplicates / new files ──
        new_lines = read_text(path)

        if new_lines:
            # ── Near-duplicate check ──
            match = find_most_similar(new_lines, repo_index)
            if match and match[1] >= SIMILARITY_THRESHOLD:
                rel_path, ratio, existing_lines = match
                lines.append(f"**STATUS: NEAR-DUPLICATE** (similarity: {ratio:.1%})  \n")
                lines.append(f"Most similar to: `{rel_path}`  \n\n")

                # Unified diff (capped)
                diff = list(difflib.unified_diff(
                    existing_lines, new_lines,
                    fromfile=f"existing: {rel_path}",
                    tofile=f"inbox: {path.name}",
                    n=3, lineterm=""
                ))
                if diff:
                    lines.append("**Diff (existing → new):**\n\n")
                    lines.append("```diff\n")
                    lines.extend(l + "\n" for l in diff[:MAX_DIFF_LINES])
                    if len(diff) > MAX_DIFF_LINES:
                        lines.append(f"... ({len(diff) - MAX_DIFF_LINES} more lines not shown)\n")
                    lines.append("```\n\n")
                else:
                    lines.append("*(No textual diff — may be whitespace-only difference)*\n\n")

                lines.append("### Claude action required:\n")
                lines.append("- Review the diff above for any **new mathematical content** not in the existing file\n")
                lines.append("- Check for **discrepancies** (conflicting claims) — if found, open a court case in `02-court/active/`\n")
                lines.append("- If the new file adds meaningful new content, integrate it; otherwise discard\n\n")

                dest = archive(path, "near-duplicate")
                new_persist[file_hash] = str(dest.relative_to(INBOX_DIR.parent))
                lines.append("---\n\n")
                continue

            # ── New file ──
            lines.append("**STATUS: NEW CONTRIBUTION 🆕**\n\n")
            lines.append("**Content preview:**\n\n")
            lines.append("```\n")
            lines.extend(new_lines[:PREVIEW_LINES])
            if len(new_lines) > PREVIEW_LINES:
                lines.append(f"\n... ({len(new_lines) - PREVIEW_LINES} more lines)\n")
            lines.append("```\n\n")

            # Suggest where this might belong based on extension/content
            suggestions = _suggest_destination(path, new_lines)
            if suggestions:
                lines.append("**Suggested integration:**\n")
                for s in suggestions:
                    lines.append(f"- {s}\n")
                lines.append("\n")

            lines.append("### Claude action required:\n")
            lines.append("- Read the full file (it's been archived; path shown below)\n")
            lines.append("- Extract theorems → `01-canon/theorems/`\n")
            lines.append("- Extract open questions → `00-navigation/OPEN-QUESTIONS.md`\n")
            lines.append("- Extract tangents → `00-navigation/TANGENTS.md`\n")
            lines.append("- Extract mistakes → `01-canon/MISTAKES.md`\n")
            lines.append("- Note any discrepancies with existing canon → open court case if needed\n")
            lines.append("- Update `00-navigation/SESSION-LOG.md`\n\n")

        else:
            # Binary or unreadable file
            lines.append("**STATUS: NEW BINARY/UNREADABLE FILE**\n\n")
            lines.append("File could not be read as text. If this is an image, move it to `03-artifacts/images/` and update `03-artifacts/images/INDEX.md`.\n\n")

        dest = archive(path, "new")
        new_persist[file_hash] = str(dest.relative_to(INBOX_DIR.parent))
        lines.append(f"*Archived to: `{dest.relative_to(REPO_ROOT)}`*\n\n")
        lines.append("---\n\n")

    # ── Summary ──
    lines.append("## Summary\n\n")
    lines.append(f"Processed {len(inbox_files)} file(s). See sections above for required actions.\n\n")
    lines.append("When done integrating, add an entry to `00-navigation/SESSION-LOG.md`.\n")

    # Save updated hash db
    with open(HASH_DB_PATH, "w") as f:
        json.dump(new_persist, f, indent=2)

    _write_report(lines)
    print(f"\nReport written to: {REPORT_PATH}")
    print("Hand Claude this report and ask it to integrate the content.")


def _write_report(lines):
    with open(REPORT_PATH, "w", encoding="utf-8") as f:
        f.writelines(lines)


def _suggest_destination(path: Path, lines: list[str]) -> list[str]:
    """Heuristic suggestions for where this file's content belongs."""
    content = "".join(lines[:100]).lower()
    suggestions = []

    if any(kw in content for kw in ["theorem", "lemma", "proof", "corollary", "claim"]):
        suggestions.append("Contains theorem/proof material → extract to `01-canon/theorems/`")
    if any(kw in content for kw in ["mistake", "bug", "error", "wrong", "incorrect"]):
        suggestions.append("Contains error/bug report → add to `01-canon/MISTAKES.md`")
    if any(kw in content for kw in ["open question", "conjecture", "unsolved", "unknown"]):
        suggestions.append("Contains open questions → add to `00-navigation/OPEN-QUESTIONS.md`")
    if any(kw in content for kw in ["tangent", "idea", "connection to", "suggests", "hypothesis"]):
        suggestions.append("Contains speculative ideas → add to `00-navigation/TANGENTS.md`")
    if path.suffix == ".py":
        suggestions.append("Python script → archive to `03-artifacts/code/` and index in CODE INDEX")
    if path.suffix == ".tex":
        suggestions.append("LaTeX source → archive to `03-artifacts/drafts/`, extract math content to Markdown")
    if path.suffix in {".png", ".jpg", ".svg", ".pdf"}:
        suggestions.append("Image/figure → move to `03-artifacts/images/` and update INDEX.md")

    return suggestions


if __name__ == "__main__":
    process_inbox()

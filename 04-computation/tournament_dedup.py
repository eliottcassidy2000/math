#!/usr/bin/env python3
"""
tournament_dedup.py -- Tournament-Powered Data Deduplicator
kind-pasteur-2026-03-24-S20cq

Fast near-duplicate detection using tournament hashing + structural fingerprinting.
Combines tournament_hash (LSH) with dataprint (structural analysis) for
multi-level deduplication:

  Level 1: Exact hash (SHA-256) -- identical files
  Level 2: Tournament hash -- structural duplicates (same structure, different encoding)
  Level 3: Fingerprint similarity -- similar structure (e.g. same template, different data)

APPLICATIONS:
  - Dataset cleaning for ML training (remove near-duplicates)
  - Document deduplication (legal, academic, corporate)
  - Image similarity detection (via byte-level structural hash)
  - Log file deduplication (different timestamps, same events)
  - Backup optimization (avoid storing near-identical versions)
  - Code clone detection (find copy-pasted code)

USAGE:
  python tournament_dedup.py scan /path/to/files     # find duplicates
  python tournament_dedup.py --threshold 0.2 scan .   # stricter matching
  python tournament_dedup.py compare file1 file2      # compare two files
  python tournament_dedup.py cluster /path/to/files   # group by similarity

LICENSE: MIT
"""

import os
import sys
import hashlib
import time
import json
import argparse
from collections import defaultdict
from typing import List, Tuple, Dict, Optional

# Add script dir to path
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
if SCRIPT_DIR not in sys.path:
    sys.path.insert(0, SCRIPT_DIR)

from dataprint import fingerprint, similarity as fp_similarity, DataPrint
from tournament_hash import TournamentHasher, TournamentHash

__version__ = "1.0.0"


class TournamentDedup:
    """Multi-level deduplication engine."""

    def __init__(self, threshold: float = 0.3, n_vertices: int = 12):
        """
        Args:
            threshold: tournament hash distance threshold for near-duplicates
            n_vertices: tournament hash size (12 = 66 bits, good balance)
        """
        self.threshold = threshold
        self.hasher = TournamentHasher(n_vertices=n_vertices)
        self._files = {}       # path -> FileInfo
        self._sha256 = {}      # sha256 -> [paths]
        self._thash = {}       # path -> TournamentHash
        self._fprint = {}      # path -> DataPrint

    def add_file(self, path: str) -> dict:
        """Add a file to the dedup database.

        Returns: dict with hash info and any detected duplicates.
        """
        with open(path, 'rb') as f:
            data = f.read()

        size = len(data)
        sha = hashlib.sha256(data).hexdigest()

        # Level 1: Exact duplicate check
        exact_dups = []
        if sha in self._sha256:
            exact_dups = list(self._sha256[sha])

        if sha not in self._sha256:
            self._sha256[sha] = []
        self._sha256[sha].append(path)

        # Level 2: Tournament hash
        thash = self.hasher.hash_bytes(data)
        self._thash[path] = thash

        # Level 3: Structural fingerprint
        fprint = fingerprint(data)
        self._fprint[path] = fprint

        # Find near-duplicates by tournament hash distance
        near_dups = []
        for other_path, other_hash in self._thash.items():
            if other_path == path:
                continue
            dist = thash.normalized_distance(other_hash)
            if dist <= self.threshold:
                near_dups.append((other_path, dist))

        near_dups.sort(key=lambda x: x[1])

        info = {
            'path': path,
            'size': size,
            'sha256': sha,
            'tournament_hash': thash.hex(),
            'class': fprint.class_label,
            'entropy': fprint.entropy,
            'exact_duplicates': exact_dups,
            'near_duplicates': near_dups[:10],
        }
        self._files[path] = info
        return info

    def scan_directory(self, directory: str, extensions: Optional[list] = None,
                       max_size: int = 100 * 1024 * 1024) -> dict:
        """Scan a directory for duplicates.

        Args:
            directory: path to scan
            extensions: list of extensions to include (None = all)
            max_size: skip files larger than this (default 100MB)

        Returns: dict with summary and duplicate groups.
        """
        t0 = time.time()
        files_scanned = 0
        files_skipped = 0
        total_bytes = 0

        for root, dirs, files in os.walk(directory):
            # Skip hidden dirs and common non-content dirs
            dirs[:] = [d for d in dirs if not d.startswith('.') and d not in
                       ('node_modules', '__pycache__', '.git', 'venv', '.venv')]

            for fname in files:
                if extensions:
                    ext = os.path.splitext(fname)[1].lower()
                    if ext not in extensions:
                        continue

                path = os.path.join(root, fname)
                try:
                    size = os.path.getsize(path)
                    if size > max_size or size == 0:
                        files_skipped += 1
                        continue
                    self.add_file(path)
                    files_scanned += 1
                    total_bytes += size
                except (PermissionError, OSError):
                    files_skipped += 1

        elapsed = time.time() - t0

        # Build duplicate groups
        exact_groups = {}
        for sha, paths in self._sha256.items():
            if len(paths) > 1:
                exact_groups[sha[:16]] = paths

        # Near-duplicate clusters (greedy clustering)
        near_groups = self._cluster_near_duplicates()

        return {
            'files_scanned': files_scanned,
            'files_skipped': files_skipped,
            'total_bytes': total_bytes,
            'elapsed_seconds': elapsed,
            'exact_duplicate_groups': exact_groups,
            'near_duplicate_clusters': near_groups,
            'unique_files': len(self._sha256),
            'wasted_bytes': sum(
                sum(os.path.getsize(p) for p in paths[1:])
                for paths in self._sha256.values()
                if len(paths) > 1
            ),
        }

    def _cluster_near_duplicates(self) -> list:
        """Cluster files by tournament hash similarity."""
        paths = list(self._thash.keys())
        visited = set()
        clusters = []

        for i, p1 in enumerate(paths):
            if p1 in visited:
                continue
            cluster = [p1]
            visited.add(p1)

            for j in range(i + 1, len(paths)):
                p2 = paths[j]
                if p2 in visited:
                    continue
                dist = self._thash[p1].normalized_distance(self._thash[p2])
                if dist <= self.threshold:
                    cluster.append(p2)
                    visited.add(p2)

            if len(cluster) > 1:
                clusters.append(cluster)

        return clusters

    def compare(self, path1: str, path2: str) -> dict:
        """Detailed comparison between two files."""
        with open(path1, 'rb') as f:
            data1 = f.read()
        with open(path2, 'rb') as f:
            data2 = f.read()

        sha1 = hashlib.sha256(data1).hexdigest()
        sha2 = hashlib.sha256(data2).hexdigest()
        th1 = self.hasher.hash_bytes(data1)
        th2 = self.hasher.hash_bytes(data2)
        fp1 = fingerprint(data1)
        fp2 = fingerprint(data2)

        return {
            'file1': path1,
            'file2': path2,
            'size1': len(data1),
            'size2': len(data2),
            'exact_match': sha1 == sha2,
            'tournament_distance': th1.normalized_distance(th2),
            'fingerprint_similarity': fp_similarity(fp1, fp2),
            'class1': fp1.class_label,
            'class2': fp2.class_label,
            'entropy1': fp1.entropy,
            'entropy2': fp2.entropy,
            'verdict': _verdict(sha1 == sha2, th1.normalized_distance(th2),
                               fp_similarity(fp1, fp2)),
        }


def _verdict(exact: bool, t_dist: float, fp_sim: float) -> str:
    if exact:
        return "EXACT DUPLICATE"
    if t_dist < 0.1:
        return "NEAR-IDENTICAL (structural)"
    if t_dist < 0.3 and fp_sim > 0.95:
        return "SIMILAR (same template)"
    if fp_sim > 0.9:
        return "RELATED (same data type)"
    return "DIFFERENT"


# ============================================================================
# CLI
# ============================================================================

def cmd_scan(args):
    dedup = TournamentDedup(threshold=args.threshold)
    result = dedup.scan_directory(args.directory)

    print(f"Tournament Dedup v{__version__}")
    print("=" * 60)
    print(f"  Scanned:     {result['files_scanned']} files ({result['total_bytes']:,} bytes)")
    print(f"  Skipped:     {result['files_skipped']} files")
    print(f"  Time:        {result['elapsed_seconds']:.1f}s")
    print(f"  Unique SHA:  {result['unique_files']}")

    if result['exact_duplicate_groups']:
        print(f"\n  EXACT DUPLICATES ({len(result['exact_duplicate_groups'])} groups):")
        for sha_prefix, paths in result['exact_duplicate_groups'].items():
            print(f"    [{sha_prefix}...] {len(paths)} copies:")
            for p in paths[:5]:
                print(f"      {p}")

    if result['near_duplicate_clusters']:
        print(f"\n  NEAR-DUPLICATE CLUSTERS ({len(result['near_duplicate_clusters'])} groups):")
        for cluster in result['near_duplicate_clusters'][:10]:
            print(f"    Group ({len(cluster)} files):")
            for p in cluster[:5]:
                print(f"      {p}")

    wasted = result['wasted_bytes']
    if wasted > 0:
        print(f"\n  WASTED SPACE: {wasted:,} bytes ({wasted/1024/1024:.1f} MB)")


def cmd_compare(args):
    dedup = TournamentDedup(threshold=args.threshold)
    result = dedup.compare(args.file1, args.file2)

    print(f"Tournament Dedup v{__version__} -- Compare")
    print("=" * 60)
    print(f"  File 1:      {result['file1']} ({result['size1']:,} bytes, {result['class1']})")
    print(f"  File 2:      {result['file2']} ({result['size2']:,} bytes, {result['class2']})")
    print(f"  Exact match: {result['exact_match']}")
    print(f"  T-distance:  {result['tournament_distance']:.3f}")
    print(f"  FP similarity: {result['fingerprint_similarity']:.3f}")
    print(f"  Verdict:     {result['verdict']}")


def cmd_demo(args):
    """Run demo on repo files."""
    import tempfile

    print(f"Tournament Dedup v{__version__} -- Demo")
    print("=" * 60)

    # Create some test files
    tmpdir = tempfile.mkdtemp(prefix='td_demo_')

    files = {}
    # Exact duplicate pair
    content = b"the quick brown fox jumps over the lazy dog\n" * 100
    for name in ['doc_a.txt', 'doc_a_copy.txt']:
        path = os.path.join(tmpdir, name)
        with open(path, 'wb') as f:
            f.write(content)
        files[name] = path

    # Near-duplicate (slightly modified)
    content2 = b"the fast brown fox leaps over the lazy dog\n" * 100
    path = os.path.join(tmpdir, 'doc_b.txt')
    with open(path, 'wb') as f:
        f.write(content2)
    files['doc_b.txt'] = path

    # Different file
    content3 = b"quantum computing uses qubits for parallel processing\n" * 100
    path = os.path.join(tmpdir, 'doc_c.txt')
    with open(path, 'wb') as f:
        f.write(content3)
    files['doc_c.txt'] = path

    # Binary file
    path = os.path.join(tmpdir, 'data.bin')
    with open(path, 'wb') as f:
        f.write(bytes(range(256)) * 40)
    files['data.bin'] = path

    # Scan
    dedup = TournamentDedup(threshold=0.3)
    result = dedup.scan_directory(tmpdir)

    print(f"\n  Created {len(files)} test files in {tmpdir}")
    print(f"  Scanned: {result['files_scanned']} files")
    print(f"  Exact duplicate groups: {len(result['exact_duplicate_groups'])}")
    print(f"  Near-duplicate clusters: {len(result['near_duplicate_clusters'])}")

    # Pairwise comparisons
    print(f"\n  Pairwise comparisons:")
    names = list(files.keys())
    for i in range(len(names)):
        for j in range(i + 1, len(names)):
            comp = dedup.compare(files[names[i]], files[names[j]])
            print(f"    {names[i]:>20} vs {names[j]:<20}: "
                  f"t_dist={comp['tournament_distance']:.3f} "
                  f"fp_sim={comp['fingerprint_similarity']:.3f} "
                  f"-> {comp['verdict']}")

    # Cleanup
    import shutil
    shutil.rmtree(tmpdir)
    print(f"\n  Cleaned up temp files")


def main():
    parser = argparse.ArgumentParser(
        description=f'Tournament Dedup v{__version__}',
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('--threshold', '-t', type=float, default=0.3,
                       help='Near-duplicate threshold (default 0.3)')

    sub = parser.add_subparsers(dest='command')

    p_scan = sub.add_parser('scan', help='Scan directory for duplicates')
    p_scan.add_argument('directory', help='Directory to scan')

    p_compare = sub.add_parser('compare', help='Compare two files')
    p_compare.add_argument('file1')
    p_compare.add_argument('file2')

    p_demo = sub.add_parser('demo', help='Run demo')

    args = parser.parse_args()

    if args.command == 'scan':
        cmd_scan(args)
    elif args.command == 'compare':
        cmd_compare(args)
    elif args.command == 'demo':
        cmd_demo(args)
    else:
        parser.print_help()


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Check packaged hashes, Python syntax, and local Markdown destinations."""
import ast
import hashlib
import json
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parent


def main():
    manifests = []
    for manifest in sorted(ROOT.rglob('SHA256SUMS')):
        count = 0
        for line in manifest.read_text().splitlines():
            if not line.strip():
                continue
            digest, name = line.split(None, 1)
            target = manifest.parent / name.lstrip('*')
            if not target.is_file():
                raise RuntimeError(f'Missing hash target: {target}')
            if hashlib.sha256(target.read_bytes()).hexdigest() != digest:
                raise RuntimeError(f'Hash mismatch: {target}')
            count += 1
        manifests.append({'path': manifest.relative_to(ROOT).as_posix(), 'entries': count})
    scripts = sorted(ROOT.rglob('*.py'))
    for script in scripts:
        ast.parse(script.read_text(), filename=script.relative_to(ROOT).as_posix())
    links = 0
    frozen_documents = 0
    for document in sorted(ROOT.rglob('*.md')):
        if 'inputs' in document.relative_to(ROOT).parts:
            frozen_documents += 1
            continue
        for target in re.findall(r'\[[^\]]*\]\(([^\s)]+)\)', document.read_text()):
            if '://' in target or target.startswith(('mailto:', '#')):
                continue
            path = target.split('#', 1)[0]
            if not (document.parent / path).exists():
                raise RuntimeError(f'Missing local link: {document.relative_to(ROOT)} -> {target}')
            links += 1
    root_manifest = ROOT / 'SHA256SUMS'
    if not root_manifest.is_file():
        raise RuntimeError('Root SHA256SUMS is missing')
    print(json.dumps({'status': 'PASS', 'manifests': manifests,
                      'python_syntax_files': len(scripts), 'authored_local_markdown_destinations': links,
                      'frozen_input_markdown_excluded_from_link_check': frozen_documents,
                      'root_manifest_sha256': hashlib.sha256(root_manifest.read_bytes()).hexdigest()},
                     indent=2, sort_keys=True))


if __name__ == '__main__':
    main()

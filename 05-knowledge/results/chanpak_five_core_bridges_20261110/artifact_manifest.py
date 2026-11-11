#!/usr/bin/env python3
"""Check the frozen research packet and its executable external dependency."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
MANIFEST = HERE / 'MANIFEST.json'
DEPENDENCIES = ['05-knowledge/results/session251_creative_completion_20261111/ternary/compiler.py']


def digest(path):
    # Git checkout newline conversion is not a mathematical change.
    return hashlib.sha256(path.read_bytes().replace(b'\r\n',b'\n')).hexdigest()


def packet_files():
    return sorted(p for p in HERE.rglob('*') if p.is_file()
                  and p.suffix in ('.md','.py','.cpp','.json','.txt')
                  and p != MANIFEST and '__pycache__' not in p.parts)


def local_links():
    checked = 0
    for source in HERE.rglob('*.md'):
        body = source.read_text(encoding='utf-8')
        for raw in re.findall(r'\]\(([^)]+)\)',body):
            target = raw.strip('<>').split('#',1)[0]
            if not target or '://' in target or target.startswith('mailto:'):
                continue
            destination = (source.parent / target).resolve()
            if os.name == 'nt' and not str(destination).startswith('\\\\?\\'):
                destination = Path('\\\\?\\'+str(destination))
            if not destination.exists():
                raise RuntimeError(f'Broken local link: {source.relative_to(HERE)} -> {raw}')
            checked += 1
    return checked


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--write',action='store_true',help='Freeze after reviewed changes')
    args = parser.parse_args()
    current = {p.relative_to(HERE).as_posix():digest(p) for p in packet_files()}
    dependencies = {p:digest(ROOT/p) for p in DEPENDENCIES}
    if args.write:
        data = {'status':'COMPLETE research packet; Chan–Pak and all-height five-core coverage OPEN',
                'hash_convention':'SHA256 after replacing CRLF with LF; all listed files are text',
                'initial_baseline':'0f96e904a570',
                'concurrent_prior_algebra_commit':'cb7a17db97',
                'files':current,'executable_external_dependencies':dependencies}
        MANIFEST.write_text(json.dumps(data,indent=2,sort_keys=True)+'\n',encoding='utf-8')
    else:
        frozen = json.loads(MANIFEST.read_text(encoding='utf-8'))
        if frozen['files'] != current:
            differences = sorted(k for k in set(current)|set(frozen['files'])
                                 if current.get(k) != frozen['files'].get(k))
            raise RuntimeError(f'Packet differs from frozen manifest: {differences}')
        if frozen['executable_external_dependencies'] != dependencies:
            raise RuntimeError('Inherited executable dependency changed')
    print(json.dumps({'status':'PASS','packet_files':len(current),
                      'external_dependencies':len(dependencies),
                      'local_links':local_links()}))


if __name__ == '__main__':
    main()

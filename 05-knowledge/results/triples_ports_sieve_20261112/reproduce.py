"""Replay this packet only; pin all required external repository inputs.

Run at repository root or from any directory. Text hashes and comparisons normalize CRLF to LF;
all subprocesses use -B, and --optimized adds -O without disabling checks.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
ENGINES = [
    ('arithmetic/wheel_refinement.py', 'arithmetic/wheel_refinement.json', []),
    ('arithmetic/joint_branch_sieve.py', 'arithmetic/joint_branch_sieve.json', []),
    ('arithmetic/audit_root.py', 'arithmetic/audit_root.json', []),
    ('dihedral/check.py', 'dihedral/check.json', []),
    ('triples/check.py', 'triples/check.json', []),
    ('squares/verify.py', 'squares/receipt.json', ['--output', 'receipt.json']),
    ('squares/audit_mutation.py', 'squares/mutation_audit.json', []),
    ('ports/verify.py', 'ports/receipt.json', []),
    ('ports/four_local.py', 'ports/four_local_receipt.json', []),
    ('ports/scout_returns.py', 'ports/scout_receipt.json', []),
]


def native(path):
    return Path('\\\\?\\'+str(path)) if os.name == 'nt' and not str(path).startswith('\\\\?\\') else path


def digest(path):
    raw = native(path).read_bytes()
    if path.suffix in (".md", ".py", ".json"):
        raw = raw.replace(b"\r\n", b"\n")
    return hashlib.sha256(raw).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--optimized', action='store_true')
    parser.add_argument('--seal-only', action='store_true')
    args = parser.parse_args()
    rows = []
    if not args.seal_only:
        for script, receipt, extra in ENGINES:
            path = HERE/receipt
            before = digest(path) if path.exists() else None
            cmd = [sys.executable]+(['-O'] if args.optimized else [])+['-B', str(HERE/script)]+extra
            run = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True, encoding='utf-8')
            if run.returncode:
                raise RuntimeError(f'{script}: {run.stdout}\n{run.stderr}')
            after = digest(path)
            if before and before != after:
                raise RuntimeError(f'{receipt} changed; inspect the new receipt before resealing')
            rows.append({'script':script,'receipt':receipt,'sha256':after,'exit_code':0})
            print(f'PASS {script}', flush=True)
        (HERE/('REPLAY_OPTIMIZED.json' if args.optimized else 'REPLAY_NORMAL.json')).write_text(
            json.dumps({'optimized':args.optimized,'engines':rows},indent=2)+'\n',encoding='utf-8')
    # Every linked local source is pinned; missing links fail instead of silently
    # treating a historical or renamed file as an available dependency.
    linked = {
        ROOT/'05-knowledge/results/session290_core_descent_20261113/descent/profiles.json.gz',
        ROOT/'05-knowledge/results/endpoint_coverage_compiler_20261111/wheel.json',
    }
    links = 0
    for note in sorted(HERE.rglob('*.md')):
        prose = re.sub(r'```.*?```', '', note.read_text(encoding='utf-8'), flags=re.S)
        prose = re.sub(r'`[^`]*`', '', prose)
        for target in re.findall(r'\]\(([^)]+)\)', prose):
            if '://' in target or target.startswith('#'):
                continue
            target = target.split('#')[0].strip('<>')
            if not target:
                continue
            path = (note.parent/target).resolve()
            if not native(path).is_file():
                raise RuntimeError(f'Broken link {note.relative_to(HERE)} -> {target}')
            links += 1
            if not path.is_relative_to(HERE):
                linked.add(path)
    (HERE/'SOURCE_PINS.json').write_text(json.dumps(
        {p.relative_to(ROOT).as_posix():digest(p) for p in sorted(linked)},indent=2)+'\n',encoding='utf-8')
    (HERE/'LINK_CHECK.json').write_text(json.dumps({'local_links':links,'status':'PASS'},indent=2)+'\n')
    manifest = {}
    for p in sorted(HERE.rglob('*')):
        if p.is_file() and '__pycache__' not in p.parts and p.name != 'MANIFEST.json':
            manifest[p.relative_to(HERE).as_posix()] = digest(p)
    (HERE/'MANIFEST.json').write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf-8')
    print(f'SEAL {len(manifest)} files, {links} local links, {len(linked)} external source pins')


if __name__ == '__main__':
    main()

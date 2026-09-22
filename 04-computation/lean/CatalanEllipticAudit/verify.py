"""Fail-fast clean build and axiom-free audit of the twelve finite witnesses."""
from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parent
SOURCES = (
    'lean-toolchain', 'lakefile.toml', 'lake-manifest.json',
    'CatalanEllipticAudit.lean', 'CatalanEllipticAudit/Certificates.lean',
    'AxiomAudit.lean', 'README.md', 'verify.py', '.gitignore',
)


def require(test, message):
    if not test:
        raise RuntimeError(message)


def write(path, text):
    path.write_text(text, encoding='utf-8', newline='\n')


def main():
    sys.stdout.reconfigure(encoding='utf-8')
    write(ROOT / 'verification.json', '{"status":"RUNNING"}\n')
    transcript = []

    def run(args):
        executable = shutil.which(args[0])
        require(executable is not None, f'Missing executable: {args[0]}')
        result = subprocess.run([executable, *args[1:]], cwd=ROOT,
                                text=True, encoding='utf-8', errors='replace',
                                capture_output=True, check=False)
        output = result.stdout + result.stderr
        transcript.append(f"$ {' '.join(args)}\n{output}exit_code={result.returncode}\n")
        write(ROOT / 'verification.log', '\n'.join(transcript))
        require(result.returncode == 0, f"Nonzero exit: {' '.join(args)}\n{output}")
        return output

    for name in SOURCES:
        require((ROOT / name).is_file(), f'Missing source: {name}')
        raw = (ROOT / name).read_bytes()
        raw.decode('utf-8')
        require(b'\r' not in raw, f'Non-LF source: {name}')
        if name.endswith('.lean'):
            source = raw.decode('utf-8')
            require(not re.search(r'\b(sorry|admit|native_decide|sorryAx)\b', source),
                    f'Forbidden proof escape: {name}')
            require(not re.search(r'^\s*axiom\s', source, re.MULTILINE),
                    f'Custom axiom: {name}')
    actual_lean = {p.relative_to(ROOT).as_posix() for p in ROOT.rglob('*.lean')
                   if '.lake' not in p.parts}
    require(actual_lean == {s for s in SOURCES if s.endswith('.lean')},
            'Unexpected or missing Lean source')
    names = re.findall(r'^theorem\s+(\w+)',
                       (ROOT / 'CatalanEllipticAudit/Certificates.lean').read_text(encoding='utf-8'),
                       flags=re.MULTILINE)
    require(len(names) == len(set(names)) == 12, 'Public theorem inventory changed')
    require((ROOT / 'CatalanEllipticAudit.lean').read_text(encoding='utf-8') ==
            'import CatalanEllipticAudit.Certificates\n', 'Public import reach changed')
    audit = (ROOT / 'AxiomAudit.lean').read_text(encoding='utf-8')
    require(audit.startswith('import CatalanEllipticAudit\n'), 'Audit bypasses public root')
    for name in names:
        require(f'#print axioms CatalanEllipticAudit.{name}\n' in audit,
                f'Missing theorem in audit: {name}')
    manifest = json.loads((ROOT / 'lake-manifest.json').read_text(encoding='utf-8'))
    require(manifest.get('packages') == [], 'External dependency introduced')

    version = run(['lean', '--version']).strip()
    require('version 4.30.0' in version, 'Incorrect Lean version')
    run(['lake', 'clean'])
    run(['lake', 'build'])
    output = run(['lake', 'env', 'lean', 'AxiomAudit.lean'])
    axiom_map = {}
    for name in names:
        qualified = 'CatalanEllipticAudit.' + name
        require(f"'{qualified}' does not depend on any axioms" in output,
                f'Theorem is not certified axiom-free: {qualified}')
        axiom_map[qualified] = []
    record = {'status': 'PASS', 'lean_version': version,
              'scope': 'Twelve concrete cleared arithmetic witnesses, not global classifications',
              'public_root': 'CatalanEllipticAudit', 'external_packages': [],
              'theorems_audited': len(names), 'axioms': axiom_map,
              'sha256': {name: sha256((ROOT / name).read_bytes()).hexdigest() for name in SOURCES}}
    write(ROOT / 'verification.json', json.dumps(record, sort_keys=True, indent=2) + '\n')
    print('PASS: clean public-root build; all 12 concrete certificates are axiom-free.')


if __name__ == '__main__':
    try:
        main()
    except Exception as error:
        write(ROOT / 'verification.json', json.dumps({'status': 'FAIL', 'error': str(error)}, indent=2) + '\n')
        raise

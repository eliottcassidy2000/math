#!/usr/bin/env python3
"""Rebuild optional detailed audit ledgers outside the frozen packet, then hash them."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parent
OMITTED = tuple('reconstructed/independent_sources_N%d.json' % n for n in (512, 1024, 4096))

def require(ok, message):
    if not ok:
        raise RuntimeError(message)

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--audit-root', type=Path, help='Optional external audit root for pre-assembly validation.')
    parser.add_argument('--wrong-expected-hash', action='store_true',
                        help='Hostile control: deliberately corrupt one expected detailed-ledger hash.')
    args = parser.parse_args()
    audit = args.audit_root.resolve() if args.audit_root else ROOT / 'session75-dyadic-audit'
    raw_manifest = (audit / 'MANIFEST.audit.json').read_bytes()
    manifest = json.loads(raw_manifest)
    expected = {name: manifest['files'][name] for name in OMITTED}
    stored_count = 0
    for name, digest in manifest['files'].items():
        if name in OMITTED:
            continue
        relative = Path(name)
        require(not relative.is_absolute() and '..' not in relative.parts and name == relative.as_posix(),
                'unsafe stored audit-manifest path')
        file = audit / relative
        require(file.is_file() and not file.is_symlink(), 'missing/nonregular stored audit entry')
        require(hashlib.sha256(file.read_bytes()).hexdigest() == digest,
                'stored audit-manifest SHA256 mismatch: ' + name)
        stored_count += 1
    require(all(len(value) == 64 and all(c in '0123456789abcdef' for c in value)
                for value in expected.values()), 'invalid detailed-ledger digest')
    if args.wrong_expected_hash:
        expected[OMITTED[0]] = '0' * 64
    rows = []
    with tempfile.TemporaryDirectory(prefix='session75-dyadic-reconstructed-') as temp:
        output = Path(temp).resolve()
        require(not output.is_relative_to(ROOT.resolve()), 'generated output must be outside frozen packet')
        flags = ['-I', '-B'] + (['-O'] if sys.flags.optimize else [])
        env = {k: v for k, v in os.environ.items() if k not in ('PYTHONHOME', 'PYTHONPATH')}
        env['PYTHONDONTWRITEBYTECODE'] = '1'
        proc = subprocess.run([sys.executable, *flags, str(audit / 'independent_verify.py'),
                               '--producer', str(audit / 'producer_revision2'),
                               '--write-sources', str(output)], cwd=audit, env=env,
                              capture_output=True, timeout=1200)
        require(proc.returncode == 0, 'independent detailed-source reconstruction failed:\n' +
                proc.stderr.decode(errors='replace'))
        expected_receipt = (audit / 'receipt.compared.json').read_bytes()
        require(proc.stdout == expected_receipt, 'independent native-source comparison receipt differs')
        require(sorted(p.name for p in output.iterdir()) == sorted(Path(p).name for p in OMITTED),
                'unexpected detailed-source output inventory')
        for name in OMITTED:
            file = output / Path(name).name
            require(file.is_file() and not file.is_symlink(), 'nonregular reconstructed source file')
            data = file.read_bytes()
            digest = hashlib.sha256(data).hexdigest()
            require(digest == expected[name], 'regenerated omitted artifact SHA256 mismatch: ' + name)
            rows.append({'path': name, 'bytes': len(data), 'sha256': digest})
        receipt_sha256 = hashlib.sha256(proc.stdout).hexdigest()
    require(not output.exists(), 'temporary detailed-source outputs not removed')
    print(json.dumps({'status': 'PASS', 'generated_outside_frozen_packet': True,
                      'temporary_outputs_removed': True, 'stored_audit_manifest_entries_verified': stored_count,
                      'unchanged_audit_manifest_sha256': hashlib.sha256(raw_manifest).hexdigest(),
                      'complete_native_source_comparison_receipt_sha256': receipt_sha256,
                      'verified_omitted_outputs': rows}, indent=2, sort_keys=True))

if __name__ == '__main__':
    main()

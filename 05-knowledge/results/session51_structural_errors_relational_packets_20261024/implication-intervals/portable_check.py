#!/usr/bin/env python3
"""Real-copy replay, optimized replay, and disposable manifest mutation controls."""
import hashlib
import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def run(args, cwd, expect_success=True):
    p = subprocess.run(args, cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if expect_success and p.returncode:
        raise RuntimeError(p.stderr.decode() or p.stdout.decode())
    if not expect_success and not p.returncode:
        raise RuntimeError('mutation was not rejected: ' + repr(args))
    return p.stdout


def main():
    reference = (ROOT / 'receipt.json').read_bytes()
    with tempfile.TemporaryDirectory(prefix='session51-implication-portable-') as tmp:
        root = Path(tmp)
        packet = root / 'packet'
        unrelated = root / 'unrelated-working-directory'
        unrelated.mkdir()
        shutil.copytree(ROOT, packet, ignore=shutil.ignore_patterns('__pycache__'))
        run([sys.executable, str(packet / 'verify_manifest.py')], unrelated)
        normal = run([sys.executable, str(packet / 'check.py')], unrelated)
        optimized = run([sys.executable, '-O', str(packet / 'check.py')], unrelated)
        if not normal == optimized == reference:
            raise RuntimeError('normal/optimized/receipt mismatch in copied packet')
        rejected = []
        for label in ('source-change', 'input-change', 'missing-theorem', 'extra-file', 'extra-nested-manifest'):
            damaged = root / label
            shutil.copytree(packet, damaged)
            if label == 'source-change':
                p = damaged / 'check.py'
                p.write_bytes(p.read_bytes() + b'\n# hostile manifest mutation\n')
            elif label == 'input-change':
                p = damaged / 'inputs' / 'SESSION49_NONRETURN.md'
                p.write_bytes(p.read_bytes() + b'\nHostile input mutation.\n')
            elif label == 'missing-theorem':
                (damaged / 'THEOREM.md').unlink()
            elif label == 'extra-file':
                (damaged / 'UNDECLARED.txt').write_text('unmanifested file')
            else:
                (damaged / 'inputs' / 'MANIFEST.json').write_text('{}')
            run([sys.executable, '-O', str(damaged / 'verify_manifest.py')], unrelated, False)
            rejected.append(label)
    print(json.dumps({'schema': 'session51-implication-portability-v1', 'real_temporary_copy': True, 'changed_working_directory': True, 'normal_optimized_receipt_byte_equal': True, 'receipt_sha256': hashlib.sha256(reference).hexdigest(), 'manifest_checked_before_replay': True, 'manifest_mutations_rejected': rejected}, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()

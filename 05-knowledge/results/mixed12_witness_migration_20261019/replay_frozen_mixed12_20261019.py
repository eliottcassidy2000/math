#!/usr/bin/env python3
"""Verify the frozen package, then replay both exact certificates under -O.
No third-party dependencies. Use --check-only for integrity/normalization checks.
Full replay uses a temporary directory and leaves the frozen files unchanged.
"""
import argparse
import ast
import gzip
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parent
SCRIPTS = ('verify_mixed12_witness_migration_20261019.py',
           'verify_mixed12_exception_independent_20261019.py')
RECEIPTS = ('mixed12_witness_migration_audit_20261019.json',
            'mixed12_exception_independent_audit_20261019.json')


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def sha(data):
    return hashlib.sha256(data).hexdigest()


class ExplicitChecks(ast.NodeTransformer):
    def __init__(self):
        self.count = 0

    def visit_Assert(self, node):
        self.count += 1
        message = node.msg or ast.Constant('Verification failed')
        return ast.copy_location(ast.If(
            test=ast.UnaryOp(op=ast.Not(), operand=node.test),
            body=[ast.Raise(exc=ast.Call(func=ast.Name(id='RuntimeError', ctx=ast.Load()),
                                       args=[message], keywords=[]), cause=None)],
            orelse=[]), node)


class CanonicalDiagnostic(ast.NodeTransformer):
    # Ignore only generated error-message wording, which ast.unparse can format
    # differently between Python releases. No condition or arithmetic is changed.
    def visit_Raise(self, node):
        call = node.exc
        if (isinstance(call, ast.Call) and isinstance(call.func, ast.Name)
                and call.func.id == 'RuntimeError' and len(call.args) == 1
                and isinstance(call.args[0], ast.Constant)
                and isinstance(call.args[0].value, str)
                and call.args[0].value.startswith('Verification failed: ')):
            call.args[0] = ast.Constant('Verification failed')
        return node


def inspect_package():
    manifest = json.loads((ROOT / 'manifest.json').read_text())
    for entry in manifest['files']:
        data = (ROOT / entry['path']).read_bytes()
        require(len(data) == entry['bytes'] and sha(data) == entry['sha256'],
                'Manifest mismatch: ' + entry['path'])
    checks = []
    for name in SCRIPTS:
        transformer = ExplicitChecks()
        expected = transformer.visit(ast.parse((ROOT / 'originals' / name).read_bytes()))
        actual = CanonicalDiagnostic().visit(ast.parse((ROOT / name).read_bytes()))
        require(ast.dump(expected, include_attributes=False) == ast.dump(actual, include_attributes=False),
                'Changes beyond assertion normalization: ' + name)
        require(not any(isinstance(node, ast.Assert) for node in ast.walk(actual)),
                'An optimization-sensitive assertion remains: ' + name)
        checks.append({'script': name, 'explicit_checks': transformer.count})
    main = gzip.decompress((ROOT / (RECEIPTS[0] + '.gz')).read_bytes())
    reference = manifest['preserved_main_receipt']
    require(len(main) == reference['uncompressed_bytes'] and
            sha(main) == reference['uncompressed_sha256'], 'Main receipt preservation mismatch')
    independent = (ROOT / RECEIPTS[1]).read_bytes()
    return manifest, checks, (main, independent)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check-only', action='store_true')
    args = parser.parse_args()
    manifest, checks, references = inspect_package()
    report = {'status': 'PASS', 'mode': 'integrity-only' if args.check_only else 'full-exact-replay',
              'checks_active_under_python_O': checks,
              'scope': 'Mixed (1,2) or dual (2,1) maximum-span critical roots with c <= 3 only.'}
    if not args.check_only:
        replay = []
        with tempfile.TemporaryDirectory(prefix='mixed12-exact-replay-') as temp:
            work = Path(temp)
            for i, name in enumerate(SCRIPTS):
                shutil.copyfile(ROOT / name, work / name)
                command = [sys.executable, '-O', str(work / name)]
                if i == 0:
                    command += ['--max-common', '3']
                completed = subprocess.run(command, check=True, text=True, capture_output=True)
                generated = (work / RECEIPTS[i]).read_bytes()
                require(generated == references[i], 'Receipt replay differs: ' + RECEIPTS[i])
                replay.append({'receipt': RECEIPTS[i], 'byte_for_byte_identical': True,
                               'bytes': len(generated), 'sha256': sha(generated),
                               'stdout_status': json.loads(completed.stdout)['status']})
        report['replayed_receipts'] = replay
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()

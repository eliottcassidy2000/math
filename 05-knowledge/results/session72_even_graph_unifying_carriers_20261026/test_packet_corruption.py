#!/usr/bin/env python3
"""Exercise frozen-packet rejection gates without modifying the real packet."""
import argparse
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

def need(ok, why):
    if not ok:
        raise RuntimeError(why)

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('packet', type=Path)
    args = parser.parse_args()
    source = args.packet.resolve()
    native = json.loads((source / 'FILE_MANIFEST.json').read_text())
    chosen = next(r['path'] for r in native['files'] if r['path'] != 'verify_packet.py')
    cases = [('baseline', None), ('append-byte', 'size mismatch'), ('bit-flip', 'SHA256 mismatch'),
             ('missing-file', 'missing/nonregular'), ('extra-file', 'unlisted or absent'),
             ('file-symlink', 'symlink path'), ('mode-change', 'mode mismatch'),
             ('duplicate-row', 'unique sorted'), ('wrong-source-pin', 'research source pin'),
             ('absolute-path', 'unsafe path'), ('parent-traversal', 'unsafe path'),
             ('wrong-schema', 'manifest schema'), ('unlisted-replay', 'unlisted replay')]
    if native.get('evidence_bindings'):
        cases.append(('wrong-evidence-binding', 'native/review evidence mismatch'))
    results = []
    for case, error in cases:
        with tempfile.TemporaryDirectory(prefix='session72-corrupt-') as temp:
            packet = Path(temp) / 'packet'
            shutil.copytree(source, packet, ignore=shutil.ignore_patterns('__pycache__'))
            mp = packet / 'FILE_MANIFEST.json'
            manifest = json.loads(mp.read_text())
            file = packet / chosen
            if case == 'append-byte':
                file.write_bytes(file.read_bytes() + b'\n')
            elif case == 'bit-flip':
                data = bytearray(file.read_bytes())
                need(data, 'Corruption target cannot be empty')
                data[0] ^= 1
                file.write_bytes(data)
            elif case == 'missing-file':
                file.unlink()
            elif case == 'extra-file':
                (packet / 'UNLISTED-CORRUPTION').write_text('unexpected\n')
            elif case == 'file-symlink':
                file.unlink()
                file.symlink_to(packet / 'FILE_MANIFEST.json')
            elif case == 'mode-change':
                file.chmod(file.stat().st_mode ^ 0o100)
            elif case == 'duplicate-row':
                manifest['files'].append(manifest['files'][0])
            elif case == 'wrong-source-pin':
                manifest['source_pin'] = '0' * 40
            elif case in ('absolute-path', 'parent-traversal'):
                manifest['files'][0]['path'] = '/not-allowed' if case == 'absolute-path' else '../not-allowed'
                manifest['files'].sort(key=lambda r: r['path'])
            elif case == 'wrong-schema':
                manifest['schema'] = 'unknown-schema'
            elif case == 'unlisted-replay':
                manifest['replays'] = [{'script': 'not-listed.py', 'normal': chosen, 'optimized': chosen}]
            elif case == 'wrong-evidence-binding':
                manifest['evidence_bindings'][0]['sha256'] = '0' * 64
            mp.write_text(json.dumps(manifest, indent=2) + '\n')
            for optimized in (False, True):
                flags = ['-I', '-B'] + (['-O'] if optimized else [])
                proc = subprocess.run([sys.executable, *flags, str(packet / 'verify_packet.py')],
                                      cwd=temp, capture_output=True)
                if case == 'baseline':
                    need(proc.returncode == 0, 'Baseline verifier failed: ' + proc.stderr.decode(errors='replace'))
                else:
                    need(proc.returncode != 0 and not proc.stdout and error.encode() in proc.stderr,
                         'Corruption missed intended gate: ' + case + '\n' + proc.stderr.decode(errors='replace'))
                results.append({'case': case, 'optimized': optimized,
                                'expected': 'PASS' if case == 'baseline' else 'REJECT', 'observed': 'PASS'})
    print(json.dumps({'status': 'PASS', 'cases_including_baseline': len(cases),
                      'normal_and_optimized_checks': len(results), 'checks': results}, indent=2, sort_keys=True))

if __name__ == '__main__':
    main()

"""Seal exact artifact bytes and inherited Git sources; no network or mutations by default."""
from pathlib import Path
import argparse
import hashlib
import json
import os
import subprocess

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
BASE = '9beb1921510af39d30ec91d31090ad90afdadb78'
INCOMING = '34b61472cc'
LATE = '3c1180551a'
SOURCES = {
    '05-knowledge/results/session112_counting_interfaces_20261030/arithmetic/MODULAR_OBSERVERS.md': BASE,
    '05-knowledge/results/scaled_count_fft_carriers_20261030/SOURCE_AUDIT.md': BASE,
    '05-knowledge/results/scaled_count_fft_carriers_20261030/source_receipt.json': BASE,
    '05-knowledge/results/scaled_count_fft_carriers_20261030/THEORY.md': BASE,
    '05-knowledge/results/scaled_count_fft_carriers_20261030/HANDOFF.md': BASE,
    '05-knowledge/results/regular_parent_symmetry_20261029/verify.py': BASE,
    '05-knowledge/results/session117_ramsey_circuit_interfaces_20261030/circuits/POINT_DEFECT_CIRCUITS.md': BASE,
    '05-knowledge/results/session117_ramsey_circuit_interfaces_20261030/challenge_audit/POTENTIAL_FRONTIER.md': INCOMING,
    '05-knowledge/results/session117_ramsey_circuit_interfaces_20261030/affine_ramsey/AFFINE_RAMSEY_GATE_GAP.md': INCOMING,
    '05-knowledge/results/session117_ramsey_circuit_interfaces_20261030/SOURCE_MODELS.md': INCOMING,
    '05-knowledge/results/ramsey_circuit_interfaces_20261030/CIRCUITS.md': INCOMING,
    '05-knowledge/results/terminal_fees_phase_covers_20261030/LRC.md': LATE,
}


def sha(data):
    return hashlib.sha256(data).hexdigest()


def git(*args):
    return subprocess.check_output(['git', *args], cwd=ROOT)


def lp(path):
    p = str(path.resolve())
    return Path('\\\\?\\'+p) if os.name == 'nt' and not p.startswith('\\\\?\\') else Path(p)


def build():
    artifacts = {}
    for p in sorted(HERE.rglob('*')):
        if not p.is_file() or p.name == 'manifest.json' or '__pycache__' in p.parts:
            continue
        data = p.read_bytes()
        if b'\r' in data:
            raise RuntimeError('Artifact must have LF bytes: '+p.name)
        artifacts[p.relative_to(HERE).as_posix()] = sha(data)
    sources = {path: {'revision': git('rev-parse', rev).decode().strip(),
                      'sha256': sha(git('show', f'{rev}:{path}'))}
               for path, rev in SOURCES.items()}
    return {'sources': sources, 'artifacts': artifacts,
            'basis': 'Exact Git source blobs and exact LF packet bytes.',
            'native_runtime': 'Python standard library; inherited native 81/27 helper.',
            'formal_scope': 'Pinned external EngineGroup imports plus packet-local root ReducedEngine.lean.',
            'status': 'VERIFIED formal semantics / FINITE-EXACT controls; new TM2 time bound OPEN.'}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--write', action='store_true')
    ap.add_argument('--staged', action='store_true')
    ap.add_argument('--allow-incoming-missing', action='store_true')
    args = ap.parse_args()
    result = build()
    target = HERE/'manifest.json'
    if args.write:
        target.write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8', newline='\n')
    elif json.loads(target.read_text()) != result:
        raise RuntimeError('Manifest mismatch')
    for path, entry in result['sources'].items():
        local = lp(ROOT/path)
        if not local.exists() and args.allow_incoming_missing and SOURCES[path] == INCOMING:
            continue
        if sha(local.read_bytes().replace(b'\r\n', b'\n')) != entry['sha256']:
            raise RuntimeError('Current dependency differs from pin: '+path)
    if args.staged:
        for name, digest in result['artifacts'].items():
            rel = (HERE/name).relative_to(ROOT).as_posix()
            if sha(git('show', ':'+rel)) != digest:
                raise RuntimeError('Staged artifact differs: '+name)
        if git('show', ':'+target.relative_to(ROOT).as_posix()) != target.read_bytes():
            raise RuntimeError('Staged manifest differs')
    print(json.dumps({'sources': len(SOURCES), 'artifacts': len(result['artifacts']),
                      'write': args.write, 'staged': args.staged,
                      'allow_incoming_missing': args.allow_incoming_missing}))


if __name__ == '__main__':
    main()

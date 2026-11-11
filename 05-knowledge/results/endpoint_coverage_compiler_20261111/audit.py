"""Replay all bounded claims in normal and optimized Python; freeze hashes."""
from pathlib import Path
import hashlib
import json
import subprocess
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
STEPS = [
    (['ray_probe.py'],'ray_probe.json'),
    (['inflation_bank.py'],'inflation_bank.json'),
    (['inflation_bank.py','--max-n','7'],'inflation_bank_7.json'),
    (['composition_probe.py'],'composition_probe.json'),
    (['joint_ports.py'],'joint_ports.json'),
    (['endpoint_checks.py'],'endpoint_checks.json'),
    (['wheel.py'],'wheel.json'),
    (['verify.py'],'verification.json'),
]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    checks = []
    for args,output in STEPS:
        result = subprocess.run([sys.executable]+args,cwd=HERE,capture_output=True,text=True)
        if result.returncode:
            raise RuntimeError(f'{args}: {result.stdout}\n{result.stderr}')
        normal = (HERE/output).read_bytes()
        result = subprocess.run([sys.executable,'-O']+args,cwd=HERE,capture_output=True,text=True)
        if result.returncode or (HERE/output).read_bytes() != normal:
            raise RuntimeError(f'Optimized mismatch: {args}\n{result.stdout}\n{result.stderr}')
        checks.append(dict(command='python3 '+' '.join(args),normal='PASS',optimized='PASS',byte_match=True))
        print('PASS normal/-O: '+' '.join(args),flush=True)
    prior = HERE.parent / 'modular_centers_collatz_pell_20261111'
    pins = {str((prior/name).relative_to(ROOT)).replace('\\','/'):sha(prior/name)
            for name in ('content_checks.py','tournament_checks.py')}
    expected = {'content_checks.py':'bf894b45e3d3f29a5e953b72bc12fec6ba3a07db7868a768758624917eacc63d',
                'tournament_checks.py':'25a193fe1d66a6c6780451f5de6c393398c1b7788e743744317c258f41aa6a7a'}
    for name,digest in expected.items():
        if sha(prior/name) != digest:
            raise RuntimeError(f'Borrowed source changed: {name}')
    report = dict(scope='same-agent review, independent counting implementation',
                  python=sys.version.split()[0],inherited_commit='e753f4be7660',
                  borrowed_source_sha256=pins,checks=checks)
    (HERE/'audit.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf8')
    manifest = {p.name:sha(p) for p in sorted(HERE.iterdir())
                if p.is_file() and p.name != 'manifest.json'}
    (HERE/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf8')


if __name__ == '__main__':
    main()

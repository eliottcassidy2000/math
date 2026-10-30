"""Check the overlay against the pinned external project; save a build receipt."""
from pathlib import Path
import argparse
import hashlib
import json
import re
import subprocess

HERE = Path(__file__).resolve().parent
PIN = 'c91cc3f4947caca8977e81a5d14244f73b5ae2b1'
MATHLIB = '56a005855ea17e0b0958f73ec801ae140b0a1679'


def sha(data):
    return hashlib.sha256(data).hexdigest()


def run(argv, cwd):
    return subprocess.run(argv, cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                          encoding='utf-8', errors='strict', check=True).stdout


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--external', type=Path, required=True)
    args = ap.parse_args()
    repo = args.external.resolve()
    if run(['git', 'rev-parse', 'HEAD'], repo).strip() != PIN:
        raise RuntimeError('Wrong external revision')
    if run(['git', 'status', '--porcelain'], repo).strip():
        raise RuntimeError('External checkout must be clean')
    if run(['git', 'rev-parse', 'HEAD'], repo/'.lake/packages/mathlib').strip() != MATHLIB:
        raise RuntimeError('Wrong Mathlib revision')
    source = HERE/'ReducedEngine.lean'
    output = run(['lake', 'env', 'lean', str(source)], repo).replace('\r\n', '\n')
    expected = re.findall(r'^#print axioms (.+)$', source.read_text(encoding='utf-8'), re.M)
    checked = re.findall(r"'([^']+)' depends on axioms: \[([^]]+)\]", output)
    if [x[0] for x in checked] != expected or any(
        set(x[1].split(', ')) != {'propext', 'Classical.choice', 'Quot.sound'} for x in checked):
        raise RuntimeError('Unexpected theorem/axiom audit: '+output)
    (HERE/'lean_output.txt').write_text(output, encoding='utf-8', newline='\n')
    receipt = {'external_commit': PIN, 'mathlib_commit': MATHLIB,
               'toolchain': (repo/'lean-toolchain').read_text().strip(),
               'lean_version': run(['lake', 'env', 'lean', '--version'], repo).strip(),
               'command': 'lake env lean <this-packet>/ReducedEngine.lean',
               'working_directory': 'clean pinned external checkout with cached imports',
               'exit_code': 0, 'source_sha256': sha(source.read_bytes()),
               'output_sha256': sha(output.encode()), 'checked_theorems': expected,
               'axioms': ['propext', 'Classical.choice', 'Quot.sound'],
               'scope': 'Packet root module imports external EngineGroup. No revised TM2 compiler or time theorem.',
               'dependencies': {p: sha((repo/p).read_bytes().replace(b'\r\n', b'\n')) for p in [
                   'FactoringAux/Level.lean', 'FactoringAux/EngineGrid.lean',
                   'FactoringAux/EngineGroup.lean', 'FactoringAux/EngineTop.lean',
                   'FactoringAux/UniParams.lean', 'FactoringAux/Transform.lean',
                   'FactoringAux/Cluster.lean', 'FactoringAux/MulNet.lean']}}
    (HERE/'lean_receipt.json').write_text(json.dumps(receipt, indent=2)+'\n', encoding='utf-8', newline='\n')
    print(output, end='')


if __name__ == '__main__':
    main()

#!/usr/bin/env python3
"""Adversarial source/input-boundary tests, not the prediction entrypoint."""
import hashlib
import importlib.util
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def main():
    here = Path(__file__).resolve().parent
    spec = importlib.util.spec_from_file_location('blind_predictor_for_guard_tests', here / 'predictor.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    source = (here / 'predictor.py').read_text()
    controls = []
    for label, suffix in [
            ('reject nonstandard generation import', '\nimport frozen_growth_input_only\n'),
            ('reject child canonicalizer call', '\ncanonical(())\n'),
            ('reject deletion lookup call', '\nlookup(())\n'),
            ('reject dynamic execution', '\nexec("pass")\n'),
            ('reject target literal as input', '\nANSWER = 6880\n'),
            ('reject trace literal as input', '\nTRACE = 55040\n')]:
        try:
            module.enforce_source_boundary(source + suffix)
        except RuntimeError:
            controls.append(label)
        else:
            raise RuntimeError('hostile source was accepted: ' + label)
    baseline = (here / 'receipt.json').read_bytes()
    with tempfile.TemporaryDirectory(prefix='blind-refinement-hostiles-') as directory:
        root = Path(directory)
        shutil.copyfile(here / 'predictor.py', root / 'predictor.py')
        shutil.copyfile(here / 'parents-through-seven.json', root / 'parents-through-seven.json')
        (root / 'receipt.json').write_text('{"fake_answer": -1}\n')
        (root / 'transfer.py').write_text('raise RuntimeError("forbidden generator imported")\n')
        (root / 'secret.txt').write_text('Forbidden external information\n')
        command = [sys.executable, str(root / 'predictor.py'), '--receipt', str(root / 'fresh.json')]
        run = subprocess.run(command, cwd=root, capture_output=True, check=False)
        require(run.returncode == 0, run.stderr.decode())
        require((root / 'fresh.json').read_bytes() == baseline, 'poisoned controls cannot influence prediction')
        controls.append('target receipts and generation code irrelevant to blind result')
        corrupt = json.loads((root / 'parents-through-seven.json').read_text())
        corrupt['levels'][-1]['parents'].pop()
        (root / 'corrupt-parents.json').write_text(json.dumps(corrupt))
        run = subprocess.run(command + ['--parents', str(root / 'corrupt-parents.json')],
                             cwd=root, capture_output=True, check=False)
        require(run.returncode != 0 and b'frozen input atlas hash' in run.stderr, 'tampered atlas rejected')
        controls.append('reject incomplete or altered atlas by pin')
        run = subprocess.run(command + ['--target-order', '9'], cwd=root, capture_output=True, check=False)
        require(run.returncode != 0, 'unsupported higher-order correction rejected')
        controls.append('reject unsupported target order nine')
        probe = '\n'.join([
            'import predictor', 'from pathlib import Path',
            'p = Path.cwd()',
            'predictor.restrict_io(p / "predictor.py", p / "parents-through-seven.json", p / "unused.json")',
            'try:', '    (p / "secret.txt").read_bytes()',
            'except RuntimeError as error:',
            '    print(str(error))', 'else:', '    raise RuntimeError("unpermitted read succeeded")'])
        run = subprocess.run([sys.executable, '-c', probe], cwd=root, capture_output=True, check=False)
        require(run.returncode == 0 and b'blind prediction file access boundary' in run.stdout,
                'runtime unexpected-file-read guard enforced')
        controls.append('reject runtime unpermitted file read')
    print(json.dumps({'status': 'PASS', 'controls': controls,
                      'baseline_receipt_sha256': hashlib.sha256(baseline).hexdigest()}, sort_keys=True))


if __name__ == '__main__':
    main()

#!/usr/bin/env python3
from pathlib import Path
from hashlib import sha256
from tempfile import TemporaryDirectory
import json
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parent


def main():
    manifest = json.loads((ROOT/'MANIFEST.json').read_text())
    expected = (ROOT/'check.normal.out').read_bytes()
    if expected != (ROOT/'check.optimized.out').read_bytes():
        raise RuntimeError('saved normal and optimized reports differ')
    commands = []
    with TemporaryDirectory(prefix='session50-boundaries-isolated-') as directory:
        target = Path(directory)/'packet'
        target.mkdir()
        for name in [r['path'] for r in manifest['files']]+['MANIFEST.json']:
            destination = target/name
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(ROOT/name, destination)
        for args in (['check.py'], ['-O', 'check.py'], ['verify_manifest.py']):
            run = subprocess.run([sys.executable, *args], cwd=target, capture_output=True, check=True)
            if args[-1] == 'check.py' and run.stdout != expected:
                raise RuntimeError('copied-directory replay differs')
            commands.append({'command':'python '+' '.join(args), 'returncode':run.returncode,
                             'stdout_sha256':sha256(run.stdout).hexdigest()})
    print(json.dumps({'status':'PASS', 'method':'Copied only manifest-listed files plus MANIFEST.json into a new temporary directory; reran both modes and manifest verification there.',
                      'normal_optimized_outputs_identical':True, 'commands':commands},indent=2,sort_keys=True))


if __name__ == '__main__':
    main()

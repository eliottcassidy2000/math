#!/usr/bin/env python3
"""Audit normal/-O parity, true relocation, and fail-closed input checking."""
from pathlib import Path
import hashlib
import json
import shutil
import subprocess
import sys
import tempfile

HERE=Path(__file__).resolve().parent

def sha(data):
    return hashlib.sha256(data).hexdigest()

def run(script, out, optimized=False):
    command=[sys.executable]+(['-O'] if optimized else [])+[str(script),'--receipt',str(out)]
    result=subprocess.run(command,cwd=script.parent,capture_output=True)
    if result.returncode:
        raise RuntimeError(result.stderr.decode())
    if result.stdout != out.read_bytes():
        raise RuntimeError('stdout does not equal receipt')
    return result.stdout

def main():
    normal=run(HERE/'check.py',HERE/'receipt.json')
    optimized=run(HERE/'check.py',HERE/'receipt-optimized.json',True)
    if normal != optimized:
        raise RuntimeError('normal/-O mismatch')
    with tempfile.TemporaryDirectory(prefix='fifteen-kernels-realcopy-') as tmp:
        dest=Path(tmp)/'packet'
        dest.mkdir()
        for name in ('check.py','source-manifest.json'):
            shutil.copy2(HERE/name,dest/name)
        shutil.copytree(HERE/'inputs',dest/'inputs')
        if any(p.is_symlink() for p in dest.rglob('*')):
            raise RuntimeError('relocation must be an actual copy')
        copied=run(dest/'check.py',dest/'receipt.json')
        copied_optimized=run(dest/'check.py',dest/'receipt-optimized.json',True)
        if not normal == copied == copied_optimized:
            raise RuntimeError('relocated output mismatch')
        source=dest/'inputs'/'bolza-kernel-action.json'
        source.write_text(source.read_text()+'\n')
        tamper=subprocess.run([sys.executable,'-O',str(dest/'check.py')],cwd=dest,capture_output=True)
        if tamper.returncode == 0 or b'input digest mismatch' not in tamper.stderr:
            raise RuntimeError('tampered input was not rejected')
    result={'status':'PASS', 'normal_optimized_byte_identical':True,
            'real_copy_normal_optimized_byte_identical':True,
            'real_copy_contained_only_declared_checker_manifest_and_inputs':True,
            'no_symlinks_in_real_copy':True,
            'tampered_input_rejected_under_optimized':True,
            'receipt_sha256':sha(normal),
            'checker_sha256':sha((HERE/'check.py').read_bytes()),
            'audit_sha256':sha(Path(__file__).read_bytes()),
            'active_gates':json.loads(normal)['checks']}
    text=json.dumps(result,indent=2,sort_keys=True)+'\n'
    (HERE/'audit-receipt.json').write_text(text)
    print(text,end='')

if __name__=='__main__':
    main()

"""Relocate only manifest-declared files, rerun both modes, and test rejection."""
from pathlib import Path
import hashlib
import json
import shutil
import subprocess
import sys
import tempfile

def run():
    src = Path(__file__).resolve().parent
    manifest = json.loads((src / 'MANIFEST.json').read_text())
    receipt = {'schema': 's7-portability-receipt-v1', 'status': 'PASS', 'copy_policy': 'MANIFEST.json and only its declared files', 'regeneration': [], 'rejection_controls': []}
    with tempfile.TemporaryDirectory(prefix='s7-manifest-relocation-') as td:
        dst = Path(td) / 'relocated'
        dst.mkdir()
        for name in ['MANIFEST.json'] + [r['path'] for r in manifest['files']]:
            target = dst / name
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(src / name, target)
        def invoke(program, optimized=False):
            args=[sys.executable, '-B'] + (['-O'] if optimized else []) + [str(dst / program)]
            return subprocess.run(args, cwd=dst, capture_output=True)
        for optimized in (False, True):
            valid=invoke('verify.py', optimized)
            if valid.returncode: raise RuntimeError(valid.stderr.decode())
            result=invoke('audit.py', optimized)
            expected=(dst / ('audit.optimized.out' if optimized else 'audit.out')).read_bytes()
            if result.returncode or result.stdout != expected or result.stderr:
                raise RuntimeError('relocated regeneration differs')
            receipt['regeneration'].append({'mode': '-O' if optimized else 'normal', 'exit_code': result.returncode, 'stdout_sha256': hashlib.sha256(result.stdout).hexdigest(), 'stderr_bytes': len(result.stderr)})
        cases = ['modified_code', 'missing_source', 'undeclared_file', 'duplicate_manifest_path', 'traversal_manifest_path']
        for case in cases:
            manifest_bytes=(dst/'MANIFEST.json').read_bytes()
            target=None; saved=None
            if case=='modified_code':
                target=dst/'audit.py'; saved=target.read_bytes();target.write_bytes(saved+b'\n# tamper\n')
            elif case=='missing_source':
                target=dst/'sources/inherited/04-computation/experiments/cooper_s7_weight1_frobenius_check_20261018.gp';saved=target.read_bytes();target.unlink()
            elif case=='undeclared_file':
                target=dst/'not-declared.txt';target.write_text('unlisted')
            else:
                altered=json.loads(manifest_bytes)
                altered['files'].append(dict(altered['files'][0]))
                if case=='traversal_manifest_path': altered['files'][-1]['path']='../outside.txt'
                (dst/'MANIFEST.json').write_text(json.dumps(altered))
            for optimized in (False, True):
                result=invoke('verify.py',optimized)
                if result.returncode==0 or b'REJECT:' not in result.stderr:
                    raise RuntimeError('rejection control accepted: '+case)
                receipt['rejection_controls'].append({'case':case,'mode':'-O' if optimized else 'normal','exit_code':result.returncode,'rejected':True})
            (dst/'MANIFEST.json').write_bytes(manifest_bytes)
            if target is not None:
                if saved is not None:target.write_bytes(saved)
                else:target.unlink()
        if invoke('verify.py').returncode:raise RuntimeError('restoration failed')
    return receipt
if __name__=='__main__':
    print(json.dumps(run(),sort_keys=True,indent=2))

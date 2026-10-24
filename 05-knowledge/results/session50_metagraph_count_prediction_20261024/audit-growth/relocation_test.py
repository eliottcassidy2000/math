#!/usr/bin/env python3
"""Copy only manifest-listed files; run normal/-O checks and closure hostiles."""
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
ROOT=Path(__file__).resolve().parent

def need(ok,message):
    if not ok: raise RuntimeError(message)

manifest=json.loads((ROOT/'manifest.json').read_text())
with tempfile.TemporaryDirectory(prefix='marked-growth-relocated-') as tmp:
    dest=Path(tmp)/'packet'; dest.mkdir()
    for name in list(manifest['files'])+['manifest.json']:
        target=dest/name; target.parent.mkdir(parents=True,exist_ok=True)
        shutil.copyfile(ROOT/name,target)
    results={}
    for optimized in (False,True):
        args=[sys.executable]+(['-O'] if optimized else [])+[str(dest/'verify.py'),'--deep']
        r=subprocess.run(args,cwd=tmp,text=True,capture_output=True)
        need(r.returncode==0,'actual-copy verification: '+r.stderr)
        results['optimized' if optimized else 'normal']=json.loads(r.stdout)
    extra=dest/'unlisted-hostile.txt'; extra.write_text('unexpected payload\n')
    r=subprocess.run([sys.executable,'-O',str(dest/'verify.py')],cwd=tmp,text=True,capture_output=True)
    need(r.returncode!=0 and 'manifest file closure' in r.stderr,'extra file was not rejected')
    extra.unlink()
    target=dest/'AUDIT.md'; original=target.read_bytes(); target.write_bytes(original+b'\nmutation\n')
    r=subprocess.run([sys.executable,'-O',str(dest/'verify.py')],cwd=tmp,text=True,capture_output=True)
    need(r.returncode!=0 and 'digest mismatch' in r.stderr,'payload mutation was not rejected')
    target.write_bytes(original)
    link=dest/'unexpected-symlink'; link.symlink_to('AUDIT.md')
    r=subprocess.run([sys.executable,'-O',str(dest/'verify.py')],cwd=tmp,text=True,capture_output=True)
    need(r.returncode!=0 and 'symlink' in r.stderr,'symlink was not rejected')
    link.unlink()
print(json.dumps({'status':'PASS','manifest_only_actual_copy':True,'changed_working_directory':True,
                  'normal_and_optimized':results,'unlisted_file_rejected':True,
                  'mutated_payload_rejected':True,'symlink_rejected':True},sort_keys=True,indent=2))

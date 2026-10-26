#!/usr/bin/env python3
"""Hostile controls for the packet seal, run without modifying the real packet."""
import json, shutil, subprocess, sys, tempfile
from pathlib import Path
root=Path(__file__).resolve().parent
rows=[]
for optimized in (False,True):
    for kind in ('changed','extra','missing','symlink'):
        with tempfile.TemporaryDirectory(prefix='session67-packet-hostile-') as tmp:
            p=Path(tmp)/'packet';shutil.copytree(root,p)
            target=p/'README.md'
            if kind=='changed':target.write_bytes(target.read_bytes()+b'corruption\n')
            elif kind=='extra':(p/'unexpected.txt').write_text('corruption\n')
            elif kind=='missing':target.unlink()
            else:
                target.unlink();target.symlink_to(root/'README.md')
            command=[sys.executable]+(['-O'] if optimized else [])+[str(root/'verify_packet.py'),'--root',str(p)]
            result=subprocess.run(command,capture_output=True,text=True)
            if result.returncode==0 or 'RuntimeError' not in result.stderr:
                raise RuntimeError('failed to reject '+kind)
            rows.append({'python_optimized':optimized,'corruption':kind,'rejected':True})
print(json.dumps({'status':'PASS','controls':rows},sort_keys=True))

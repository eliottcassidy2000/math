"""Reproduce this focused reader packet without touching inherited files."""
from pathlib import Path
import hashlib
import json
import subprocess
import sys

HERE=Path(__file__).resolve().parent


def need(ok,message):
    if not ok:
        raise RuntimeError(message)


def main():
    outputs=[]
    receipts=[]
    generated=['witnesses.json','tail.json','sidecars.json','ENTRY_DEFECTS.json','receipt.json']
    for optimized in (False,True):
        command=[sys.executable]+(['-O'] if optimized else [])
        output=b''
        for name in ('export_reader.py','probe_tail.py','check_reader.py'):
            result=subprocess.run(command+[str(HERE/name)],capture_output=True,check=True)
            need(not result.stderr,'unexpected checker stderr')
            output+=result.stdout.replace(b'\r\n',b'\n')
        outputs.append(output)
        receipts.append([(HERE/name).read_bytes() for name in generated])
    need(outputs[0]==outputs[1] and receipts[0]==receipts[1],'normal/-O mismatch')
    for name,data in zip(('normal.out','optimized.out'),outputs):
        (HERE/name).write_bytes(data)
    names=['PROOF.md','export_reader.py','probe_tail.py','check_reader.py','endpoint_profile.cpp',
           'reproduce.py','normal.out','optimized.out']+generated
    pins={name:hashlib.sha256((HERE/name).read_bytes().replace(b'\r\n',b'\n')).hexdigest()
          for name in names}
    (HERE/'manifest.json').write_text(json.dumps(dict(normal_optimized_identical=True,files=pins,
        scope='32 finite arithmetic targets,20 actual retained-name targets,40 fresh named graph profiles; no all-t or regenerated-entry claim.'),indent=2)+'\n',encoding='utf-8',newline='\n')
    print('PASS: focused reader packet normal/-O bytes agree; all owned files pinned')


if __name__=='__main__':
    main()

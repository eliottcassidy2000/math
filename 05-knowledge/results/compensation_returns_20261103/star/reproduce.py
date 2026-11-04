"""Reproduce only the exact fixed-star packet; never rewrite old inputs."""
from pathlib import Path
import hashlib
import json
import subprocess
import sys
HERE=Path(__file__).resolve().parent
def need(ok,message):
    if not ok:raise RuntimeError(message)
def digest(path):return hashlib.sha256(path.read_bytes().replace(b'\r\n',b'\n')).hexdigest()
def main():
    outputs=[];receipts=[]
    artifacts=['census.json','typed_profiles.json','nearest.json','controls.json','receipt.json']
    for optimized in (False,True):
        r=subprocess.run([sys.executable]+(['-O'] if optimized else [])+[str(HERE/'check_star.py')],capture_output=True,check=True)
        need(not r.stderr,'unexpected reproduction stderr')
        outputs.append(r.stdout.replace(b'\r\n',b'\n'))
        receipts.append([(HERE/p).read_bytes() for p in artifacts])
    need(outputs[0]==outputs[1] and receipts[0]==receipts[1],'normal/-O output mismatch')
    for name,value in zip(['normal.out','optimized.out'],outputs):(HERE/name).write_bytes(value)
    files=['PROOF.md','check_star.py','seed_fields.cpp','reproduce.py','normal.out','optimized.out']+artifacts
    summary=dict(normal_optimized_identical=True,files={f:digest(HERE/f) for f in files},
        inherited_dependencies=json.loads((HERE/'receipt.json').read_text())['dependencies'])
    (HERE/'manifest.json').write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8',newline='\n')
    print('PASS: complete fixed-star packet normal/-O byte-identical')
if __name__=='__main__':main()

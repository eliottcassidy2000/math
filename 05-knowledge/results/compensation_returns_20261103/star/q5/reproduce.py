"""Pin only the bounded q5 packet and verify prior packet pins unchanged."""
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
    artifacts=['census.json','witnesses.json','controls.json','receipt.json'];outputs=[];data=[]
    parents=[HERE.parent,HERE.parent/'q3',HERE.parent/'q3'/'initial']
    pins=[(d,json.loads((d/'manifest.json').read_text())['files']) for d in parents]
    need(all(digest(d/f)==h for d,files in pins for f,h in files.items()),'prior packet pins before q5')
    for optimized in (False,True):
        r=subprocess.run([sys.executable]+(['-O'] if optimized else [])+[str(HERE/'check_q5.py')],capture_output=True,check=True)
        need(not r.stderr,'unexpected q5 reproduction stderr')
        outputs.append(r.stdout.replace(b'\r\n',b'\n'));data.append([(HERE/f).read_bytes() for f in artifacts])
    need(outputs[0]==outputs[1] and data[0]==data[1],'normal/-O q5 mismatch')
    need(all(digest(d/f)==h for d,files in pins for f,h in files.items()),'prior packet pins after q5')
    for name,value in zip(['normal.out','optimized.out'],outputs):(HERE/name).write_bytes(value)
    files=['PROOF.md','check_q5.py','reproduce.py','normal.out','optimized.out']+artifacts
    result=dict(normal_optimized_identical=True,prior_packets_unchanged=True,files={f:digest(HERE/f) for f in files},
        inherited_dependencies=json.loads((HERE/'receipt.json').read_text())['dependencies'])
    (HERE/'manifest.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
    print('PASS: bounded q5 complete star census normal/-O identical; prior packet pins unchanged')
if __name__=='__main__':main()

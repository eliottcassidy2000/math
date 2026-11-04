"""Reproduce and pin the final initial-cap packet without modifying parents."""
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
    artifacts=['census.json','witnesses.json','receipt.json'];outputs=[];data=[]
    pins=json.loads((HERE.parent/'manifest.json').read_text())['files']
    need(all(digest(HERE.parent/f)==h for f,h in pins.items()),'parent q3 pins before run')
    for optimized in (False,True):
        r=subprocess.run([sys.executable]+(['-O'] if optimized else [])+[str(HERE/'check_initial.py')],capture_output=True,check=True)
        need(not r.stderr,'unexpected initial reproduction stderr')
        outputs.append(r.stdout.replace(b'\r\n',b'\n'));data.append([(HERE/f).read_bytes() for f in artifacts])
    need(outputs[0]==outputs[1] and data[0]==data[1],'normal/-O initial mismatch')
    need(all(digest(HERE.parent/f)==h for f,h in pins.items()),'parent q3 pins after run')
    for name,value in zip(['normal.out','optimized.out'],outputs):(HERE/name).write_bytes(value)
    files=['PROOF.md','check_initial.py','reproduce.py','normal.out','optimized.out']+artifacts
    result=dict(normal_optimized_identical=True,parent_packet_unchanged=True,files={f:digest(HERE/f) for f in files},
        inherited_dependencies=json.loads((HERE/'receipt.json').read_text())['dependencies'])
    (HERE/'manifest.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
    print('PASS: initial-cap complete star census normal/-O identical; parent pins unchanged')
if __name__=='__main__':main()

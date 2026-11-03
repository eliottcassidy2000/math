"""Replay the complete finite calculations normally and with -O."""
from pathlib import Path
import hashlib,json,subprocess,sys

HERE=Path(__file__).resolve().parent
FILES=['birth.json','growth.json','LEDGER.md','path56.json']
def main():
    receipt={};outputs={};before=None
    for flags,label in [([], 'normal'),(['-O'], 'optimized')]:
        for name in ['birth.py','growth.py']:
            r=subprocess.run([sys.executable,'-B',*flags,str(HERE/name)],
                             capture_output=True,text=True,check=True)
            outputs[label+':'+name]=r.stdout.strip()
        hashes={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in FILES}
        receipt[label]=hashes
        if before is not None and hashes!=before:raise RuntimeError('artifact disagreement')
        before=hashes
    for name in ['birth.py','growth.py']:
        if outputs['normal:'+name]!=outputs['optimized:'+name]:
            raise RuntimeError('stdout disagreement')
    receipt['normal_optimized_equal']=True;receipt['stdout']=outputs
    (HERE/'checks.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n',encoding='utf-8')
    print(json.dumps({'normal_optimized_equal':True,'artifacts':before},sort_keys=True))
if __name__=='__main__':main()

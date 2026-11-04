#!/usr/bin/env python3
from pathlib import Path
import json,hashlib,subprocess,tempfile,sys,shutil
P=Path(__file__).resolve().parent

def need(ok,msg):
    if not ok:raise RuntimeError(msg)
def main():
    roots=[Path(sys.argv[1]).resolve()if len(sys.argv)>1 else P.parent/'session189-odd-transitive',Path(sys.argv[2]).resolve()if len(sys.argv)>2 else P.parent/'session189-collision-next'];pins={};outputs={}
    for label,src in zip(['odd','anchor'],roots):
        manifest=json.loads((src/'RELEASE_FILES.json').read_text());pins[label]={f:hashlib.sha256((src/f).read_bytes()).hexdigest()for f in manifest};need(all(pins[label][f]==v['sha256']for f,v in manifest.items()),'frozen allowlist '+label)
        with tempfile.TemporaryDirectory(prefix='odd-anchor-audit-')as td:
            D=Path(td)
            for f in manifest:shutil.copyfile(src/f,D/f)
            for opt,mode in [([], 'normal'),(['-O'],'optimized')]:
                if label=='odd':
                    subprocess.check_output([sys.executable,*opt,str(D/'prepare.py')]);out=subprocess.check_output([sys.executable,*opt,str(D/'verify.py')]);need(out==(src/'receipt.json').read_bytes(),'odd producer replay');q=subprocess.check_output([sys.executable,*opt,str(D/'query_controls.py')]);need(q==(src/'query_receipt.json').read_bytes(),'actual grouped replay');outputs['query.'+mode+'.json']=q
                else:
                    subprocess.check_output([sys.executable,*opt,str(D/'controls.py'),'prepare']);out=subprocess.check_output([sys.executable,*opt,str(D/'controls.py')]);need(out==(src/'receipt.json').read_bytes(),'anchor producer replay')
                outputs[label+'.'+mode+'.json']=out
                for f in manifest:need(hashlib.sha256((D/f).read_bytes()).hexdigest()==pins[label][f],'regenerated frozen file '+label+'/'+f)
        need(all(hashlib.sha256((src/f).read_bytes()).hexdigest()==h for f,h in pins[label].items()),'source moved during review')
    for f,data in outputs.items():(P/f).write_bytes(data)
    receipt={'status':'PASS','producer_packets':2,'normal_optimized_exact':True,'source_pins':pins,'sources':{'replay_audit.py':hashlib.sha256((P/'replay_audit.py').read_bytes()).hexdigest()},'artifacts':{f:hashlib.sha256((P/f).read_bytes()).hexdigest()for f in sorted(outputs)}}
    (P/'replay_receipt.json').write_text(json.dumps(receipt,sort_keys=True,indent=2)+'\n');print(json.dumps(receipt,sort_keys=True))
if __name__=='__main__':main()

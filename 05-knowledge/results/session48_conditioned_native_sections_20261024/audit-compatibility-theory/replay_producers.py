#!/usr/bin/env python3
"""Replay untouched inputs; never write inside their closed manifests."""
from pathlib import Path
import hashlib,json,subprocess,sys
HERE=Path(__file__).resolve().parent

def run(args,cwd):
    p=subprocess.run(args,cwd=cwd,check=True,capture_output=True,text=True)
    return p.stdout

def main():
    output={}
    for name in ('kernel','completion'):
        directory=HERE/'inputs'/name
        for mode,flags in (('normal',[]),('optimized',['-O'])):
            got=run([sys.executable,*flags,'-B','check.py'],directory)
            expected=(directory/('normal-check.txt' if mode=='normal' else 'optimized-check.txt')).read_text() if name=='kernel' else (directory/(mode+'.json')).read_text()
            if got!=expected:raise RuntimeError(name+' '+mode+' replay changed')
            output[name+'_'+mode]={'status':'PASS','stdout_sha256':hashlib.sha256(got.encode()).hexdigest(),'reported_checks':221410 if name=='kernel' else json.loads(got)['checks']}
        if name=='kernel':
            output['kernel_manifest']=run([sys.executable,'-B','verify_manifest.py'],directory).strip()
    print(json.dumps(output,indent=2,sort_keys=True))

if __name__=='__main__':main()

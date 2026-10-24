#!/usr/bin/env python3
"""Copy only manifest-listed inputs to fresh /tmp and reproduce in isolation."""
from pathlib import Path
import json, shutil, subprocess, sys, tempfile
HERE=Path(__file__).resolve().parent

def main():
    subprocess.run([sys.executable,'-I','-B',str(HERE/'verify_manifest.py')],check=True,cwd='/tmp',capture_output=True,text=True)
    manifest=json.loads((HERE/'MANIFEST.json').read_text())
    with tempfile.TemporaryDirectory(prefix='session48-affine-phases-',dir='/tmp') as tmp:
        root=Path(tmp)
        for name in ['MANIFEST.json']+sorted(manifest['files']):
            dest=root/name;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(HERE/name,dest)
        outputs=[]
        for flags in [('-I','-B'),('-I','-O','-B')]:
            result=subprocess.run([sys.executable,*flags,str(root/'check.py')],check=True,cwd='/tmp',capture_output=True,text=True)
            outputs.append(result.stdout)
        if outputs[0]!=outputs[1]:raise RuntimeError('Normal and optimized isolated checks disagree')
        verified=subprocess.run([sys.executable,'-I','-B',str(root/'verify_manifest.py')],check=True,cwd='/tmp',capture_output=True,text=True)
        print('PASS fresh /tmp closure: copied only manifest-listed files and MANIFEST.json')
        print('PASS normal and optimized isolated runs: identical exact output')
        print(outputs[0],end='')
        print(verified.stdout,end='')

if __name__=='__main__':main()

#!/usr/bin/env python3
"""Ensure malformed producer metadata cannot satisfy independent audit."""
from pathlib import Path
import argparse,copy,hashlib,json,shutil,subprocess,sys,tempfile
H=Path(__file__).resolve().parent

def need(c,m):
    if not c:raise RuntimeError(m)

def run():
    results=[]
    for label,kind in [('engine byte tamper','byte'),('forged internal scalar with refreshed checksum','internal'),('forged blue boundary with refreshed checksum','blue')]:
        with tempfile.TemporaryDirectory(prefix='balanced-hostile-') as d:
            root=Path(d)/'producer';shutil.copytree(H/'producer',root)
            if kind=='byte':
                with (root/'engine.py').open('a') as f:f.write('\n# deliberate hostile mutation\n')
                command=[sys.executable]+(['-O'] if not __debug__ else [])+[str(root/'verify_manifest.py')]
            else:
                receipt=json.loads((root/'normal.json').read_text())
                row=next(x for x in receipt['counts'] if x['n']==10)
                row['internal_ports' if kind=='internal' else 'boundary_blue_ports']+=1
                (root/'normal.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n')
                manifest=json.loads((root/'MANIFEST.json').read_text())
                rec=next(x for x in manifest['files'] if x['path']=='normal.json');rec['sha256']=hashlib.sha256((root/'normal.json').read_bytes()).hexdigest()
                (root/'MANIFEST.json').write_text(json.dumps(manifest,indent=2,sort_keys=True)+'\n')
                command=[sys.executable]+(['-O'] if not __debug__ else [])+[str(H/'compare.py'),'--producer',str(root),'--audit',str(H),'--output',str(Path(d)/'out.json')]
            completed=subprocess.run(command,capture_output=True,text=True)
            need(completed.returncode!=0,'accepted hostile '+label)
            results.append({'test':label,'status':'REJECTED','last_error_line':completed.stderr.strip().splitlines()[-1]})
    return {'status':'PASS','tests':results}
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',default='hostile-packet.json');a=p.parse_args();r=run();Path(a.output).write_text(json.dumps(r,indent=2,sort_keys=True)+'\n');print(json.dumps(r,indent=2,sort_keys=True))

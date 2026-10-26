#!/usr/bin/env python3
"""Exercise active provenance rejection independently of the outer packet hashes."""
import argparse,json,shutil,subprocess,sys,tempfile
from pathlib import Path
def need(ok,why):
    if not ok:raise RuntimeError(why)
def main():
    a=argparse.ArgumentParser();a.add_argument('packet',type=Path);args=a.parse_args()
    source=args.packet.resolve();base=json.loads((source/'PROVENANCE_RULES.json').read_text())
    cases=[('baseline',None),('schema','provenance schema')]
    if base.get('repository_origins'):
        cases += [('origin-sha256','provenance SHA256 mismatch'),('origin-size','provenance size mismatch'),('origin-blob','source Git blob mismatch'),('origin-pin','unapproved source pin'),('origin-path','unsafe provenance path')]
    if base.get('reviewed_bindings'):cases.append(('reviewed-binding','provenance SHA256 mismatch'))
    rows=[]
    for case,error in cases:
        with tempfile.TemporaryDirectory(prefix='session72-provenance-mutation-') as t:
            p=Path(t)/'packet';shutil.copytree(source,p,ignore=shutil.ignore_patterns('__pycache__'))
            path=p/'PROVENANCE_RULES.json';r=json.loads(path.read_text())
            if case=='schema':r['schema']='wrong'
            elif case=='origin-sha256':r['repository_origins'][0]['sha256']='0'*64
            elif case=='origin-size':r['repository_origins'][0]['bytes']+=1
            elif case=='origin-blob':r['repository_origins'][0]['git_blob']='0'*40
            elif case=='origin-pin':r['repository_origins'][0]['pin']='0'*40
            elif case=='origin-path':r['repository_origins'][0]['path']='../escape'
            elif case=='reviewed-binding':r['reviewed_bindings'][0]['sha256']='0'*64
            path.write_text(json.dumps(r))
            for optimized in (False,True):
                q=subprocess.run([sys.executable,'-I','-B',*(['-O'] if optimized else []),str(p/'verify_provenance.py')],capture_output=True,cwd=t)
                if error:need(q.returncode!=0 and not q.stdout and error.encode() in q.stderr,'Wrong rejection: '+case+' '+q.stderr.decode())
                else:need(q.returncode==0,'Baseline failed: '+q.stderr.decode())
                rows.append({'case':case,'optimized':optimized,'expected':'REJECT' if error else 'PASS','observed':'PASS'})
    print(json.dumps({'status':'PASS','normal_and_optimized_checks':len(rows),'checks':rows},indent=2,sort_keys=True))
if __name__=='__main__':main()

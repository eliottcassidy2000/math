#!/usr/bin/env python3
"""Replay independent audit normally/-O, reject explicit corruptions; optionally replay producers."""
import argparse,hashlib,json,pathlib,subprocess,sys
ROOT=pathlib.Path(__file__).resolve().parent

def need(ok,msg):
    if not ok:raise RuntimeError(msg)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def execute(cmd):
    r=subprocess.run(cmd,capture_output=True)
    return r

def main():
    p=argparse.ArgumentParser();p.add_argument('--producers',action='store_true');args=p.parse_args()
    pins=json.loads((ROOT/'SOURCE_PINS.json').read_text())
    for row in pins['files']:
        data=(ROOT/row['snapshot']).read_bytes()
        need(hashlib.sha256(data).hexdigest()==row['sha256'],'inherited SHA256 pin')
        need(hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest()==row['blob'],'inherited Git blob pin')
    rows=[]
    for optimized in (False,True):
        prefix=[sys.executable,'-B']+(['-O'] if optimized else [])
        cmd=prefix+[str(ROOT/'check.py')];r=execute(cmd)
        expected=ROOT/('optimized.json' if optimized else 'normal.json')
        need(r.returncode==0 and r.stdout==expected.read_bytes(),'independent replay mismatch')
        rows.append({'command':cmd,'returncode':r.returncode,'receipt_sha256':sha(expected),'byte_match':True})
        for control in ('field','trajectory','sparse'):
            cmd=prefix+[str(ROOT/'check.py'),'--inject',control];r=execute(cmd)
            expected_error={'field':b'RuntimeError: deliberate field-rank corruption','trajectory':b'RuntimeError: deliberate trajectory conflation','sparse':b'RuntimeError: deliberate pairwise gluing corruption'}[control]
            need(r.returncode!=0 and r.stderr.strip().endswith(expected_error),'corruption passed or failed unexpectedly')
            rows.append({'command':cmd,'returncode':r.returncode,'expected_failure':True,'stderr_sha256':hashlib.sha256(r.stderr).hexdigest()})
    producers=[]
    if args.producers:
        for name,receipt in [('session72-response','normal.json'),('session72-cross-domain','check.normal.json')]:
            directory=ROOT.parent/name
            before={p.name:sha(p) for p in directory.glob('*.md')}
            before['check.py']=sha(directory/'check.py')
            for optimized in (False,True):
                cmd=[sys.executable,'-B']+(['-O'] if optimized else [])+[str(directory/'check.py')]
                r=execute(cmd)
                need(r.returncode==0 and r.stdout==(directory/receipt).read_bytes(),'producer exact replay mismatch '+name)
                producers.append({'name':name,'optimized':optimized,'command':cmd,'returncode':r.returncode,'byte_match':True,'receipt_sha256':sha(directory/receipt),'source_hashes':before})
            need(all(sha(directory/n)==s for n,s in before.items()),'producer source changed during replay')
    print(json.dumps({'status':'PASS','inherited_source_pins_verified':len(pins['files']),'independent_replays':rows,'producer_replays':producers},sort_keys=True,indent=2))
if __name__=='__main__':main()

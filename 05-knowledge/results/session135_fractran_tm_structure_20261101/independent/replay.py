#!/usr/bin/env python3
"""Authenticate inputs and rerun all seven independent audits in both Python modes."""
import hashlib,json,subprocess,sys,argparse
from pathlib import Path
HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
SCRIPTS={'early':'audit_early.py','density':'audit_density.py','square_tm':'audit_square_tm.py','translation':'audit_translation.py','lrc':'audit_lrc.py','acceleration':'audit_acceleration.py','periodic_points':'audit_periodic_points.py'}
def require(ok,msg):
    if not ok:raise RuntimeError(msg)
def main():
    pins=json.loads((HERE/'INPUT_PINS.json').read_text())
    for name,wanted in pins.items():
        require(hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==wanted,('input hash',name))
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--only',nargs='+',choices=list(SCRIPTS));options=parser.parse_args()
    selected=[label for label in SCRIPTS if not options.only or label in options.only]
    records=[]
    for label in selected:
        script=SCRIPTS[label]
        for mode,mode_flags in [('normal',[]),('optimized',['-O'])]:
            result=subprocess.run([sys.executable,'-B']+mode_flags+[str(HERE/script)],cwd=HERE,capture_output=True)
            require(result.returncode==0,(script,mode,result.stderr.decode()))
            expected=(HERE/(label+'.'+mode+'.json')).read_bytes()
            require(result.stdout==expected,('receipt mismatch',script,mode))
            records.append({'audit':label,'mode':mode,'sha256':hashlib.sha256(result.stdout).hexdigest()})
    require(len(records)==2*len(selected),'incomplete selected-audit replay')
    print(json.dumps({'status':'PASS','authenticated_input_files':len(pins),'replays':records},indent=2))
if __name__=='__main__':main()

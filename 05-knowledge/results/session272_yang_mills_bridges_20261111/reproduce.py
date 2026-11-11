#!/usr/bin/env python3
"""Replay new controls only. Does not build/reverify trusted YM or inherited work."""
import argparse,json,subprocess,sys
from pathlib import Path
HERE=Path(__file__).resolve().parent

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--with-uv',action='store_true',help='also run slower mpmath precision comparison')
    args=parser.parse_args()
    records=[]
    for path,output in [('limits_measure/probe.py','limits_measure/probe.out'),('spectral_loops/check.py','spectral_loops/check.out')]:
        proc=subprocess.run([sys.executable,str(HERE/path)],capture_output=True,text=True,check=True)
        expected=(HERE/output).read_text()
        assert proc.stdout==expected,(path,'output changed')
        records.append({'script':path,'status':'PASS','matches_committed_output':True})
    if args.with_uv:
        before=json.loads((HERE/'uv_ir/OUTPUT.json').read_text())
        proc=subprocess.run([sys.executable,str(HERE/'uv_ir/check_uv_ir.py')],capture_output=True,text=True,check=True)
        after=json.loads((HERE/'uv_ir/OUTPUT.json').read_text())
        assert before==after,'UV diagnostic output changed'
        records.append({'script':'uv_ir/check_uv_ir.py','status':'PASS_NUMERICAL_NOT_INTERVAL','matches_committed_output':True})
    else:
        records.append({'script':'uv_ir/check_uv_ir.py','status':'NOT_RUN_IN_THIS_INVOCATION'})
    result={'status':'PASS','scope':'new finite controls and optional UV diagnostic, not source proof verification','checks':records}
    print(json.dumps(result,indent=2))
if __name__=='__main__':main()

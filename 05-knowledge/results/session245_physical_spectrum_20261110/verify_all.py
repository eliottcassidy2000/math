#!/usr/bin/env python3
"""Replay every sealed proof-facing receipt; no external repo dependency."""
import json,pathlib,subprocess,sys
R=pathlib.Path(__file__).resolve().parent
JOBS={'interfaces/verify.py':'interfaces/CHECK.json','interfaces/verify_chain.py':'interfaces/CHAIN_CHECK.json','forest/verify.py':'forest/RECEIPT.json','forest/verify_run_budget.py':'forest/RUN_BUDGET_CHECK.json','descent/verify.py':'descent/CHECK.json','audit/independent.py':'audit/independent.json'}
FAIL={'interfaces/verify.py','forest/verify_run_budget.py','descent/verify.py'}
def need(ok,msg):
 if not ok:raise RuntimeError(msg)
def main():
 rows=[]
 for script,receipt in JOBS.items():
  expected=json.loads((R/receipt).read_text())
  for opt in (False,True):
   cmd=[sys.executable,'-I','-B']+(['-O'] if opt else [])+[str(R/script)]
   p=subprocess.run(cmd,capture_output=True)
   need(p.returncode==0,('execution failed',script,opt,p.stderr.decode()[-2000:]))
   need(json.loads(p.stdout)==expected,('receipt mismatch',script,opt))
   if script in FAIL:
    bad=subprocess.run(cmd+['--fail'],capture_output=True)
    need(bad.returncode!=0,('hostile sentinel accepted',script,opt))
  rows.append({'script':script,'normal':'PASS','optimized':'PASS','hostile':'PASS' if script in FAIL else 'not applicable'})
 print(json.dumps({'status':'PASS','checks':rows,'scope':'exact finite replays; all-target Chan-Pak remains OPEN'},indent=2,sort_keys=True))
if __name__=='__main__':main()

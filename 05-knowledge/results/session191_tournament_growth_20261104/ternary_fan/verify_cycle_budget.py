#!/usr/bin/env python3
"""Regenerate all inputs and certify the sharp second-level cycle budget."""
import hashlib,json,subprocess
from pathlib import Path
ROOT=Path(__file__).resolve().parent
PIN='55409706be4f8c71f5e2b2cd0777c9c0f0f8de97eac4b174d648f50b638dc903'
def need(x,m):
 if not x:raise ValueError(m)
def main():
 generator=ROOT/'generate_bank8.bin';scan=ROOT/'verify_cycle_budget.bin'
 for src,out in [('generate_bank8.cpp',generator),('verify_cycle_budget.cpp',scan)]:subprocess.run(['g++','-std=c++17','-O3',str(ROOT/src),'-o',str(out)],check=True)
 gen=subprocess.run([str(generator)],capture_output=True,check=True);need(hashlib.sha256(gen.stdout).hexdigest()==PIN,'exact regenerated bank')
 bank=json.loads(gen.stdout);need(len(bank)==6880,'bank count');tours=set()
 for r in bank:
  for mask in range(1,256):
   names=[v for v in range(8) if mask>>v&1]
   if len(names)<3:continue
   tours.add(tuple(sum(1<<j for j,v in enumerate(names) if r[u]>>v&1) for u in names))
 ordered=sorted(tours,key=lambda r:(len(r),r));data=''.join(str(len(r))+' '+' '.join(map(str,r))+'\n' for r in ordered)
 result=json.loads(subprocess.run([str(scan)],input=data,text=True,capture_output=True,check=True).stdout)
 need(result['input_rows']==33278 and result['order8_parents']==701 and result['order9_extensions']==179456,'exact universe sizes')
 need(result['cycle_sets']=={'89':[2,3,6],'91':[2,6],'95':[7],'97':[2]},'complete cycle sets')
 need(result['target_H_order9_records']=={'89':18,'91':18,'95':18,'97':75},'nonempty order9 target scan')
 need(5**3>97,'strong order10 cutoff')
 result.update({'status':'PASS','bank_sha256':PIN,'cutoff':{'minimum_at_order10':125,'maximum_target':97},'minimum_uniform_budget':7})
 raw=(json.dumps(result,indent=2,sort_keys=True)+'\n').encode();(ROOT/'cycle_budget_receipt.json').write_bytes(raw)
 print(json.dumps({'status':'PASS','complete_cycle_sets':result['cycle_sets'],'order9_extensions':result['order9_extensions'],'minimum_uniform_budget':7,'receipt_sha256':hashlib.sha256(raw).hexdigest()},sort_keys=True))
if __name__=='__main__':main()

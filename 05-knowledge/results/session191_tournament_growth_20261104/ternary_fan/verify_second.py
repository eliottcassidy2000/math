#!/usr/bin/env python3
"""Regenerate the exact bank and reproduce fixed-core continuation boundaries."""
import collections, hashlib, json, subprocess
from pathlib import Path
from verify_fan import require,profile,brute
ROOT=Path(__file__).resolve().parent
PIN='55409706be4f8c71f5e2b2cd0777c9c0f0f8de97eac4b174d648f50b638dc903'
def oddpart(n):
 while n%2==0:n//=2
 return n

def main():
 generator=ROOT/'generate_bank8.bin'
 subprocess.run(['g++','-std=c++17','-O3',str(ROOT/'generate_bank8.cpp'),'-o',str(generator)],check=True)
 generated=subprocess.run([str(generator)],capture_output=True,check=True)
 data=generated.stdout;require(hashlib.sha256(data).hexdigest()==PIN,'regenerated bank hash')
 bank=json.loads(data);require(len(bank)==6880,'bank class count')
 tours=set()
 for rows in bank:
  for mask in range(1,256):
   names=[i for i in range(8) if mask>>i&1]
   if len(names)<3:continue
   tours.add(tuple(sum(1<<j for j,v in enumerate(names) if rows[u]>>v&1) for u in names))
 ordered=sorted(tours,key=lambda r:(len(r),r));require(len(ordered)==33278,'induced input count')
 text=''.join(str(len(r))+' '+' '.join(map(str,r))+'\n' for r in ordered)
 binary=ROOT/'verify_second.bin';subprocess.run(['g++','-std=c++17','-O2',str(ROOT/'verify_second.cpp'),'-o',str(binary)],check=True)
 result=json.loads(subprocess.run([str(binary)],input=text,text=True,capture_output=True,check=True).stdout)
 require(result['total']==33278,'processed total')
 require(5**8>61**3,'Busch cutoff')
 grouped={h:[] for h in [53,55,59,61]}
 for row in result['target_profiles']:grouped[row['H']].append(row)
 require([len(grouped[h]) for h in grouped]==[39,22,65,0],'target record counts')
 units={h:sum(1 for x in grouped[h] if 1 in x['J']) for h in grouped}
 require(units=={53:2,55:0,59:0,61:0},'unit-port records')
 require(all(oddpart(f)!=63 for f in result['F61']+result['P2_61']),'dyadic obstruction')
 # Independent Python recurrence rechecks every H61 record and every candidate target record.
 checked={}
 for rows in result['rows61']:
  rows=tuple(rows)
  p=profile(rows,0);require(p['H']==61,'H61 independent')
  checked[rows]=True
 for row in result['target_profiles']:
  p=profile(row['rows'],row['root']);require(p['H']==row['H'] and p['F']==63 and p['J']==row['J'],'target independent')
 # Explicit existing second-level witnesses checked additionally by literal permutations.
 witnesses=[(53,5,0,[2,16,3,7,45,79,159,63]),(55,6,None,[2,24,11,33,13,87,31]),(59,5,None,[6,12,32,37,15,83,31])]
 selected=[]
 for h,r,v,rows in witnesses:
  p=profile(rows,r);require(p==brute(rows,r),'independent witness permutations')
  require(p['H']==h and p['F']==63 and p['strong'],'witness')
  if v is not None:require(p['J'][v]==1,'unit witness')
  selected.append({'root':r,'unit_port':v,'rows':rows,**p})
 embedded=[250,244,249,226,8,16,48,112]
 p=profile(embedded,1);require(p['H']==53 and p['F']==63 and p['J'][3]==1,'physical lower continuation')
 from verify_fan import CORES
 require([r&63 for r in embedded[:6]]==CORES['A'],'physical lower induced parent')
 result['physical_lower_continuation']={'rows':embedded,'root':1,'unit_port':3,**p}
 result.update({'bank_sha256':PIN,'cutoff':{'five_power_eight':5**8,'sixtyone_cubed':61**3},'unit_record_counts':units,'witnesses':selected,'status':'PASS'})
 raw=(json.dumps(result,indent=2,sort_keys=True)+'\n').encode();(ROOT/'second_receipt.json').write_bytes(raw)
 print(json.dumps({'status':'PASS','input_rows':len(ordered),'target_records':{h:len(grouped[h]) for h in grouped},'unit_records':units,'F61_count':len(result['F61']),'P2_61':result['P2_61'],'receipt_sha256':hashlib.sha256(raw).hexdigest()},sort_keys=True))
if __name__=='__main__':main()

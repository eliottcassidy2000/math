"""Reproduce and check the finite corroboration of the all-order transfer proof."""
from pathlib import Path
from itertools import product
from collections import Counter
import hashlib,json,subprocess,sys
D=Path(__file__).resolve().parent

def need(ok,why):
 if not ok:raise RuntimeError(why)
def sha(name):return hashlib.sha256((D/name).read_bytes()).hexdigest()
def main():
 subprocess.run(['g++','-std=c++17','-O3','-Wall','-Wextra','-Werror',str(D/'kernels.cpp'),'-o',str(D/'kernels')],check=True)
 raw=subprocess.check_output([str(D/'kernels')]);need(raw==(D/'kernels.txt').read_bytes(),'retained complete exact cut kernels reproduce')
 seen=set();by_order=Counter();scalar_even=0
 for line in raw.decode().splitlines():
  n,code,i,k,*counts=map(int,line.split());key=(n,code,i,k)
  need(key not in seen,'distinct kernel record');seen.add(key);by_order[n]+=1
  need(len(counts)==9,'full cut-kernel vector')
  need(n in range(2,6) and code in range(1<<(n*(n-1)//2)) and i in range(n) and k in (2,3,5),'kernel domain')
  for h,value in enumerate(counts,1):
   if h>=2:need(value%2==int(h=={2:3,3:3,5:5}[k]),'universal higher-cover parity table')
   if h>n+k-1:need(value==0,'complete zero tail')
  if k in (3,5):need(counts[0]%2==0,'scalar return-segment parity');scalar_even+=1
 expected={(n,c,i,k) for n in range(2,6) for c in range(1<<(n*(n-1)//2)) for i in range(n) for k in (2,3,5)}
 need(seen==expected,'every labelled tournament and distinguished label')
 # Recompute actual physical counts, then audit retained rows/masks and arithmetic.
 result=subprocess.check_output([sys.executable,str(D/'physical_checks.py')]);log=json.loads(result)
 from physical_checks import build,valid
 data=json.loads((D/'physical_records.json').read_text());rows=(D/'physical_values.txt').read_text().splitlines();at=0;contrasts=Counter()
 for r in data['records']+data['complete_outside_domain_bank']:
  blocks=r.get('child_rows') or [data['child_blocks'][j] for j in r['child_type_indices']]
  for side,s in enumerate(r['sides']):
   expected_side=build(r['Q'],blocks,*r['block_query'],r.get('ambient',False),side)
   need({k:v for k,v in s.items() if k!='ordered_covers'}==expected_side,'physical labels and actual augmented graph')
   need(valid(s['augmented_rows']),'actual physical tournament')
   values=list(map(int,rows[at].split()));at+=1;need(values==s['ordered_covers'],'literal cover data aligned')
   for h,value in enumerate(values,1):
    import math
    need(value%math.factorial(h)==0,'ordered covers have free path-slot action on physical vertices')
  contrast=[x-y for x,y in zip(r['sides'][0]['ordered_covers'],r['sides'][1]['ordered_covers'])]
  need(contrast==r['contrast'],'exact contrast alignment')
  if not r.get('ambient',False):contrasts[contrast[0]%16]+=1
 need(at==len(rows)==944,'all physical cases retained')
 bad=data['outside_domain_counterexample'];need((bad['contrast'][1]-2*bad['contrast'][0])%16!=0 or any(x%16 for x in bad['contrast'][2:]),'out-of-domain hostile failure')
 files=['PROOF.md','kernels.cpp','kernels.txt','physical.cpp','physical_checks.py','physical_records.json','physical_values.txt']
 out=dict(status='PASS',scope='All-order proof with finite corroboration; no local nonvanishing claim.',kernel_records=len(seen),kernel_records_by_order=dict(sorted(by_order.items())),scalar_return_segment_checks=scalar_even,complete_kernel_coordinate_count=len(seen)*9,physical_control_summary=log,outside_domain_first_counterexample=dict(quotient=bad['Q'],child_sizes=list(map(len,bad['child_rows'])),block_query=bad['block_query'],exact_contrast=bad['contrast'],residues=[v%16 for v in bad['contrast']]),sha256={n:sha(n) for n in files})
 print(json.dumps(out,sort_keys=True,separators=(',',':')))
if __name__=='__main__':main()

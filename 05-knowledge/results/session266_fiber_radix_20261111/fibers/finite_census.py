"""Replay complete inherited ALL-H census through order 8, no enumeration.
Run: python finite_census.py. All generated CSVs use LF serialization.
"""
import csv,hashlib,json,math
from pathlib import Path
P=Path(__file__).parent
D=json.loads((P/'finite_input.json').read_text())
counts={int(n):{int(h):r for h,r in byh.items()} for n,byh in D['counts_by_order'].items()}
checks=0
expected=[1,1,2,4,12,56,456,6880]
max_h=max(h for byh in counts.values() for h in byh)
fields=['n','H','ordinary','self_converse','merged','labeled','cumulative_ordinary','cumulative_self_converse','cumulative_merged','cumulative_labeled']
channels=['ordinary','self_converse','merged','labeled'];totals=[]
with open(P/'ALL_H_THROUGH_8.csv','w',newline='') as f:
 writer=csv.DictWriter(f,fieldnames=fields,lineterminator='\n');writer.writeheader()
 for n in range(1,9):
  census=counts[n]
  if sum(r['ordinary'] for r in census.values())!=expected[n-1]:raise RuntimeError(('ordinary total',n))
  if sum(r['labeled'] for r in census.values())!=2**math.comb(n,2):raise RuntimeError(('labeled total',n))
  checks+=2
  total={key:sum(r[key] for r in census.values()) for key in channels}
  total.update(n=n,cumulative_ordinary=sum(expected[:n]));totals.append(total)
  for h in range(1,max_h+1,2):
   r=census.get(h,{key:0 for key in channels})
   if 2*r['merged']!=r['ordinary']+r['self_converse']:raise RuntimeError(('converse',n,h))
   if h in (7,21) and any(r.values()):raise RuntimeError(('gap',n,h))
   row=dict(n=n,H=h,**r)
   for key in channels:row['cumulative_'+key]=sum(counts[m].get(h,{}).get(key,0) for m in range(1,n+1))
   writer.writerow(row);checks+=2
receipt={'scope':'Every Hamiltonian-path count, every tournament isomorphism class, sizes1..8. Even H and odd H above max_h are uniformly zero in this size range. No census was rerun.','max_h':max_h,'active_checks':checks,'totals':totals,'source':D['source'],'csv_sha256':hashlib.sha256((P/'ALL_H_THROUGH_8.csv').read_bytes()).hexdigest()}
(P/'finite_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(receipt))

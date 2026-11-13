"""Exact newborn-cut kernels; inherited labelled row recovery and a scalar-state hostile."""
from collections import Counter,defaultdict
from io import StringIO
from math import factorial
from pathlib import Path
import csv,json,subprocess
from check import BANK,SOURCE_REF,load,tables,opposite,need
HERE=Path(__file__).resolve().parent
CSV='05-knowledge/results/spectrum_multiplicity_signed_tables_20261111/multiplicity/spectrum_n1_n8.csv'

def kernel(T):
 n=len(T);full=(1<<n)-1;hp,end=tables(T);_,start=tables(opposite(T))
 W=[[sum(end[m][i]*start[full^m][j] for m in range(1,full)) for j in range(n)] for i in range(n)]
 return hp[-1],start[-1],end[-1],W

def main():
 rows=[];groups=defaultdict(list);labelled=Counter();checks=0
 for rec in load(BANK)['records']:
  if rec['n']>5:continue
  T=tuple(rec['adjacency_outmasks'])
  for flip in range(rec['ordinary_multiplicity']):
   U=T if flip==0 else opposite(T);n=len(U);h,S,E,W=kernel(U);hist=Counter()
   for x in range(1<<n):
    Q=tuple(row|((1<<n) if not(x>>i&1) else 0) for i,row in enumerate(U))+(x,)
    actual=tables(Q)[0][-1]
    predicted=h+sum((S[i]-E[i]) for i in range(n) if x>>i&1)+sum(W[i][j] for i in range(n) if not(x>>i&1) for j in range(n) if x>>j&1)
    need(predicted==actual,'quadratic cut kernel');hist[actual]+=1;checks+=1
   row=dict(n=n,H=h,aut=rec['aut'],rows=U,histogram=dict(sorted(hist.items())))
   rows.append(row);groups[n,h].append(row)
   if n==5:
    for q,count in hist.items():labelled[q]+=factorial(n)//rec['aut']*count
 csvdata=list(csv.DictReader(StringIO(subprocess.check_output(['git','show',SOURCE_REF+':'+CSV],text=True))))
 target={int(r['H']):int(r['labelled_tournaments']) for r in csvdata if r['n']=='6' and int(r['labelled_tournaments'])}
 need(dict(labelled)==target,'exact inherited labelled row atn6')
 need(sum(labelled.values())==2**15,'complete labelled total')
 need(sum(h*c for h,c in labelled.items())==factorial(6)*2**10,'Hamilton incidence first moment')
 hostiles=[(n,h,rs) for (n,h),rs in sorted(groups.items()) if len({str(r['histogram']) for r in rs})>1]
 need(hostiles[0][:2]==(5,3),'first scalar-source collision in complete orders3..5')
 end=(26,28,25,16,0);prepared=(18,20,17,23,0);x=9;preparation=[]
 for T,expected in ((end,17),(prepared,29)):
  U=tuple(row|((1<<5) if not(x>>i&1) else 0) for i,row in enumerate(T))+(x,)
  need(tables(T)[0][-1]==3 and tables(U)[0][-1]==expected,'literal samecut preparation')
  preparation.append(dict(source=T,cut=x,child=U,H_child=expected))
 out=dict(status='PASS',source_ref=SOURCE_REF,preparation=preparation,scope='All ordinary inherited tournament classes atorders3..5, all newborn cuts; exactweightedn5ton6labelledflow',
  quadratic_controls=checks,ordinary_source_classes=len(rows),first_hostile=dict(n=5,H=3,sources=hostiles[0][2]),labelled_n6=dict(sorted(labelled.items())))
 (HERE/'birth_flow_receipt.json').write_text(json.dumps(out,indent=2)+'\n')
 print(json.dumps({'status':'PASS','quadratic_controls':checks,'ordinary_source_classes':len(rows),'first_hostile':out['first_hostile']},indent=2))

if __name__=='__main__':main()

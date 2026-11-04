#!/usr/bin/env python3
from pathlib import Path
from itertools import permutations
import json,hashlib
from basic import *
P=Path(__file__).resolve().parent

def cycle_count(r):
 n=len(r);D=[[0]*n for _ in range(1<<n)];D[1][0]=1
 for mask in range(1,1<<n):
  if not mask&1:continue
  for i,x in enumerate(D[mask]):
   if not x:continue
   opts=r[i]&~mask
   while opts:
    b=opts&-opts;opts-=b;D[mask|b][b.bit_length()-1]+=x
 return sum(D[-1][i]for i in range(1,n)if r[i]&1)
def literal_cycles(r):
 return[(0,)+p for p in permutations(range(1,len(r)))if all(r[a]>>b&1 for a,b in zip((0,)+p,p+(0,)))]
def one_exit(r,t):
 n=len(r);out=list(r)+[1<<t]
 for i in range(n):
  if i!=t:out[i]|=1<<n
 return tuple(out)
def main():
 pins=json.load(open(P/'INCOMING.json'))
 for f,x in pins.items():
  b=(P/f).read_bytes();need(hashlib.sha256(b).hexdigest()==x['sha256'],'pinned input hash');need(hashlib.sha1(b'blob '+str(len(b)).encode()+b'\0'+b).hexdigest()==x['git_blob'],'pinned Git blob')
 T=(128,9,3,85,39,15,567,4478,639,1215,511,2047,3967);need(valid(T)and strong(T)and T[0]==1<<7,'new strong named witness');U,ids=delete(T,0);H=sum(endpoints(T));digit=endpoints(U)[ids.index(7)];need((H,digit)==(12875,2),'actual new H and digit');Y,_=phi(T,0);need(strong(U)and strong(Y)and cycle_count(Y)==2,'correct forward marked image')
 child=C2(T,0);need(strong(child)and sum(endpoints(child))==38629,'actual strong +4 child')
 pad=tuple(sum(1<<j for j in range(i+1,15))for i in range(2))+tuple(x<<2 for x in T);pu,pids=delete(pad,2);need(valid(pad)and not strong(pad)and sum(endpoints(pad))==12875 and endpoints(pu)[pids.index(9)]==2,'weak order15 padding boundary')
 RT=(6,12,24,17,3);cycs=literal_cycles(RT);need(cycle_count(RT)==len(cycs)==2 and all(x.bit_count()==2 for x in RT),'regular five exception');edges=lambda p:set(zip(p,p[1:]+p[:1]));diff=len(edges(cycs[0])-edges(cycs[1]));need(diff==5,'regular five cycles not three-edge swap')
 core=(2,4,9,3);rev=tuple(((1<<4)-1)^x^(1<<i)for i,x in enumerate(core));starts=endpoints(rev);r=next(i for i,s in enumerate(starts)if s==2);S=one_exit(core,r);a=4;need(strong(S)and cycle_count(S)==2,'five-name sole-exit two-cycle class');cyc2=literal_cycles(S);d2=len(edges(cyc2[0])-edges(cyc2[1]));need(d2==3,'three-edge difference boundary example');ready,root=phi(S,a);du,di=delete(ready,a);need(strong(ready)and endpoints(du)[di.index(root)]==2,'canonical admissible inverse');admissible=[i for i,x in enumerate(S)if x.bit_count()==1 and strong(delete(S,i)[0])];need(admissible==[a],'unique admissible mark')
 wrong=(2,4,9,3);wronginv,_=phi(wrong,1);need(strong(wrong)and cycle_count(wrong)==1 and not strong(wronginv),'earlier marked counterexample persists')
 control={'new_12875':{'rows':T,'H':H,'strong':True,'a':0,'root':7,'digit':digit,'Phi_rows':Y,'Phi_cycles':2,'deleted_core_strong':True,'C2_rows':child,'C2_H':38629,'C2_strong':True,'order15_prefix_rows':pad,'order15_prefix_strong':False},'RT5':{'rows':RT,'cycles':cycs,'one_sided_arc_difference':diff,'sole_exit_count':0},'sole_exit_order5':{'rows':S,'cycles':cyc2,'one_sided_arc_difference':d2,'admissible_exit':a,'ready_inverse_rows':ready,'ready_inverse_H':sum(endpoints(ready))},'older_marked_counterexample':{'Y':wrong,'mark':1,'inverse':wronginv,'inverse_strong':False}}
 (P/'controls.json').write_text(json.dumps(control,sort_keys=True,indent=2)+'\n');receipt={'status':'PASS','scope':'Targeted source comparison and explicit witnesses only; no rerun of the incoming order9 cycle census, higher ready generators, or all-order segment-swap/count theorem.','source_sha256':{f:hashlib.sha256((P/f).read_bytes()).hexdigest()for f in ['verify.py','basic.py','INCOMING.json']},'controls_sha256':hashlib.sha256((P/'controls.json').read_bytes()).hexdigest()};(P/'receipt.json').write_text(json.dumps(receipt,sort_keys=True,indent=2)+'\n');print(json.dumps(receipt,sort_keys=True))
if __name__=='__main__':main()

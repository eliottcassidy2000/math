#!/usr/bin/env python3
from pathlib import Path
from itertools import combinations
import json,hashlib,subprocess,tempfile
P=Path(__file__).resolve().parent

def need(b,s):
 if not b:raise RuntimeError(s)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def literal(a):
 n=len(a);full=(1<<n)-1;H=C=0
 def go(first,last,S):
  nonlocal H,C
  if S==full:
   H+=1
   if first==0 and a[last]&1:C+=1
   return
  q=a[last]&~S
  while q:
   bit=q&-q;q-=bit;go(first,bit.bit_length()-1,S|bit)
 for v in range(n):go(v,v,1<<v)
 return H,C

def triangle(a,V):
 for x,y,z in combinations(sorted(V),3):
  if a[x]>>y&1 and a[y]>>z&1 and a[z]>>x&1:return[x,y,z]
  if a[x]>>z&1 and a[z]>>y&1 and a[y]>>x&1:return[x,z,y]
 return None

def transitive_order(a,V):
 V=set(V);out=[]
 while V:
  x=next((x for x in V if all(x==y or a[x]>>y&1 for y in V)),None)
  if x is None:return None
  out.append(x);V.remove(x)
 return out

def main():
 expected=json.loads((P/'c3_complete_certificate.json').read_text())
 need(sha(P/'verify_c3_complete.cpp')=='e49f0b9ff6e29f7b7260d6516ccb0fafe03ea4bbcf06dccc0664b7fddead806f','frozen complete census source')
 need(5**12>533**3,'strong-order cutoff')
 with tempfile.TemporaryDirectory(prefix='feedback533-')as d:
  d=Path(d);census=d/'census';filter=d/'filter';csv=d/'all.csv'
  subprocess.run(['g++','-O2','-std=c++17',str(P/'verify_c3_complete.cpp'),'-o',str(census)],check=True)
  gen=subprocess.run([str(census),str(csv)],check=True,text=True,capture_output=True)
  need(not gen.stderr and sha(csv)==expected['all_records_sha256'],'complete canonical bank')
  subprocess.run(['g++','-O2','-std=c++17',str(P/'verify_feedback533.cpp'),'-o',str(filter)],check=True)
  selected=subprocess.run([str(filter),str(csv)],check=True,text=True,capture_output=True)
  need(int(selected.stderr)==12077,'all C1/C2 bank graphs recounted')
  records=[list(map(int,line.split(',')))for line in selected.stdout.splitlines()]
  need(len(records)==6 and all(r[:2]==[11,2]for r in records),'exact six H533 cores')
  controls=[]
  for n,c,*a in records:
   need(literal(a)==(533,2),'independent literal H/C')
   lower=[]
   for deleted in combinations(range(n),2):
    tri=triangle(a,set(range(n))-set(deleted));need(tri is not None,'feedback pair exists')
    need(all(a[x]>>y&1 for x,y in zip(tri,tri[1:]+tri[:1])),'triangle certificate')
    lower.append(dict(deleted=list(deleted),cycle=tri))
   upper=None
   for deleted in combinations(range(n),3):
    order=transitive_order(a,set(range(n))-set(deleted))
    if order is not None:upper=dict(deleted=list(deleted),transitive_order=order);break
   need(upper is not None,'feedback three witness')
   controls.append(dict(rows=a,H=533,C=2,feedback_number=3,pair_deletion_cycle_witnesses=lower,three_deletion_transitive_witness=upper))
 birth=[8,9,10,7,6,5,4,3,2,1,0]
 moves=[[1,1],[1,1],[2,1],[3,0],[2,0],[2,0],[2,0],[2,0],[2,0]]
 a=controls[0]['rows'];r=[2,0];ports=[0,1]
 for mask,drop in moves:
  n=len(r);r.append(0)
  for u in range(n):
   if u in ports and mask>>ports.index(u)&1:r[n]|=1<<u
   else:r[u]|=1<<n
  ports.append(n);del ports[drop]
 expected_rows=[sum(1<<j for j,v in enumerate(birth)if a[u]>>v&1)for u in birth]
 need(r==expected_rows and len(ports)==2,'literal changing-port construction')
 changing=dict(birth_order=birth,moves=moves,final_ports=[birth[i]for i in ports])
 out=dict(status='PASS',changing_two_port_witness=changing,strong_H=533,maximum_strong_order=12,cutoff={'H_cubed':533**3,'five_to_twelve':5**12},complete_bank_sha256=expected['all_records_sha256'],recounted_C1_C2_graphs=12077,C1_hits=0,C2_hits=6,core_order=11,all_feedback_numbers=3,cores=controls,dependencies={name:sha(P/name)for name in ['C3_COMPLETE.md','verify_c3_complete.cpp','c3_complete_certificate.json','PROOF.md']})
 text=json.dumps(out,sort_keys=True,indent=2)+'\n';(P/'feedback533_certificate.json').write_text(text);print(text,end='')
if __name__=='__main__':main()

from pathlib import Path
from itertools import product
import json,subprocess
D=Path(__file__).resolve().parent

def need(c,m):
 if not c:raise RuntimeError(m)
def op(Q):return [((1<<len(Q))-1)^(1<<i)^r for i,r in enumerate(Q)]
def tt(n):return [sum(1<<j for j in range(i+1,n)) for i in range(n)]
BLOCKS=[[0,1,3,7],[0,9,3,5],[2,4,9,3],[2,12,9,1],tt(6)]
# Each list is checked as an actual tournament, not accepted by its intended type.
def valid(T):return all(not(T[i]>>i&1) and all(((T[i]>>j)&1)+((T[j]>>i)&1)==1 for j in range(len(T)) if j!=i) for i in range(len(T)))
def build(Q,blocks,a,c,ambient,side):
 target=[op(B) if side and i%2==0 else B for i,B in enumerate(blocks)];R=op(Q) if side else Q;T=[];sets=[]
 for B in target:sets.append(list(range(len(T),len(T)+len(B))));T.extend([0]*len(B))
 for i,B in enumerate(target):
  off=sets[i][0]
  for j,row in enumerate(B):T[off+j]=(row<<off)|sum(1<<v for k,S in enumerate(sets) if R[i]>>k&1 for v in S)
 aa=sum(1<<v for i,S in enumerate(sets) if a>>i&1 for v in S);cc=sum(1<<v for i,S in enumerate(sets) if c>>i&1 for v in S);exterior=[]
 if ambient:
  n=len(T);T=[row|(1<<n) for row in T]+[1<<(n+1),(1<<n)-1];aa|=1<<(n+1);cc|=1<<(n+1);exterior=[n,n+1]
 n=len(T);T=[row|(((aa>>i)&1)<<n)|(((cc>>i)&1)<<(n+1)) for i,row in enumerate(T)]+[(((1<<n)-1)^aa)|(1<<(n+1)),((1<<n)-1)^cc]
 return dict(augmented_rows=T,physical_query=[aa,cc],child_sets=sets,port_names=[n,n+1],exterior_names=exterior,target_child_rows=target)

def main():
 need(all(valid(B) for B in BLOCKS),'actual child tournaments');queries3=[(0,0),(1,0),(0,1),(1,2),(1,3),(3,5),(5,6),(7,7)];configs=[]
 for ids in list(product(range(4),repeat=2))+[(4,0),(1,4),(4,4)]:configs.append(([2,0],ids,list(product(range(4),repeat=2))))
 for Q in [[6,4,0],[2,4,1]]:
  for ids in [(0,0,0),(1,2,3),(3,1,2),(4,0,1),(2,4,0),(1,3,4)]:configs.append((Q,ids,queries3))
 records=[];cases=[]
 for Q,ids,queries in configs:
  for a,c in queries:
   record=dict(Q=Q,child_type_indices=ids,block_query=[a,c],ambient=False,sides=[])
   for side in range(2):p=build(Q,[BLOCKS[i] for i in ids],a,c,False,side);record['sides'].append(p);cases.append(p['augmented_rows'])
   records.append(record)
 # A selected returned-root ambient witnesses the odd multiplier5, not a scalar-only claim.
 ambient_base=[0,1,6,20,250,304,320,336]
 for base in ambient_base:
  original=records[base];record=dict(Q=original['Q'],child_type_indices=original['child_type_indices'],block_query=original['block_query'],ambient=True,local_record=base,sides=[])
  for side in range(2):p=build(record['Q'],[BLOCKS[i] for i in record['child_type_indices']],*record['block_query'],True,side);record['sides'].append(p);cases.append(p['augmented_rows'])
  records.append(record)
 hostile=[]
 for a,c in product(range(8),repeat=2):
  record=dict(Q=[2,4,1],child_rows=[tt(2),tt(2),tt(3)],block_query=[a,c],sides=[])
  for side in range(2):p=build(record['Q'],record['child_rows'],a,c,False,side);record['sides'].append(p);cases.append(p['augmented_rows'])
  hostile.append(record)
 need(max(map(len,cases))<=18,'exact integer bound18')
 subprocess.run(['g++','-std=c++17','-O3','-Wall','-Wextra','-Werror',str(D/'physical.cpp'),'-o',str(D/'physical')],check=True);inp=str(len(cases))+'\n'+''.join(str(len(T))+'\n'+' '.join(map(str,T))+'\n' for T in cases);raw=subprocess.check_output([str(D/'physical')],input=inp.encode()).decode();values=[list(map(int,l.split())) for l in raw.splitlines()];need(len(values)==len(cases),'all literal ordered-cover counts');(D/'physical_values.txt').write_text(raw)
 at=0
 for r in records+hostile:
  for s in r['sides']:s['ordered_covers']=values[at];at+=1
  r['contrast']=[x-y for x,y in zip(r['sides'][0]['ordered_covers'],r['sides'][1]['ordered_covers'])]
 for r in records:
  delta=r['contrast']
  if not r['ambient']:need((delta[1]-2*delta[0])%16==0 and all(x%16==0 for x in delta[2:]),'full all-even-child contrast law')
  else:need((delta[0]-5*records[r['local_record']]['contrast'][0])%16==0,'actual cyclic ambient multiplier5')
 bad=[r for r in hostile if (r['contrast'][1]-2*r['contrast'][0])%16 or any(x%16 for x in r['contrast'][2:])];need(bad,'outside-domain boundary fails the asserted law')
 data=dict(status='PASS',child_blocks=BLOCKS,records=records,complete_outside_domain_bank=hostile,outside_domain_counterexample=bad[0],physical_cases=len(cases),maximum_augmented_order=max(map(len,cases)));(D/'physical_records.json').write_text(json.dumps(data,separators=(',',':'))+'\n');print(json.dumps(dict(status='PASS',all_even_local_queries=sum(not r['ambient'] for r in records),cyclic_ambient_queries=sum(r['ambient'] for r in records),hostile_queries=len(hostile),physical_cases=len(cases),maximum_augmented_order=max(map(len,cases)))))
if __name__=='__main__':main()

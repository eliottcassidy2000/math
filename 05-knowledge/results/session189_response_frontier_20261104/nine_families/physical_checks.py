from pathlib import Path
import json,subprocess
D=Path(__file__).resolve().parent

def need(ok,msg):
 if not ok:raise RuntimeError(msg)
def op(Q):return [((1<<len(Q))-1)^(1<<i)^r for i,r in enumerate(Q)]
def tt(n):return [sum(1<<j for j in range(i+1,n)) for i in range(n)]
def build(Q,blocks,a,c,side,ambient=False):
 R=op(Q) if side else Q;Bs=[op(B) if side and i%2==0 else B for i,B in enumerate(blocks)];sets=[];T=[]
 for B in Bs:sets.append(list(range(len(T),len(T)+len(B))));T.extend([0]*len(B))
 for i,B in enumerate(Bs):
  for j,row in enumerate(B):T[sets[i][j]]=(row<<sets[i][0])|sum(1<<v for k,S in enumerate(sets) if R[i]>>k&1 for v in S)
 aa=sum(1<<v for i,S in enumerate(sets) if a>>i&1 for v in S);cc=sum(1<<v for i,S in enumerate(sets) if c>>i&1 for v in S);ext=[]
 if ambient:
  n=len(T);T=[r|(1<<n) for r in T]+[1<<(n+1),(1<<n)-1];aa|=1<<(n+1);cc|=1<<(n+1);ext=[n,n+1]
 n=len(T);T=[r|(((aa>>i)&1)<<n)|(((cc>>i)&1)<<(n+1)) for i,r in enumerate(T)]+[(((1<<n)-1)^aa)|(1<<(n+1)),((1<<n)-1)^cc]
 return dict(augmented_rows=T,physical_query=[aa,cc],child_sets=sets,target_child_rows=Bs,port_names=[n,n+1],exterior_names=ext)
def main():
 families={q['name']:q for q in json.loads((D/'families.json').read_text())};records=[];cases=[]
 configs=[]
 for name in ('critical_U9','critical_W9'):
  for V in (0,511,1,128,129,73,438):
   blocks=[tt(2) if V>>i&1 else tt(1) for i in range(9)]
   if V==73:blocks[3]=[2,4,1] # actual nontransitive child, still nonsingleton
   if name=='critical_U9':queries=[(16,0,4,2)]
   else:queries=[(0,256,4,2*((V&1)^((V>>7)&1))),(4,0,8,4*(1^(V&1)^((V>>7)&1)))]
   for a,c,mod,residue in queries:configs.append((name,V,blocks,a,c,mod,residue,False))
 for name in ('circulant_00','circulant_01'):
  for V in (0,1,2,3,73,255,438):
   blocks=[tt(2) if V>>i&1 else tt(1) for i in range(9)]
   if V==3:blocks[0]=[2,4,1]
   if V==0:a,c,mod,residue=1,2,4,2
   elif name=='circulant_00':
    i=next(i for i in range(9) if ((V>>((i+4)%9))^(V>>((i+5)%9)))&1);a,c,mod,residue=0,1<<i,4,2
   else:
    i=next(i for i in range(9) if ((V>>((i-1)%9))^(V>>((i+1)%9)))&1);a,c,mod,residue=1<<i,0,8,4
   configs.append((name,V,blocks,a,c,mod,residue,False))
 for idx in [0,7,21,28]:
  old=configs[idx];configs.append(old[:-1]+(True,))
 for name,V,blocks,a,c,mod,residue,ambient in configs:
  r=dict(name=name,quotient_rows=families[name]['rows'],nontrivial_child_mask=V,source_child_rows=blocks,block_query=[a,c],modulus=mod,expected_residue=residue,ambient=ambient,sides=[])
  for side in range(2):
   s=build(r['quotient_rows'],blocks,a,c,side,ambient);r['sides'].append(s);cases.append(s['augmented_rows'])
  records.append(r)
 need(max(map(len,cases))<=20,'Hamilton integer bound20')
 subprocess.run(['g++','-std=c++17','-O3','-Wall','-Wextra','-Werror',str(D/'physical_h.cpp'),'-o',str(D/'physical_h')],check=True)
 inp=str(len(cases))+'\n'+''.join(str(len(T))+'\n'+' '.join(map(str,T))+'\n' for T in cases);raw=subprocess.check_output([str(D/'physical_h')],input=inp.encode()).decode();values=list(map(int,raw.split()));need(len(values)==len(cases),'every literal physical side');(D/'physical_values.txt').write_text(raw)
 for j,r in enumerate(records):
  for side in range(2):r['sides'][side]['H']=values[2*j+side]
  r['contrast']=values[2*j]-values[2*j+1];need(r['contrast']%r['modulus']==r['expected_residue'],'physical affine or singleton query including fixed cyclic ambient')
 (D/'physical_records.json').write_text(json.dumps(records,separators=(',',':'))+'\n');print(json.dumps(dict(status='PASS',physical_query_pairs=len(records),physical_sides=len(cases),maximum_augmented_order=max(map(len,cases)))))
if __name__=='__main__':main()

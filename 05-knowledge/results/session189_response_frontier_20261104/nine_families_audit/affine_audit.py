"""Independent bounded-family construction, modules, path completion and affine anchors."""
from pathlib import Path
from functools import lru_cache
import json,hashlib
D=Path(__file__).resolve().parent;P=D.parent/'nine_families'
if not P.exists():P=D.parent/'session189-nine-families'

def need(c,s):
 if not c:raise RuntimeError(s)
def op(R):return [((1<<len(R))-1)^(1<<i)^r for i,r in enumerate(R)]
def mods(R):
 n=len(R);full=(1<<n)-1;out=[]
 for S in range(1,full):
  if S.bit_count()<2:continue
  i=(S&-S).bit_length()-1
  if all(not(S>>j&1) or ((R[i]^R[j])&(full^S))==0 for j in range(n)):out.append(S)
 return out

def endpoints(R):
 n=len(R);full=(1<<n)-1;J=[[0]*n for _ in R]
 for first in range(n):
  def walk(S,last):
   if not S:J[first][last]+=1;return
   choices=S&R[last]
   while choices:
    b=choices&-choices;choices-=b;walk(S^b,b.bit_length()-1)
  walk(full^(1<<first),first)
 return J

def starts(R,repeated=-1):
 n=len(R);full=(1<<n)-1
 if repeated<0:
  @lru_cache(None)
  def suffix(S,last):
   if not S:return 1
   out=0;choices=S&R[last]
   while choices:
    b=choices&-choices;choices-=b;out+=suffix(S^b,b.bit_length()-1)
   return out
  return [suffix(full^(1<<i),i) for i in range(n)]
 other=full^(1<<repeated)
 @lru_cache(None)
 def suffix(S,left,last):
  if not S and not left:return 1
  out=0
  if left and R[last]>>repeated&1:out+=suffix(S,left-1,repeated)
  choices=S&R[last]
  while choices:
   b=choices&-choices;choices-=b;out+=suffix(S^b,left,b.bit_length()-1)
  return out
 return [suffix(other,1,i) if i==repeated else suffix(other^(1<<i),2,i) for i in range(n)]

def oneport(R,c):
 n=len(R);return [r|(((c>>i)&1)<<n) for i,r in enumerate(R)]+[((1<<n)-1)^c]
def rank(rows):
 basis={}
 for row in rows:
  while row:
   k=row.bit_length()-1
   if k not in basis:basis[k]=row;break
   row^=basis[k]
 return len(basis)
def main():
 data=json.loads((P/'families.json').read_text());need(len(data)==18,'exact18-model bounded family')
 expected=[]
 for code in range(16):
  steps=[j if code>>(j-1)&1 else 9-j for j in range(1,5)];expected.append([sum(1<<((i+j)%9) for j in steps) for i in range(9)])
 U=expected[15].copy()
 for i in range(4):
  for j in range(i+1,4):U[i]^=1<<j;U[j]^=1<<i
 W=[sum(1<<j for j in range(i+1,8))+((i%2)<<8) for i in range(8)]+[85];expected.extend([U,W]);records=[];primitive_count=0
 for index,(q,R) in enumerate(zip(data,expected)):
  need(q['rows']==R,'exact labeled construction');need(all(not(R[i]>>i&1) and all(((R[i]>>j)&1)+((R[j]>>i)&1)==1 for j in range(9) if i!=j) for i in range(9)),'valid quotient tournament');M=mods(R);need(M==q['nontrivial_modules'] and bool(not M)==q['prime'],'independent complete module test')
  J=endpoints(R);S=list(map(sum,J));E=[sum(row[j] for row in J) for j in range(9)];p=[(s+e)%2 for s,e in zip(S,E)];need(J==q['joint_endpoints'] and S==q['S'] and E==q['E'] and p==q['p'],'independent literal Hamilton endpoint enumeration')
  B=op(R);rsk=[[0]*9 for _ in range(9)];base=[];kmat=[]
  for j in range(9):
   a=starts(R,j);b=starts(B,j)
   for i in range(9):rsk[i][j]=a[i]-b[i]
  for i in range(9):
   A=oneport(R,1<<i);C=oneport(B,1<<i);base.append(sum(starts(A))-sum(starts(C)));kmat.append([sum(starts(A,j))-sum(starts(C,j)) for j in range(9)])
  need(rsk==q['start_run_skews'] and base==q['one_port_base'] and kmat==q['one_port_double'],'all independent exact run primitives and orientation signs');primitive_count+=171
  patterns=[];blind=[]
  for V in range(512):
   d=[(base[i]+2*sum(kmat[i][j] for j in range(9) if V>>j&1))%4 for i in range(9)];t=[(S[i]-E[i]+2*sum(rsk[i][j] for j in range(9) if V>>j&1))%4 for i in range(9)];a=[(2*d[i]+2*t[i])%8 if not p[i] else None for i in range(9)];row=dict(nontrivial_child_mask=V,first=d,second=a);need(row==q['patterns'][V],'every singleton-position anchor vector')
   if not any(p) and not any(d) and not any(a):blind.append(V)
   if index==0:need(all(d[i]==2*(((V>>((i+4)%9))^(V>>((i+5)%9)))&1) for i in range(9)),'first connected-cycle identities')
   if index==1:need(all(a[i]==4*(((V>>((i-1)%9))^(V>>((i+1)%9)))&1) for i in range(9)),'second connected-cycle identities')
   if index==17:need(d[8]==2*((V&1)^((V>>7)&1)) and a[2]==4*(1^(V&1)^((V>>7)&1)),'W9 two fixed bits sum to one')
   patterns.append(row)
  need(blind==q['anchor_blind_patterns'],'exact residual list');first=[sum((x%2)<<j for j,x in enumerate(row)) for row in kmat];second=[sum(((x+y)%2)<<j for j,(x,y) in enumerate(zip(a,b))) for a,b in zip(kmat,rsk)];records.append(dict(index=index,name=q['name'],rows=R,proper_modules=M,joint_endpoints=J,S=S,E=E,p=p,start_run_skews=rsk,one_port_base=base,one_port_double=kmat,patterns=patterns,anchor_blind_patterns=blind,first_rank=rank(first),second_rank=rank(second),combined_rank=rank(first+second)))
 need(records[16]['p']==[0,0,0,0,1,0,0,0,1],'U9 endpoint certificate');need(sum(records[0]['S'])==3267 and sum(records[1]['S'])==3357,'distinct representatives by Hamilton count')
 transports=[];classes=[[0,3,5,10,12,15],[1,4,7,8,11,14]]
 for inds in classes:
  A=expected[inds[0]]
  for dest in inds:
   candidates=[]
   for u in [1,2,4,5,7,8]:
    perm=[i*u%9 for i in range(9)]
    if all(((A[i]>>j)&1)==((expected[dest][perm[i]]>>perm[j])&1) for i in range(9) for j in range(9)):candidates.append(perm)
   need(candidates,'every claimed unit transport');perm=candidates[0];inverse=[perm.index(i) for i in range(9)];masks=[sum(1<<perm[j] for j in range(9) if mask>>j&1) for mask in range(512)];need(len(set(masks))==512,'complete physical-mask transport bijective');transports.append(dict(source=inds[0],target=dest,source_to_target=perm,target_to_source=inverse,all_mask_transports=masks))
 for j in range(16):
  prime=j in classes[0]+classes[1];need((not records[j]['proper_modules'])==prime,'exact12 prime circulants')
  need(records[j]['anchor_blind_patterns']==([0,511] if prime else [0,73,146,219,292,365,438,511]),'all prime/composite residual-domain patterns')
  if not prime:need(records[j]['proper_modules']==[73,146,292],'exact three composite physical modules')
 outputs={'affine_recount.json':records,'physical_transports.json':transports}
 for name,obj in outputs.items():
  file=D/name;txt=json.dumps(obj,separators=(',',':'))+'\n'
  if file.exists():need(file.read_text()==txt,'independent data stable')
  else:file.write_text(txt)
 return dict(status='PASS',bounded_models=18,complete_joint_endpoint_coordinates=18*81,independent_exact_anchor_primitives=primitive_count,singleton_position_patterns=18*512,prime_circulants=12,unit_transport_records=len(transports),physical_mask_transports=512*len(transports),source_sha256=hashlib.sha256((D/'affine_audit.py').read_bytes()).hexdigest(),output_pins={name:hashlib.sha256((D/name).read_bytes()).hexdigest() for name in outputs},producer_proof_sha256=hashlib.sha256((P/'PROOF.md').read_bytes()).hexdigest())
if __name__=='__main__':print(json.dumps(main(),sort_keys=True,indent=2))

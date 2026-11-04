"""Modulo-eight profile and ambient transfer. Standard library only."""
from pathlib import Path
import importlib.util,itertools,json,random,hashlib
ROOT=Path(__file__).resolve().parent
s=importlib.util.spec_from_file_location('covered',ROOT.parent/'covered_response'/'verify.py');v=importlib.util.module_from_spec(s);s.loader.exec_module(v)
def hc(a):
 n=len(a)
 if n<3:return 0
 D=[[0]*n for _ in range(1<<n)];D[1][0]=1
 for mask in range(1,1<<n):
  if not mask&1:continue
  for u in range(n):
   count=D[mask][u]
   if not count:continue
   dest=a[u]&~mask
   while dest:
    bit=dest&-dest;dest-=bit;D[mask|bit][bit.bit_length()-1]+=count
 return sum(D[-1][u] for u in range(n) if a[u]&1)
def stirling(n,k):
 a=[1]+[0]*k
 for _ in range(n):a=[0]+[j*a[j]+a[j-1] for j in range(1,k+1)]
 return a[k]
def check(a):
 n=len(a);q=v.profile(a)+[0,0,0];H=q[1];z=hc(a) if n%2 else 0
 v.need((q[2]-2*H-(1<<n)+4+4*z)%8==0,'q2 exact mod8 signature')
 expected=6*stirling(n,3)+(2*H-2-4*z if n%2 else 0)
 v.need((q[3]-expected)%8==0,'q3 exact mod8 signature');v.need(all(x%8==0 for x in q[4:]),'factorial higher cover divisibility')
 return {'order':n,'rows':a,'q':q[:n+1],'hamilton_cycles_if_odd':z}
def substitute(Q,r,A):
 outside=[u for u in range(len(Q)) if u!=r];n=len(A);rows=list(A)+[0]*len(outside)
 for i,u in enumerate(outside):
  for j,w in enumerate(outside):
   if Q[u]>>w&1:rows[n+i]|=1<<(n+j)
  if Q[u]>>r&1:rows[n+i]|=(1<<n)-1
  else:
   for j in range(n):rows[j]|=1<<(n+i)
 return rows
records=[]
for n in range(1,6):
 for c in range(1<<(n*(n-1)//2)):records.append(check(v.tournament(n,c)))
rng=random.Random(188800)
for n,count in [(6,512),(7,32),(8,32),(9,8)]:
 for _ in range(count):records.append(check(v.tournament(n,rng.randrange(1<<(n*(n-1)//2)))))
A=v.tt(4);three=next(v.tournament(4,c) for c in range(64) if v.hp(v.tournament(4,c))==3);five=[8,1,3,6];pairs=[(three,A),(five,A)];ambient=[]
for n in range(1,5):
 for code in range(1<<(n*(n-1)//2)):
  Q=v.tournament(n,code)
  for r in range(n):
   factor=v.hp(substitute(Q,r,[2,0]));v.need(factor%2==1,'odd ambient multiplier')
   for B,A in pairs:
    delta=v.hp(B)-v.hp(A);actual=v.hp(substitute(Q,r,B))-v.hp(substitute(Q,r,A));v.need((actual-factor*delta)%8==0 and actual%8!=0,'literal even ambient transfer');ambient.append([n,code,r,delta,factor,actual])
def cyclic(blocks,reverse=False):
 sizes=list(map(len,blocks));off=[sum(sizes[:i]) for i in range(3)];rows=[0]*sum(sizes)
 for i,a in enumerate(blocks):
  j=(i+(2 if reverse else 1))%3
  for u,row in enumerate(a):rows[off[i]+u]=(row<<off[i])|(((1<<sizes[j])-1)<<off[j])
 return rows
controls=[]
for C in ([2,0],[6,4,0],[2,4,1],[14,12,8,0]):
 T=cyclic([[2,0],[2,0],C]);U=cyclic([[2,0],[2,0],C],True);X=v.extend(T,3,12);Y=v.extend(U,3,12);p=check(X);q=check(Y);delta=[a-b for a,b in zip(p['q'],q['q'])];cx=substitute([2,4,1],0,X);cy=substitute([2,4,1],0,Y);ambient_delta=v.hp(cx)-v.hp(cy);v.need(ambient_delta==3*delta[1]+delta[2],'literal cyclic ambient coefficients')
 if len(C)%2==0:v.need(ambient_delta%8==4,'even port-core arbitrary transfer control')
 else:v.need(ambient_delta%8==0 and ambient_delta!=0 and (hc(X)-hc(Y))%2==1,'odd parity obstruction is not exact collision')
 controls.append({'child_C':C,'old_T':T,'old_U':U,'query':[3,12],'port_core_T':p,'port_core_U':q,'profile_difference':delta,'ambient_T_rows':cx,'ambient_U_rows':cy,'ambient_difference':ambient_delta})
print(json.dumps({'status':'PASS','scope':'All-order even-module odd-unit mod8 transfer, and exact odd-order Hamilton-cycle correction. No exact collision or unrestricted reconstruction claim.','tournaments_checked':len(records),'profile_records':records,'literal_ambient_substitutions':len(ambient),'ambient_records':ambient,'cyclic_port_core_controls':controls,'source_sha256':{str(p.relative_to(ROOT.parent)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [ROOT/'PROOF.md',ROOT/'verify.py',ROOT.parent/'covered_response'/'verify.py']}},sort_keys=True,indent=2))

#!/usr/bin/env python3
import json,hashlib,math
from pathlib import Path
Q=(30,64,10,66,174,335,21,367,95)
TAU=(0,7,8,5,6,3,4,1,2)
P=Path(__file__).resolve().parent

def need(p,s):
 if not p:raise RuntimeError(s)
def tr(n):return tuple(sum(1<<j for j in range(i+1,n))for i in range(n))
def converse(a):return tuple(sum(1<<j for j in range(len(a))if a[j]>>i&1)for i in range(len(a)))
def subset_h(a):
 n=len(a);N=1<<n;D=[[0]*n for _ in range(N)]
 for v in range(n):D[1<<v][v]=1
 for S in range(1,N):
  for v,c in enumerate(D[S]):
   if c:
    choices=a[v]&~S
    while choices:
     bit=choices&-choices;choices-=bit;D[S|bit][bit.bit_length()-1]+=c
 H=[sum(row)for row in D];H[0]=1;return H
def profile(a):
 H=subset_h(a);full=(1<<len(a))-1
 return H[-1],sum(H[S]*H[full^S]for S in range(full+1))
def build(m,U):
 blocks=[(0,)]*9;blocks[0]=tr(m);blocks[2]=U;blocks[8]=converse(U)
 names=[(i,j)for i,B in enumerate(blocks)for j in range(len(B))];index={p:v for v,p in enumerate(names)};a=[0]*len(names);tau=[]
 for u,(i,j)in enumerate(names):
  tau.append(index[(TAU[i],m-1-j if i==0 else j)])
  for v,(k,l)in enumerate(names):
   if i==k:arc=blocks[i][j]>>l&1
   else:arc=Q[i]>>k&1
   if arc:a[u]|=1<<v
 n=len(a)
 need(sorted(tau)==list(range(n))and all(tau[tau[u]]==u for u in range(n)),'anti permutation')
 for u in range(n):
  need(not(a[u]>>u&1),'loop')
  for v in range(u+1,n):need((a[u]>>v&1)+(a[v]>>u&1)==1,'tournament')
  for v in range(n):need((a[u]>>v&1)==(a[tau[v]]>>tau[u]&1),'physical anti')
 for u in range(n):
  S=1<<u
  while True:
   R=S
   for v in range(n):
    if S>>v&1:R|=a[v]
   if R==S:break
   S=R
  need(S==(1<<n)-1,'physical strongness')
 return tuple(a),tau,names

def main():
 controls=[]
 for U in [tr(1),tr(2),tr(3),(2,4,1)]:
  h,b=profile(U);need(h%2==1 and b%2==0 and b>=2,'physical profile parity')
  for m in range(1,5):
   a,tau,_=build(m,U);H=subset_h(a)[-1];want=h*h+18*b*b*3**m
   need(H==want and H%8==1,'literal paired identity');controls.append([list(U),m,len(a),H,h,b])
 arithmetic=0
 for h in range(1,32,2):
  for b in range(2,33,2):
   for m in range(1,9):
    H=h*h+18*b*b*3**m;y=3**(m//2+1)*b
    need(H==h*h+(6 if m%2 else 2)*y*y and H%8==1,'norm identity')
    if H%3:need(H%24==1,'prime-eligible residue')
    if H%7==0:
     if m%2==0:need(h%7==y%7==0 and H%49==0 and H>=649,'even-m obstruction')
     else:
      if (h-y)%7==0:x,z=(h+6*y)//7,(y-h)//7
      else:need((h+y)%7==0,'odd norm sign');x,z=(h-6*y)//7,(h+y)//7
      need(H//7==x*x+6*z*z and x%2==z%2==1,'quotient norm certificate')
    arithmetic+=1
 a,tau,names=build(9,tr(5));h,b=profile(tr(5));H=h*h+18*b*b*3**9;p=H//7
 need(len(a)==25 and H==6**11+1==362797057==7*51828151,'large exact identity value')
 bound=math.isqrt(p);need(bound==7199,'exact trial bound');remainders=[p%d for d in range(2,bound+1)]
 need(all(remainders),'prime quotient exhaustive trial division')
 need(p==6665**2+6*1111**2,'explicit quotient norm')
 need(6**5+1==7*11*101 and 6**7+1==7**2*29*197,'composite prime-exponent hostiles')
 witness={'rows':a,'anti_involution':tau,'vertex_block_and_local_index':names,'root_chain_length':9,'paired_module_order':5,'H':H,'prime_quotient':p,'H_count_method':'proved three-module identity, not a25-vertex endpoint DP','norm_coordinates':{'h':1,'y':7776,'quotient_x':-6665,'quotient_z':1111},'prime_certificate':{'method':'all integer trial divisors','first':2,'last':bound,'count':len(remainders),'all_remainders_nonzero':True,'remainders_sha256':hashlib.sha256(json.dumps(remainders,separators=(',',':')).encode()).hexdigest()}}
 data=json.dumps(witness,sort_keys=True,indent=2)+'\n';(P/'witness.json').write_text(data)
 print(json.dumps({'status':'PASS','literal_physical_controls':controls,'literal_graphs':len(controls),'integer_norm_controls':arithmetic,'large_witness_order':25,'large_H':H,'large_prime_quotient':p,'large_count_method':'proved three-module identity','witness_sha256':hashlib.sha256(data.encode()).hexdigest()},sort_keys=True,indent=2))
if __name__=='__main__':main()

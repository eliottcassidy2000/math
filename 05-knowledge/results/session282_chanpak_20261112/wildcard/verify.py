#!/usr/bin/env python3
"""Standalone exact receiver-carry and local-action controls. No project imports."""
from functools import lru_cache
from itertools import permutations,product
from math import gcd
import json,random

def check(ok,msg):
 if not ok:raise RuntimeError(msg)
def reverse(T):return tuple(sum(1<<u for u in range(len(T)) if T[u]>>v&1) for v in range(len(T)))
def hp_subsets(T):
 n=len(T);incoming=reverse(T)
 @lru_cache(None)
 def end(S,v):
  R=S^(1<<v)
  if not R:return 1
  return sum(end(R,u) for u in range(n) if R>>u&1 and incoming[v]>>u&1)
 H=[1]+[sum(end(S,v) for v in range(n) if S>>v&1) for S in range(1,1<<n)]
 return H,end
def packet(T):
 H,end=hp_subsets(T);RH,rend=hp_subsets(reverse(T));F=(1<<len(T))-1;D=F^1;vals=[]
 for HH,ee in ((RH,rend),(H,end)):
  mass=signed=0
  for S in range(1,F+1,2):
   x=ee(S,0)*HH[F^S];mass+=x;signed+=x*(1 if S.bit_count()%2 else -1)
  vals.extend((mass,signed))
 U,d,V,e=vals;Q=sum(H[S]*H[D^S] for S in range(F+1) if not S&1)
 return dict(H=H[F],h=H[D],Q=Q,U=U,V=V,d=d,e=e)
def permutation_packet(T):
 n=len(T);F=(1<<n)-1;paths={0:[()]}
 for S in range(1,F+1):paths[S]=[p for p in permutations([i for i in range(n) if S>>i&1]) if all(T[a]>>b&1 for a,b in zip(p,p[1:]))]
 U=V=d=e=0
 for S in range(1,F+1,2):
  sg=1 if S.bit_count()%2 else -1;o=len(paths[F^S]);a=sum(p[0]==0 for p in paths[S])*o;b=sum(p[-1]==0 for p in paths[S])*o
  U+=a;V+=b;d+=sg*a;e+=sg*b
 return dict(H=len(paths[F]),h=len(paths[F^1]),Q=sum(len(paths[S])*len(paths[(F^1)^S]) for S in range(F+1) if not S&1),U=U,V=V,d=d,e=e)
def mat(P):
 h=P['h'];return [[h,P['U']-h,P['d']-h,P['H']-h],[0,P['Q'],0,P['V']-h],[0,0,0,P['e']-h],[0,0,0,h]]
def mul(A,B):return [[sum(A[i][r]*B[r][j] for r in range(4)) for j in range(4)] for i in range(4)]
def decode(M):
 h=M[0][0];return dict(H=h+M[0][3],h=h,Q=M[1][1],U=h+M[0][1],V=h+M[1][3],d=h+M[0][2],e=h+M[2][3])
def word(Ps):
 M=[[int(i==j) for j in range(4)] for i in range(4)]
 for P in Ps:M=mul(M,mat(P))
 return decode(M)
def attach(D,cut):return (cut<<1,)+tuple((r<<1)|(1-((cut>>i)&1)) for i,r in enumerate(D))
def fuse(A,B):
 n=len(A);m=len(B);T=[0]*(n+m-1)
 for C,ids in ((A,list(range(n))),(B,[0]+list(range(n,n+m-1)))):
  for i,u in enumerate(ids):
   for j,v in enumerate(ids):
    if C[i]>>j&1:T[u]|=1<<v
 for i in range(1,n):
  for j in range(n,n+m-1):T[i]|=1<<j
 return tuple(T)
def carries(L,R):
 h=L[0]['h'];Q=L[0]['Q'];g=gcd(h,Q);p=h//g;q=Q//g;x=[0];y=[0]
 check(len(L)==len(R) and all((P['h'],P['Q'])==(h,Q) for P in L+R),'common base')
 for a,b in zip(L,R):
  t=p*x[-1]+b['V']-a['V']
  if t%q:return None
  x.append(t//q)
 for a,b in zip(reversed(L),reversed(R)):
  t=p*y[-1]+b['U']-a['U']
  if t%q:return None
  y.append(t//q)
 y.reverse()
 if x[-1] or y[0]:return None
 check(all(abs(z)<g for z in x+y),'strict bound');return x,y
def local_action(L,R,xy):
 h=L[0]['h'];p=h//gcd(h,L[0]['Q']);x,y=xy;J=h*sum(b['H']-a['H'] for a,b in zip(L,R))
 for i in range(len(L)-1):
  J-=p*((L[i]['U']-h)*x[i+1]+(R[i+1]['V']-h)*y[i+1])
  J+=(R[i]['d']-h)*(R[i+1]['e']-h)-(L[i]['d']-h)*(L[i+1]['e']-h)
 return J

def main():
 D=(18,8,3,5,14);Ts=[attach(D,c) for c in range(32)];Ps=[packet(T) for T in Ts]
 for T,P in zip(Ts,Ps):check(P==permutation_packet(T),'independent permutation packet')
 A=fuse(Ts[9],Ts[9]);B=fuse(Ts[3],Ts[24]);AP=packet(A);BP=packet(B);L=[Ps[9],Ps[9]];R=[Ps[3],Ps[24]]
 check(AP==word(L) and BP==word(R),'literal11 profile')
 check([AP['H'],BP['H']]==[2747,2749] and all(AP[k]==BP[k] for k in AP if k!='H'),'fullcentral+2')
 xy=carries(L,R);check(xy==([0,-1,0],[0,1,0]) and local_action(L,R,xy)==2,'nonzero carries and action')
 words=list(product(range(32),repeat=2));packets=[word([Ps[i] for i in w]) for w in words];accepted=central=0
 for i,w in enumerate(words):
  left=[Ps[t] for t in w];lp=packets[i]
  for j,z in enumerate(words):
   right=[Ps[t] for t in z];rp=packets[j];xy=carries(left,right);same=(lp['U'],lp['V'])==(rp['U'],rp['V'])
   check((xy is not None)==same,'alltwo carry iff')
   if same:
    accepted+=1;check(local_action(left,right,xy)==rp['H']-lp['H'],'alltwo action')
    central+=int((lp['d'],lp['e'])==(rp['d'],rp['e']))
 rng=random.Random(282);alphabets=[[packet(attach(C,c)) for c in range(1<<len(C))] for C in (D,(14,28,24,16,1),(32,17,3,7,13,30))]
 for alphabet in alphabets:
  h=alphabet[0]['h'];Q=alphabet[0]['Q'];g=gcd(h,Q)
  for _ in range(1000):
   k=rng.randrange(1,10);left=[rng.choice(alphabet) for _ in range(k)];right=[rng.choice(alphabet) for _ in range(k)]
   if rng.randrange(2):right=left[:]
   lp=word(left);rp=word(right);xy=carries(left,right)
   check((xy is not None)==((lp['U'],lp['V'])==(rp['U'],rp['V'])),'long carry iff')
   if xy is not None:
    J=local_action(left,right,xy);check(J==(rp['H']-lp['H'])*h if k==1 else J*h**(k-2)==rp['H']-lp['H'],'long action')
   check(lp['H']%g**max(0,k-2)==0,'content budget')
 # Closedform and residue theorem; first extra13divisor afterk2 occursk34.
 extra=[]
 for k in range(2,105):
  a=word([Ps[9]]*k);b=word([Ps[3],Ps[24]]+[Ps[9]]*(k-2));check(b['H']-a['H']==2*13**(k-2),'scaled rewriting')
  numer=8**k-7*k-1;check(numer%49==0,'integral geometric sum')
  val=169+312*k+1914*(numer//49)+40*(k-1);check(val*13**(k-2)==a['H'],'closedform')
  check(val%13==(8*k-pow(8,k,13))%13,'residue formula')
  check((val%13==0)==(k%52 in (1,34,44,51)),'fourresidue criterion')
  if val%13==0:extra.append(k)
 contexts=[(0,),(2,0),(0,1),(4,1,2)];ct=0
 for C in contexts:
  for E in contexts:
   if len(C)+len(E)>5:continue
   aa=fuse(fuse(C,A),E);bb=fuse(fuse(C,B),E);ch=packet(C)['h'] if len(C)>1 else 1;eh=packet(E)['h'] if len(E)>1 else 1
   check(hp_subsets(bb)[0][-1]-hp_subsets(aa)[0][-1]==2*ch*eh,'literal context transport');ct+=1
 result=dict(status='PASS',core=list(D),letter_profiles={str(i):Ps[i] for i in (9,3,24)},graphs=[list(A),list(B)],packets=[AP,BP],carry_x=[0,-1,0],carry_y=[0,1,0],local_action_terms=[26,-26,2],all_two_letter_comparisons=len(words)**2,matched_UV_pairs=accepted,full_central_pairs=central,random_length_controls=3000,context_controls=ct,closedform_maxk=104,extra13_indices=extra,tested_bases=[dict(h=a[0]['h'],Q=a[0]['Q'],g=gcd(a[0]['h'],a[0]['Q'])) for a in alphabets],scope='Exact carries/localaction/fullcentral identity; no absolute oddtarget coverage')
 print(json.dumps(result,sort_keys=True,indent=2))
if __name__=='__main__':main()

#!/usr/bin/env python3
import itertools,math,json,functools

def req(p,s):
 if not p:raise RuntimeError(s)

def relabel(T,p):return tuple(sum(1<<j for j,v in enumerate(p) if T[u]>>v&1) for u in p)

def canonical(T):
 groups=[[i for i in range(len(T)) if T[i].bit_count()==d] for d in sorted({r.bit_count() for r in T})]
 return min(relabel(T,sum(blocks,())) for blocks in itertools.product(*(itertools.permutations(g) for g in groups)))

def rooted(T,r):return min(relabel(T,(r,)+p) for p in itertools.permutations([v for v in range(len(T)) if v!=r]))

def extend(T,mask):
 n=len(T);return tuple(row|((not(mask>>i&1))<<n) for i,row in enumerate(T))+(mask,)

@functools.lru_cache(None)
def literal(T):return sum(all(T[u]>>v&1 for u,v in zip(w,w[1:])) for w in itertools.permutations(range(len(T))))

def induce(T,V):return tuple(sum(1<<j for j,v in enumerate(V) if T[u]>>v&1) for u in V)

def response(S,a):
 q=0;n=len(S)
 for mask in range(1<<n):
  if not(mask>>a&1):continue
  V=[v for v in range(n) if mask>>v&1];outside=tuple(v for v in range(n) if not(mask>>v&1))
  ending=sum(w[-1]==a and all(S[u]>>v&1 for u,v in zip(w,w[1:])) for w in itertools.permutations(V))
  q+=ending*literal(induce(S,outside))
 return q

def build(S,a):
 t=len(S);r=2*t;T=[0]*(r+1);T[r]=(1<<t)-1
 for i in range(t):
  T[t+i]|=1<<r
  for j in range(t):
   if S[i]>>j&1:T[i]|=1<<j;T[t+j]|=1<<(t+i)
   if i==j==a:T[i]|=1<<(t+j)
   else:T[t+j]|=1<<i
 return tuple(T)

def maps(T,anti=False):
 n=len(T);deg=[r.bit_count() for r in T]
 candidates=[[v for v in range(n) if (deg[u]+deg[v]==n-1 if anti else deg[u]==deg[v])] for u in range(n)]
 order=sorted(range(n),key=lambda u:len(candidates[u]));p=[-1]*n;answer=[]
 def go(depth,used):
  if depth==n:answer.append(tuple(p));return
  u=order[depth]
  for v in candidates[u]:
   if used>>v&1:continue
   if all(bool(T[u]>>w&1)==bool(T[p[w]]>>v&1 if anti else T[v]>>p[w]&1) for w in order[:depth]):
    p[u]=v;go(depth+1,used|(1<<v));p[u]=-1
 go(0,0);return answer

def reaches(T,root,reverse=False):
 reached={root};todo=[root]
 while todo:
  u=todo.pop()
  for v in range(len(T)):
   arc=(T[v]>>u&1) if reverse else(T[u]>>v&1)
   if arc and v not in reached:reached.add(v);todo.append(v)
 return reached

def hp(T):
 n=len(T);dp=[[0]*n for _ in range(1<<n)]
 for v in range(n):dp[1<<v][v]=1
 for mask in range(1,1<<n):
  for v in range(n):
   if not(mask>>v&1) or mask==(1<<v):continue
   rest=mask^(1<<v);dp[mask][v]=sum(dp[rest][u] for u in range(n) if rest>>u&1 and T[u]>>v&1)
 return sum(dp[-1])

def main():
 bank={(0,)};A=[1];R=[0];C=[0];L=[0];CL=[0];classes=halfwords=automaps=antimaps=0;fullgroups=[]
 for t in range(1,6):
  if t>1:bank={canonical(extend(T,m)) for T in bank for m in range(1<<(t-1))}
  sides={rooted(S,a) for S in bank for a in range(t)};A.append(len(bank));R.append(len(sides));C.append(0);L.append(0);CL.append(0);outputs=set()
  for S in sorted(sides):
   a=0;G=[p for p in maps(S) if p[a]==a];g=len(G);hs=literal(S);q=response(S,a)
   T=build(S,a);r=2*t;tau=tuple(list(range(t,2*t))+list(range(t))+[r]);auts=maps(T);anti=[p for p in maps(T,True) if all(p[p[i]]==i for i in range(2*t+1))]
   rootauts=[p for p in auts if p[r]==r];rootanti=[p for p in anti if p[r]==r];cent=[p for p in auts if all(p[tau[i]]==tau[p[i]] for i in range(2*t+1))]
   req(len(rootauts)==g*g and len(rootanti)==g and len(cent)==g,'root group, anti and full centralizer')
   req(len(anti)==len(auts)//g,'full anti-involution conjugacy orbit')
   req({p[r] for p in auts}=={next(i for i in range(2*t+1) if p[i]==i) for p in anti},'unmarked fixed-root orbit')
   req(hs%g==0 and q%g==0 and(hs//g)%2==1,'rooted cover orbit divisibility')
   actual=hp(T);req(actual==hs*hs+2*q*q,'independent full graph count')
   req(actual%(g*g)==0 and actual//(g*g)==(hs//g)**2+2*(q//g)**2,'root-stabilizer path orbits')
   req(actual%len(auts)==0,'full unmarked path orbits')
   fixed=0
   for bits in range(1<<t):
    selected=[i+t if bits>>i&1 else i for i in range(t)]
    for first in itertools.permutations(selected):
     w=first+(r,)+tuple(tau[v] for v in reversed(first));fixed+=all(T[u]>>v&1 for u,v in zip(w,w[1:]));halfwords+=1
   req(fixed==hs,'full half-word fixed set')
   strong=len(reaches(T,r))==len(T) and len(reaches(T,r,True))==len(T)
   req(strong==(len(reaches(S,a,True))==t),'terminal side criterion')
   mass=hs//g;L[t]+=mass
   if strong:C[t]+=1;CL[t]+=mass
   outputs.add(canonical(T));classes+=1;automaps+=len(auts);antimaps+=len(anti);fullgroups.append((t,g,len(auts)))
  req(len(outputs)==len(sides),'unmarked physical injection')
  req(L[t]==t*2**((t-1)*(t-2)//2),'total rank-one mass')
  req(R[t]==sum(A[i]*C[t-i] for i in range(t)),'ordinary strong convolution')
  M=[1]+[2**((j-1)*(j-2)//2) for j in range(1,t+1)]
  req(t*M[t]==sum(M[i]*CL[t-i] for i in range(t)),'weighted strong convolution')
 req(A==[1,1,1,2,4,12] and R==[0,1,2,4,12,48],'complete side/root inventories')
 req(C==[0,1,1,2,7,33] and CL==[0,1,1,4,25,281],'strong extension through side order5')
 req((1,1,3) in fullgroups,'C3 full group differs from stabilizer')
 M=[1]+[2**((j-1)*(j-2)//2) for j in range(1,51)];CM=[0]*51
 for n in range(1,51):
  CM[n]=n*M[n]-sum(M[i]*CM[n-i] for i in range(1,n));req(CM[n]%2==int(n%3!=0),'optional all-order parity consequence control')
 print(json.dumps({'status':'PASS','unlabelled_tournaments':A,'rooted_side_classes':R,'strong_rank_one_classes':C,'rank_one_mass':L,'strong_rank_one_mass':CL,'physical_source_classes':classes,'full_automorphisms_enumerated':automaps,'full_anti_involutions_enumerated':antimaps,'fixed_half_words_tested':halfwords,'unmarked_output_canonicalization':'complete through side order5','parity_coefficient_controls_through':50},sort_keys=True,indent=2))
if __name__=='__main__':main()

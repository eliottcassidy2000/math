"""Exact marked-pair quotient-converse controls, standard library only."""
from itertools import combinations,product
from functools import lru_cache
from pathlib import Path
import hashlib,json,random


def need(ok,msg):
 if not ok:raise RuntimeError(msg)
def rows(n,code):
 T=[0]*n
 for k,(u,v) in enumerate(combinations(range(n),2)):
  if code>>k&1:T[u]|=1<<v
  else:T[v]|=1<<u
 return T
def opposite(T):return [((1<<len(T))-1)^r^(1<<i) for i,r in enumerate(T)]
def extend(T,a,c):
 n=len(T);R=list(T)+[1<<(n+1),0]
 for u in range(n):
  if a>>u&1:R[u]|=1<<n
  else:R[n]|=1<<u
  if c>>u&1:R[u]|=1<<(n+1)
  else:R[n+1]|=1<<u
 return R
def hp(T):
 n=len(T);N=(1<<n)-1
 @lru_cache(None)
 def f(S,u):
  if S==N:return 1
  return sum(f(S|1<<v,v) for v in range(n) if not S>>v&1 and T[u]>>v&1)
 return sum(f(1<<u,u) for u in range(n)) if n else 1
def endpoints(T):
 n=len(T);D=[[0]*n for _ in range(1<<n)]
 for u in range(n):D[1<<u][u]=1
 for S in range(1,1<<n):
  for u,z in enumerate(D[S]):
   if not z:continue
   rem=T[u]&~S
   while rem:
    bit=rem&-rem;rem-=bit;D[S|bit][bit.bit_length()-1]+=z
 return D
def covers(T):
 n=len(T);N=(1<<n)-1;E=endpoints(T);S=endpoints(opposite(T));C=[[0]*n for _ in T]
 for M in range(1,N):
  for u in range(n):
   for v in range(n):C[u][v]+=E[M][u]*S[N^M][v]
 return C
def q2(T):return sum(map(sum,covers(T)))
def induced(T,vs):return [sum(1<<j for j,v in enumerate(vs) if T[u]>>v&1) for u in vs]
def blow(Q,blocks):
 off=[];n=0
 for A in blocks:off.append(n);n+=len(A)
 T=[0]*n;sets=[]
 for i,A in enumerate(blocks):
  sets.append(list(range(off[i],off[i]+len(A))))
  for u,r in enumerate(A):T[off[i]+u]|=r<<off[i]
  for j,B in enumerate(blocks):
   if Q[i]>>j&1:
    for u in sets[-1]:T[u]|=((1<<len(B))-1)<<off[j]
 return T,sets
def strong(T):
 N=(1<<len(T))-1
 for u in range(len(T)):
  S=1<<u
  while True:
   M=S
   for v in range(len(T)):
    if S>>v&1:M|=T[v]
   if M==S:break
   S=M
  if S!=N:return False
 return True
def prime(T):
 n=len(T);N=(1<<n)-1
 for M in range(1,N):
  if M.bit_count()<2:continue
  if all(T[u]&M in (0,M) for u in range(n) if not M>>u&1):return False
 return True

def battery(Q,blocks,z,I,J):
 need(blocks[z]==[2,0],'pair orientation p to q')
 T,sets=blow(Q,blocks);U,_=blow(opposite(Q),blocks);p,q=sets[z]
 A=sum((sets[i] for i in I),[]);B=sum((sets[i] for i in J),[])
 need(not(set(A)&set(B)) and z not in I and z not in J,'disjoint exterior marked sets')
 am=sum(1<<i for i in A);bm=sum(1<<i for i in B);data=[];coef=0
 for bits in product((0,1),repeat=4):
  a=(bits[0]<<q)|(am if bits[1] else 0);c=(bits[2]<<p)|(bm if bits[3] else 0)
  f=hp(extend(T,a,c));g=hp(extend(U,a,c));sign=(-1)**(4-sum(bits));coef+=sign*(f-g)
  data.append({'bits':bits,'a':a,'c':c,'T':f,'U':g,'sign':sign})
 keep=[v for v in range(len(T)) if v not in (p,q)];R=induced(T,keep);index={v:i for i,v in enumerate(keep)};C=covers(R)
 want=sum(C[index[j]][index[i]]-C[index[i]][index[j]] for i in A for j in B)
 need(coef==want,'grouped pair coefficient')
 return {'quotient':Q,'children':blocks,'pair_index':z,'marked_child_sets':[I,J],'rows_T':T,'rows_U':U,'quartic_coefficient':coef,'cover_skew':want,'literal_battery':data}

def main():
 rng=random.Random(188821);general=[]
 for k in (3,4,5):
  for trial in range(5):
   Q=rows(k,rng.randrange(1<<(k*(k-1)//2)));z=trial%k
   blocks=[[0] for _ in Q];blocks[z]=[2,0]
   for j in range(k):
    if j!=z and (j+trial)%3==0:blocks[j]=[2,0] if trial%2 else [0,1]
   remaining=[j for j in range(k) if j!=z];I=remaining[:max(1,len(remaining)//2)];J=remaining[len(I):]
   general.append(battery(Q,blocks,z,I,J))
 # All small arbitrary ordered-cut profiles, not just transitive examples.
 cutchecks=0
 small=[rows(n,c) for n in range(1,4) for c in range(1<<(n*(n-1)//2))]
 for A,B in product(small,repeat=2):
  T,sets=blow([2,0],[A,B]);C=covers(T);AB=sum(C[u][v] for u in sets[0] for v in sets[1]);BA=sum(C[v][u] for u in sets[0] for v in sets[1]);ha,hb=hp(A),hp(B);qa,qb=q2(A),q2(B)
  need(AB==ha*hb and BA==(ha+qa)*(hb+qb),'ordered-cut exact skew')
  need((BA>AB)==(len(A)+len(B)>=3),'sharp two-singleton boundary');cutchecks+=1
 Q=[14,28,8,16,5];need(strong(Q) and prime(Q),'explicit strong prime quotient')
 need(not strong(induced(Q,[0,1,2,3])),'deletion is not strong')
 witnesses=[]
 for blocks in ([[2,0]]*5,[[2,4,1],[2,0],[0],[0],[2,0]],[[0],[0],[0],[0],[2,0]]):
  w=battery(Q,blocks,4,[0,1],[2,3]);need(w['quartic_coefficient']>0,'prime cut separator');witnesses.append(w)
 need(witnesses[0]['quartic_coefficient']==224,'ten-name positive quartic')
 regular=[6,12,24,17,3];need(prime(regular) and all(strong(induced(regular,[v for v in range(5) if v!=u])) for u in range(5)),'regular-five residual control')
 root=Path(__file__).resolve().parent
 print(json.dumps({'status':'PASS','general_grouped_coefficients':general,'ordered_cut_checks':cutchecks,'strong_prime_controls':witnesses,'regular_five_outside_cut_corollary':regular,'source_sha256':{str(p.relative_to(root.parent)):hashlib.sha256(p.read_bytes()).hexdigest() for p in (root/'PROOF.md',root/'verify.py',root.parent/'response'/'PROOF.md')}},indent=2,sort_keys=True))
if __name__=='__main__':main()

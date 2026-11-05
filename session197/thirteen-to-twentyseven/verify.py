from itertools import permutations,product,combinations
import json

def need(b,s):
 if not b:raise RuntimeError(s)
def hp(T):
 n=len(T);N=1<<n;D=[[0]*n for _ in range(N)]
 for v in range(n):D[1<<v][v]=1
 for M in range(1,N):
  for u in range(n):
   c=D[M][u]
   if not c:continue
   bits=T[u]&~M
   while bits:
    b=bits&-bits;bits-=b;D[M|b][b.bit_length()-1]+=c
 return sum(D[-1])
def canonical(T):
 classes={}
 for v,x in enumerate(T):classes.setdefault(x.bit_count(),[]).append(v)
 choices=[list(permutations(classes[d]))for d in sorted(classes)];best=None
 for parts in product(*choices):
  p=sum(parts,());rows=tuple(sum(((T[a]>>b)&1)<<j for j,b in enumerate(p))for a in p)
  if best is None or rows<best:best=rows
 return best

def extend(T,M):
 n=len(T);U=list(T)+[M]
 for a in range(n):
  if not(M>>a&1):U[a]|=1<<n
 return U

def induced(T,vs):return [sum(((T[a]>>b)&1)<<j for j,b in enumerate(vs))for a in vs]
def strong(T):
 n=len(T);full=(1<<n)-1
 for v in range(n):
  bits=1<<v;old=0
  while old!=bits:
   old=bits
   for a in range(n):
    if bits>>a&1:bits|=T[a]
  if bits!=full:return False
 return True

def main():
 bank={3:{canonical([2,4,1]):3}};levels=[dict(order=3,candidates=0,accepted=1,classes=1)]
 for n in range(4,9):
  nxt={};candidates=accepted=0
  for T in sorted(bank[n-1]):
   for M in range(1,(1<<(n-1))-1):
    candidates+=1;U=extend(T,M);h=hp(U)
    if h>27:continue
    need(strong(U),'mixed extension strong');accepted+=1;C=canonical(U)
    if C in nxt:need(nxt[C]==h,'canonical count consistency')
    nxt[C]=h
  bank[n]=nxt;levels.append(dict(order=n,candidates=candidates,accepted=accepted,classes=len(nxt)))
 need(not bank[8],'empty order-eight cap frontier')
 thirteen=[];twentyseven=[];small=[]
 for n,B in bank.items():
  for T,h in sorted(B.items()):
   if h==13:
    need(n==5,'all strong H13 cores have order5');thirteen.append(dict(rows=T,H=h))
   if h==27:
    counts=sorted(hp(induced(T,vs))for vs in combinations(range(n),5))
    need(13 not in counts,'no induced H13 core');twentyseven.append(dict(rows=T,H=h,five_subset_H=counts))
   if n in [3,4]:small.append([n,h])
 need(len(thirteen)==1 and len(twentyseven)==5,'exact core counts')
 # Check the H27 hostile family directly rather than treating noncontainment as nonexistence.
 need(hp([2,32,11,19,7,29])==27,'H27 positive control')
 need(hp([34,12,41,49,39,2])==25,'opposite binary branch positive control')
 W=[2,4,73,83,103,79,3];V=induced(W,[0,1,2,3,4,6])
 need(hp(W)==29 and hp(V)==15,'H15-to-H29 induced hostile')
 # Full labelled <=5 controls compare exactly with the generated small banks.
 label_counts=[]
 for n in range(3,6):
  E=list(combinations(range(n),2));got=set();labelled=0
  for code in range(1<<len(E)):
   T=[0]*n
   for z,(a,b)in enumerate(E):T[a if code>>z&1 else b]|=1<<(b if code>>z&1 else a)
   if strong(T) and hp(T)<=27:got.add(canonical(T));labelled+=1
  need(got==set(bank[n]),'independent labelled small bank');label_counts.append([n,labelled,len(got)])
 out=dict(status='PASS',empty_frontier_order=8,H15_to_H29_hostile=dict(target=W,subset=[0,1,2,3,4,6],source=V),strong_cap=27,levels=levels,labelled_small_controls=label_counts,H13_cores=thirteen,H27_cores=twentyseven,bank=[dict(order=n,rows=T,H=h)for n,B in bank.items()for T,h in sorted(B.items())])
 print(json.dumps(out,sort_keys=True,indent=2))
if __name__=='__main__':main()

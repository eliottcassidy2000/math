from pathlib import Path
from itertools import combinations
from collections import Counter
import json,random,hashlib
D=Path(__file__).resolve().parent

def need(ok,msg):
 if not ok:raise RuntimeError(msg)
def edge(u,v):return min(u,v),max(u,v)
def edges(q):return {edge(u,v) for u,v in zip(q,q[1:])}
def blue(T,q):
 r=dict(zip(q,q[::-1]));return all((r[v],r[u]) in T for u,v in T)
def path(T,q):return len(set(q))==len(q) and all(e in T for e in zip(q,q[1:]))
def move(T,q):
 E=edges(q);return {(u,v) if edge(u,v) in E else (v,u) for u,v in T}
def source(n,h,rng):
 I={h,n-3,n-2};T=set();seen=set()
 for u,v in combinations(range(1,n-1),2):
  if (u,v) in seen:continue
  orb={(u,v),(n-1-v,n-1-u)};seen|=orb;bit=v==u+1 or bool(rng.getrandbits(1));T.update((a,b) if bit else (b,a) for a,b in orb)
 J={n-1-v for v in I}
 for v in range(1,n-1):T.add((v,0) if v in I else (0,v));T.add((n-1,v) if v in J else (v,n-1))
 T.add((n-1,0));return T

def odd_words(n,l):
 d,u=tuple(range(3,l+3)),tuple(range(l+3,n-1));return ((0,2)+d+u+(n-1,1),d+(0,2)+u+(n-1,1),d+(0,1,2)+u+(n-1,))
def x_words(n):
 d=tuple(range(3,n-3));return ((0,1,2)+d+(n-1,n-3,n-2),(0,2)+d+(n-1,n-3,1,n-2),(0,2)+d+(n-3,1,n-2,n-1))
def replay(T,n,h,words):
 P=tuple(range(n));need(path(T,P) and blue(T,P),'marked source');need({u for u,v in T if v==0}=={h,n-3,n-2,n-1},'actual ports');need(sum(u==n-1 for u,v in T)==4,'score')
 now=T;xor=set()
 for q in words:need(set(q)==set(P) and path(now,q),'directed');need(not blue(now,q),'black');xor^=edges(q);now=move(now,q)
 need(xor==edges(P) and now==move(T,P),'literal target')

def cyc(p):
 seen=set();out=[]
 for i in range(len(p)):
  if i in seen:continue
  v=i;c=[]
  while v not in seen:seen.add(v);c.append(v);v=p[v]
  out.append(c)
 return out

def main():
 rng=random.Random(202610210019);counts=Counter();digest=hashlib.sha256()
 for n in range(15,61):
  for h in [2,3,4,5,n-4]+[rng.randrange(2,n-3) for _ in range(25)]:
   for _ in range(2):
    T=source(n,h,rng)
    for l in range(1,n//2-2,2):
     if n<2*(l+3) or h in {2,l+2} or (l>=3 and h in {4,l-1}) or (l+3,2) not in T:continue
     words=odd_words(n,l);replay(T,n,h,words);counts[f'odd_prefix_{l}']+=1;digest.update(json.dumps([n,h,l,sorted(T)],separators=(',',':')).encode())
    if h==2 and (n-2,1) in T and (2,n-2) in T:
     replay(T,n,h,x_words(n));counts['second_port_case']+=1
 perm_checks=0
 for n in range(15,201):
  for l in range(1,n//2-2,2):
   if n<2*(l+3):continue
   q=odd_words(n,l)[2];r=dict(zip(q,q[::-1]));g=[r[n-1-i] for i in range(n)];lens=[len(c) for c in cyc(g) if len(c)>1];L=l+3
   from math import gcd
   need(lens and set(lens)=={L//gcd(L,3)} and all(k%2==0 for k in lens),'even block cycle formula');need(sum(lens)==2*L,'two disjoint full blocks');perm_checks+=1
 hostiles=[]
 for w in json.loads((D/'old-template-hostiles.json').read_text()):
  n,h=w['n'],w['h'];T={(u,v) for u,row in enumerate(w['rows']) for v in range(n) if row>>v&1};need(len(T)==n*(n-1)//2 and all((v,u) not in T for u,v in T),'complete tournament');now=T;flags=[]
  for q in w['words']:need(path(now,q),'old words legal');flags.append(blue(now,q));now=move(now,q)
  need(flags==[False,False,True] and now==move(T,tuple(range(n))),'genuine old final-word hostile');need((6,2) in T and h not in {2,4,5},'uniform cut-six guards');words=odd_words(n,3);replay(T,n,h,words);hostiles.append({'n':n,'h':h,'old_blue_steps':flags,'repair_words':words})
 print(json.dumps({'status':'PASS','checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'hostile_source_sha256':hashlib.sha256((D/'old-template-hostiles.json').read_bytes()).hexdigest(),'controls':dict(counts),'permutation_checks':perm_checks,'repaired_same_labelled_sources':hostiles,'replay_digest_sha256':digest.hexdigest()},sort_keys=True,indent=2))
if __name__=='__main__':main()

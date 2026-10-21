from pathlib import Path
from itertools import combinations
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
def source(n,I,rng):
 T=set();seen=set()
 for u,v in combinations(range(1,n-1),2):
  if (u,v) in seen:continue
  orb={(u,v),(n-1-v,n-1-u)};seen|=orb;bit=v==u+1 or bool(rng.getrandbits(1))
  T.update((a,b) if bit else (b,a) for a,b in orb)
 J={n-1-v for v in I}
 for v in range(1,n-1):T.add((v,0) if v in I else (0,v));T.add((n-1,v) if v in J else (v,n-1))
 T.add((n-1,0));return T

def replay(T,n,I,words,required):
 P=tuple(range(n));need(path(T,P) and blue(T,P),'marked source');need(sum(v==0 for u,v in T)==4 and sum(u==n-1 for u,v in T)==4,'degrees')
 need({u for u,v in T if v==0}==I|{n-1},'actual ports');now=T;xor=set();flags=[]
 for q in words:need(path(now,q) and set(q)==set(P),'directed');flags.append(blue(now,q));xor^=edges(q);now=move(now,q)
 need(xor==edges(P) and now==move(T,P),'literal target')
 if required:need(not any(flags),'three black')
 return flags

def chain(n):
 y,c,z=n-3,n-2,n-1;W=tuple(range(3,y));q1=(0,2)+W+(y,c,z,1);q2=(z,)+W+(y,2,1,0,c);r=dict(zip(q2,q2[::-1]));E=edges(q1)
 pairs=[(0,4),(1,3),(n-4,c),(7,z),(0,n-8)];total=0
 for k,((u,v),target) in enumerate(zip(pairs,pairs[1:])):
  if k%2:need(target==(n-1-v,n-1-u),'P reflection link')
  else:
   a,b=r[v],r[u];need(edge(a,b)==target,'Q reflection link');total^=int(a>b)^int((u,v) in E)^int(target in E)
 need(total==1 and all(v not in {5,y,c,z} for v in (4,n-8)),'two-true-port contradiction')

def main():
 rng=random.Random(2026102108);count=0;digest=hashlib.sha256()
 for n in range(15,41):
  y,c,z=n-3,n-2,n-1;W=tuple(range(3,y));words=((0,2)+W+(y,c,z,1),(z,)+W+(y,2,1,0,c),(1,z)+W+(y,2,0,c))
  for h in range(3,y):
   if h==n-4:continue
   I={h,y,c}
   for _ in range(16):
    T=source(n,I,rng)
    if not {(c,1),(c,2),(2,y)}<=T:continue
    flags=replay(T,n,I,words,True);count+=1;digest.update(json.dumps([n,h,sorted(T),flags],separators=(',',':')).encode())
 for n in range(15,201):chain(n)
 hostiles=[]
 for w in json.loads((D/'witnesses.json').read_text()):
  n,h=w['n'],w['h'];T={(u,v) for u,row in enumerate(w['rows']) for v in range(n) if row>>v&1};need(len(T)==n*(n-1)//2 and all((v,u) not in T for u,v in T),'complete tournament');need({(n-2,1),(n-2,2),(n-3,2)}<=T,'hostile branch')
  flags=replay(T,n,{h,n-3,n-2},w['words'],False);need(flags[2] and blue(move(T,tuple(range(n))),w['words'][2]),'actual blue final word');hostiles.append({'n':n,'h':h,'blue_steps':flags})
 print(json.dumps({'status':'PASS','checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'witnesses_sha256':hashlib.sha256((D/'witnesses.json').read_bytes()).hexdigest(),'three_black_controls':count,'chain_orders':[15,200],'hostiles':hostiles,'replay_digest_sha256':digest.hexdigest()},sort_keys=True,indent=2))
if __name__=='__main__':main()

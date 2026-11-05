import itertools,json,random,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'ribbon'))
from transfer import graph,held_karp

def components(a):
 n=len(a);r=[{i}|{j for j in range(n) if a[i][j]} for i in range(n)]
 for k in range(n):
  for i in range(n):
   if k in r[i]:r[i]|=r[k]
 seen=set();cs=[]
 for i in range(n):
  if i not in seen:
   c={j for j in r[i] if i in r[j]};seen|=c;cs.append(c)
 return cs

def triangles(a):
 return [set(t) for t in itertools.combinations(range(len(a)),3) if all(sum(a[i][j] for j in t)==1 for i in t)]

def check(word,w):
 a=graph(word,w);h=held_karp(a)
 keep=sorted(set().union(*(c for c in components(a) if len(c)>1)))
 if not keep:assert h==1;return
 b=[[a[i][j] for j in keep] for i in keep];n=len(b)
 assert held_karp(b)==h
 assert all(not b[j][i] or j-i<=w for j in range(n) for i in range(j))
 ts=triangles(b)
 assert all(max(t)-min(t)<=2*w for t in ts)
 assert set().union(*ts)==set(range(n))
 selected=[next(t for t in ts if i in t) for i in range(0,n,4*w+1)]
 assert len(set().union(*selected))==3*len(selected)
 assert h>=3**len(selected)
 k=0;p=1
 while p*3<=h:p*=3;k+=1
 assert n<=(4*w+1)*k

def paths(a,vs):
 return [p for p in itertools.permutations(vs) if all(a[x][y] for x,y in zip(p,p[1:]))]

def merge(a,p,q):
 p=list(p);q=list(q);out=[]
 while p and q:
  out.append(p.pop(0) if a[p[0]][q[0]] else q.pop(0))
 return tuple(out+p+q)

if __name__=='__main__':
 results={}
 for w,n in [(1,7),(2,7),(3,6)]:
  alphas=[range(0,1<<w,1<<max(0,w-j)) for j in range(n)]
  total=0
  for word in itertools.product(*alphas):check(word,w);total+=1
  results[f'compression_width{w}_order{n}']=total
 rng=random.Random(197)
 cases=0
 for n in range(2,9):
  for _ in range(12):
   a=graph([rng.randrange(8) for _ in range(n)],3)
   for k in range(1,n):
    pp=paths(a,range(k));qq=paths(a,range(k,n));outs=set()
    for p in pp:
     for q in qq:
      z=merge(a,p,q)
      assert all(a[x][y] for x,y in zip(z,z[1:]))
      assert tuple(x for x in z if x<k)==p
      assert tuple(x for x in z if x>=k)==q
      outs.add(z)
    assert len(outs)==len(pp)*len(qq);cases+=1
 results['merge_split_checks']=cases
 print(json.dumps(results,indent=2))
 Path(__file__).with_name('VALIDATION.json').write_text(json.dumps(results,indent=2)+'\n')

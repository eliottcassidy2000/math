import json,itertools
from pathlib import Path
from collections import Counter,defaultdict

def matrix(n,w):
 A=[[False]*n for _ in range(n)]
 for i in range(n-1):
  for j in range(i+1,n-1): A[i][j]=True
 for i,b in enumerate(w): A[n-1][i]=b=='1'; A[i][n-1]=b=='0'
 for i in range(n-1):
  for j in range(i): A[i][j]=not A[j][i]
 return A

def canonical(A):
 n=len(A);words=[]
 for a in range(n):
  B=[i for i in range(n) if i!=a]
  scores={i:sum(A[i][j] for j in B) for i in B}
  B.sort(key=lambda i:-scores[i])
  if all(A[B[i]][B[j]] for i in range(n-1) for j in range(i+1,n-1)):
   words.append(''.join(str(int(A[a][b])) for b in B))
 if not words:return None
 ordinary=min(words); merged=min(words+[w[::-1].translate(str.maketrans('01','10')) for w in words])
 return ordinary,merged

def run(n):
 pairs=[(i,j) for i in range(n) for j in range(i+2,n)]; pindex={p:k for k,p in enumerate(pairs)}; full=(1<<len(pairs))-1
 def reflect(c):return sum(((c>>pindex[n-1-j,n-1-i])&1)<<k for k,(i,j) in enumerate(pairs))
 tilings={}; ordinary=set(); nodes=set()
 for bits in itertools.product('01',repeat=n-1):
  A=matrix(n,bits);o,m=canonical(A);ordinary.add(o);nodes.add(m)
  def visit(path,left):
   if not left:
    c=sum(int(A[path[i]][path[j]])<<k for k,(i,j) in enumerate(pairs))
    assert c not in tilings or tilings[c]==(o,m)
    tilings[c]=(o,m);return
   for v in left:
    if not path or A[path[-1]][v]:visit(path+[v],[x for x in left if x!=v])
  visit([],list(range(n)))
 orbits={};internal=[]
 for c,(o,m) in tilings.items():
  if c^full not in tilings:continue
  r=reflect(c);orb=min(c,c^full,r,r^full)
  if orb not in orbits:
   e={'id':orb,'endpoints':sorted([m,tilings[c^full][1]]),'color':'blue' if r==c else 'black','tilings':sorted(set([c,c^full,r,r^full]))};orbits[orb]=e;internal.append(e)
 bypair=Counter(tuple(e['endpoints']) for e in internal)
 touched=set(x for e in internal for x in e['endpoints']);H=Counter(m for o,m in tilings.values());h_ord=Counter(o for o,m in tilings.values())
 blue=sum(reflect(c)==c for c in tilings)
 result={'n':n,'ordinary_classes':len(ordinary),'merged_nodes':len(nodes),'selected_tilings':len(tilings),'blue_incidences':blue,'internal_tilings':sum(c^full in tilings for c in tilings),'internal_edges':len(internal),'internal_blue_edges':sum(e['color']=='blue' for e in internal),'internal_black_edges':sum(e['color']=='black' for e in internal),'internal_loops':sum(e['endpoints'][0]==e['endpoints'][1] for e in internal),'isolated_nodes':len(nodes-touched),'edges':sorted(internal,key=lambda e:(e['endpoints'],e['id'])),'nodes':[{'word':w,'ordinary_tiling_mass':h_ord.get(w),'merged_tiling_mass':H[w]} for w in sorted(nodes)]}
 print({k:v for k,v in result.items() if k not in ['edges','nodes']},flush=True)
 return result

if __name__=='__main__':
 results=[run(n) for n in range(3,9)]
 with open(Path(__file__).with_name('catalogue.json'),'w') as f:json.dump(results,f,indent=2)
 print(json.dumps(results[-1]['edges'],indent=2))

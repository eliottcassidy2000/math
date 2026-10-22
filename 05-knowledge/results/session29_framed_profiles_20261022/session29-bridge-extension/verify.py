"""Independent literal verification. No assertions; works unchanged under -O."""
from itertools import combinations,permutations
from pathlib import Path
import json
D=Path(__file__).parent
R=json.loads((D/'results.json').read_text());Z=json.loads((D/'detours-results.json').read_text())
def need(x,s):
 if not x:raise RuntimeError(s)
def pairs(n):return list(combinations(range(n),2))
def mat(m,n):
 a=[[0]*n for _ in range(n)]
 for k,(i,j) in enumerate(pairs(n)):a[i][j]=(m>>k)&1;a[j][i]=1-a[i][j]
 return a
def code(a):return sum(a[i][j]<<k for k,(i,j) in enumerate(pairs(len(a))))
def reorder(a,p,op=False):return [[(a[j][i] if op else a[i][j]) for j in p] for i in p]
def canon(a):return min(code(reorder(a,p)) for p in permutations(range(len(a))))
def rooted(a,root):return min(code(reorder(a,p+(root,))) for p in permutations(tuple(i for i in range(len(a)) if i!=root)))
def target(e,f):
 a=[[0]*10 for _ in range(10)];x=mat(e,5);y=mat(f,5)
 for i in range(5):
  for j in range(5):a[i][j]=x[i][j];a[i+5][j+5]=y[i][j];a[i][j+5]=1
 return flip(a,list(range(10)))
def flip(a,p):
 need(sorted(p)==list(range(len(a))),'path vertices')
 need(all(a[i][j] for i,j in zip(p,p[1:])),'path edges')
 b=[[a[j][i] for j in range(len(a))] for i in range(len(a))]
 for i,j in zip(p,p[1:]):b[i][j]=1;b[j][i]=0
 return b
def ep(word):
 e,f=word;x=mat(e,5);y=mat(f,5)
 s=tuple(sorted((canon(x),canon(reorder(y,range(5),True)))))
 v=target(e,f);a=[row[:5] for row in v[:5]];b=[row[5:] for row in v[5:]]
 q=tuple(sorted((rooted(a,4),rooted(reorder(b,range(5),True),0))))
 return ('S',)+s,('Q',)+q
def canonical_target(q,r):
 a=[[0]*10 for _ in range(10)];x=mat(q,5);y=mat(r,5)
 for i in range(5):
  for j in range(5):a[i][j]=x[i][j];a[i+5][j+5]=y[j][i];a[j+5][i]=1
 a[4][9]=1;a[9][4]=0
 return code(a)
count=0
for rec in R['original_eight']:
 if rec['alternate_actual_words'] is None:continue
 s,t=ep(rec['actual_original_word']);cur=s;visited={s}
 for word in rec['alternate_actual_words']:
  u,v=ep(word);need(cur in (u,v),'alternate incidence');cur=v if cur==u else u
  need(cur not in visited,'simple alternate path');visited.add(cur)
 need(cur==t,'alternate terminal');count+=1
for rec in Z:
 T=mat(rec['target_mask'],10);U=flip(T,rec['target_to_U_path']);need(code(U)==rec['U_mask'],'first literal flip')
 V=flip(U,rec['U_to_V_path']);need(code(V)==rec['V_mask'],'second literal flip')
 w=rec['V_to_canonical_permutation'];p=sorted(range(10),key=lambda i:w[i]);m=code(reorder(V,p,rec['V_converse']))
 need(m==canonical_target(*rec['V_address']),'literal rooted isomorphism')
 old=next(r for r in R['original_eight'] if r['tile']==rec['tile']);need(code(target(*old['actual_original_word']))==rec['target_mask'],'original target mask')
 need(ep(rec['boundary_path_actual'][0])==(('S',41,41),('Q',)+tuple(rec['V_address'])),'mixed boundary edge')
 # Certificate U is a 4+6 join via first-four path vertices; canonicalize blocks.
 A=rec['target_to_U_path'][:4];B=rec['target_to_U_path'][4:]
 need(all(U[i][j] for i in A for j in B),'U ordered join cut')
 need(canon(reorder(U,A))==4 and canon(reorder(U,B))==81,'U component classes')
 need(4!=5 and rec['V_address'][0]!=rec['V_address'][1],'four distinct vertices')
 count+=1
base=R['base_incidence'];rd={};cd={}
for e in base:rd[e['source_class']]=rd.get(e['source_class'],0)+1;cd[e['root_class']]=cd.get(e['root_class'],0)+1
pred=[e['tile'] for e in base if rd[e['source_class']]==1 or cd[e['root_class']]==1]
need(pred==R['all_boundary_bridge_tiles'],'general bridge criterion matches BFS')
print(json.dumps({'verified_original_global_nonbridges':count,'expanded_boundary_bridges':len(pred),'source_leaf_bridge_tiles':[e['tile'] for e in base if rd[e['source_class']]==1],'target_only_boundary_leaf_tiles':[e['tile'] for e in base if cd[e['root_class']]==1 and rd[e['source_class']]>1]},sort_keys=True))

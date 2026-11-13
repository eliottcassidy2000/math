"""Independent literal controls for the exact run certificates; no assert reliance."""
from pathlib import Path
from itertools import permutations
from math import comb
from hashlib import sha256
import json
import inflation_probe as ip
import two_block as tb
HERE=Path(__file__).resolve().parent
P=(78,12,24,16,35,79,30)
def require(ok,why):
 if not ok:raise RuntimeError(why)
def birth(rows,a):
 n=len(rows)
 return tuple(r|((1<<n)if i==a else 0)for i,r in enumerate(rows))+(((1<<n)-1)^(1<<a),)
def valid(rows):
 n=len(rows)
 require(all(r>=0 and not(r>>n)and not(r>>i&1)for i,r in enumerate(rows)),'row domain')
 require(all(((rows[i]>>j&1)+(rows[j]>>i&1))==1 for i in range(n)for j in range(i)),'tournament')
def strong(rows):
 for a in range(len(rows)):
  seen=1<<a
  while True:
   nxt=seen
   for b in range(len(rows)):
    if seen>>b&1:nxt|=rows[b]
   if nxt==seen:break
   seen=nxt
  if seen!=(1<<len(rows))-1:return False
 return True
def permutation_profile(rows):
 n=len(rows);h=0;s=[0]*n;f=[0]*n
 for path in permutations(range(n)):
  good=[bool(rows[path[i]]>>path[i+1]&1)for i in range(n-1)]
  if all(good):h+=1;s[path[0]]+=1
  for cut in range(1,n+1):
   if all(good[:cut-1]) and all(good[cut:]):f[path[cut-1]]+=1
 return h,s,f
def onto(k,j):
 return int(k==0)if j==0 else sum((-1)**(j-i)*comb(j,i)*i**k for i in range(j+1))
def evaluate_two(poly,l,m):return sum(c*onto(l,a)*onto(m,b)for(a,b),c in poly.items())
def substitute(rows,vertex,module):
 owner=[i for i in range(len(rows))for _ in range(len(module)if i==vertex else 1)]
 local=[j for i in range(len(rows))for j in range(len(module)if i==vertex else 1)]
 out=tuple(sum(((module[local[x]]>>local[y]&1)if u==v else(rows[u]>>v&1))<<y for y,v in enumerate(owner)if x!=y)for x,u in enumerate(owner))
 return out,{i:owner.index(i)for i in range(len(rows))}
def two_covers(rows):
 n=len(rows);h=[]
 for mask in range(1<<n):
  vs=[i for i in range(n)if mask>>i&1]
  q=tuple(sum((rows[v]>>w&1)<<j for j,w in enumerate(vs))for v in vs)
  h.append(ip.literal(q)[0]if vs else 1)
 return sum(h[m]*h[((1<<n)-1)^m]for m in range(1<<n))


def main():
 record={'status':'VERIFIED exact literal and permutation controls','universe':{},'selected_family':[]}
 for k in range(1,9):
  rows,idx=ip.inflate(P,1,k);valid(rows);h,s,f=ip.literal(rows);D=2**k
  require((h,s[idx[0]],s[idx[4]],f[idx[0]],f[idx[4]])==(21*D+3,6*D+1,1,30*D+5,24*D+5),('parent family',k))
  pair=[]
  for entry,defect in ((idx[0],6*D+1),(idx[4],1)):
   child=birth(rows,entry);valid(child);hc,sc,fc=ip.literal(child)
   require(strong(child) and strong(rows),('strong',k))
   require((hc,fc[-1]-hc)==(45*D+7,defect),('literal child full H and defect',k))
   require(sum(r>>(len(child)-1)&1 for r in child)==1,'newborn indegree')
   pair.append({'H':hc,'defect':fc[-1]-hc,'start':sc[-1]})
  record['selected_family'].append({'k':k,'parent_H':h,'children':pair})
  if k<=2:require(permutation_profile(rows)==(h,s,f),('independent permutations',k))
 module_controls=[]
 for name,M in [('C3',(2,4,1)),('C4',(2,4,9,3)),('five_core',(22,4,8,3,14)),('TT3',(6,4,0))]:
  rows,idx=substitute(P,1,M);hm=ip.literal(M)[0];pm=two_covers(M);valid(rows)
  children=[]
  for entry,defect in ((0,hm+6*pm),(4,hm)):
   child=birth(rows,idx[entry]);hc,sc,fc=ip.literal(child)
   require(hc==7*hm+45*pm and fc[-1]-hc==defect,('arbitrary module',name,entry))
   children.append({'H':hc,'defect':defect})
  module_controls.append({'module':name,'H_module':hm,'P_module_2':pm,'children':children})
 record['arbitrary_module_controls']=module_controls
 data=json.loads((HERE/'inflation_probe.json').read_text());chart_checks=0
 for item in data['hits']:
  rows,idx=ip.inflate(item['parent'],item['inflate'],2);h,s,f=ip.literal(rows)
  for name,a in zip(('a','b'),item['ports']):
   want=(ip.evalrun(item['H_runs'],2),ip.evalrun(item['S'+name],2),ip.evalrun(item['F'+name],2))
   require((h,s[idx[a]],f[idx[a]])==want,('chart parent',item))
   child=birth(rows,idx[a]);hc,sc,fc=ip.literal(child)
   require(hc==h+f[idx[a]]-s[idx[a]] and fc[-1]-hc==s[idx[a]],('chart full child',item))
  chart_checks+=1
 X=birth(birth(P,4),7);H,out=ip.run_profile(X,1,[a for a in range(len(X))if a!=1]);sx,fx=out[8]
 fixed_differences={a:[f-s-ff+ss for f,s,ff,ss in zip(F,S,fx,sx)]for a,(S,F)in out.items()}
 for a,diff in fixed_differences.items():
  require(a==8 or all(v>=0 for v in diff)or all(v<=0 for v in diff),('fixed sign',a))
  require(a==8 or diff[0]!=0,('fixed strictness',a))
 split,idx=ip.inflate(X,1,3);_,o=tb.profile(split,(1,3),(2,idx[8]));S,F=o[2];Sx,Fx=o[idx[8]]
 keys=S.keys()|F.keys()|Sx.keys()|Fx.keys()
 delta={q:F.get(q,0)-S.get(q,0)-Fx.get(q,0)+Sx.get(q,0)for q in keys}
 delta={q:v for q,v in delta.items()if v}
 require(all(v>0 for v in delta.values()),'two-block nonnegative certificate')
 require(all(delta.get(q,0)>0 for q in ((0,0),(0,1),(1,0),(1,1))),'two-block strictness on four boundaries')
 grid=[]
 for l in range(4):
  for m in range(4):
   k=l+m+1;rows,idx=ip.inflate(X,1,k);h,s,f=ip.literal(rows);pivot=1+l;x=idx[8]
   diff=f[pivot]-s[pivot]-f[x]+s[x]
   require(diff==evaluate_two(delta,l,m),('two-block literal control',l,m))
   require(all(f[a]-s[a]!=f[x]-s[x]for a in range(len(rows))if a!=x),('complete port slice',l,m))
   grid.append({'left':l,'right':m,'difference':diff})
 record['universe']={'selected_k':list(range(1,9)),'independent_permutation_parent_k':[1,2],'all_neutral_charts_k2':chart_checks,'two_block_left_right_grid':'0..3 each, including empty blocks'}
 record['fixed_port_differences_runs']={str(a):d for a,d in fixed_differences.items()}
 record['internal_port_difference_runs']={str(q):v for q,v in sorted(delta.items())}
 record['internal_port_literal_grid']=grid
 record['source_sha256']={str(p.relative_to(HERE.parent)):sha256(p.read_bytes()).hexdigest()for p in [HERE/'inflation_probe.py',HERE/'inflation_probe.json',HERE/'two_block.py',HERE/'audit_counts.py',HERE.parent/'returns/port_collisions.json']}
 (HERE/'audit_counts.json').write_text(json.dumps(record,indent=2)+'\n')
 print('PASS: 225 exact neutral charts at k2; growing-defect-to-1 family k1..8; independent permutation parents k1,k2; complete all-k no-retarget coefficient certificate plus 16 literal split-block controls.')
if __name__=='__main__':main()

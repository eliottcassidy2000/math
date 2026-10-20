"""Independent verification of the four seeds and four-vertex extension."""
from pathlib import Path
import json,math
ROOT=Path(__file__).resolve().parent.parent

def from_mask(n,t):
 a=[[False]*n for _ in range(n)]
 for i in range(n-1):a[i][i+1]=True
 k=0
 for j in range(n-1,1,-1):
  for i in range(j-1):
   a[j][i]=bool(t>>k&1);a[i][j]=not a[j][i];k+=1
 for i in range(n-1):a[i+1][i]=False
 return a

def path_edges(p):return set(zip(p,p[1:]))
def reverse(a,p):
 b=[r[:] for r in a]
 for x,y in zip(p,p[1:]):b[x][y],b[y][x]=b[y][x],b[x][y]
 return b

def c3(a):return math.comb(len(a),3)-sum(math.comb(sum(row),2) for row in a)
def blue(a,p):
 # Reflection-fixed tiling is equivalent to path-reversal anti-isomorphism.
 n=len(p)
 return all(a[p[i]][p[j]]==a[p[n-1-j]][p[n-1-i]] for i in range(n) for j in range(i+1,n))

def verify(a,p,q):
 n=len(a);assert sorted(p)==sorted(q)==list(range(n))
 assert all(not a[i][i] for i in range(n))
 assert all(a[i][j]!=a[j][i] for i in range(n) for j in range(i+1,n))
 ep,eq=path_edges(p),path_edges(q)
 assert ep.isdisjoint(eq)
 assert all(a[x][y] for x,y in ep|eq)
 assert len({p[0],p[-1],q[0],q[-1]})==4
 ds=list(map(sum,a));alpha=ds[p[0]]-ds[p[-1]]-1;beta=ds[q[0]]-ds[q[-1]]-1
 assert alpha>0 and beta>0 and alpha!=beta and alpha%2==beta%2==(n-1)%2
 ap,aq=reverse(a,p),reverse(a,q);apq=reverse(ap,q)
 assert apq==reverse(aq,p)
 assert all(ap[x][y] for x,y in eq) and all(aq[x][y] for x,y in ep)
 vals=list(map(c3,[a,ap,apq,aq]));c=vals[0]
 assert vals==[c,c+alpha,c+alpha+beta,c+beta]
 assert len(set(vals))==4
 assert not any([blue(a,p),blue(a,q),blue(ap,q),blue(aq,p)])
 return {'n':n,'alpha':alpha,'beta':beta,'four_c3_values':vals,'four_moves_black':True}

def choose_attachment(p,q):
 forbidden={p[0],p[-1],q[0],q[-1]}
 for a,b in zip(p,p[1:]):
  for c,d in zip(q,q[1:]):
   if len({a,b,c,d})==4 and not (forbidden&{a,b,c,d}):return a,b,c,d
 raise AssertionError('No legal attachment arcs')

def extend(a,p,q,attachment):
 n=len(a);A,B,C,D=attachment
 assert (A,B) in path_edges(p) and (C,D) in path_edges(q)
 assert len(set(attachment))==4
 assert not set(attachment)&{p[0],p[-1],q[0],q[-1]}
 # Every otherwise unspecified old-new edge initially points old -> new.
 # Explicit requirements are installed through a checked assignment map.
 required={}
 def require(x,y):
  assert x!=y and (y,x) not in required
  required[x,y]=True
 x0,x1,x2,x3=range(n,n+4)
 for x,y in [(x0,x1),(x1,x2),(x2,x3),(x2,x0),(x0,x3),(x3,x1),(A,x0),(x3,B),(C,x2),(x1,D)]:require(x,y)
 for start in [p[0],q[0]]:
  for x in range(n,n+4):require(start,x)
 for end in [p[-1],q[-1]]:
  for x in range(n,n+4):require(x,end)
 z=[[False]*(n+4) for _ in range(n+4)]
 for i in range(n+4):
  for j in range(i+1,n+4):
   value=a[i][j] if j<n else True
   if (i,j) in required:value=True
   if (j,i) in required:value=False
   z[i][j]=value;z[j][i]=not value
 pp=p[:p.index(A)+1]+[x0,x1,x2,x3]+p[p.index(A)+1:]
 qq=q[:q.index(C)+1]+[x2,x0,x3,x1]+q[q.index(C)+1:]
 return z,pp,qq,(x1,x2,x0,x3)

if __name__ == '__main__':
 report=[]
 for base in [8,9,10,11]:
  seed=json.loads((ROOT/f'score_square_seed_n{base}.json').read_text())
  a=from_mask(base,seed['t']);p=seed['P'];q=seed['Q']
  row=verify(a,p,q);assert row['alpha']==seed['alpha'] and row['beta']==seed['beta'];assert list(map(sum,a))==seed['scores']
  att=choose_attachment(p,q)
  rec={'base_order':base,'initial_attachment':list(att),'checks':[row]}
  for step in range(25):
   a,p,q,att=extend(a,p,q,att);next_row=verify(a,p,q)
   assert next_row['alpha']==row['alpha']+4 and next_row['beta']==row['beta']+4
   rec['checks'].append(next_row);row=next_row
  report.append(rec)
 (ROOT/'audit/score_construction_audit.json').write_text(json.dumps({'seed_orders':[8,9,10,11],'extensions_per_seed':25,'verified_orders':sum(len(r['checks']) for r in report),'all_checks':'passed','families':report},indent=2)+'\n')
 print('All four seeds and 25 successive extensions per seed passed, covering every order 8 through 111.')

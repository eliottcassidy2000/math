import itertools,json,pathlib,sys
from fractions import Fraction as F
P=pathlib.Path(__file__).resolve().parents[1];sys.path.insert(0,str(P/'session203-positive-profile-defects'))
from counting import labelled,nine,relabel,literal

def glue(a,g,r):
 n=len(a);rows=list(a)+[0]*(len(g)-1)
 def node(i):return 0 if i==0 else n+i-1
 for i in range(len(g)):
  for j in range(len(g)):
   if g[i]>>j&1:rows[node(i)]|=1<<node(j)
 for i in range(1,n):
  for j in range(1,len(g)):rows[i]|=1<<node(j)
 return relabel(tuple(rows),[node(r)]+[i for i in range(len(rows))if i!=node(r)])
def solve(M,y):
 a=[list(map(F,row))+[F(v)] for row,v in zip(M,y)];n=len(a)
 for k in range(n):
  t=next(t for t in range(k,n)if a[t][k]);a[k],a[t]=a[t],a[k];z=a[k][k];a[k]=[v/z for v in a[k]]
  for i in range(n):
   if i!=k:
    z=a[i][k];a[i]=[u-z*v for u,v in zip(a[i],a[k])]
 return [row[-1] for row in a]
def det(M):
 a=[list(map(F,r))for r in M];d=F(1)
 for k in range(len(a)):
  t=next(t for t in range(k,len(a)) if a[t][k]);a[t],a[k]=a[k],a[t]
  if t!=k:d=-d
  x=a[k][k];d*=x
  for i in range(k+1,len(a)):
   z=a[i][k]/x;a[i]=[u-z*v for u,v in zip(a[i],a[k])]
 return d
bases=[(0,),(2,0),(4,1,2),(6,0,2)];ps=[nine(a)for a in bases];c1=[0,2,4,5];c2=[6,3,8,7];M1=[[p[i]for i in c1]for p in ps];M2=[[p[i]for i in c2]for p in ps]
contexts=checks=anchormarks=0
for m in range(1,4):
 for g in labelled(m):
  hh=1 if m==1 else nine(tuple(x>>1 for x in g[1:]))[0]
  for r in range(m):
   contexts+=1;ys=[nine(glue(a,g,r))for a in bases];H=solve(M1,[p[0]for p in ys]);L=solve(M2,[p[3]for p in ys]);c=ys[1][1]
   assert H[0]==hh,(g,r,H,hh)
   if r==0:assert L[0]==0 and L[2]==0;anchormarks+=1
   for n in range(1,5):
    for a in labelled(n):
     p=nine(a);q=nine(glue(a,g,r));checks+=1
     assert q[0]==sum(H[j]*p[i]for j,i in enumerate(c1))
     assert q[3]==sum(L[j]*p[i]for j,i in enumerate(c2))
     if n>1:assert q[1]==c*p[1]
span=json.loads((P/'session203-digital-context-search/SPAN_WITNESS.json').read_text());rows=[]
for word,p in span['word_states']:
 q=nine(literal(word));assert tuple(p)==q and q[1]==q[5]==1;rows.append([q[i]for i in span['columns']])
D=det(rows);assert D==-1152
for n in range(9):
 vals=[]
 for word in itertools.product(('O','Z'),repeat=n):
  from compiler import move,SEED
  p=SEED
  for op in word:p=move(p,op)
  assert p[1]==p[5]==1 and p[7]==2**(n+1) and p[2]==p[4]%2;vals.append(p[4])
 assert sorted(vals)==list(range(2**n+1,2**(n+1)+1))
print(json.dumps(dict(status='PASS',contexts=contexts,anchor_mark_contexts=anchormarks,literal_input_context_checks=checks,span_determinant=str(D),scope='Every labeled context through m3 and each mark; all labeled old inputs through n4; literal eight-state span; binary-register words through length8'),indent=2))

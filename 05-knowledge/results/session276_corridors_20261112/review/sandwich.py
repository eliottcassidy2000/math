import contextlib,io,json
from pathlib import Path
with contextlib.redirect_stdout(io.StringIO()):
 from independent import graph,subset_count

def one(A):
 b=list(map(int,bin(A)[2:])) if A else [];n=len(b)+1;G=[[False]*n for _ in range(n)]
 for i,x in enumerate(b,1):
  G[0][i]=bool(x);G[i][0]=not x
  for j in range(i+1,n):G[i][j]=True
 return G

def converse(G):return [list(x) for x in zip(*G)]
def fuse(L,R):
 n=len(L);m=len(R);V=[[False]*(n+m-1) for _ in range(n+m-1)]
 for G,ids in ((L,list(range(n))),(R,[0]+list(range(n,n+m-1)))):
  for i,a in enumerate(ids):
   for j,b in enumerate(ids):V[a][b]=G[i][j]
 for i in range(1,n):
  for j in range(n,n+m-1):V[i][j]=True
 return V

def f(A):
 if A==0:return 1
 a,r=divmod(A,2);return f(a)+(a+(a&1) if r==0 else 0)
count=0
for A in range(64):assert subset_count(one(A))==f(A);count+=1
for A in range(4):
 for B in range(4):
  for j,C in enumerate([graph(5,{0,1}),graph(5,{0,2})]):
   value=subset_count(fuse(fuse(one(A),C),converse(one(B))))
   expected=41+2*j+15*(f(A)+f(B)-2)+63*B+42*A+120*A*B+9*(B&1)+6*(A&1)
   assert value==expected,(A,B,j,value,expected);count+=1
R=json.loads((Path(__file__).resolve().parents[1]/'digits/REPLAY.json').read_text())
for name in ['C9_15381','C9_15383']:
 rows=R[name];G=[[bool(r>>j&1) for j in range(len(rows))] for r in rows]
 expected=int(name.split('_')[1]);assert subset_count(G)==expected;count+=1
 for A,B in [(1,0),(0,1),(1,1)]:
  value=subset_count(fuse(fuse(one(A),G),converse(one(B))))
  pred=expected+3357*(f(A)+f(B)-2)+20732*B+25915*A+50004*A*B+1492*(B&1)+1865*(A&1)
  assert value==pred,(name,A,B,value,pred);count+=1
print(json.dumps({'status':'PASS','independent_literal_checks':count,'C9_source':'15381/15383 outgoing-size-four pair'},indent=2))

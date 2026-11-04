import itertools,json
from pathlib import Path

def tables(rows):
 n=len(rows);N=1<<n;D=[[0]*n for _ in range(N)];S=[[0]*n for _ in range(N)];H=[0]*N;H[0]=1
 for mask in range(1,N):
  for v in range(n):
   if mask>>v&1:
    sub=mask^(1<<v)
    D[mask][v]=sum(D[sub][u] for u in range(n) if sub>>u&1 and rows[u]>>v&1) if sub else 1
    S[mask][v]=sum(S[sub][u] for u in range(n) if sub>>u&1 and rows[v]>>u&1) if sub else 1
  H[mask]=sum(D[mask])
 return D,S,H

def fields(rows,W):
 D,S,H=tables(rows);V=(1<<len(rows))-1
 F=sum(sum(D[m][u] for u in W)*H[V^m] for m in range(1,V+1))
 J=sum(sum(D[m][u] for u in W)*sum(S[V^m][v] for v in W) for m in range(1,V+1))
 return H[V],F,J

def attach(a,b,W):
 na,nb=len(a),len(b);r=na+nb;out=[0]*(r+1)
 for i in range(na):out[i]=a[i]|(((1<<nb)-1)<<na)
 for j in range(nb):out[na+j]=b[j]<<na
 for i in range(r):
  if i>=na and i-na in W:out[i]|=1<<r
  else:out[r]|=1<<i
 return out,r

def check(ok,msg):
 if not ok:raise RuntimeError(msg)

count=0
for n in range(1,5):
 edges=list(itertools.combinations(range(n),2))
 for bits in range(1<<len(edges)):
  b=[0]*n
  for k,(i,j) in enumerate(edges):
   if bits>>k&1:b[i]|=1<<j
   else:b[j]|=1<<i
  for mask in range(1,1<<n):
   W=[v for v in range(n) if mask>>v&1];h,F,J=fields(b,W)
   for a in ([0],[2,0],[2,4,1]):
    _,_,ha=tables(a);full=(1<<len(a))-1;p2=sum(ha[m]*ha[full^m] for m in range(full+1))
    rows,r=attach(a,b,W);D,S,H=tables(rows);v=(1<<len(rows))-1
    pred=ha[-1]*(h-J)+p2*F
    check(H[-1]==pred,('count',b,W,a,H[-1],pred))
    er=sum(D[m][r]*H[v^m] for m in range(v+1))
    check(er-H[-1]==ha[-1]*J,('defect',b,W,a))
    if n<=3 and len(a)<=2:
     brute=sum(all(rows[p[i]]>>p[i+1]&1 for i in range(len(p)-1)) for p in itertools.permutations(range(len(rows))))
     check(brute==pred,'permutation mismatch')
    count+=1
print(json.dumps({'attachment_and_defect_cases':count,'status':'PASS','counts':'exact Python integers','permutation_controls':'all suffixes n<=3 and prefix TT1/TT2'}))

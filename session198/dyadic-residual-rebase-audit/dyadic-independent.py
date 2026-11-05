"""Raw-coordinate and independently reconstructed all-edge dyadic audit."""
from pathlib import Path
import json
HERE=Path(__file__).resolve().parent

def step(s,l):
 D,E,U,d=s
 if l=='O':return D,0,2*U,d
 if l=='Z':return D-E+U,d,2*U-d,d
 return D-E+U,D-E,D+U,D

bank={(1,0,2,1)};queue=list(bank);cap=2000
for s in queue:
 for l in 'OZJ':
  t=step(s,l)
  if t[2]<=cap and t not in bank:bank.add(t);queue.append(t)
assert len(bank)==32548
p2=lambda n:n>0 and not n&(n-1)
dyadic={s for s in bank if s[2]>s[0]and s[2]%s[0]==0 and p2(s[2]//s[0])}
expected=set()
for D,E,U,d in bank|{(1,1,1,1)}:
 if D!=U:continue
 u=2*U
 while u<=cap:expected.add((D,0,u,d));u*=2
assert dyadic==expected

def append(T,c,g):
 n=len(T);A=T.copy()+[0]
 for u in range(n):
  if g!='O'and u==c:A[n]|=1<<u
  else:A[u]|=1<<n
 return A,n if g=='J'else c

def word(T,c,w):
 for g in w:T,c=append(T,c,g)
 return T,c

def inventory(T):
 n=len(T);size=1<<n;F=size-1;reverse=[F^row^(1<<i)for i,row in enumerate(T)]
 P=[[0]*n for _ in range(size)];Q=[[0]*n for _ in range(size)];C=[[0]*n for _ in range(size)]
 for i in range(n):P[1<<i][i]=Q[1<<i][i]=1
 C[1][0]=1
 for S in range(1,size):
  for u in range(n):
   for table,rows in [(P,T),(Q,reverse)]:
    if table[S][u]:
     rem=rows[u]&~S
     while rem:
      bit=rem&-rem;rem-=bit;v=bit.bit_length()-1;table[S|bit][v]+=table[S][u]
   if C[S][u]:
    rem=T[u]&~S
    while rem:
     bit=rem&-rem;rem-=bit;v=bit.bit_length()-1;C[S|bit][v]+=C[S][u]
 H=[sum(row)for row in P];H[0]=1
 U=[sum(H[F^S]*Q[S][r]for S in range(1,size))for r in range(n)]
 return H[-1],P[-1],Q[-1],U,[H[F^(1<<r)]for r in range(n)],sum(C[-1][v]for v in range(n)if T[v]&1)

def induced(T,V):return[sum(((T[u]>>v)&1)<<j for j,v in enumerate(V))for u in V]
def check_map(T,X,V):assert all(bool(X[i]>>j&1)==bool(T[V[i]]>>V[j]&1)for i in range(len(X))for j in range(len(X)))
def remove_sink(T,kept):
 mask=sum(1<<v for v in kept);sinks=[v for v in kept if not T[v]&mask];assert len(sinks)==1
 z=sinks[0];kept.remove(z);return z

small_tail=small_ratio=0
for n in range(1,6):
 edges=[(i,j)for i in range(n)for j in range(i+1,n)]
 for mask in range(1<<len(edges)):
  T=[0]*n
  for k,(i,j)in enumerate(edges):
   if mask>>k&1:T[i]|=1<<j
   else:T[j]|=1<<i
  D,E,S,U,d,C=inventory(T)
  for r in range(n):
   if S[r]!=1:continue
   if U[r]>D and U[r]%D==0 and p2(U[r]//D):
    k=(U[r]//D).bit_length()-1;kept=list(range(n));tail=[remove_sink(T,kept)for _ in range(k)];assert r in kept
    W=induced(T,kept);c=kept.index(r);dd,ee,ss,uu,_,_=inventory(W);assert dd==uu[c]and ss[c]==1
    X,c=word(W,c,'O'*k);check_map(T,X,kept+tail[::-1]);small_tail+=1
   if n>1 and U[r]==D and D%E[r]==0 and p2(D//E[r]-1)and D//E[r]-1>=4:
    k=(D//E[r]-1).bit_length()-2;ports=[v for v in range(n)if T[r]>>v&1];assert len(ports)==1;c=ports[0]
    kept=[v for v in range(n)if v!=r];z=remove_sink(T,kept);tail=[remove_sink(T,kept)for _ in range(k)];assert c in kept
    W=induced(T,kept);cc=kept.index(c);dd,ee,ss,uu,_,_=inventory(W);assert dd==E[r]and dd==uu[cc]and ss[cc]==1
    X,cc=word(W,cc,'O'*(k+1)+'J');check_map(T,X,kept+tail[::-1]+[z,r]);small_ratio+=1

pairs=[];cores=[([0],0)]+[word([2,0],0,w+'OJ')for w in['','O','Z','J']]
for T,c in cores:
 h,E,S,U,d,C=inventory(T);assert h==U[c]and S[c]==1
 for k in range(1,4):
  A,a=word(T,c,'O'*(k+1)+'JO'+'Z'*k+'J');B,b=word(T,c,'O'*(2*k+3)+'J')
  ha,ea,sa,ua,da,ca=inventory(A);hb,eb,sb,ub,db,cb=inventory(B)
  expect=((1<<(2*k+3))+1)*h
  assert ha==hb==expect and ca==cb==1 and len(A)==len(B)==len(T)+2*k+4
  assert not any(sa[r]==1 and ua[r]==ha for r in range(len(A)))and sb[b]==1 and ub[b]==hb
  pairs.append({'core_H':h,'k':k,'H':ha,'order':len(A)})
out={'status':'PASS','raw_source_profiles':len(bank),'exact_dyadic_profiles':len(dyadic),'all_edge_tail_reconstructions':small_tail,'all_edge_balanced_ratio_reconstructions':small_ratio,'independent_graph_pairs':pairs}
(HERE/'DYADIC_INDEPENDENT.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n');print(json.dumps(out,sort_keys=True))

import itertools,json
from pathlib import Path

def ck(x,msg):
 if not x: raise RuntimeError(msg)
def tt(n):return [sum(1<<j for j in range(i+1,n)) for i in range(n)]
def join(a,b):return [x|(((1<<len(b))-1)<<len(a)) for x in a]+[x<<len(a) for x in b]
def tables(r):
 n=len(r);end=[[0]*n for _ in range(1<<n)];H=[0]*(1<<n);H[0]=1
 for v in range(n):end[1<<v][v]=1
 for S in range(1,1<<n):
  for v in range(n):
   if S>>v&1:
    T=S^(1<<v)
    if T:end[S][v]=sum(end[T][u] for u in range(n) if r[u]>>v&1)
  H[S]=sum(end[S])
 return H,end

def sig(r,W):
 n=len(r);H,E=tables(r);rev=[sum(1<<j for j in range(n) if r[j]>>i&1) for i in range(n)];_,S=tables(rev);U=(1<<n)-1
 e=[sum(z[i] for i in W) for z in E];s=[sum(z[i] for i in W) for z in S]
 return [H[U],sum(H[x]*H[U^x] for x in range(U+1)),sum(e[x]*H[U^x] for x in range(U+1)),sum(H[x]*s[U^x] for x in range(U+1)),sum(e[x]*s[U^x] for x in range(U+1)),e[U],s[U]]
def attach(a,b,W):
 r=join(a,b);n=len(r);r.append((1<<len(a))-1)
 for v in range(len(b)):
  i=len(a)+v
  if v in W:r[i]|=1<<n
  else:r[n]|=1<<i
 return r

def strong(r):
 for v in range(len(r)):
  seen=1<<v
  while True:
   nex=seen
   for i in range(len(r)):
    if seen>>i&1:nex|=r[i]
   if nex==seen:break
   seen=nex
  if seen!=(1<<len(r))-1:return False
 return True
b0=tt(3);b2=join(tt(2),[2,4,1]);b4=tt(5);blocks=[(b0,{0,1,2}),(b2,{0,1,2}),(b4,{0,2,4})]
ck(sig(tt(2),{0,1})==[1,4,3,3,2,1,1] and sig(join(tt(1),[2,4,1]),{0,1})==[3,24,9,15,4,1,3],'short base')
ck([sig(*x) for x in blocks]==[[1,8,7,7,6,1,1],[3,48,21,39,16,1,3],[1,32,21,21,12,1,1]],'base')
count=0;full=0
for n in range(4):
 pairs=list(itertools.combinations(range(n),2))
 for code in range(1<<len(pairs)):
  L=[0]*n
  for k,(i,j) in enumerate(pairs):
   if code>>k&1:L[i]|=1<<j
   else:L[j]|=1<<i
  for mask in range(1<<n):
   W={i for i in range(n) if mask>>i&1}
   h,p,f,g,j,e,s=sig(L,W) if n else [1,1,0,0,0,0,0]
   out=[(join(L,b),W|{n+x for x in V}) for b,V in blocks]
   ss=[sig(*x) for x in out];c=[x[0]-x[4] for x in ss]
   ck(ss[1][2]-3*ss[0][2]==0 and c[1]-3*c[0]==2*h,'lift2')
   ck(ss[2][2]-3*ss[0][2]==-2*f and c[2]-3*c[0]==4*h+2*j+2*e,'lift4')
   for A in [tt(1),tt(2),[2,4,1]]:
    HA=tables(A)[0][-1];PA=sig(A,set())[1]
    rr=[attach(A,*x) for x in out];hh=[tables(x)[0][-1] for x in rr]
    ck(all(strong(x) for x in rr),'strong')
    ck(hh[1]-3*hh[0]==2*HA*h,'full2')
    ck(hh[2]-3*hh[0]==HA*(4*h+2*j+2*e)-2*PA*f,'full4')
    full+=3
   count+=1
result={'context_signatures':count,'literal_full_graphs':full,'outer_types':['TT1','TT2','C3'],'max_context_order':3,'base_signatures':[sig(*x) for x in blocks],'status':'PASS'}
Path(__file__).with_name('receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))

import itertools,json,math,pathlib
OUT=pathlib.Path(__file__).resolve().parent

def need(x,msg):
 if not x:raise RuntimeError(msg)
def tournaments(n):
 for mask in range(1<<(n*(n-1)//2)):
  rows=[0]*n;k=0
  for i in range(n):
   for j in range(i+1,n):
    if mask>>k&1:rows[i]|=1<<j
    else:rows[j]|=1<<i
    k+=1
  yield rows

def h(rows):
 n=len(rows);d=[[0]*n for _ in range(1<<n)]
 for v in range(n):d[1<<v][v]=1
 for mask in range(1,1<<n):
  for v,c in enumerate(d[mask]):
   if c:
    opts=rows[v]&~mask
    while opts:
     b=opts&-opts;opts-=b;d[mask|b][b.bit_length()-1]+=c
 return sum(d[-1])
def clone(a,x,t):
 old=[v for v in range(len(a)) if v!=x];labels=old+[x]*t;n=len(labels);b=[0]*n
 for i in range(n):
  for j in range(i+1,n):
   forward=(labels[i]==labels[j]) or bool(a[labels[i]]>>labels[j]&1)
   if forward:b[i]|=1<<j
   else:b[j]|=1<<i
 return b

def coeffs(a,x):
 U=[v for v in range(len(a)) if v!=x];n=len(U);A=[]
 for r in range(1,n//2+2):
  count=0
  for lengths in itertools.product(range(n+1),repeat=r+1):
   if sum(lengths)!=n or any(k<2 for k in lengths[1:-1]):continue
   for perm in itertools.permutations(U):
    segs=[];p=0
    for size in lengths:segs.append(perm[p:p+size]);p+=size
    if any(any(not(a[u]>>v&1) for u,v in zip(seg,seg[1:])) for seg in segs):continue
    if segs[0] and not(a[segs[0][-1]]>>x&1):continue
    if segs[-1] and not(a[x]>>segs[-1][0]&1):continue
    if any(not(a[x]>>seg[0]&1) or not(a[seg[-1]]>>x&1) for seg in segs[1:-1]):continue
    count+=1
  A.append(count)
 return A

def surj(t,r):return sum((-1)**(r-i)*math.comb(r,i)*i**t for i in range(r+1))
def main():
 cases=0;clones=0
 for n in range(1,5):
  for a in tournaments(n):
   for x in range(n):
    A=coeffs(a,x);need(A[0]==h(a),'A1')
    for r,z in enumerate(A,1):need(z%math.factorial(r-1)==0,'factorial divisibility')
    for t in range(1,6):
     count=h(clone(a,x,t));need(count==sum(z*surj(t,r) for r,z in enumerate(A,1)),'independent clone spectrum');need(count%2==1,'parity');clones+=1
    cases+=1
 a=[20,1,2,7,14];A=coeffs(a,0);need(A==[9,10,2],'mixed coefficient')
 for t in range(1,8):need(h(clone(a,0,t))==2*3**t+4*2**t-5,'mixed count')
 pure=[26,24,67,116,36,71,19];PA=coeffs(pure,6);need(PA==[145,150,50,0],'pure ternary coefficients')
 pure_counts=[]
 for t in range(1,7):
  count=h(clone(pure,6,t));need(count==50*3**t-5,'pure ternary count');pure_counts.append(count)
 out={'pure_ternary_A':PA,'pure_ternary_counts_t1_to6':pure_counts,'status':'PASS','method':'Independent outside-permutation and segment-length enumeration, no imports from author','coefficient_cases':cases,'clone_count_checks':clones,'mixed_A':A,'mixed_t':list(range(1,8))}
 (OUT/'independent.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()

import json

def need(x,msg):
 if not x:raise RuntimeError(msg)
def f(word):
 n=0;h=1
 for bit in word:
  if bit=='0':h+=n+(n&1)
  n=2*n+(bit=='1')
 return n,h

def roots(a,k):
 # Retain all four roots and extend by one top bit; no single-root lift formula.
 out=[1,3,5,7]
 for j in range(3,k):out=[x+t*(1<<j)for x in out for t in(0,1)if(x+t*(1<<j))**2%(1<<(j+1))==a%(1<<(j+1))]
 need(len(out)==4 and all(x*x%(1<<k)==a%(1<<k)for x in out),'four exact square roots')
 return out

def word_for(k,D):
 M=1<<k;w='1'*k+'0'*k
 _,h=f(w);j=((h-D)//2)%(M//2)
 # Reverse the loop order, and reverse each two-loop pair, relative to the source.
 for i in reversed(range(k-1)):
  if j>>i&1:
   w+='11'+'0'*k
   if i:w+='1'*(i+1)+'0'*k
 return w+'0'*k

def physical_profile(bits):
 # Partition the transitive deletion into the words before and after its apex.
 a={(0,0):1}
 for incoming in bits:
  b={}
  for (last,first),v in a.items():
   key=(1 if incoming else 2,first);b[key]=b.get(key,0)+v
   key=(last,first or(2 if incoming else 1));b[key]=b.get(key,0)+v
  a=b
 h=sum(v for(last,first),v in a.items()if last in(0,1)and first in(0,1))
 R=1+sum(1<<(len(bits)-1-i)for i,x in enumerate(bits)if not x)
 return h,int(bool(bits[-1])),int(not bits[0]),R,1

def source_graph(bits):
 n=len(bits);A=[0]*(n+1)
 for i,x in enumerate(bits):
  if x:A[i]|=1<<n
  else:A[n]|=1<<i
  for j in range(i+1,n):A[i]|=1<<j
 return A,n

def diagonal(V,c):
 n=len(V);A=[0]*(2*n+1);b=n;start=n+1
 for i in range(n):
  for j in range(n):
   if V[i]>>j&1:A[i]|=1<<j;A[start+j]|=1<<(start+i)
 for i in range(n):
  if i==c:A[b]|=1<<i
  else:A[i]|=1<<b
 for i in range(n+1):
  for j in range(n):
   u,v=(start+j,i)if j==c and i in(c,b)else(i,start+j);A[u]|=1<<v
 anti=list(range(start,start+n))+[b]+list(range(n))
 need(all(((A[i]>>j)&1)==((A[anti[j]]>>anti[i])&1)for i in range(len(A))for j in range(len(A))),'all anti edges')
 return A

def literal(T):
 n=len(T);full=(1<<n)-1;h=cyc=0
 def dfs(first,last,S):
  nonlocal h,cyc
  if S==full:h+=1;cyc+=first==0 and bool(T[last]&1);return
  rem=T[last]&~S
  while rem:
   b=rem&-rem;rem-=b;dfs(first,b.bit_length()-1,S|b)
 for x in range(n):dfs(x,x,1<<x)
 return h,cyc

small=large=profiles=edges=physical=0
for k in range(3,11):
 M=1<<k
 for r in range(1,M,4):
  t=1 if r%8==5 else 2;u=1<<t;x=min(roots((r-2*u*(u-1))%M,k));D=(x-u)%M
  word=word_for(k,D);n,g=f(word);bits=[0]+[int(b)for b in word[::-1]]+[0]*t
  h,e,s,U,d=physical_profile(bits);profiles+=1
  need((h,e,s,d)==(g+n,0,1,1),'full physical weak source profile')
  need(h%M==D and U%M==u,'physical source target residues')
  H=(h+U)**2+2*U*(U-1);need(H%M==r and H%4==1,'physical diagonal residue')
  need(2*(len(bits)+1)+1<=8*k*k+6*k+9,'actual graph order bound')
  if k<=5:
   V,c=source_graph(bits);T=diagonal(V,c);edges+=len(T)*(len(T)-1)//2
  small+=1
for k in(16,32,64):
 for r in(1,5,(1<<k)-3,(1<<k)-7):
  M=1<<k;t=1 if r%8==5 else 2;u=1<<t;x=min(roots((r-2*u*(u-1))%M,k));D=(x-u)%M
  word=word_for(k,D);h,e,s,U,d=physical_profile([0]+list(map(int,word[::-1]))+[0]*t)
  need(((h+U)**2+2*U*(U-1))%M==r,'large physical assignment count');large+=1
for L in range(2,10):
 word=bin(L-1)[2:]
 for t in(1,2):
  bits=[0]+list(map(int,word[::-1]))+[0]*t;h,e,s,U,d=physical_profile(bits);V,c=source_graph(bits);T=diagonal(V,c)
  H,C=literal(T);need((H,C)==((h+U)**2+2*U*(U-1),2),'literal final Hamilton words and cycles');physical+=1
print(json.dumps(dict(status='PASS',all_residues_through_2pow10=small,independent_large_residues=large,physical_assignment_profiles=profiles,literal_anti_edges=edges,literal_final_graph_counts=physical),sort_keys=True))

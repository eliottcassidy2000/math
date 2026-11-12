"""Finite strong cores and exact all-length generating functions, two apexes."""
import itertools,collections,json,sys
from capped import count

def graph(w):
 m=len(w);n=m+2;g=[0]*n;u=m;v=m+1;g[u]|=1<<v
 for i,l in enumerate(w):
  for j in range(i+1,m):g[i]|=1<<j
  for a,b in ((u,l&1),(v,l>>1)):
   if b:g[a]|=1<<i
   else:g[i]|=1<<a
 return g

def components(g):
 n=len(g);r=[g[i]|1<<i for i in range(n)]
 for k in range(n):
  for i in range(n):
   if r[i]>>k&1:r[i]|=r[k]
 groups=[];unused=set(range(n))
 while unused:
  i=min(unused);C={j for j in unused if r[i]>>j&1 and r[j]>>i&1};groups.append(C);unused-=C
 groups.sort(key=lambda C:sum(bool(g[next(iter(C))]>>next(iter(D))&1) for D in groups if C!=D),reverse=True)
 return groups

def ham(g):
 n=len(g);d=[[0]*n for _ in range(1<<n)]
 for i in range(n):d[1<<i][i]=1
 for S in range(1,1<<n):
  for i,c in enumerate(d[S]):
   q=g[i]&~S
   while c and q:
    b=q&-q;j=b.bit_length()-1;d[S|b][j]+=c;q-=b
 return sum(d[-1])

def one_graph(w):
 m=len(w);g=[0]*(m+1)
 for i,l in enumerate(w):
  g[i]=sum(1<<j for j in range(i+1,m))
  if l:g[m]|=1<<i
  else:g[i]|=1<<m
 return g

def cutoff(N):
 # Every strong n-tournament has 5^(n-1)<=H^3 (inherited Busch).
 b=0
 while 5**(b+1)<=N**3:b+=1
 return b

def core_polynomials(N):
 B=cutoff(N);a={1:collections.Counter({0:1})};c={};samples={}
 for m in range(2,B+1):
  for mid in itertools.product(range(2),repeat=m-2):
   w=(1,)+mid+(0,);h=ham(one_graph(w))
   if h<=N:a.setdefault(h,collections.Counter())[m]+=1
 for m in range(1,B):
  for w in itertools.product(range(4),repeat=m):
   g=graph(w)
   if len(components(g))!=1:continue
   h=count(w)
   if h<=N:
    c.setdefault(h,collections.Counter())[m]+=1
    samples.setdefault(h,[]).append(''.join(map(str,w)))
 return a,c,samples

def product_numerator(N,a):
 p=collections.Counter()
 for h,ap in a.items():
  if N%h or N//h not in a:continue
  for i,x in ap.items():
   for j,y in a[N//h].items():p[i+j]+=x*y
 return p

def coeff(m,p,den):
 from math import comb
 return sum(x*comb(m-i+den-1,den-1) for i,x in p.items() if m>=i)

if __name__=='__main__':
 N=int(sys.argv[1]) if len(sys.argv)>1 else 63
 a,c,words=core_polynomials(N)
 out={}
 for h in range(1,N+1,2):
  p=product_numerator(h,a);q=c.get(h,{})
  out[h]={'one_apex':dict(a.get(h,{})),'two_apex_strong':dict(q),'separate_apex_numerator':dict(p)}
 # independent entire word count to length7, directly DP up to length4
 for m in range(8):
  hist=collections.Counter()
  for w in itertools.product(range(4),repeat=m):
   h=count(w)
   if m<=4 and h!=ham(graph(w)):raise RuntimeError(('DP',w,h))
   if h<=N:hist[h]+=1
  for h,rec in out.items():
   predicted=coeff(m,rec['two_apex_strong'],2)+coeff(m,rec['separate_apex_numerator'],3)
   if predicted!=hist[h]:raise RuntimeError(('GF',m,h,predicted,hist[h]))
 print(json.dumps({'bound':N,'core_chain_bound':cutoff(N),'fibres':out,'controls':{'all_words_lengths':'0..7','independent_subset_DP':'0..4'},'H13strongcores':words.get(13,[])},indent=2))

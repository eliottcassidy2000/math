#!/usr/bin/env python3
import itertools,json,hashlib

def require(c,msg):
 if not c: raise RuntimeError(msg)
def graph(word):
 n=len(word);out=[0]*(n+1)
 for i in range(n):
  for j in range(i+1,n):out[i]|=1<<j
  if word[i]=='1':out[n]|=1<<i
  else:out[i]|=1<<n
 return out

def counts(out):
 n=len(out);dp=[[0]*n for _ in range(1<<n)];h=[0]*(1<<n);h[0]=1
 for v in range(n):dp[1<<v][v]=1
 for mask in range(1,1<<n):
  for v,val in enumerate(dp[mask]):
   if not val:continue
   avail=out[v]&~mask
   while avail:
    bit=avail&-avail;avail-=bit;dp[mask|bit][bit.bit_length()-1]+=val
  h[mask]=sum(dp[mask])
 return h

def stream(word):
 p=l=0;h=1
 for c in word:
  if c=='0':h+=p+l;p*=2;l=0
  else:p=2*p+1;l=1
 return h

def literal(word):
 if not word:return 1
 return int(word[0])+1-int(word[-1])+sum(word[i:i+2]=='01' for i in range(len(word)-1))+sum(2**(i-j-1) for i in range(len(word)) for j in range(i) if word[j]=='1' and word[i]=='0')

def census(B):
 maxn=2+(B-2).bit_length()-1 if B>=3 else 0
 out={1:''}
 for n in range(2,maxn+1):
  for bits in itertools.product('01',repeat=n-2):
   word='1'+''.join(bits)+'0';h=stream(word)
   require(h>=2**(n-2)+2,'core bound')
   if h<=B:out.setdefault(h,word)
 return out,maxn

def labelled(n,code):
 out=[0]*n;k=0
 for i in range(n):
  for j in range(i+1,n):
   if code>>k&1:out[i]|=1<<j
   else:out[j]|=1<<i
   k+=1
 return out

def wrapper(a,b):
 m=len(a);n=len(b);out=a[:]+[x<<m for x in b]+[0]
 for i in range(m):out[i]|=((1<<n)-1)<<m
 for i in range(m,m+n):out[i]|=1<<(m+n)
 out[-1]=(1<<m)-1
 return out

def main():
 checks=0
 for n in range(9):
  for bits in itertools.product('01',repeat=n):
   w=''.join(bits);h=counts(graph(w))[-1]
   require(h==stream(w)==literal(w),'binary formula '+w);checks+=1
 for a in range(1,12):
  for b in range(1,12):require(stream('1'+'0'*a+'1'*b+'0')==(2**a+1)*(2**b+1),'four runs')
 bank,maxn=census(10000)
 profiles={};witness13=None
 for n in range(1,6):
  for code in range(1<<(n*(n-1)//2)):
   g=labelled(n,code);h=counts(g);p=(h[-1],sum(h[s]*h[(1<<n)-1-s] for s in range(1<<n)))
   profiles.setdefault(p,g)
   if p[0]==13:witness13=g
 require(witness13 is not None,'13 witness')
 wrappers=0;seven=None
 for (a,f),A in profiles.items():
  for (b,g),B in profiles.items():
   expected=(f-a)*(g-b)+2*a*b
   actual=counts(wrapper(A,B))[-1]
   require(actual==expected,'asymmetric wrapper');wrappers+=1
   if actual%7==0 and actual%49!=0 and seven is None:seven={'A':A,'B':B,'count':actual,'profiles':[[a,f],[b,g]]}
 require(seven is not None,'odd 7 valuation example')
 insertion_checks=0
 for n in range(1,5):
  for code in range(1<<(n*(n-1)//2)):
   A=labelled(n,code)
   paths=[p for p in itertools.permutations(range(n)) if all(A[p[i]]>>p[i+1]&1 for i in range(n-1))]
   endpoint={}
   for mask in range(1,1<<n):
    vs=[v for v in range(n) if mask>>v&1];beg=[0]*n;end=[0]*n
    for p in itertools.permutations(vs):
     if all(A[p[i]]>>p[i+1]&1 for i in range(len(p)-1)):beg[p[0]]+=1;end[p[-1]]+=1
    endpoint[mask]=(beg,end)
   for b in range(1<<n):
    G=A[:]+[b]
    for v in range(n):
     if not b>>v&1:G[v]|=1<<n
    d=sum(sum(bool(b>>p[i]&1) and not b>>p[i+1]&1 for i in range(n-1)) for p in paths)
    orphans=0
    for p in itertools.permutations(range(n+1)):
     if not all(G[p[i]]>>p[i+1]&1 for i in range(n)):continue
     old=tuple(x for x in p if x!=n)
     if not all(A[old[i]]>>old[i+1]&1 for i in range(n-1)):orphans+=1
    formula=0
    for mask in range(1,(1<<n)-1):
     comp=((1<<n)-1)^mask
     for x in range(n):
      for y in range(n):
       if mask>>x&1 and comp>>y&1 and not b>>x&1 and b>>y&1 and A[y]>>x&1:
        formula+=endpoint[mask][1][x]*endpoint[comp][0][y]
    require(orphans==formula,'orphan subset formula')
    require(counts(G)[-1]==len(paths)+d+orphans,'insertion partition');insertion_checks+=1
 for r in range(1,8):
  A=labelled(r,(1<<(r*(r-1)//2))-1)
  require(counts(wrapper(A,[8,1,3,6]))[-1]==27*2**r-17,'positive lane')
 for k in range(100):
  H=108*8**k-17
  require(H%7==0 and ((H//7)%7==0)==(k%7==5),'7 valuation')
 result={'status':'PASS','insertion_orphan_checks':insertion_checks,'literal_subset_dp_binary_checks':checks,'complete_class_bound':10000,'max_core_length':maxn,'attained_counts_through_10000':len(bank),'missing_through_199':[x for x in range(1,200,2) if x not in bank],'five_vertex_H13_outmasks':witness13,'distinct_small_profiles':len(profiles),'asymmetric_wrapper_checks':wrappers,'odd_7_valuation_asymmetric_example':seven,'bank_sha256':hashlib.sha256(json.dumps(bank,sort_keys=True).encode()).hexdigest()}
 print(json.dumps(result,indent=2))
if __name__=='__main__':main()

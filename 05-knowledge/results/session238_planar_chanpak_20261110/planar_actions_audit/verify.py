import sys,json,itertools
from fractions import Fraction as F
from functools import lru_cache
checks=0
def check(ok,msg):
 global checks
 checks+=1
 if not ok: raise RuntimeError(msg)
def maps(k,d):
 return {'A':lambda h,r:(h,2*r),'a':lambda h,r:(h,r/2),'C':lambda h,r:(h+k*r,4*r-2*d),'c':lambda h,r:(h-k*(r+2*d)/4,(r+2*d)/4)}
def apply(word,p,ops):
 for w in word:p=ops[w](*p)
 return p
def inverse(w):return ''.join(c.swapcase() for c in w[::-1])
# Chronological words, independently expanded from named compositions.
D='aaC'; K=inverse(D)+'a'+D+'A'; J=inverse(K)+inverse(D)+K+D
for x,d,h,r in itertools.product([F(1,3),F(2),F(3),F(7,2)],[F(1),F(3,2)],[F(-2),F(1)],[F(0),F(2),F(6)]):
 k=x/2;ops=maps(k,d)
 check(apply(D,(h,r),ops)==(h+k*r/4,r-2*d),'D')
 check(apply(K,(h,r),ops)==(h-k*r/8-k*d/4,r-2*d),'K')
 check(apply(J,(h,r),ops)==(h-3*k*d/4,r),'commutator')
 check(apply(inverse(J),(h,r),ops)==(h+3*k*d/4,r),'inverse sign')
 # Integrate proposed Hamiltonian vector fields with R(s)=r-2ds.
 check(h+k*r/4-k*d/4+k*d/4==h+k*r/4,'Hamilton D')
 check(h-k*r/8+k*d/8-3*k*d/8==h-k*r/8-k*d/4,'Hamilton K')
for x in [F(2),F(3)]:
 ops=maps(x/2,F(1))
 check(apply(inverse(J)*8,(F(1),F(2)),ops)==(1+3*x,2),'gap one')
 check(apply(inverse(J)*24,(1+x,F(6)),ops)==(1+10*x,6),'gap two')
 for seed,word,power,target in [((F(1),F(2)),'AAA',8,(1+3*x,16)),((1+x,F(6)),'AA',24,(1+10*x,24))]:
  padded=apply(word,seed,ops)
  check(apply(inverse(J)*power,padded,ops)==target,'padded gap')
  check(target[0]<=x*target[1]/2,'padded physical cone')
# Construct all positive words of total weight <=14 independently by recursion.
levels={0:[('',0)]}; seen={}
for n in range(15):
 if n:
  levels[n]=[(w+'A',2*m) for w,m in levels[n-1]]
  if n>=2:levels[n]+=[(w+'C',4*m+1) for w,m in levels[n-2]]
 for w,m in levels[n]:
  check(m & (m<<1)==0,'separated bits')
  check(m==0 or m.bit_length()<=n-1,'top bit')
  r=apply(w,(F(1),F(2)),maps(F(1),F(1)))[1]
  check(r==2**(n+1)-2*m,'normal form')
  check(r not in seen,'seed injectivity');seen[r]=w
  check(F(4,3)*2**n+F(2,3)<=r<=2**(n+1),'interval')
 admissible=[m for m in range(2**max(0,n-1)) if not(m&(m<<1))]
 check(sorted(m for w,m in levels[n])==admissible,'surjective word code')
# Actual tournaments, odd-cycle packing counts and rooted path complements.
def add(g,bridge):
 n=len(g); out=[v|(0 if bridge and i==0 else 1<<n) for i,v in enumerate(g)]
 return tuple(out+[1 if bridge else 0])
def stats(g,x):
 n=len(g); full=(1<<n)-1; cycles=[]
 for root in range(n):
  def visit(path,used):
   v=path[-1]
   if len(path)>=3 and len(path)%2 and g[v]>>root&1:cycles.append(used)
   for nxt in range(root+1,n):
    if not used>>nxt&1 and g[v]>>nxt&1:visit(path+[nxt],used|1<<nxt)
  visit([root],1<<root)
 @lru_cache(None)
 def pack(mask):
  if not mask:return 1
  b=mask&-mask
  return pack(mask^b)+sum(x*pack(mask^c) for c in cycles if c&b and c&mask==c)
 R=0;delta=0
 def paths(v,mask,length):
  nonlocal R,delta
  z=pack(full^mask);R+=z;delta+=z*(-1 if length%2 else 1)
  for u in range(n):
   if not mask>>u&1 and g[v]>>u&1:paths(u,mask|1<<u,length+1)
 paths(0,1,1)
 return pack(full),R,delta,pack(full^1)
physical=0
for seed in [(0,),(2,4,1)]:
 start=add(seed,False)
 for x in [2,3]:
  initial=stats(start,x)
  for n in range(5):
   for word,m in levels[n]:
    g=start
    for c in word:
     if c=='C':g=add(g,True)
     g=add(g,False)
    actual=stats(g,x); pred=apply(word,tuple(map(F,initial[:2])),maps(F(x,2),F(initial[3])))
    check(actual==(pred[0],pred[1],0,initial[3]),'physical reconstruction')
    physical+=1
if '--inject-failure' in sys.argv:check(False,'intentional audit failure')
print(json.dumps({'status':'PASS','checks':checks,'physical_instances':physical,'positive_words':len(seen),'chronological_D':D,'chronological_K':K,'chronological_J':J},sort_keys=True))

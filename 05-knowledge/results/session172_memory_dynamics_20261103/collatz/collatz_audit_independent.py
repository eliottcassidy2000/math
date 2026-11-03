import itertools,json

def req(b,s):
 if not b:raise RuntimeError(s)
def step(x):return (3*x+1)//2 if x%2 else x//2
stats=[];tested=0;exception_tests=0
for k in range(1,13):
 survivors=0;maps=set();byweight={}
 for w in itertools.product((0,1),repeat=k):
  prefix=[]
  for j in range(1,k+1):
   bits=w[:j];s=sum(bits)
   a=sum(2**i*3**sum(bits[i+1:]) for i,bit in enumerate(bits) if bit)
   prefix.append((s,a))
  s,a=prefix[-1];r=(-a*pow(3**s,-1,2**k))%(2**k);rp=r or 2**k
  maps.add((s,a));byweight.setdefault(s,set()).add(a)
  hs=[a//(2**j-3**s) for j,(s,a) in enumerate(prefix,1) if 2**j>3**s]
  if not hs:survivors+=1;continue
  h=min(hs)
  # Inspect representatives immediately around the predicted threshold, including exceptions far beyond the first six.
  boundary=(h-rp)//(2**k)
  for q in sorted(set([0,1,max(0,boundary-1),max(0,boundary),max(0,boundary+1),max(0,boundary+2)])):
   n=rp+q*2**k;x=n;bits=[];desc=False
   for j in range(k):bits.append(x%2);x=step(x);desc|=x<n
   req(tuple(bits)==w,'cylinder inverse');req(desc==(n>h),'boundary threshold')
   tested+=1;exception_tests+=not desc
 req(len(maps)==2**k,'affine injectivity')
 stats.append([k,survivors])
# Independent survivor recursion propagates actual accepted words rather than weight buckets.
words=[()]
for k in range(1,19):
 words=[w+(b,) for w in words for b in (0,1) if 3**(sum(w)+b)>2**k]
 if k<=12:req(len(words)==stats[k-1][1],'survival enumeration')
# Unbounded-horizon witnesses much deeper than main check.
for k in [65,128,257,1024]:
 n=3*2**k-1;x=n
 for j in range(1,k+1):
  req(x%2==1,'all-one branch');x=step(x);req(x==3**j*3*2**(k-j)-1 and x>n,'long ascent')
print(json.dumps({'status':'PASS','threshold_boundary_representatives':tested,'non_descent_boundary_cases':exception_tests,'survivor_counts':stats,'depth18_survivors':len(words),'long_ascent_depths':[65,128,257,1024]},indent=2))

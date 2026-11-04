import itertools,json

def tournament(n,bits):
    a=[0]*n;k=0
    for i in range(n):
      for j in range(i+1,n):
        if bits>>k&1:a[i]|=1<<j
        else:a[j]|=1<<i
        k+=1
    return a

def counts(a):
    n=len(a);d=[[0]*n for _ in range(1<<n)];h=[0]*(1<<n);h[0]=1
    for v in range(n):d[1<<v][v]=1
    for s in range(1,1<<n):
      h[s]=sum(d[s])
      for v,x in enumerate(d[s]):
        if not x:continue
        r=a[v]&~s
        while r:
          b=r&-r; r-=b; w=b.bit_length()-1;d[s|b][w]+=x
    full=(1<<n)-1
    return h[full],sum(h[s]*h[full^s] for s in range(1,full))

def cyc(parts):
    sizes=list(map(len,parts));offset=[0,sizes[0],sizes[0]+sizes[1]];out=[]
    for i,p in enumerate(parts):
      j=(i+1)%3;cross=((1<<sizes[j])-1)<<offset[j]
      out.extend((x<<offset[i])|cross for x in p)
    return out

def require(ok,msg):
    if not ok:raise RuntimeError(msg)

Ts={n:[tournament(n,b) for b in range(1<<(n*(n-1)//2))] for n in range(1,5)}
profiles={tuple(a):counts(a) for arr in Ts.values() for a in arr}
checks=0; values=set();cases=[]
for na in range(1,5):
 for nb in range(1,5):
  for a in Ts[na]:
   for b in Ts[nb]:
    x,u=profiles[tuple(a)];y,v=profiles[tuple(b)]
    got=counts(cyc([a,b,[0]]))[0];pred=3*x*y+x*v+y*u+u*v
    require(got==pred,('formula',a,b,got,pred));checks+=1
    if got<=45:values.add(got)
for sizes in [(2,2,2),(3,2,2),(3,3,2)]:
 vals=set()
 for parts in itertools.product(*(Ts[n] for n in sizes)):
  vals.add(counts(cyc(parts))[0])
 cases.append({'sizes':sizes,'spectrum':sorted(vals)})
 if sizes==(2,2,2):require(vals=={45},'base45')
 else:require(min(vals)>45,'strict45')
 values.update(v for v in vals if v<=45)
# Add transitive size-five one-side example, the only remaining low-count size.
t5=tournament(5,(1<<10)-1)
require(counts(cyc([t5,[0],[0]]))[0]==33,'five transitive')
values.add(33)
expected={3,5,9,11,15,17,23,27,33,37,45}
require(values==expected,('low spectrum',values))
# Independent ordered-cover enumeration and full kernel controls.
def full_profile(a):
 n=len(a);h=[0]*(1<<n);h[0]=1
 for mask in range(1,1<<n):
  vertices=[v for v in range(n) if mask>>v&1]
  h[mask]=sum(all(a[p[i]]>>p[i+1]&1 for i in range(len(p)-1)) for p in itertools.permutations(vertices))
 q=[[0]*(n+1) for _ in range(1<<n)];q[0][0]=1
 for mask in range(1,1<<n):
  sub=mask
  while sub:
   for k in range(1,n+1):q[mask][k]+=h[sub]*q[mask^sub][k-1]
   sub=(sub-1)&mask
 return q[-1]+[0,0]
ps={tuple(a):full_profile(a) for arr in Ts.values() for a in arr}
full_checks=0
small=[a for n in range(1,4) for a in Ts[n]]
for parts in itertools.product(small,repeat=3):
 profiles3=[ps[tuple(a)] for a in parts];pred=0
 def get(i,k):return profiles3[i][k] if k<len(profiles3[i]) else 0
 for r in range(1,max(map(len,parts))+1):
  pred+=3*get(0,r)*get(1,r)*get(2,r)
  for high in [(0,),(1,),(2,),(0,1),(0,2),(1,2)]:
   z=1
   for i in range(3):z*=get(i,r+(i in high))
   pred+=z
 require(pred==counts(cyc(parts))[0],('full kernel',parts,pred));full_checks+=1
# Direct permutation implementation checks selected literal wrappers independently of DP.
permutation_checks=[]
for sizes in [(1,1,1),(2,2,1),(3,2,1),(2,2,2),(3,2,2)]:
 parts=[tournament(n,(1<<(n*(n-1)//2))-1) for n in sizes];a=cyc(parts)
 literal=sum(all(a[p[i]]>>p[i+1]&1 for i in range(len(p)-1)) for p in itertools.permutations(range(len(a))))
 require(literal==counts(a)[0],('permutation',sizes));permutation_checks.append([sizes,literal])
print(json.dumps({'formula_cases':checks,'full_kernel_cases':full_checks,'small_all_n_cyclic_module_spectrum_le45':sorted(values),'multi_module_tests':cases,'permutation_controls':permutation_checks,'status':'PASS'},indent=2))

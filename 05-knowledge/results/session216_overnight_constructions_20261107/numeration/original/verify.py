from pathlib import Path
import itertools,json,random
D=Path(__file__).parent
C=[tuple(map(int,l.split())) for l in (D/'CORE_STATES.txt').read_text().splitlines()]
CORES=[(8,1,3,6),(2,4,9,19,7),(2,8,3,21,7)]
def need(x,m):
 if not x:raise RuntimeError(m)
def pair_f(n):
 p=[i for i in range(n.bit_length()) if n>>i&1]
 return sum((1<<(j-i-1))+(j==i+1) for i,j in itertools.combinations(p,2))
def formula(c,b,p,r):
 h,q,A,G,T,S,E=C[c];rev=sum(((r>>i)&1)<<(b-i-1) for i in range(b));L=(1<<(b+1))-rev
 return (q*p+A)*L+h*(2*r+1-pair_f(r))-p*G-h*pair_f(p)-T-(p&1)*S-E*((r>>(b-1))&1 if b else 0)
def graph(c,b,p,r):
 ci=0 if c<16 else 1+(c-16)//32;cut=c if c<16 else (c-16)%32;core=CORES[ci];a=p.bit_length();n=len(core);l=a+n+b;v=l;x=l+1;out=[0]*(l+2)
 for i in range(l):
  for j in range(i+1,l):out[i]|=1<<j
 for i in range(n):out[a+i]=(out[a+i]&~(((1<<n)-1)<<a))|(core[i]<<a)
 selected={i for i in range(a) if p>>(a-i-1)&1}|{a+i for i in range(n) if cut>>i&1}|{a+n+i for i in range(b) if r>>(b-i-1)&1}
 for i in range(l):
  out[i]|=1<<v
  if i in selected:out[x]|=1<<i
  else:out[i]|=1<<x
 out[v]|=1<<x
 need(a>0 or cut>0,'strong cut')
 # All retained large examples have a positive prefix.
 hp=next(w for w in itertools.permutations(range(n)) if all(core[w[i]]>>w[i+1]&1 for i in range(n-1)))
 cycle=list(range(a))+[a+i for i in hp]+list(range(a+n,l))+[v,x]
 if a:need(all(out[cycle[i]]>>cycle[(i+1)%len(cycle)]&1 for i in range(len(cycle))),'Hamilton cycle')
 for i in range(l+2):
  for j in range(i+1,l+2):need(bool(out[i]>>j&1)!=bool(out[j]>>i&1),'tournament')
 return out,cycle

def hp(out):
 n=len(out);dp=[{} for _ in range(1<<n)]
 for i in range(n):dp[1<<i][i]=1
 for mask in range(1,1<<n):
  for end,c in dp[mask].items():
   rest=out[end]&~mask
   while rest:
    bit=rest&-rest;j=bit.bit_length()-1;dest=dp[mask|bit];dest[j]=dest.get(j,0)+c;rest-=bit
 return sum(dp[-1].values())
records=[]
for file in ('large.out','large_repair.out','short_side_repairs.out'):
 for line in (D/file).read_text().splitlines():
  z=json.loads(line)
  if z['core']<0:continue
  need(formula(z['core'],z['b'],z['P'],z['R'])==z['N'],'large pair-sum formula')
  out,cycle=graph(z['core'],z['b'],z['P'],z['R']);z.update(outgoing_masks=out,Hamilton_cycle=cycle,order=len(out));records.append(z)
(D/'LARGE_GRAPHS.json').write_text(json.dumps(records,indent=2)+'\n')
physical=0
for c in range(80):
 vals=[]
 for u,v in itertools.product(range(2),repeat=2):
  p=4+u+2*v;out,_=graph(c,0,p,0);H=hp(out);need(H==formula(c,0,p,0),'small physical formula');vals.append(H);physical+=1
 need(vals[3]-vals[2]-vals[1]+vals[0]==-2*C[c][0],'mixed curvature')
rng=random.Random(216);unique=0
for _ in range(20000):
 c=rng.randrange(80);b=rng.randrange(10);r=rng.randrange(1<<b);p=rng.randrange(1,2000);h,q,A,G,T,S,E=C[c];N=formula(c,b,p,r);rev=sum(((r>>i)&1)<<(b-i-1) for i in range(b));L=(1<<(b+1))-rev;K=q*L-G;const=formula(c,b,0,r);x=N-const
 lo=max(1,(x+K-1)//K);hi=(x-h+S)//(K-h);need(lo<=p<=hi,'candidate bound')
 if q*q*(1<<(2*b))>h*N:need(lo==hi==p,'unique balanced');unique+=1
print(json.dumps({'large_literal_graphs':len(records),'physical_subset_counts':physical,'interval_controls':20000,'unique_interval_controls':unique,'all_passed':True},indent=2))

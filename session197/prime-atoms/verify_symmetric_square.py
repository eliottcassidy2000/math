#!/usr/bin/env python3
"""Standalone support census, literal hostile control, and compiled exact source scan."""
import hashlib,itertools,json,math,pathlib,subprocess,tempfile
ROOT=pathlib.Path(__file__).resolve().parent

def require(p,why):
 if not p:raise RuntimeError(why)

def matching(M):
 t=len(M);mate=[-1]*t
 def grow(i,seen):
  for j in range(t):
   if M[i]>>j&1 and j not in seen:
    seen.add(j)
    if mate[j]<0 or grow(mate[j],seen):mate[j]=i;return True
  return False
 return sum(grow(i,set()) for i in range(t))

def classify(M):
 t=len(M);loops={i for i in range(t) if M[i]>>i&1};edges=[{i,j} for i in range(t) for j in range(i+1,t) if M[i]>>j&1]
 if not edges:
  require(len(loops)==2,'diagonal pair');return 'two_loops'
 if len(edges)==1 and loops==edges[0]:return 'full_block'
 centers=set.intersection(*edges)
 require(any(loops<={c} for c in centers),'star with at most its center loop')
 return 'star'

def support_census():
 out=[]
 for t in range(1,6):
  positions=list(itertools.combinations_with_replacement(range(t),2));hist={'star':0,'two_loops':0,'full_block':0}
  for bits in range(1<<len(positions)):
   M=[0]*t
   for k,(i,j) in enumerate(positions):
    if bits>>k&1:M[i]|=1<<j;M[j]|=1<<i
   if matching(M)==2:hist[classify(M)]+=1
  expected=5*math.comb(t,2)+2*t*(2**(t-1)-t)
  require(sum(hist.values())==expected,'complete symmetric rank-two matrix count')
  out.append({'half_order':t,'all_symmetric_matrices':1<<len(positions),'rank2_types':hist})
 return out

def hp(rows):
 n=len(rows);dp=[[0]*n for _ in range(1<<n)]
 for v in range(n):dp[1<<v][v]=1
 for mask in range(1,1<<n):
  for v in range(n):
   if not(mask>>v&1) or mask==(1<<v):continue
   rest=mask^(1<<v)
   dp[mask][v]=sum(dp[rest][u] for u in range(n) if rest>>u&1 and rows[u]>>v&1)
 return sum(dp[-1])

def clone(T,r,m):
 old=[v for v in range(len(T)) if v!=r];labels=old+[r]*m;n=len(labels);out=[0]*n
 for i in range(n):
  for j in range(i+1,n):
   forward=(i<j) if labels[i]==labels[j] else bool(T[labels[i]]>>labels[j]&1)
   if forward:out[i]|=1<<j
   else:out[j]|=1<<i
 return out

def literal_runs(T,r):
 n=len(T);outside=[v for v in range(n) if v!=r];ans=[0,0,0]
 for word in itertools.permutations(outside):
  for k in range(1,4):
   for gaps in itertools.combinations(range(n),k):
    full=[]
    for j in range(n):
     if j in gaps:full.append(r)
     if j<n-1:full.append(word[j])
    if all(T[u]>>v&1 for u,v in zip(full,full[1:])):ans[k-1]+=1
 return ans

def anti(T,p):
 return all(p[p[i]]==i for i in range(len(T))) and all(bool(T[i]>>j&1)==bool(T[p[j]]>>p[i]&1) for i,j in itertools.combinations(range(len(T)),2))

def controls():
 T=[2,4,9,19,39,79,31];r=2;n=len(T);aa=[list(p) for p in itertools.permutations(range(n)) if anti(T,p)];a=literal_runs(T,r)
 require(aa==[list(reversed(range(n)))],'hostile unique anti-involution')
 require(a==[31,30,6],'literal hostile run coefficients')
 require(aa[0][3]==3 and aa[0][r]!=r,'exterior fixed point, root not fixed')
 require(not any(a[2]%c==0 and math.isqrt(a[2]//c)**2==a[2]//c for c in (2,4)),'hostile top coefficient is not twice/four times a square')
 counts=[hp(clone(T,r,m)) for m in (1,2,3)]
 require(counts==[31,91,247],'independent hostile clone counts')
 P=[414,260,8,50,294,327,415,318,12];p=4
 require(anti(P,list(reversed(range(9)))),'marked positive-offset anti-involution')
 positive=[hp(clone(P,p,m)) for m in (1,2,3,4)]
 require(positive==[217,649,1945,5833],'positive-offset clock control')
 require((positive[1]-positive[0])//2==216 and (positive[2]-positive[0]-6*216)//6==72,'positive-offset runs')
 return {'unmarked_hostile':{'rows':T,'root':r,'anti_involutions':aa,'run_coefficients':a,'clone_counts':counts},'mirror_star_scope_control':{'rows':P,'root':p,'run_coefficients':[217,216,72],'clone_counts':positive,'claim':'This marked self-converse +1 clock is not in the ordered mirror-star subclass.'}}

def main():
 support=support_census();hostiles=controls()
 with tempfile.TemporaryDirectory(prefix='session197-square-') as folder:
  exe=pathlib.Path(folder)/'check'
  subprocess.run(['g++','-O3','-std=c++17',str(ROOT/'verify_symmetric_square.cpp'),'-o',str(exe)],check=True)
  census=json.loads(subprocess.check_output([str(exe)],text=True))
 require(census['status']=='PASS' and census['first_equation_survivors']==0,'compiled census')
 require(5**12>485**3,'Busch bound excludes strong order13 or greater for H<=485')
 allowed=[2*q*q for q in range(1,12,2)]
 require([A for A in allowed if 78<=A<242]==[98,162],'square amplitude gap')
 print(json.dumps({'status':'PASS','support_census':support,'exact_census':census,'controls':hostiles,'all_order_amplitude_lower_bound':242,'attainment_claimed':False,'dependencies':['SYMMETRIC_TERNARY.md: marked symmetric amplitude>=78 and complete exclusion through order9','MINIMUM_TERNARY.md: root-run transfer and root-SCC reduction','EVEN_CONVERSE.md: even self-converse H=1 modulo4','Busch Theorem1: 5^(n-1)<=H^3 for strong tournaments'],'source_sha256':hashlib.sha256((ROOT/'verify_symmetric_square.cpp').read_bytes()).hexdigest()},sort_keys=True,indent=2))
if __name__=='__main__':main()

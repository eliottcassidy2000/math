# Full ternary recursive compiler, inherited Newton transform THM-3134.
import importlib.util
from pathlib import Path
_spec=importlib.util.spec_from_file_location("ternary_verified_core",Path(__file__).with_name("verify.py"))
v=importlib.util.module_from_spec(_spec);_spec.loader.exec_module(v)
from math import comb,factorial
from functools import lru_cache
import json
from hashlib import sha256
@lru_cache(None)
def F(p):return tuple(v.cover(p))
@lru_cache(None)
def Bs(p,s,M):
 f=F(p)
 return tuple(sum(comb(k,r)*f[s+r] for r in range(min(k,len(f)-s-1)+1)) for k in range(M+1))
@lru_cache(None)
def transform(a,b,c):
 if len(a)==2:return v.cyc(b,c)
 if len(b)==2:return v.cyc(a,c)
 if len(c)==2:return v.cyc(a,b)
 ps=(a,b,c);M=sum(len(p)-1 for p in ps);q=min(len(p)-1 for p in ps)
 D=[]
 for s in range(q+1):
  arr=[1]*(M+1)
  for p in ps:arr=[x*y for x,y in zip(arr,Bs(p,s,M))]
  ds=[]
  while arr:
   ds.append(arr[0]);arr=[y-x for x,y in zip(arr,arr[1:])]
  D.append(ds)
 ff=[0]+[sum(2**j*comb(d,j)*sum(comb(d+t-1,t)*D[t+j][d-j] for t in range(q-j+1)) for j in range(min(d,q)+1)) for d in range(1,M+1)]
 out=[0]*(M+1);fall=(1,)
 for d in range(1,M+1):
  v.need(ff[d]%factorial(d)==0,'Newton cover divisibility');fall=v.mul(fall,(-(d-1),1))
  for i,x in enumerate(fall):out[i]+=x*(ff[d]//factorial(d))
 result=tuple(out);v.need(v.val(result,1)==v.cscalar(*(F(p) for p in ps)),'Newton scalar gate')
 return result

def compile(seeds,X=4095):
 N=2+(X-1).bit_length()-1;bank=[{} for _ in range(N+1)]
 def put(n,p,d,e):
  if v.val(p,1)<=X and (p not in bank[n] or (d,e)<bank[n][p]):bank[n][p]=(d,e)
 for name,rows,p in seeds:put(len(rows),p,0,name)
 for n in range(2,N+1):
  for i in range(1,n//2+1):
   for a,(da,ea) in sorted(bank[i].items()):
    for b,(db,eb) in sorted(bank[n-i].items()):put(n,v.mul(a,b),max(da,db),'J('+ea+','+eb+')')
  for i in range(1,n//3+1):
   for j in range(i,(n-i)//2+1):
    k=n-i-j
    for a,(da,ea) in sorted(bank[i].items()):
     for b,(db,eb) in sorted(bank[j].items()):
      if i==j and b<a:continue
      for c,(dc,ec) in sorted(bank[k].items()):
       if j==k and c<b:continue
       h=v.cscalar(F(a),F(b),F(c))
       if h>X:continue
       p=transform(a,b,c);put(n,p,1+max(da,db,dc),'T('+ea+','+eb+','+ec+')')
  print('order',n,'profiles',len(bank[n]),flush=True,file=__import__('sys').stderr)
 generators={}
 for row in bank:
  for p,(d,e) in row.items():
   h=v.val(p,1)
   if h not in generators or (d,e)<generators[h]:generators[h]=(d,e)
 reachable={1:(0,'s0')}
 for h in range(3,X+1,2):
  choices=[generators[h]] if h in generators else []
  choices += [(max(d,reachable[h//g][0]),'J('+e+','+reachable[h//g][1]+')') for g,(d,e) in generators.items() if 1<g<h and h%g==0 and h//g in reachable]
  if choices:reachable[h]=min(choices)
 return bank,reachable

def graph(expr,seedrows):
 if expr in seedrows:return seedrows[expr]
 body=expr[2:-1];args=[];level=0;start=0
 for i,c in enumerate(body):
  if c=='(':level+=1
  elif c==')':level-=1
  elif c==',' and level==0:args.append(body[start:i]);start=i+1
 args.append(body[start:]);blocks=[graph(a,seedrows) for a in args]
 return v.join(*blocks) if expr[0]=='J' else v.subst((2,4,1),blocks)

def main():
 distinct={}
 for n in range(1,6):
  for rows in v.tournaments(n):distinct.setdefault((n,v.permutation_poly(rows)),rows)
 seeds=[('s'+str(i),r,p) for i,((n,p),r) in enumerate(sorted(distinct.items()))]
 # independent whole-Q literal controls for the inherited Newton transform
 checks=0
 for ns in __import__('itertools').product(range(1,4),repeat=3):
  ps=[(0,)*n+(1,) for n in ns];p=transform(*ps);rows=v.subst((2,4,1),[v.tt(n) for n in ns])
  v.need(v.qvalues(rows,len(rows))==[v.val(p,k) for k in range(len(p))],'literal Newton full profile');checks+=1
 old=v.compile_profiles(seeds,13);oldgen={v.val(p,1) for row in old for p in row};oldvals={1}
 for h in range(3,4096,2):
  if h in oldgen or any(h%g==0 and h//g in oldvals for g in oldgen if 1<g<h):oldvals.add(h)
 output={'status':'PASS','full_profile_literal_controls':checks,'cutoff':4095,'old_scalar_count':len(oldvals)}
 for name,ss in [('seed5',seeds),('seed5_plus_H19',seeds+[('H19',(40,1,3,6,15,30),v.permutation_poly((40,1,3,6,15,30)))])]:
  bank,reach=compile(ss);added=sorted(set(reach)-oldvals)
  seedrows={tag:rows for tag,rows,p in ss};literal=[]
  for h in added:
   rows=graph(reach[h][1],seedrows)
   if len(rows)<=14:
    v.need(v.hsub(rows)[-1]==h,'literal full-grammar new recipe');literal.append(h)
  output[name]={'literal_new_recipe_checks':len(literal),'scalar_count' :len(reach),'depth_counts':{str(d):sum(x[0]==d for x in reach.values()) for d in sorted(set(x[0] for x in reach.values()))},'profile_counts':[len(row) for row in bank[1:]],'added_over_old':added,'first_missing':sorted(set(range(1,4096,2))-reach.keys())[:40],'scalar_depth_sha256':sha256(json.dumps(sorted((h,d) for h,(d,e) in reach.items()),separators=(',',':')).encode()).hexdigest(),'recipes':{str(h):e for h,(d,e) in sorted(reach.items()) if h in added[:20]}}
 print(json.dumps(output,indent=2,sort_keys=True))
if __name__=='__main__':main()

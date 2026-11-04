"""Independent unordered-cover slot assignment and predecessor run-count audit."""
from pathlib import Path
from itertools import product,permutations,combinations
from collections import Counter
from functools import lru_cache
import hashlib,json,math,subprocess,tempfile
HERE=Path(__file__).resolve().parent;BASE=HERE.parent;SRC=BASE/'named_marked_reduction'
def need(x,msg):
 if not x:raise RuntimeError(msg)
def tt(n):return [((1<<n)-1)^((1<<(u+1))-1) for u in range(n)]
def blow(sizes,rev):
 offs=[sum(sizes[:i]) for i in range(3)];T=[0]*sum(sizes)
 for i,n in enumerate(sizes):
  for u,row in enumerate(tt(n)):T[offs[i]+u]=(row<<offs[i])|sum(((1<<sizes[j])-1)<<offs[j] for j in range(3) if i!=j and((i-j if rev else j-i)%3)==1)
 return T,offs

def marked(T):
 n=len(T);N=1<<n;paths={}
 for M in range(1,N):paths[M]=[p for p in permutations(v for v in range(n) if M>>v&1) if all(T[u]>>v&1 for u,v in zip(p,p[1:]))]
 @lru_cache(None)
 def parts(M):
  if not M:return ((),)
  b=M&-M;R=M^b;S=R;out=[]
  while True:
   P=S|b
   for path in paths[P]:
    for tail in parts(M^P):out.append((path,)+tail)
   if not S:break
   S=(S-1)&R
  return tuple(out)
 patterns={s:[p for p in product((1,2,3),repeat=s) if sum(v.bit_count() for v in p)<=4] for s in range(5)};tables=Counter();covers=0
 for cover in parts(N-1):
  h=len(cover);covers+=math.factorial(h)
  for s in range(min(h,4)+1):
   factor=math.factorial(h-s)
   for chosen in permutations(range(h),s):
    for pattern in patterns[s]:
     marks=tuple(tuple(([cover[j][0]] if role&1 else[])+([cover[j][-1]] if role&2 else[])) for j,role in zip(chosen,pattern));tables[h,s,pattern,marks]+=factor
 need(all(v%math.factorial(h-s)==0 for(h,s,p,m),v in tables.items()),'exact complete unmarked-slot divisibility');high=[v for(h,s,p,m),v in tables.items() if h>=s+4];need(all(v%8==0 for v in high),'every forbidden entry vanishes modulo8')
 zero=single=0
 for u in range(n):
  cardinal=Counter(len(z) for z in parts((N-1)^(1<<u)))
  for h in range(1,n+1):
   actual=tables[h,1,(3,),((u,u),)];expected=cardinal[h-1]*math.factorial(h-1);need(actual==expected,'same-name double endpoint forces exactly a singleton path');single+=1
   if h>=2:need(tables[h,2,(1,1),((u,),(u,))]==0,'repeated name in distinct slots is impossible');zero+=1
 # Counter lookups above return zero without adding keys.
 return {'rows':T,'ordered_covers':covers,'nonzero_marked_entries':len(tables),'high_entries':len(high),'largest_high_count':max(high,default=0),'tables_sha256':hashlib.sha256(repr(sorted(tables.items())).encode()).hexdigest()},zero,single

def main():
 paths=[SRC/x for x in('PROOF.md','verify.py','runs.cpp')];pins={str(p.relative_to(BASE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths};data=json.loads((SRC/'normal.json').read_text());need((SRC/'normal.json').read_bytes()==(SRC/'optimized.json').read_bytes(),'source modes agree');records=[];zeros=singletons=0
 for expected in data['marked_profile_controls']:
  got,z,s=marked(expected['rows']);need(got==expected,'every complete marked table hash and statistic matches');records.append(got);zeros+=z;singletons+=s
 need(len(records)==78,'all78 actual child controls');inputs=[];expected=[];faces=[];positive=0
 for r in data['physical_four_variable_faces']:
  sizes=r['child_sizes'];m=sizes[0];s=r['exceptional_names_in_first_child'];T,off=blow(sizes,r['reverse']);need(T==r['old_rows'],'physical cyclic substitution names');roles=r['roles'];base=r['base_masks'];qs=[]
  for mask in range(16):
   a,c=base
   for j,(port,u) in enumerate(roles):
    if mask>>j&1:
     if port:c^=1<<u
     else:a^=1<<u
   qs.append([a,c]);inputs.append((T,m,a,c))
   ordinary=[(a>>u&1,c>>u&1) for u in range(s,m)];need(len(set(ordinary))==1,'joint port pattern constant outside exceptional names')
  need(qs==r['queries'],'complete four-variable actual faces');vv=r['count_by_first_child_runs'];need(len(vv)==16 and all(len(v)==m+1 for v in vv),'all run-number coordinates');expected.extend(vv);need([sum(v) for v in vv]==r['total_values'],'complete literal totals');coeff=[]
  for mask in range(16):coeff.append([sum((-1)**((mask^sub).bit_count())*vv[sub][h] for sub in range(16) if sub&mask==sub) for h in range(m+1)])
  need(coeff==r['complete_face_coefficients_by_runs'],'every face coefficient including shared-name terms');need(all(v[h]%8==0 for v in vv for h in range(s+4,m+1)),'each high-run literal term divisible8');need(all(v[h]%8==0 for v in coeff for h in range(s+4,m+1)),'all high-run face coefficients');positive+=sum(sum(v[s+4:])>0 for v in vv);faces.append({'sizes':sizes,'reverse':r['reverse'],'exceptional_names':s,'positive_high_run_queries':sum(sum(v[s+4:])>0 for v in vv),'complete_run_vectors_sha256':hashlib.sha256(json.dumps(vv,separators=(',',':')).encode()).hexdigest()})
 with tempfile.TemporaryDirectory() as td:
  exe=Path(td)/'runs';subprocess.run(['g++','-std=c++17','-O2','-Wall','-Wextra','-Werror',str(HERE/'runs.cpp'),'-o',str(exe)],check=True,capture_output=True);wire=str(len(inputs))+'\n'+'\n'.join(f'{len(T)} {m} '+' '.join(map(str,T))+f' {a} {c}' for T,m,a,c in inputs)+'\n';out=subprocess.check_output([str(exe)],input=wire.encode()).decode().splitlines();need([list(map(int,z.split())) for z in out]==expected,'independent predecessor counts by exact first-child run number')
 need(positive==data['positive_high_run_queries'] and positive>0,'nonvacuous modular cancellation');need(pins=={str(p.relative_to(BASE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},'source freeze')
 print(json.dumps({'status':'PASS','complete_marked_table_controls':len(records),'complete_table_records_sha256':hashlib.sha256(json.dumps(records,separators=(',',':')).encode()).hexdigest(),'repeated_name_zero_controls':zeros,'same_name_singleton_controls':singletons,'physical_faces':faces,'complete_actual_run_vectors':len(inputs),'positive_high_run_queries':positive,'source_sha256':pins,'audit_sha256':{x:hashlib.sha256((HERE/x).read_bytes()).hexdigest() for x in('audit.py','runs.cpp')},'scope':'All-order finite marked-profile sufficiency and coefficient/literal cutoffs; no nonannihilation or collision theorem.'},sort_keys=True,indent=2))
if __name__=='__main__':main()

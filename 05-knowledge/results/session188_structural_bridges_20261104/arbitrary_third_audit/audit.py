"""Independent finite-run coefficient and named-family audits, no producer imports."""
from pathlib import Path
from functools import lru_cache
from itertools import product,combinations
import json,hashlib,subprocess,sys
HERE=Path(__file__).resolve().parent
SOURCES=[HERE.parent/'arbitrary_third_block',HERE.parent/'arbitrary_third_root_family']
def need(x,msg):
 if not x:raise RuntimeError(msg)
def extend(R,a,c):
 n=len(R);U=list(R)+[1<<(n+1),0]
 for u in range(n):
  for mask,v in ((a,n),(c,n+1)):
   if mask>>u&1:U[u]|=1<<v
   else:U[v]|=1<<u
 return tuple(U)
def opposite(R):return tuple(((1<<len(R))-1)^r^(1<<i) for i,r in enumerate(R))
def TT(n):return tuple(sum(1<<j for j in range(i+1,n)) for i in range(n))
def blow(blocks,theta=0,bits=0):
 active=[i for i,B in enumerate(blocks) if len(B)>1];offs=[0]
 for B in blocks:offs.append(offs[-1]+len(B))
 R=[0]*offs[-1]
 for i,B in enumerate(blocks):
  if i in active and bits>>active.index(i)&1:B=opposite(B)
  j=(i+(-1 if theta else 1))%3
  for u,row in enumerate(B):R[offs[i]+u]=(row<<offs[i])|(((1<<len(blocks[j]))-1)<<offs[j])
 return tuple(R)
@lru_cache(None)
def H(R):
 n=len(R);dp=[{} for _ in range(1<<n)]
 for u in range(n):dp[1<<u][u]=1
 for used in range(1,1<<n):
  for u,z in dp[used].items():
   rem=R[u]&~used
   while rem:
    bit=rem&-rem;rem-=bit;v=bit.bit_length()-1;dest=dp[used|bit];dest[v]=dest.get(v,0)+z
 return sum(dp[-1].values())
def words(Q,a,c,r):
 R=extend(Q,a,c)
 @lru_cache(None)
 def visit(left,last):
  if not any(left):return 1
  out=0
  for v,k in enumerate(left):
   if k and (last==5 or R[last]>>v&1):
    new=list(left);new[v]-=1;out+=visit(tuple(new),v)
  return out
 return visit(tuple(r)+(1,1),5)
def modules(R):
 return [M for M in range(1,(1<<len(R))-1) if M.bit_count()>1 and all(row&M in(0,M) for v,row in enumerate(R) if not M>>v&1)]
def xor(R,sets):
 U=list(R)
 for u,v in combinations(range(len(R)),2):
  if sum(bool(M>>u&1 and M>>v&1) for M in sets)%2:U[u]^=1<<v;U[v]^=1<<u
 return tuple(U)
def main():
 pins={str(S.name+'/'+p):hashlib.sha256((S/p).read_bytes()).hexdigest() for S in SOURCES for p in ('PROOF.md','verify.py')}
 data=[json.loads((S/'normal.json').read_text()) for S in SOURCES]
 forms={};wordcoeffs=0
 for m,q in [(2,[0,1,2]),(3,[0,1,6,6])]:
  forms[str(m)]=[]
  for a,c in [(1,2),(2,4),(4,1)]:
   vec=[]
   for k in range(1,m+4):
    value=0
    for i,j in product(range(1,m+1),repeat=2):
     value+=q[i]*q[j]*(words((2,4,1),a,c,(i,j,k))-words((4,1,2),a,c,(i,j,k)));wordcoeffs+=1
    vec.append(value)
   forms[str(m)].append(vec)
 need(forms==data[0]['forms'],'independent literal finite-run coefficient forms')
 literals=cubestates=cubequeries=0
 for case in data[0]['literal_cases']:
  blocks=(TT(case['m']),TT(case['m']),tuple(case['C']));T,U=blow(blocks),blow(blocks,1)
  need([H(T),H(U)]==case['scalar'],'literal scalar controls')
  for rec in case['queries']:need(H(extend(T,rec['a'],rec['c']))==rec['T'] and H(extend(U,rec['a'],rec['c']))==rec['U'],'literal arbitrary third queries');literals+=1
 for case in data[0]['cubes']:
  blocks=(TT(case['m']),TT(case['m']),tuple(case['C']));r=sum(len(B)>1 for B in blocks);rootcount=1 if case['m']==2 else 2
  need(case['states']==2**(r+1) and case['adaptive_bound']==rootcount+r and case['fixed_count']==rootcount+2*r,'cube sizes and queries')
  seen=set()
  for rec in case['all_fixed_signatures']:
   R=blow(blocks,rec['theta'],rec['bits']);vals=tuple(H(extend(R,*q)) for q in case['queries']);need(vals==tuple(rec['values']) and vals not in seen,'independent cube signature');seen.add(vals);cubestates+=1;cubequeries+=len(vals)
 families=familyqueries=0;domains=[]
 for case in data[1]['cases']:
  T=tuple(case['T']);need(T==blow((TT(case['m']),TT(case['m']),tuple(case['C']))),'root source physical names');need(modules(T)==case['proper_modules'],'complete proper original module list')
  need(case['families']==1<<len(case['proper_modules']),'complete root family universe')
  base=tuple(H(extend(T,*q)) for q in case['battery']);need(base==tuple(case['baseline_values']),'root baseline')
  equal=strict=0
  for rec in case['records']:
   expected=[(1<<len(T))-1]+[M for j,M in enumerate(case['proper_modules']) if rec['family_bits']>>j&1];need(rec['modules']==expected,'exact original XOR family')
   U=xor(T,expected);vals=tuple(H(extend(U,*q)) for q in case['battery']);need(vals==tuple(rec['values']) and vals!=base and vals[0]>=base[0],'independent root family query battery');equal+=vals[0]==base[0];strict+=vals[0]>base[0];families+=1;familyqueries+=len(vals)
  need(equal==case['scalar_equal'] and strict==case['scalar_strict'],'equality strict census');domains.append({'label':case['label'],'families':equal+strict,'equal':equal,'strict':strict})
 for S,record in zip(SOURCES,data):
  normal=subprocess.check_output([sys.executable,str(S/'verify.py')]);optimized=subprocess.check_output([sys.executable,'-O',str(S/'verify.py')]);need(normal==optimized and json.loads(normal)==record,'unchanged producer reproducible')
  (HERE/(S.name+'.normal.stdout')).write_bytes(normal);(HERE/(S.name+'.optimized.stdout')).write_bytes(optimized)
 need(pins=={str(S.name+'/'+p):hashlib.sha256((S/p).read_bytes()).hexdigest() for S in SOURCES for p in ('PROOF.md','verify.py')},'frozen sources')
 out={'status':'PASS','independent_run_word_terms':wordcoeffs,'six_exact_coefficient_forms':forms,'independent_literal_queries':literals,'independent_cube_states':cubestates,'independent_cube_query_values':cubequeries,'independent_root_families':families,'independent_root_query_values':familyqueries,'root_domains':domains,'source_sha256':pins,'audit_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
 (HERE/'receipt.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,sort_keys=True))
if __name__=='__main__':main()

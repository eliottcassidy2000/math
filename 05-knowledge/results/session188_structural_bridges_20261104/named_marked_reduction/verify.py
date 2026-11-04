"""Exact slot-factorial controls and full four-variable physical run banks."""
from pathlib import Path
from itertools import permutations,product,combinations
from collections import Counter
import hashlib,json,math,subprocess,tempfile
ROOT=Path(__file__).resolve().parent
def need(x,msg):
 if not x:raise RuntimeError(msg)
def tt(n):return [((1<<n)-1)^((1<<(i+1))-1) for i in range(n)]
def rows(n,code):
 T=[0]*n
 for b,(i,j) in enumerate(combinations(range(n),2)):
  u,v=(i,j) if code>>b&1 else(j,i);T[u]|=1<<v
 return T
def marked(T):
 n=len(T);tables=Counter();covers=0
 patterns={s:[p for p in product([1,2,3],repeat=s) if sum(v.bit_count() for v in p)<=4] for s in range(5)}
 for word in permutations(range(n)):
  bad=sum((not(T[u]>>v&1))<<i for i,(u,v) in enumerate(zip(word,word[1:])))
  for cuts in range(1<<(n-1)):
   if cuts&bad!=bad:continue
   paths=[];start=0
   for i in range(n):
    if i==n-1 or cuts>>i&1:paths.append(word[start:i+1]);start=i+1
   h=len(paths);covers+=1
   for s in range(min(4,h)+1):
    for pattern in patterns[s]:
     marks=tuple(tuple(([paths[j][0]] if p&1 else [])+([paths[j][-1]] if p&2 else [])) for j,p in enumerate(pattern));tables[h,s,pattern,marks]+=1
 for (h,s,p,m),v in tables.items():need(v%math.factorial(h-s)==0,'exact free unmarked-slot factorial')
 high=[(k,v) for k,v in tables.items() if k[0]>=k[1]+4];need(all(v%8==0 for k,v in high),'all nonzero forbidden marked tables divisible8')
 return {'rows':T,'ordered_covers':covers,'nonzero_marked_entries':len(tables),'high_entries':len(high),'largest_high_count':max([v for k,v in high],default=0),'tables_sha256':hashlib.sha256(repr(sorted(tables.items())).encode()).hexdigest()}
def blow(sizes,reverse):
 offs=[];T=[]
 for n in sizes:offs.append(len(T));T.extend([0]*n)
 for i,n in enumerate(sizes):
  for u,r in enumerate(tt(n)):
   T[offs[i]+u]|=r<<offs[i]
   for j,m in enumerate(sizes):
    if i!=j and ((i-j if reverse else j-i)%3)==1:T[offs[i]+u]|=((1<<m)-1)<<offs[j]
 return T,offs

def main():
 ts=[rows(n,k) for n in range(1,5) for k in range(1<<(n*(n-1)//2))]+[tt(5),tt(6),[58,60,57,16,32,8]]
 controls=[marked(T) for T in ts];need(any(r['high_entries'] for r in controls),'nonvacuous factorial tests')
 cases=[];meta=[]
 for s,sizes in [(0,[4,2,2]),(1,[5,2,2]),(2,[6,3,3])]:
  m=sizes[0]
  for rev in [False,True]:
   T,offs=blow(sizes,rev)
   if s==0:roles=[(0,offs[1]),(0,offs[2]),(1,offs[1]),(1,offs[2])];base=[(1<<m)-1,0]
   elif s==1:roles=[(0,0),(0,offs[1]),(1,0),(1,offs[2])];base=[(1<<m)-1,0]
   else:roles=[(0,0),(0,1),(1,0),(1,1)];base=[((1<<sizes[1])-1)<<offs[1],((1<<sizes[2])-1)<<offs[2]]
   queries=[]
   for mask in range(16):
    v=base.copy()
    for j,(port,name) in enumerate(roles):
     if mask>>j&1:v[port]^=1<<name
    queries.append(v)
   cases.append((T,m,queries));meta.append({'exceptional_names_in_first_child':s,'child_sizes':sizes,'reverse':rev,'old_rows':T,'roles':roles,'base_masks':base,'queries':queries})
 lines=[]
 for T,m,qq in cases:lines.extend([f'{len(T)} {m} {len(qq)}',' '.join(map(str,T))]+[f'{a} {c}' for a,c in qq])
 with tempfile.TemporaryDirectory() as td:
  exe=Path(td)/'runs';subprocess.run(['g++','-std=c++17','-O2','-Wall','-Wextra','-Werror',str(ROOT/'runs.cpp'),'-o',str(exe)],check=True);output=subprocess.check_output([str(exe)],input=('\n'.join(lines)+'\n').encode()).decode().splitlines()
 at=0;positive=0
 for r in meta:
  m=r['child_sizes'][0];s=r['exceptional_names_in_first_child'];vv=[list(map(int,line.split())) for line in output[at:at+16]];at+=16;need(all(len(v)==m+1 for v in vv),'all run coordinates consumed');r['count_by_first_child_runs']=vv;r['total_values']=[sum(v) for v in vv]
  coeff=[]
  for mask in range(16):
   cf=[sum((-1)**((mask^sub).bit_count())*vv[sub][h] for sub in range(16) if sub&mask==sub) for h in range(m+1)];coeff.append(cf)
   need(all(cf[h]%8==0 for h in range(s+4,m+1)),'complete four-variable coefficient high-run cancellation')
  r['complete_face_coefficients_by_runs']=coeff
  need(all(sum(v[s+4:])%8==0 for v in vv),'literal query truncation modulo8');positive+=sum(sum(v[s+4:])>0 for v in vv)
 need(at==len(output) and positive>0,'positive high-run counts, not vacuous bounds')
 print(json.dumps({'status':'PASS','marked_profile_controls':controls,'physical_four_variable_faces':meta,'physical_counts':len(output),'positive_high_run_queries':positive,'source_pins':{f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in ['PROOF.md','verify.py','runs.cpp']},'scope':'Factorial marked-profile sufficiency for complete named responses modulo8, not reconstruction or a collision.'},sort_keys=True,indent=2))
if __name__=='__main__':main()

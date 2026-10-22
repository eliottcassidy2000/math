import itertools as it, math, json, csv, hashlib, pathlib
from collections import Counter
from fractions import Fraction
import sys, subprocess
O=pathlib.Path(__file__).resolve().parent
R=pathlib.Path(sys.argv[1]).resolve() if len(sys.argv)>1 else pathlib.Path(subprocess.check_output(['git','-C',str(O),'rev-parse','--show-toplevel'],text=True).strip())

def arcs(n,code):
 e=list(it.combinations(range(n),2));a=[[False]*n for _ in range(n)]
 for b,(i,j) in enumerate(e):a[i][j]=bool(code>>b&1);a[j][i]=not a[i][j]
 return a

def hp(a):return sum(all(a[p[i]][p[i+1]] for i in range(len(a)-1)) for p in it.permutations(range(len(a))))
def aut(a):return sum(all(a[i][j]==a[p[i]][p[j]] for i in range(len(a)) for j in range(len(a))) for p in it.permutations(range(len(a))))
def strong(a):
 r=[row[:] for row in a];n=len(a)
 for i in range(n):r[i][i]=True
 for k in range(n):
  for i in range(n):
   for j in range(n):r[i][j]|=r[i][k] and r[k][j]
 return all(map(all,r))
# Independent Hamiltonian-permutation enumeration, no subset DP imported.
counts={};atoms={};weighted=[Counter() for _ in range(6)];weighted[0][1]=Fraction(1)
for n in range(1,6):
 counts[n]=Counter();atoms[n]=Counter()
 for c in range(2**(n*(n-1)//2)):
  a=arcs(n,c);h=hp(a);counts[n][h]+=1
  if strong(a):atoms[n][h]+=1
 for m in range(1,n+1):
  for h,v in atoms[m].items():
   for hh,vv in weighted[n-m].items():weighted[n][h*hh]+=Fraction(v,math.factorial(m))*vv
 assert weighted[n]=={h:Fraction(v,math.factorial(n)) for h,v in counts[n].items()}
rows=list(csv.DictReader((R/'05-knowledge/results/actual_hamilton_multiplicities_20261021/session24-multiplicities/n8_classes.csv').open()))
h63=[r for r in rows if int(r['H'])==63];assert len(rows)==6880 and len(h63)==2
for r in h63:
 code=int(r['code']);a=[[bool(code>>(8*i+j)&1) for j in range(8)] for i in range(8)]
 assert hp(a)==63 and aut(a)==1 and strong(a)
assert 5**8>63**3
# Independently enumerate E8 roots in doubled coordinates.
roots=[]
for x in it.product(range(-2,3),repeat=8):
 if sum(t*t for t in x)!=8:continue
 if (all(t%2==0 for t in x) or all(t%2 for t in x)) and sum(x)%4==0:roots.append(x)
assert len(roots)==240
assert sum(all(t%2 for t in x) and sum(x)==0 for x in roots)==70
# Colored partition recurrence, different from inverse-product implementation.
N=12
p=[1]+[0]*N
for n in range(1,N+1):
 value=24*sum(sum(d for d in range(1,k+1) if k%d==0)*p[n-k] for k in range(1,n+1));assert value%n==0;p[n]=value//n
assert p[:8]==[1,24,324,3200,25650,176256,1073720,5930496]
# Exact strong counts to 40, separate math.comb recurrence rather than mod Pascal.
s=[0]*41;t=[2**(n*(n-1)//2) for n in range(41)]
for n in range(1,41):s[n]=t[n]-sum(math.comb(n,k)*s[k]*t[n-k] for k in range(1,n))
assert s[21]%21==pow(2,10,21)==16 and pow(2,20,21)==4
# Catalogue structural reconciliation plus every present match hash and line preview.
c=json.loads((O.parent/'session27-repunits/catalogue.json').read_text());changed=[];checked=0
for e in c['entries']:
 if not e['sha256']:continue
 path=R/e['path']
 if not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest()!=e['sha256']:changed.append(e['path']);continue
 checked+=1;lines=path.read_text().splitlines()
 assert all(lines[m['line']-1][:700]==m['text'] for m in e['matches'])
assert c['tracked_paths']==c['text_files_scanned']+c['missing_worktree_files']+len(c['large_files_skipped'])+c['nontext_or_read_error_count']
result={'status':'PASS','independent_permutation_tournaments':sum(sum(v.values()) for v in counts.values()),'n_le_5_histograms':counts,'n8_inventory_rows':len(rows),'n8_H63_witnesses_independently_HP_aut_strong_checked':h63,'E8_coordinate_roots':len(roots),'balanced_half_roots':70,'oscillator24':p,'Camion21':s[21]%21,'catalogue_entries':len(c['entries']),'catalogue_lines':sum(len(e['matches']) for e in c['entries']),'catalogue_hashes_and_previews_checked':checked,'catalogue_paths_changed_since_capture':changed}
print(json.dumps(result,indent=2))

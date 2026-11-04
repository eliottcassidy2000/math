"""Independent root-skeleton/inverse-Stirling audit, no source helper imports."""
import itertools,math,json,hashlib,subprocess,sys
from pathlib import Path
SRC=Path(__file__).resolve().parent.parent/'skeleton'
def need(x,msg):
 if not x:raise RuntimeError(msg)
def hp(rows):
 n=len(rows)
 if n==0:return 1
 D=[{} for _ in range(1<<n)]
 for v in range(n):D[1<<v][v]=1
 for mask in range(1,1<<n):
  for last,count in D[mask].items():
   targets=rows[last]&~mask
   while targets:
    bit=targets&-targets;targets-=bit;nxt=bit.bit_length()-1
    E=D[mask|bit];E[nxt]=E.get(nxt,0)+count
 return sum(D[-1].values())
def inflate(rows,root,r):
 ext=[v for v in range(len(rows)) if v!=root];n=len(ext);R=[0]*(n+r)
 for i,u in enumerate(ext):
  for j,v in enumerate(ext):
   if rows[u]>>v&1:R[i]|=1<<j
  for j in range(n,n+r):
   if rows[u]>>root&1:R[i]|=1<<j
   else:R[j]|=1<<i
 for i in range(n,n+r):
  for j in range(i+1,n+r):R[i]|=1<<j
 return R

def skeleton(rows,root):
 ext=[v for v in range(len(rows)) if v!=root];ans=[0]*(len(ext)+2)
 for p in itertools.permutations(ext):
  for mask in range(1<<(len(ext)+1)):
   word=[]
   for k in range(len(ext)+1):
    if mask>>k&1:word.append(root)
    if k<len(ext):word.append(p[k])
   if all(rows[u]>>v&1 for u,v in zip(word,word[1:])):ans[mask.bit_count()]+=1
 while ans[-1]==0:ans.pop()
 return ans

def rows_from_mask(n,m):
 r=[0]*n
 for bit,(u,v) in enumerate(itertools.combinations(range(n),2)):
  if m>>bit&1:r[u]|=1<<v
  else:r[v]|=1<<u
 return r
S=[[0]*8 for _ in range(8)];s=[[0]*8 for _ in range(8)];S[0][0]=s[0][0]=1
for n in range(1,8):
 for k in range(n+1):
  S[n][k]=(S[n-1][k-1] if k else 0)+k*S[n-1][k]
  s[n][k]=(s[n-1][k-1] if k else 0)-(n-1)*s[n-1][k]
records=[]
counts={'tournaments':0,'rooted_sources':0,'literal_inflations':0,'inverse_coefficients':0,'power_expansions':0}
for n in range(1,5):
 for mask in range(1<<(n*(n-1)//2)):
  rows=rows_from_mask(n,mask);counts['tournaments']+=1
  for root in range(n):
   B=skeleton(rows,root);d=len(B)-1;counts['rooted_sources']+=1
   H=[hp(inflate(rows,root,r)) for r in range(7)]
   records.append({'rows':rows,'root':root,'B':B,'H':H})
   for r in range(7):
    need(H[r]==sum(B[k]*math.factorial(k)*S[r][k] for k in range(min(d,r)+1)),'all-order formula sample')
    counts['literal_inflations']+=1
   for k in range(7):
    T=sum(s[k][r]*H[r] for r in range(k+1))
    need(T%math.factorial(k)==0 and T>=0,'inverse integrality/nonnegativity')
    need(T//math.factorial(k)==(B[k] if k<=d else 0),'inverse exact coefficient')
    counts['inverse_coefficients']+=1
   C=[sum(B[k]*(-1)**(k-j)*math.comb(k,j) for k in range(j,d+1)) for j in range(d+1)]
   need(C[d]==B[d] and B[d]>0,'top exponential coefficient')
   for r in range(7):
    need(H[r]==sum(C[j]*j**r for j in range(d+1)),'power expansion including zero mode')
    counts['power_expansions']+=1
   need(d<=min(rows[root].bit_count(),n-1-rows[root].bit_count())+1,'run-degree bound')
def substitute(rows,root,A):
 ext=[v for v in range(len(rows)) if v!=root];n=len(ext);R=inflate(rows,root,len(A))
 for i,row in enumerate(A):
  R[n+i]&=(1<<n)-1
  R[n+i]|=row<<n
 return R

def covers(A):
 n=len(A);N=(1<<n)-1;h=[]
 for mask in range(N+1):
  p=[v for v in range(n) if mask>>v&1]
  h.append(hp([sum(1<<j for j,v in enumerate(p) if A[u]>>v&1) for u in p]))
 q=[[0]*(N+1) for _ in range(n+1)];q[0][0]=1
 for k in range(1,n+1):
  for mask in range(1,N+1):
   sub=mask
   while sub:
    q[k][mask]+=h[sub]*q[k-1][mask^sub];sub=(sub-1)&mask
 return [q[k][N] for k in range(n+1)]
modules=[]
for n in range(4):
 for mask in range(1<<(n*(n-1)//2)):
  A=rows_from_mask(n,mask);modules.append((A,covers(A)))
general_substitutions=0
for X in records:
 for A,q in modules:
  actual=hp(substitute(X['rows'],X['root'],A))
  need(actual==sum(X['B'][k]*q[k] for k in range(min(len(X['B']),len(q)))),'arbitrary-module path-cover substitution')
  general_substitutions+=1
spacing={'pairs':0,'passing_tail_tests':0,'separate_r0_failures':0};r0_hostile=None
for A in records:
 for B in records:
  spacing['pairs']+=1;D=max(len(A['B']),len(B['B']))-1;c=A['H'][1]-B['H'][1]
  finite=all(A['H'][r]-B['H'][r]==c for r in range(1,D+1))
  coeff=all((A['B'][k] if k<len(A['B']) else 0)==(B['B'][k] if k<len(B['B']) else 0) for k in range(2,D+1)) and A['B'][1]-B['B'][1]==c
  need(finite==coeff,'finite spacing equivalent to coefficient criterion')
  if finite:
   spacing['passing_tail_tests']+=1
   need(all(A['H'][r]-B['H'][r]==c for r in range(1,7)),'finite spacing extends on additional controls')
   need((A['H'][0]-B['H'][0]==c)==(A['B'][0]-B['B'][0]==c),'separate r0 condition')
   if A['H'][0]-B['H'][0]!=c:
    spacing['separate_r0_failures']+=1
    if r0_hostile is None:r0_hostile={'A':A,'B':B,'tail_spacing':c}
seeds=((16,5,33,39,46,3),(16,5,49,39,42,3),(18,4,49,39,42,3))
data=[]
for j,rows in enumerate(seeds):
 B=skeleton(rows,5);need(B==[9+2*j,27+2*j,21,4],'seed skeleton')
 H=[hp(inflate(rows,5,r)) for r in range(8)]
 need(H[0]==9+2*j,'deleted-root endpoint')
 for r in range(1,8):need(H[r]==4*3**r+9*2**r-3+2*j,'gate formula')
 C=[sum(B[k]*(-1)**(k-i)*math.comb(k,i) for k in range(i,4)) for i in range(4)]
 need(C==[-1,-3+2*j,9,4],'EGF polynomial')
 data.append({'j':j,'B':B,'power_coefficients':C,'H_r0_to_7':H})
# Execute source exactly as published, without changing its packet.
normal=subprocess.check_output([sys.executable,str(SRC/'verify.py')]);optimized=subprocess.check_output([sys.executable,'-O',str(SRC/'verify.py')])
need(normal==optimized,'source normal/optimized equality')
Path(__file__).with_name('source.normal.stdout').write_bytes(normal);Path(__file__).with_name('source.optimized.stdout').write_bytes(optimized)
result={'status':'PASS','counts':counts,'general_module_substitutions':general_substitutions,'spacing':spacing,'r0_hostile':r0_hostile,'seed_controls':data,'source_normal_optimized_byte_match':True,'source_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (SRC/'PROOF.md',SRC/'verify.py')},'audit_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
Path(__file__).with_name('receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,sort_keys=True))

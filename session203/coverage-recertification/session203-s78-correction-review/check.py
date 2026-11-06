from pathlib import Path
from functools import lru_cache
import hashlib,json,re,math
P=Path(__file__).resolve().parent.parent/'session203-independent-audit'
def req(x,m):
 if not x:raise RuntimeError(m)
checked=0
for line in (P/'COVERAGE_MANIFEST.sha256').read_text().splitlines():
 h,name=line.split(None,1);req(hashlib.sha256((P/name.strip()).read_bytes()).hexdigest()==h,'manifest '+name);checked+=1
@lru_cache(None)
def g(v):
 if v==0:return 1
 u=v//2
 return g(u) if v%2 else g(u)+u+u%2

def direct(v):
 if not v:return 1
 return 1+(v%2==0)+sum(not(v>>(i+1)&1) and bool(v>>i&1) for i in range(v.bit_length()-1))+sum(v>>(i+1) for i in range(v.bit_length()) if not(v>>i&1))

def split(v):
 k,r=divmod(v,1<<16)
 if not k:return direct(r)
 if k>=1<<16:return direct(v)
 Z=sum(1<<(15-i) for i in range(16) if not(r>>i&1))
 Q=(r%2==0)+sum(not(r>>(i+1)&1) and bool(r>>i&1) for i in range(15))+sum(r>>(i+1) for i in range(16) if not(r>>i&1))
 return direct(k)+k*Z+Q-(k%2==0 and not(r>>15&1))
probes=[0,1,2,3,(1<<32)-2,(1<<32)-1,1<<32]
probes += [(k<<16)|r for k in (1,2,3,32767,32768,65535,65536) for r in (0,1,2,3,32767,32768,65534,65535)]
for v in probes:req(split(v)==direct(v)==g(v),'endpoint/split '+str(v))
for k in range(1,33):req(g(1<<k)==(1<<k)+1,'powers of two')
bank=[tuple(map(int,s.split())) for s in (P/'aud_bank7_inherited.txt').read_text().splitlines()]
req(len(bank)==2732 and len(set(bank))==2732,'bank size')
req(all(h>=e>=0 and R>=1 and h>=1 for h,e,R in bank),'physical bank necessities')
X=1<<32;CH=1<<22
f1pairs=sum(X//(u+1)-u for u in range(math.isqrt(X)))
f2eval=sum(X//R for h,e,R in bank)
f1items=sum((X//(u+1)-u+CH-1)//CH for u in range(math.isqrt(X)))
f2items=sum((X//R+CH-1)//CH for h,e,R in bank)
req((f1pairs,f2eval,f1items,f2items)==(47964383242,134010742236,72787,33303),'full domain totals')
req(sum(1<<(n*(n-1)//2) for n in range(1,8))==2131019,'labelled count')
out=(P/'CORRECTED_NOF0_32.out').read_text();old=(P/'S78_FULL32_FROZEN.out').read_text()
pat=r'k=(\d+) (\d+) '
req(re.findall(pat,out)==re.findall(pat,old),'all dyadic block counts')
req(len(re.findall(pat,out))==22,'number of dyadic block counts')
req('F0 SKIPPED (NOF0 build)' in out,'NOF0 evidence')
req('covered by none of F0, F1, F2: 7 21  [total 2]' in out,'full odd scan evidence')
# Independently rerun a small prefix using recursive g and ordinary Python sets.
small=1<<16
vals1=set();vals2=set()
for u in range(math.isqrt(small)):
 for v in range(u,small//(u+1)):
  n=(u+1)*(v+1)+g(u)*g(v)-(u%2==0 and v%2==0)
  if n<=small:vals1.add(n)
for h,e,R in bank:
 for v in range(small//R):
  n=g(v)*h-(e if v%2==0 else 0)+(v+1)*R
  if n<=small:vals2.add(n)
miss=sorted(set(range(1,small+1,2))-vals1-vals2)
req(miss==[7,21],'independent smaller prefix')
# Small physical join controls, explicitly building each graph and counting paths.
def apex(v):
 w=bin(v)[2:] if v else ''
 n=len(w);out=[0]*(n+1)
 for i in range(n):
  for j in range(i+1,n):out[i]|=1<<j
  if w[i]=='1':out[n]|=1<<i
  else:out[i]|=1<<n
 return out,n

def hp(out):
 n=len(out);full=(1<<n)-1;dp=[{} for _ in range(1<<n)]
 for i in range(n):dp[1<<i][i]=1
 for mask in range(1,1<<n):
  for i,val in dp[mask].items():
   av=out[i]&~mask
   while av:
    b=av&-av;av-=b;j=b.bit_length()-1
    dp[mask|b][j]=dp[mask|b].get(j,0)+val
 return sum(dp[full].values())

def join(A,a,B,b):
 na,nb=len(A),len(B);out=list(A)+[0]*nb
 for i in range(nb):
  for j in range(nb):
   if B[j]>>i&1:out[na+i]|=1<<(na+j) # converse of B
 for i in range(na):
  for j in range(nb):
   if i==a and j==b:out[na+j]|=1<<i
   else:out[i]|=1<<(na+j)
 return out
ngraph=0
for u in range(16):
 A,a=apex(u)
 req(hp(A)==g(u),'apex physical')
 for v in range(16):
  B,b=apex(v)
  req(hp(join(A,a,B,b))==(u+1)*(v+1)+g(u)*g(v)-(u%2==0 and v%2==0),'join physical')
  ngraph+=1
report={'status':'PASS','frozen_manifest_files':checked,'endpoint_and_split_controls':len(probes),'F1_pairs_full':f1pairs,'F2_evaluations_full':f2eval,'F1_work_items_full':f1items,'F2_work_items_full':f2items,'bank_triples':len(bank),'bank_max_coordinates':[max(x[i] for x in bank) for i in range(3)],'matching_dyadic_blocks':22,'independent_prefix_cap':small,'independent_prefix_misses':miss,'literal_join_graphs':ngraph,'full_pass_independently_rerun_by_this_reviewer':False,'scope':'review of corrected S78 implementation and supplied full-run receipts, with independent smaller checks; not a third independently authored full-range sieve'}
print(json.dumps(report,indent=2,sort_keys=True))

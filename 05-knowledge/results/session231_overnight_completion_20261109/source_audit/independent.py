import hashlib,json,pathlib,subprocess,random,os
A=pathlib.Path(__file__).resolve().parent
P=A.parent/'source'
if not P.is_dir():P=A.parent/'session231-source-coverage'
if os.environ.get('MATH_REPO'):
 R=pathlib.Path(os.environ['MATH_REPO']).resolve()
else:
 R=next((p for p in A.parents if (p/'.git').exists()), A.parent/'math-wt-session231-overnight')
BASE='9c48b9ac1e2a9c22b0814686c76e3358ad56051d'
def ck(x,msg):
 if not x: raise RuntimeError(msg)
def H(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def g(n):
 # Independent left-to-right state machine, not producer's right-to-left sum.
 q,h=0,1
 for bit in bin(n)[2:]:
  if bit=='0':h+=q+q%2
  q=2*q+int(bit)
 return h
def F(y,v):
 a,b,c=y
 return a*g(v)-b*(1-v%2)+c*(v+1)
def modulus(a,c,N):
 M=2
 while c*c*M*M<4*a*N or c*M<4*a:M*=2
 return M
def solve(y,N):
 a,b,c=y;M=modulus(a,c,N);out=[]
 for r in range(M):
  C=M-1-int(format(r,'0%db'%(M.bit_length()-1))[::-1],2)
  E=int(r<M//2);L=c*M+a*C;Q=a*(g(r)-1)-b*(1-r%2)+c*(r+1);x=N-Q
  ck(Q>=0,'Q positivity')
  if x<a:continue
  lo=max(0,-((2*a-x)//(L+a)));hi=(x-a)//L
  ck(hi-lo<=0,'interval unique')
  if hi>=lo and L*lo+a*g(lo)+a*E*(lo%2)==x:out.append(M*lo+r)
 return set(out)
Y=[tuple(map(int,s.split())) for s in (P/'bank.txt').read_text().splitlines()]
ck(len(Y)==len(set(Y))==2732,'bank cardinality')
ck(all(c>=a>=b>=0 and a>0 and a%2 for a,b,c in Y),'bank cone')
head=subprocess.check_output(['git','rev-parse',BASE],cwd=R,text=True).strip()
bankpath='04-computation/experiments/rooted_triple_bank_order7_20261104.txt'
blob=subprocess.check_output(['git','show',BASE+':'+bankpath],cwd=R)
ck(blob==(P/'bank.txt').read_bytes(),'bank provenance')
for s in (P/'SHA256SUMS').read_text().splitlines():
 h,name=s.split();ck(H(P/name)==h,'producer digest '+name)
# Full low range per-label witness sets, not only union membership.
label_targets=0
for y in Y:
 brute={}
 for v in range(97//y[2]+1):brute.setdefault(F(y,v),set()).add(v)
 for n in range(1,97):
  ck(solve(y,n)==brute.get(n,set()),'per-label inverse')
  label_targets+=1
# Random wide suffix law with independently calculated reversal.
rng=random.Random(231);identities=0
for k in range(1,41):
 M=1<<k
 for j in range(100):
  r=rng.randrange(M);u=rng.randrange(1<<48)
  C=sum(1<<(k-1-i) for i in range(k) if not (r>>i)&1)
  ck(g(M*u+r)==g(u)+C*u+int(r<M//2)*(u%2)+g(r)-1,'wide suffix')
  identities+=1
# Max-domain symbolic envelopes for every bank entry.
cap=10**15;maxM=maxL=maxQ=maxV=maxTerms=0
for a,b,c in Y:
 M=modulus(a,c,cap);L=(a+c)*M-a
 Q=a*(M-1)+c*M # safe envelope g(r)-1<=r; drop -b e
 V=cap//c+M
 terms=(a+c)*(V+1)+b
 maxM=max(maxM,M);maxL=max(maxL,L);maxQ=max(maxQ,Q);maxV=max(maxV,V);maxTerms=max(maxTerms,terms)
 ck(max(L,Q,V,terms,2732*M)<2**63,'64-bit envelope')
 ck(c*c*M*M<2**127,'128-bit threshold')
# Rebuild executable with undefined behavior sanitizer, without modifying producer.
binary=A/'inverse-ubsan'
subprocess.run(['g++','-std=c++17','-O1','-g','-fsanitize=undefined','-fno-sanitize-recover=undefined',str(P/'inverse.cpp'),'-o',str(binary)],check=True)
targets=list(range(1,4097))+[2**49+1,2**35+1]
lines=subprocess.check_output([str(binary),str(P/'bank.txt')]+list(map(str,targets)),text=True).splitlines()
ck(len(lines)==len(targets),'complete executable response')
image=set()
for y in Y:
 for v in range(4096//y[2]+1):image.add(F(y,v))
for n,line in zip(targets,lines):
 z=line.split();ck(int(z[0])==n,'target ordering')
 if n<=4096:ck((z[1]!='NONE')==(n in image),'union exhaustive')
 if z[1]!='NONE':
  a,b,c,v=map(int,z[1:5]);ck((a,b,c) in Y and v>=0 and F((a,b,c),v)==n,'executable witness')
for s in (P/'CASES.txt').read_text().splitlines():
 z=s.split()
 if z[1]!='NONE':
  n,a,b,c,v=map(int,z[:5]);ck((a,b,c) in Y and F((a,b,c),v)==n,'case exact witness')
 else:ck(int(z[0]) in [7,21] and int(z[0]) not in image,'case absence')
print(json.dumps(dict(status='PASS',repository_head=head,bank_sha256=H(P/'bank.txt'),bank_size=len(Y),per_label_target_checks=label_targets,wide_suffix_checks=identities,executable_targets=len(targets),small_odd_misses=[n for n in range(1,4097,2) if n not in image],max_domain_envelopes=dict(M=maxM,L=maxL,Q=maxQ,V=maxV,terms=maxTerms),scope='Finite checks and symbolic inverse proof only; no all-target occupancy theorem.'),indent=2))

from pathlib import Path
import csv,json,gzip,hashlib,random,subprocess,tempfile
P=Path(__file__).resolve().parent
B=P.parents[1]/'session189-c3-frontier'
def need(x,msg):
 if not x: raise ValueError(msg)
for name,digest in json.loads((P/'DEPENDENCIES.json').read_text()).items():
 need(hashlib.sha256((P.parents[1]/name).read_bytes()).hexdigest()==digest,'dependency changed: '+name)
with tempfile.TemporaryDirectory() as td:
 exe=Path(td)/'cycle_parity'
 subprocess.run(['g++','-O3',str(P/'cycle_parity.cpp'),'-o',str(exe)],check=True)
 for inp,outp in ((B/'bank_input.txt',P/'cycle_parity.tsv'),(P/'extra_input.txt',P/'extra_cycles.tsv')):
  actual=subprocess.check_output([str(exe)],input=inp.read_bytes())
  need(actual==outp.read_bytes(),'cycle bank mismatch: '+str(outp))
with tempfile.TemporaryDirectory() as td:
 exe=Path(td)/'core_bank'
 subprocess.run(['g++','-O3',str(B/'bank.cpp'),'-o',str(exe)],check=True)
 actual=subprocess.check_output([str(exe)],input=(P/'extra_input.txt').read_bytes())
 need(actual==(P/'extra_residues.tsv').read_bytes(),'extra core bank mismatch')
def rows(n,code):
 r=[0]*n; k=n*(n-1)//2-1
 for j in range(1,n):
  for i in range(j):
   a,b=(i,j) if code>>k&1 else (j,i);r[a]|=1<<b;k-=1
 return r
def covers(r):
 n=len(r);N=1<<n;d=[[0]*n for _ in range(N)];h=[0]*N;h[0]=1
 for i in range(n):d[1<<i][i]=1
 for s in range(1,N):
  for i in range(n):
   v=d[s][i]
   if not v:continue
   h[s]+=v;t=r[i]&~s
   while t:
    b=t&-t;t-=b;d[s|b][b.bit_length()-1]+=v
 q2=sum(h[s]*h[(N-1)^s] for s in range(1,N-1))
 cyc=0
 # Each Hamilton cycle has exactly one traversal beginning at vertex zero.
 e={(1,0):1}
 for s in range(1,N,2):
  for i in range(n):
   v=e.get((s,i),0)
   if s==N-1 and r[i]&1:cyc+=v
   t=r[i]&~s
   while t:
    b=t&-t;t-=b;k=(s|b,b.bit_length()-1);e[k]=e.get(k,0)+v
 return h[-1],q2,cyc
small=[]
for n in range(3,6):
 for code in range(1<<(n*(n-1)//2)):
  h,q,c=covers(rows(n,code));need((q-2*h-4+(4*c if n%2 else 0))%8==0,'small identity');small.append((n,code,h,q,c))
rng=random.Random(29)
for n in (6,7):
 for j in range(100):
  code=rng.randrange(1<<(n*(n-1)//2));h,q,c=covers(rows(n,code));need((q-2*h-4+(4*c if n%2 else 0))%8==0,'sample identity');small.append((n,code,h,q,c))
(P/'small_checks.json').write_text(json.dumps(small))
hs=list(csv.DictReader((B/'residues.tsv').open(),delimiter='\t'));cs=list(csv.DictReader((P/'cycle_parity.tsv').open(),delimiter='\t'))
extra=next(csv.DictReader((P/'extra_residues.tsv').open(),delimiter='\t'));ec=next(csv.DictReader((P/'extra_cycles.tsv').open(),delimiter='\t'))
queries=[tuple(map(int,l.split()))for l in (B/'bank_input.txt').read_text().splitlines()[1:13]]
codes={int(r['node']):int(r['code'])for r in csv.DictReader((B/'census.tsv').open(),delimiter='\t')}
out=[];physical=[];miss=[]
for h,c in zip(hs,cs):
 node=int(h['node']);need(node==int(c['node']),'row order');chosen=None
 for i in range(12):
  gap=(5*(int(h[f'S{i}'])-int(h[f'T{i}']))-4*int(c[f'D{i}']))%8
  if gap:chosen=(i,queries[i],h,c,gap);break
 if chosen is None:
  miss.append(node);need(node==56186,'unexpected miss');chosen=(20,(2,4),extra,ec,4)
 i,(a,b),h,c,gap=chosen
 sh,th,sc,tc=[int(t[f'{tag}{i}'])for t,tag in [(h,'S'),(h,'T'),(c,'S'),(c,'T')]]
 sa=(5*sh+4-4*sc)%8;ta=(5*th+4-4*tc)%8;need((sa-ta)%8==gap and gap,'bad gap')
 out.append([node,a,b,sh,th,sc,tc,sa,ta,gap]);q=rows(9,codes[node]);graphs=[]
 for op in (False,True):
  r=[0]*31
  for z in range(27):
   child,local=divmod(z,3);r[z]|=1<<(3*child+(local+1)%3)
   for other in range(9):
    if other!=child and bool(q[child]>>other&1)!=op:r[z]|=7<<(3*other)
   for port,mask in ((29,a),(30,b)):
    if mask>>child&1:r[z]|=1<<port
    else:r[port]|=1<<z
   r[z]|=1<<27;r[28]|=1<<z
  r[27]|=1<<28;r[29]|=(1<<30)|(1<<27);r[30]|=1<<27;r[28]|=(1<<29)|(1<<30)
  for z in range(31):
   need(not(r[z]>>z&1),'loop')
   for w in range(z):need(((r[z]>>w)^(r[w]>>z))&1,'not tournament')
  graphs.append(r)
 masks=[sum(7<<(3*j)for j in range(9)if m>>j&1)|(1<<28)for m in (a,b)]
 physical.append(dict(node=node,coarse_a=a,coarse_c=b,physical_a=masks[0],physical_c=masks[1],source_rows=graphs[0],target_rows=graphs[1],source_residue=sa,target_residue=ta))
need(len(out)==2348 and miss==[56186],'coverage')
with (P/'witnesses.tsv').open('w')as f:
 w=csv.writer(f,delimiter='\t');w.writerow(['node','a','c','core_source','core_target','cycle_source','cycle_target','ambient_source','ambient_target','gap']);w.writerows(out)
with (P/'physical_witnesses.jsonl.gz').open('wb')as f:
 with gzip.GzipFile(fileobj=f,mode='wb',mtime=0)as g:
  for r in physical:g.write((json.dumps(r,separators=(',',':'))+'\n').encode())
imports=[B/'residues.tsv',B/'bank_input.txt',B/'census.tsv',B/'PROOF.md',P.parents[1]/'session189-incoming-anchor-audit/INCOMING_PROOF.md']

print(json.dumps(dict(small_checks=len(small),cyclic_cases=len(out),initial_misses=miss,extra_query=[2,4],all_gaps_four=all(x[-1]==4 for x in out))))

from pathlib import Path
import hashlib,json
ROOT=Path(__file__).resolve().parents[1]
def need(c,m):
 if not c:raise RuntimeError(m)
def step(v,l):
 D,E,U,d=v
 return (D,0,2*U,d) if l=='O' else ((D-E+U,d,2*U-d,d) if l=='Z' else (D-E+U,D-E,D+U,D))
def run(w,v=(1,0,2,1)):
 for l in w:v=step(v,l)
 return v
def N(v,w):
 D,E,U,d=v;F,G,V,e=w
 return (D-E+U)*F-(D-E+d)*G+(D+3*U-d)*V-(U+E)*e

def decode(n,x):
 q0=n%x;q1=(n//x)%x;q2=n//(x*x)
 need(q1%2==0 and q2%4==0,'digit divisibility')
 A=3*q1-9*q0-2*q2;B=(3*q1-10*q0-2*q2)//2;C=3*q2//4-q1+3*q0
 if B%2:d=B if A>=3*B else A-B-C
 else:d=B//2 if 2*A<=3*B+2*C else A-B-C
 return(A-2*d,B-d,C+d,d)

def graph(w):
 rows=[2,0];root=0
 for l in w:
  n=len(rows);rows.append(0)
  for i in range(n):
   if l!='O' and i==root:rows[n]|=1<<i
   else:rows[i]|=1<<n
  if l=='J':root=n
 return rows,root

def graft(w,v):
 left,c=graph(w);m=len(left);left.append(1<<c)
 for i in range(m):
  if i!=c:left[i]|=1<<m
 n=len(left);right,r=graph(v);R=len(right)
 rows=left+[sum(((right[j]>>i)&1)<<(n+j) for j in range(R)) for i in range(R)]
 for i in range(n):
  for j in range(R):
   if j==r and i in (c,m):rows[n+j]|=1<<i
   else:rows[i]|=1<<(n+j)
 return rows

def literal(rows,only_cycle=False):
 n=len(rows);full=(1<<n)-1
 current={(1<<i,i):1 for i in ([0] if only_cycle else range(n))}
 for length in range(1,n):
  new={}
  for (mask,last),count in current.items():
   available=rows[last]&full^0
   available &= full^mask
   while available:
    bit=available&-available;available-=bit;j=bit.bit_length()-1
    key=(mask|bit,j);new[key]=new.get(key,0)+count
  current=new
 return sum(count for (mask,last),count in current.items() if not only_cycle or rows[last]&1)

for k in range(1,20):
 x=2**k;probe=run('O'*k+'JJ'+'O'*k)
 need(probe==(1+4*x,0,2*x+4*x*x,1+2*x),'physical probe formula')
 for i in range(4):
  v=tuple(int(i==j) for j in range(4));D,E,U,d=v
  need(N(v,probe)==4*(D+3*U-d)*x*x+(6*D-6*E+8*U-2*d)*x+D-2*E,'linear identity on basis')
bank={(1,0,2,1):''};todo=list(bank)
for v in todo:
 for l in 'OZJ':
  w=step(v,l)
  if w[2]<=512 and w not in bank:bank[w]=bank[v]+l;todo.append(w)
x=1<<((14*512-1).bit_length());probe=run('O'*(x.bit_length()-1)+'JJ'+'O'*(x.bit_length()-1))
image={}
for v in bank:
 n=N(v,probe);need(decode(n,x)==v,'single-context decoding')
 need(n not in image,'context injectivity');image[n]=v
 r0,r1,r2=[N(v,run(w)) for w in ('','J','JJ')]
 D,E,U,d=v
 need((3*r0+2*r1-2*r2,(-r0+6*r1-3*r2)//2,(-3*r0+r2)//2)==(D+2*d,E+d,U-d),'three contexts')
# Independent literal graph validation, including the full bounded selector.
controls=[]
for v,w in [('', 'O'*5+'JJ'+'O'*5),('OZJ','OJJO'),('ZOZ','OJJO'),('OJ',''),('ZJ',''),('OJ','J'),('ZJ','J'),('OOJJOOZ',''),('ZOJOZOZ',''),('OOJJOOZ','JJ'),('ZOJOZOZ','JJ'),('OJOZ','J'),('ZJOO','J'),('OJOZ','JJ'),('ZJOO','JJ')]:
 rows=graft(v,w);h=literal(rows);c=literal(rows,True)
 need(h==N(run(v),run(w)) and c==2,'literal C2 graph')
 controls.append({'left':v,'right':w,'order':len(rows),'H':h,'C':c})
for words,partners in [(('OJ','ZJ'),('','J')),(('OOJJOOZ','ZOJOZOZ'),('','JJ')),(('OJOZ','ZJOO'),('J','JJ'))]:
 need(run(words[0])!=run(words[1]),'irredundant pair distinct')
 need(all(N(run(words[0]),run(w))==N(run(words[1]),run(w)) for w in partners),'irredundant pair collision')
report={'status':'PASS','proof_sha256':hashlib.sha256((ROOT/'PROOF.md').read_bytes()).hexdigest(),'script_sha256':hashlib.sha256((ROOT/'verify.py').read_bytes()).hexdigest(),'complete_bank_U_cap':512,'source_profiles':len(bank),'single_context_radix':x,'basis_identity_k1_to19':True,'literal_graph_controls':controls,'manual_audit':'All-order probe law, no-carry bounds, exact signature inverse, graph order, cubic size bound and natural-bank irredundancy are valid. The actual-source assumption is essential and correctly retained.','scope':'Response distinguishability, not graph isomorphism or target coverage.'}
Path(__file__).with_name('INDEPENDENT.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))

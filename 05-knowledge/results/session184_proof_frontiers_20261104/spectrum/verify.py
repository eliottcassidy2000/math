"""Dependency-free exact singleton/C3 grammar return audit; no assert statements."""
import itertools, json, hashlib
from pathlib import Path

def need(b,m):
 if not b: raise RuntimeError(m)
def signature(word):
 out=None
 for kind,a in word:
  z=(1,2,a,a,0,a,a) if kind==1 else (3,12,3*a,3*a,a*(a-1),a,a)
  if out is None:out=z;continue
  h,p,f,g,j,e,s=out;H,P,F,G,J,E,S=z
  out=(h*H,p*P,f*H+p*F,g*P+h*G,j*H+h*J+e*S+g*F,h*E,s*H)
 return out

def words_with_F(F):
 """All strong-attachment component words at exact positive F, all orders."""
 k=0
 while F%3**k==0:
  u=F//3**k
  def rec(pos,left,w):
   if pos>=u.bit_length():
    if left==0 and w[-1][1]>0:yield tuple(w)
    return
   yield from rec(pos+1,left,w+[(1,(u>>pos)&1)])
   if left:yield from rec(pos+2,left-1,w+[(3,(u>>pos)&3)])
  yield from rec(0,k,[])
  k+=1

def normalized_candidates(F):
 while F%2==0:F//=2
 for G in (F,2*F):
  for w in words_with_F(G):
   if w[0][1]:yield w

def graph(word):
 rows=[];ports=set()
 for kind,a in word:
  n=len(rows);block=[0] if kind==1 else [2,4,1]
  rows=[x|(((1<<kind)-1)<<n) for x in rows]+[x<<n for x in block]
  ports.update(range(n,n+a))
 return rows,ports

def tables(rows):
 n=len(rows);N=1<<n;h=[0]*N;h[0]=1;e=[[0]*n for _ in range(N)]
 for v in range(n):e[1<<v][v]=1
 for S in range(1,N):
  for v in range(n):
   if S>>v&1:
    T=S^(1<<v)
    if T:e[S][v]=sum(e[T][u] for u in range(n) if rows[u]>>v&1)
  h[S]=sum(e[S])
 return h,e

def direct_signature(rows,ports):
 n=len(rows);h,e=tables(rows);reverse=[sum(1<<j for j in range(n) if rows[j]>>i&1) for i in range(n)];_,s=tables(reverse);N=(1<<n)-1
 ep=[sum(v[i] for i in ports) for v in e];sp=[sum(v[i] for i in ports) for v in s]
 return (h[N],sum(h[x]*h[N^x] for x in range(N+1)),sum(ep[x]*h[N^x] for x in range(N+1)),sum(h[x]*sp[N^x] for x in range(N+1)),sum(ep[x]*sp[N^x] for x in range(N+1)),ep[N],sp[N])

def attachment(word,t):
 rows,W=graph(word);b=len(rows);pref=[sum(1<<j for j in range(i+1,t))|(((1<<b)-1)<<t) for i in range(t)]
 rows=pref+[v<<t for v in rows];root=t+b;rows.append((1<<t)-1)
 for i in range(b):
  if i in W:rows[t+i]|=1<<root
  else:rows[root]|=1<<(t+i)
 return rows

def strong(rows):
 N=(1<<len(rows))-1
 for v in range(len(rows)):
  S=1<<v
  while True:
   T=S
   for i in range(len(rows)):
    if S>>i&1:T|=rows[i]
   if S==T:break
   S=T
  if S!=N:return False
 return True

def all_words(limit,prefix=()):
 if prefix:yield prefix
 n=sum(x for x,_ in prefix)
 for kind in (1,3):
  if n+kind<=limit:
   for a in range(kind+1):yield from all_words(limit,prefix+((kind,a),))

def run():
 # Independent subset DP checks every component word through order seven.
 total=0;eligible=0;inverse=0;vals=0
 for w in all_words(7):
  rows,W=graph(w);z=signature(w);need(z==direct_signature(rows,W),'direct seven-coordinate mismatch')
  total+=1;k=sum(kind==3 for kind,a in w)
  if k>=2:
   A=sum(a*b for (typ,a),(typ2,b) in zip(w,w[1:]) if typ==typ2==3)
   c=z[0]-z[4];need(c%3**(k-2)==0,'valuation lower bound');need((c//3**(k-2)+A)%3==0,'leading coefficient');vals+=1
  if k==1:
   i=next(i for i,x in enumerate(w) if x[0]==3);a=w[i][1];b=w[i-1][1] if i else 0;c=w[i+1][1] if i+1<len(w) else 0
   need((z[0]-z[4]+a*(a-1+b+c))%3==0,'single triangle coefficient')
  if w[-1][1]:
   need(w in set(words_with_F(z[2])),'inverse omission');inverse+=1
   # All actual root attachments at TT1 are tested through suffix order five.
   if len(rows)<=5:
    r=attachment(w,1);need(strong(r),'terminal-port strongness');need(tables(r)[0][-1]==z[0]-z[4]+2*z[2],'full count');eligible+=1
 # Audit full return tree with every normalized target, not finite graph order cutoff.
 root=((1,1),(1,1),(1,1));queue=[('',root)];nodes=[];checks=0; audited={}
 expected={'','2','4','22','42','422','4222'}
 while queue:
  address,w=queue.pop(0);z=signature(w);c=z[0]-z[4];children={};candidate_counts={}
  for d in (2,4):
   candidates=list(normalized_candidates(3*z[2]));candidate_counts[str(d)]=len(candidates)
   found=[]
   for v in candidates:
    zz=signature(v)
    if v not in audited:
     rows,W=graph(v);need(zz==direct_signature(rows,W),'candidate literal signature')
     audited[v]={'word':[list(x) for x in v],'signature':list(zz),'c':zz[0]-zz[4]}
    need(zz[2] in (3*z[2],6*z[2]),'unexpected normalized F in tree');checks+=1
    if zz[0]-zz[4]==3*c+d:found.append(v)
   children[str(d)]=[address+str(d) for v in found]
   for v in found:queue.append((address+str(d),v))
  for t in (1,2):
   r=attachment(w,t);need(strong(r),'node strong');need(tables(r)[0][-1]==c+(1<<t)*z[2],'node count')
  nodes.append({'address':address,'word':[list(x) for x in w],'signature':list(z),'c':c,'H_TT1':c+2*z[2],'children':children,'normalized_candidate_counts':candidate_counts})
 need({x['address'] for x in nodes}==expected,'tree differs')
 # Large arbitrary words validate closed-form leading term without graph extrapolation.
 large=0
 for k in range(2,31):
  for a in (1,2,3):
   w=tuple((3,a) for _ in range(k));z=signature(w);c=z[0]-z[4];A=(k-1)*a*a
   need((c//3**(k-2)+A)%3==0 and c%3**(k-2)==0,'large valuation');large+=1
 result={'status':'PASS','direct_signature_words_order_at_most_7':total,'direct_full_TT1_attachments_suffix_order_at_most_5':eligible,'inverse_roundtrips':inverse,'small_valuation_controls':vals,'large_valuation_controls':large,'return_candidate_tests':checks,'independent_literal_candidate_audits':len(audited),'return_nodes':nodes,'scope':'Complete all-length, arbitrary fixed clock-shift return tree in singleton/C3 ordered-component grammar only.'}
 Path(__file__).with_name('candidates.json').write_text(json.dumps(list(audited.values()),indent=2)+'\n')
 Path(__file__).with_name('receipt.json').write_text(json.dumps(result,indent=2)+'\n')
 print(json.dumps(result,sort_keys=True))
if __name__=='__main__':run()

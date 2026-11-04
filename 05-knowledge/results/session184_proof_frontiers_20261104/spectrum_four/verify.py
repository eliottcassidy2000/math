"""All-order SCC<=4 return certificate. Standard library; explicit checks."""
import itertools,json
from pathlib import Path

def need(b,m):
 if not b:raise RuntimeError(m)
def mul(x,y):
 if x is None:return y
 h,p,f,g,j,e,s=x;H,P,F,G,J,E,S=y
 return(h*H,p*P,f*H+p*F,g*P+h*G,j*H+h*J+e*S+g*F,h*E,s*H)
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
def direct(rows,W):
 n=len(rows);h,e=tables(rows);rev=[sum(1<<j for j in range(n) if rows[j]>>i&1) for i in range(n)];_,s=tables(rev);N=(1<<n)-1
 ep=[sum(x[i] for i in W) for x in e];sp=[sum(x[i] for i in W) for x in s]
 return(h[N],sum(h[x]*h[N^x] for x in range(N+1)),sum(ep[x]*h[N^x] for x in range(N+1)),sum(h[x]*sp[N^x] for x in range(N+1)),sum(ep[x]*sp[N^x] for x in range(N+1)),ep[N],sp[N])
def strong(rows):
 N=(1<<len(rows))-1
 for v in range(len(rows)):
  S=1<<v
  while True:
   T=S
   for i in range(len(rows)):
    if S>>i&1:T|=rows[i]
   if T==S:break
   S=T
  if S!=N:return False
 return True
Q=(6,12,8,1)
Q_SIG=[(5,32,0,0,0,0,0),(5,32,5,7,0,1,1),(5,32,6,9,0,1,2),(5,32,11,16,4,2,3),(5,32,9,6,0,2,1),(5,32,14,13,2,3,2),(5,32,15,15,6,3,3),(5,32,20,22,12,4,4),(5,32,7,5,0,1,1),(5,32,12,12,4,2,2),(5,32,13,14,2,2,3),(5,32,18,21,10,3,4),(5,32,16,11,4,3,2),(5,32,21,18,10,4,3),(5,32,22,20,12,4,4),(5,32,27,27,22,5,5)]
BANK={**{('S',a):(1,2,a,a,0,a,a) for a in (0,1)},**{('C',a):(3,12,3*a,3*a,a*(a-1),a,a) for a in range(4)},**{('Q',a):v for a,v in enumerate(Q_SIG)}}
SIZE={'S':1,'C':3,'Q':4}
def signature(w):
 z=None
 for b in w:z=mul(z,BANK[b])
 return z
def graph(w):
 rows=[];W=set()
 for typ,a in w:
  r=(0,) if typ=='S' else (2,4,1) if typ=='C' else Q;n=len(rows);m=len(r)
  rows=[x|(((1<<m)-1)<<n) for x in rows]+[x<<n for x in r]
  U=range(a) if typ!='Q' else (i for i in range(4) if a>>i&1)
  W.update(n+i for i in U)
 return rows,W
def v2(x):
 need(x>0,'positive valuation input');return (x&-x).bit_length()-1

def search(F,c):
 hits=[];visited=0;terminal_candidates=0;maxblocks=0
 def rec(w,z):
  nonlocal visited,terminal_candidates,maxblocks
  visited+=1;maxblocks=max(maxblocks,len(w))
  for b,v in BANK.items():
   if not w and b[1]==0:continue
   zz=mul(z,v)
   if zz[2]>F:continue
   ww=w+(b,)
   if zz[2]==F and b[1]:
    terminal_candidates+=1
    if zz[0]-zz[4]==c:hits.append((ww,zz))
   if zz[1]<=F:rec(ww,zz)
 rec((),None)
 return hits,{'visited_prefixes':visited,'terminal_candidates':terminal_candidates,'maximum_recursed_blocks':maxblocks}
def all_words(n,w=()):
 if w:yield w
 used=sum(SIZE[t] for t,a in w)
 for b in BANK:
  if used+SIZE[b[0]]<=n:yield from all_words(n,w+(b,))

def run():
 for mask,z in enumerate(Q_SIG):need(direct(Q,{i for i in range(4) if mask>>i&1})==z,'Q port table')
 count4=0
 for bits in range(64):
  r=[0]*4
  for j,(u,v) in enumerate(itertools.combinations(range(4),2)):
   if bits>>j&1:r[u]|=1<<v
   else:r[v]|=1<<u
  if strong(r):
   count4+=1;need(any(all(bool(r[i]>>j&1)==bool(Q[p[i]]>>p[j]&1) for i in range(4) for j in range(4)) for p in itertools.permutations(range(4))),'strong4 type missing')
 need(count4==24,'strong4 count')
 controls=0;normalized=0
 for w in all_words(7):
  z=signature(w);r,W=graph(w);need(z==direct(r,W),'literal mixed-word signature');controls+=1
  if w[0][1]:need(v2(z[2])<=4 and v2(z[2])<v2(z[1]),'normalized valuation bound');normalized+=1
 # Seven audited source representatives are included literally, independent of other files.
 sources=[('',(('S',1),('S',1),('S',1))),('2',(('S',1),('S',1),('C',1))),('4',(('S',1),('S',0),('S',1),('S',0),('S',1))),('22',(('S',1),('S',0),('S',1),('C',2))),('42',(('S',1),('S',0),('C',1),('S',1))),('422',(('S',1),('S',0),('C',1),('C',1))),('4222',(('C',1),('S',1),('S',1),('S',1),('S',1),('S',0),('S',1)))]
 source_map=dict(sources);transitions=[];calls=0;visited=0
 for address,w in sources:
  z=signature(w)
  for d in (2,4):
   c=3*(z[0]-z[4])+d;F=3*z[2]
   while F%2==0:F//=2
   hits=[];scans=[];r=0
   while c%5**r==0:
    for a in range(5):
     hh,stats=search(F*2**a,c//5**r);calls+=1;visited+=stats['visited_prefixes']
     scans.append({'unmarked_Q_prefixes':r,'two_adic_scale':a,'F':F*2**a,'c':c//5**r,**stats})
     for word,zz in hh:hits.append({'unmarked_Q_prefixes':r,'word':[list(b) for b in word],'signature':list(zz)})
    r+=1
   next_address=address+str(d)
   if next_address in source_map:
    need(len(hits)==1 and hits[0]['unmarked_Q_prefixes']==0 and hits[0]['word']==[list(b) for b in source_map[next_address]],'expected unique return')
   else:need(not hits,'unexpected additional return')
   transitions.append({'source_address':address,'digit':d,'target_constant':c,'hits':hits,'complete_scans':scans})
 result={'status':'PASS','strong_labelled_order4_tournaments':count4,'literal_signature_words_order_at_most7':controls,'normalized_word_valuation_controls':normalized,'complete_target_scans':calls,'visited_prefixes':visited,'address_language':[a for a,w in sources],'transitions':transitions,'scope':'All suffixes whose SCCs have order at most4; all lengths and arbitrary fixed clock shifts.'}
 Path(__file__).with_name('receipt.json').write_text(json.dumps(result,indent=2)+'\n')
 print(json.dumps(result,sort_keys=True))
if __name__=='__main__':run()

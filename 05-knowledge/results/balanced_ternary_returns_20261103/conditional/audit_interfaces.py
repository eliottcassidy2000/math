#!/usr/bin/env python3
"""Independent graph construction and subset/permutation audit of I3--I4."""
from pathlib import Path
from itertools import permutations
import hashlib,json
HERE=Path(__file__).resolve().parent

def need(x,msg):
 if not x:raise RuntimeError(msg)

def source(m):
 # Independent construction in chronological physical order x,y,c1,...,cm,a,b,p,q.
 n=m+6;o=[0]*n
 def edge(a,b):o[a]|=1<<b
 x,y=0,1;cs=list(range(2,m+2));a,b,p,q=m+2,m+3,m+4,m+5
 edge(x,y)
 for c in cs:edge(y,c);edge(c,x)
 for i,c in enumerate(cs):
  for d in cs[i+1:]:edge(c,d)
 # Re-mark cm,x and append letters3,2 after y,c1,...,c(m-1).
 for c in [y]+cs[:-1]:edge(c,a);edge(c,b)
 edge(cs[-1],a);edge(x,a);edge(x,b);edge(b,cs[-1]);edge(a,b)
 for v in range(p):edge(x,p) if v==x else edge(p,v)
 for v in range(q):edge(p,q) if v==p else edge(q,v)
 for u in range(n):
  for v in range(u+1,n):need(((o[u]>>v&1)+(o[v]>>u&1))==1,'invalid tournament')
 return o,cs[-1],p,q

def data(o):
 n=len(o);full=(1<<n)-1;end=[[0]*n for _ in range(full+1)];H=[0]*(full+1);H[0]=1
 rev=[sum(1<<u for u in range(n) if o[u]>>v&1) for v in range(n)]
 for S in range(1,full+1):
  for v in range(n):
   if S>>v&1:
    old=S^(1<<v)
    if not old:end[S][v]=1
    else:
     pre=old&rev[v];total=0
     while pre:
      bit=pre&-pre;pre-=bit;total+=end[old][bit.bit_length()-1]
     end[S][v]=total
  H[S]=sum(end[S])
 F=[sum(end[S][v]*H[full^S] for S in range(1,full+1) if S>>v&1) for v in range(n)]
 # Endpoint start counts from full first-vertex recursion.
 back=[[0]*n for _ in range(full+1)]
 for S in range(1,full+1):
  for v in range(n):
   if S>>v&1:
    old=S^(1<<v)
    if not old:back[S][v]=1
    else:
     suc=old&o[v];total=0
     while suc:
      bit=suc&-suc;suc-=bit;total+=back[old][bit.bit_length()-1]
     back[S][v]=total
 start=back[full]
 return H[full],start,F

def literal(o,r):
 n=len(o);H=0;F=0;start=0
 for p in permutations(range(n)):
  good=[bool(o[p[i]]>>p[i+1]&1) for i in range(n-1)]
  ham=all(good);H+=ham;start+=ham and p[0]==r
  k=p.index(r)+1
  F+=all(good[:k-1]) and all(good[k:])
 return H,start,F

def surj(k,r):
 return sum((-1)**(r-i)*__import__('math').comb(r,i)*i**k for i in range(r+1))

def main():
 rows=[];fs=[];starts=[]
 for m in range(2,7):
  o,r,p,q=source(m);h,s,F=data(o)
  need(h==15*2**m-1,'old H formula')
  need(s[r]==1,'unique-start old clone')
  need(F[r]==26*3**m+8*2**m-25,'new cover formula')
  need(s[0]==3,'old physical x start count')
  # Independently recount the one-ready newborn using the full enlarged graph.
  new=o+[0];z=len(o)
  for v in range(z):
   if v==r:new[v]|=1<<z
   else:new[z]|=1<<v
  nh,ns,nF=data(new)
  need(nh==26*3**m+23*2**m-27 and nF[z]-nh==1,'fresh ready child')
  need(nh>3*h+4,'fresh repair overshoots')
  if m==2:need(literal(o,r)==(h,s[r],F[r]),'literal permutation crosscheck')
  indeg=[sum(o[u]>>v&1 for u in range(len(o))) for v in range(len(o))]
  need(all(indeg[v]>=3 for v in range(len(o)-2)),'all old exterior vertices degree3')
  need(indeg[p]==indeg[q]==1,'reset root indegrees')
  rows.append({'m':m,'H':h,'old_clone_start':s[r],'old_clone_F':F[r],
               'fresh_H':nh,'fresh_F':nF[z],'overshoot_plus4':nh-(3*h+4),
               'rows':o,'clone_root':r,'reset_root':p,'reset_child':q})
  fs.append(F[r]);starts.append(s[r])
 coeff=[]
 for k,val in enumerate(fs,1):
  rem=val-sum(coeff[r-1]*surj(k,r) for r in range(1,k))
  d=surj(k,k);need(rem%d==0,'nonintegral skeleton coefficient');coeff.append(rem//d)
 need(coeff==[241,250,78,0,0],'five-run skeleton coefficients')
 note=HERE.parent/'interfaces'/'NOTE.md'
 print(json.dumps({'status':'PASS','independent_scope':'No producer imports; full labelled graph construction; subset DP m2..6 and all-permutation m2; five-run coefficients reconstructed. All-order scope needs the reviewed skeleton bound.',
 'graphs':rows,'five_run_coefficients':coeff,'start_samples':starts,
 'producer_note_sha256':hashlib.sha256(note.read_bytes()).hexdigest(),
 'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},indent=2,sort_keys=True))
if __name__=='__main__':main()

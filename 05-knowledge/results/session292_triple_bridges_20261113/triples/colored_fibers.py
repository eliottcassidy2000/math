"""Relative cyclic-triple fibres: direct labelled grouping vs modular-tree formula."""
from itertools import combinations,product
from collections import Counter
from math import factorial
from pathlib import Path
import json

def need(c,m):
 if not c: raise RuntimeError(m)
def rows(n,code):
 a=[0]*n
 for k,(i,j) in enumerate(combinations(range(n),2)):
  if code>>k&1:a[i]|=1<<j
  else:a[j]|=1<<i
 return a

def tree(a,mask):
 if mask.bit_count()==1:return ('leaf',mask,[])
 modules=[]
 sub=mask
 while sub:
  if all((a[w]&sub) in (0,sub) for w in range(len(a)) if (mask>>w&1) and not(sub>>w&1)):modules.append(sub)
  sub=(sub-1)&mask
 strong=[s for s in modules if all(not(s&t) or s&t in (s,t) for t in modules)]
 children=[s for s in strong if s!=mask and not any(s!=t and s&t==s and t!=mask for t in strong)]
 need(sum(s.bit_count() for s in children)==mask.bit_count(),'partition')
 rep=lambda s:(s&-s).bit_length()-1
 scores={s:sum(bool(a[rep(s)]>>rep(t)&1) for t in children if t!=s) for s in children}
 linear=sorted(scores.values())==list(range(len(children)))
 if linear: children.sort(key=lambda s:-scores[s])
 return ('linear' if linear else 'prime',mask,[tree(a,s) for s in children])

def color(mask,cols):
 vals={cols[v] for v in range(len(cols)) if mask>>v&1}
 return next(iter(vals)) if len(vals)==1 else None

def predicted(t,cols):
 kind,mask,kids=t
 val=1
 for c in kids:val*=predicted(c,cols)
 if kind=='prime':return val*(2 if color(mask,cols) is not None else 1)
 if kind=='linear':
  prev=None;run=0
  for child in kids:
   cur=color(child[1],cols)
   if cur is not None and cur==prev:run+=1
   else:
    val*=factorial(run);run=1 if cur is not None else 0;prev=cur
  val*=factorial(run)
 return val

def engineered(t,cols):
 kind,mask,kids=t
 if kind=='leaf':cols[(mask&-mask).bit_length()-1]=0;return
 for child in kids:
  if child[0]!='leaf':engineered(child,cols)
 if kind=='linear':
  bit=0
  for child in kids:
   if child[0]=='leaf':cols[(child[1]&-child[1]).bit_length()-1]=bit;bit^=1
   else:bit=0
 else:
  for child in kids:
   if child[0]=='leaf':cols[(child[1]&-child[1]).bit_length()-1]=0
  if all(child[0]=='leaf' for child in kids):
   child=kids[0];cols[(child[1]&-child[1]).bit_length()-1]=1

def events(t):
 kind,mask,kids=t
 out=[]
 for child in kids:out.extend(events(child))
 if kind=='prime':out.append(mask)
 if kind=='linear':out.extend(a[1]|b[1] for a,b in zip(kids,kids[1:]))
 return out

def main():
 totals=[]
 for n in range(1,6):
  pairs=list(combinations(range(n),2));triples=list(combinations(range(n),3));data=[]
  for code in range(1<<len(pairs)):
   a=rows(n,code);cyc=0
   for q,(i,j,k) in enumerate(triples):
    if sum(bool(a[u]>>v&1) for u,v in ((i,j),(j,k),(k,i))) in (0,3):cyc|=1<<q
   t=tree(a,(1<<n)-1)
   cols=[-1]*n;engineered(t,cols)
   need(all(c in (0,1) for c in cols) and predicted(t,cols)==1,'engineered discriminator')
   need(len(events(t))<=n-1,'union-bound events')
   data.append((code,cyc,t))
  tested=0;maxfiber=0;random_counts={}
  palettes=(2,) if n==5 else (2,3)
  for q in palettes:
   unique=[0]*len(data)
   for cols in product(range(q),repeat=n):
    observed=Counter((cyc,code&sum(1<<j for j,(u,v) in enumerate(pairs) if cols[u]!=cols[v])) for code,cyc,t in data)
    cmask=sum(1<<j for j,(u,v) in enumerate(pairs) if cols[u]!=cols[v])
    for index,(code,cyc,t) in enumerate(data):
     v=predicted(t,cols);need(v==observed[cyc,code&cmask],('fiber',n,code,cols,v,observed[cyc,code&cmask]))
     tested+=1;maxfiber=max(maxfiber,v);unique[index]+=int(v==1)
     need((v>1)==any(color(mask,cols) is not None for mask in events(t)),'event iff')
   for index,(code,cyc,t) in enumerate(data):
    need(q**n-unique[index]<=(n-1)*q**(n-1),'random union bound')
    if cyc==0:need(unique[index]==q*(q-1)**(n-1),'transitive exact probability')
   random_counts[q]=dict(min_unique=min(unique),denominator=q**n)
  totals.append(dict(n=n,tournaments=len(data),palettes=palettes,tested=tested,max_fiber=maxfiber,random_counts=random_counts,engineered_ports=len(data)))
 out={'status':'PASS','counts':totals,'total':sum(x['tested'] for x in totals)}
 Path(__file__).with_name('colored_fibers.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()

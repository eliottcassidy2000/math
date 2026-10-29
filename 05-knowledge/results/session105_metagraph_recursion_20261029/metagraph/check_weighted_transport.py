#!/usr/bin/env python3
"""Exact native metagraph weights and strip-extension obstruction; stdlib only."""
import itertools as it
import json
import math
from collections import Counter,defaultdict
from fractions import Fraction

checks=0
def require(ok,msg):
 global checks
 checks+=1
 if not ok: raise RuntimeError(msg)

def pairs(n): return list(it.combinations(range(n),2))
def arc(x,u,v,ix): return u<v and not (x>>ix[u,v]&1) or u>v and bool(x>>ix[v,u]&1)
def relabel(x,p,ps,ix):
 return sum((not arc(x,p[u],p[v],ix))<<k for k,(u,v) in enumerate(ps))
def bank(n):
 ps=pairs(n); ix={p:i for i,p in enumerate(ps)}; perms=list(it.permutations(range(n)))
 canon={}; aut={}
 for x in range(1<<len(ps)):
  if x in canon: continue
  orbit={relabel(x,p,ps,ix) for p in perms}
  c=min(orbit); aut[c]=math.factorial(n)//len(orbit)
  for y in orbit: canon[y]=c
 return ps,ix,perms,canon,aut

def fixed_words(n):
 ps=pairs(n); off=[k for k,(u,v) in enumerate(ps) if v!=u+1]
 for t in range(1<<len(off)):
  yield sum(((t>>j)&1)<<k for j,k in enumerate(off))
def toggle(n): return sum(1<<k for k,(u,v) in enumerate(pairs(n)) if v!=u+1)
def blue(x,n,ix): return all(arc(x,u,v,ix)==arc(x,n-1-v,n-1-u,ix) for u,v in pairs(n))
def extend(x,n,bits):
 old=pairs(n); nxt={p:k for k,p in enumerate(pairs(n+1))}
 y=sum(((x>>k)&1)<<nxt[p] for k,p in enumerate(old))
 return y|sum(((bits>>i)&1)<<nxt[i,n] for i in range(n-1))
def keycount(c): return [[*k,v] for k,v in sorted(c.items())]

def main():
 banks={n:bank(n) for n in range(2,7)}; result={}; allM={}
 for n,(ps,ix,perms,canon,aut) in banks.items():
  M=Counter(); B=Counter(); fibres=Counter(); k=toggle(n)
  for x in fixed_words(n):
   a,b=canon[x],canon[x^k]; M[a,b]+=1; fibres[a]+=1
   if blue(x,n,ix): B[a,b]+=1
   require(blue(x,n,ix)==blue(x^k,n,ix),'native blue preservation')
  W=Counter(); WB=Counter(); H={}
  for a in aut:
   h=0
   for p in perms:
    if not all(arc(a,p[i],p[i+1],ix) for i in range(n-1)): continue
    path=sum(1<<ix[min(u,v),max(u,v)] for u,v in zip(p,p[1:]))
    y=a^(((1<<len(ps))-1)^path); b=canon[y]
    W[a,b]+=1; h+=1
    xp=relabel(a,p,ps,ix)
    if blue(xp,n,ix): WB[a,b]+=1
   H[a]=h
   require(h==aut[a]*fibres[a],'fibre mass')
  for a in aut:
   for b in aut:
    require(M[a,b]==M[b,a],'symmetry')
    require(B[a,b]==B[b,a],'blue symmetry')
    require(W[a,b]==aut[a]*M[a,b],'labelled witness weight')
    require(WB[a,b]==aut[a]*B[a,b],'blue labelled witness weight')
    require(Fraction(W[a,b],aut[a])==Fraction(W[b,a],aut[b]),'detailed balance')
  require(sum(M.values())==1<<((n-1)*(n-2)//2),'total mass')
  require(sum(Fraction(1,aut[a]) for a in aut)*math.factorial(n)==1<<len(ps),'labelled orbit mass')
  # Exact Markov detailed balance, pi proportional H/Aut.
  for a,b in M:
   require(Fraction(fibres[a]*W[a,b],H[a])==Fraction(fibres[b]*W[b,a],H[b]),'Markov detailed balance')
  # Fixed proposal clock, not conditioning on legality.
  proposal={}
  for a in aut:
   row={b:Fraction(W[a,b],math.factorial(n))+(Fraction(1)-Fraction(H[a],math.factorial(n)) if a==b else 0) for b in aut}
   proposal[a]=row
   require(sum(row.values())==1,'proposal row stochastic')
   require(row[a]>=Fraction(1,2),'proposal at least half lazy')
   for b in aut:
    require(Fraction(row[b],aut[a])==Fraction(W[b,a],math.factorial(n)*aut[b])+(Fraction(1,aut[b])-Fraction(H[b],math.factorial(n)*aut[b]) if a==b else 0),'proposal balance')
  if n==3:
   require(proposal[0][2]==Fraction(1,6) and proposal[2][0]==Fraction(1,2),'n3 proposal clock hostile')
  result[str(n)]={'classes':len(aut),'word_mass':sum(M.values()),'blue_word_mass':sum(B.values()),'support_ordered_entries':len(M),'diagonal_mass':sum(M[a,a] for a in aut),'class_data':[[a,aut[a],H[a],fibres[a]] for a in sorted(aut)],'M':keycount(+M),'B':keycount(+B)}
  allM[n]=M
 obstructions=[]; stronger=[]
 for n in range(2,6):
  cn=banks[n][3]; cp=banks[n+1][3]; k=toggle(n); kp=toggle(n+1)
  grouped=defaultdict(list); total=Counter()
  for x in fixed_words(n):
   sig=Counter()
   for bits in range(1<<(n-1)):
    y=extend(x,n,bits); z=extend(x^k,n,bits^((1<<(n-1))-1))
    require(z==y^kp,'native strip intertwining')
    sig[cp[y],cp[z]]+=1; total[cp[y],cp[z]]+=1
   grouped[cn[x],cn[x^k]].append((x,sig))
  require(total==allM[n+1],'exact transport recurrence')
  found=None
  for key,rows in sorted(grouped.items()):
   for i,(x,s) in enumerate(rows):
    for y,t in rows[i+1:]:
     if s!=t:
      found={'n':n,'pair_class':list(key),'word_masks':[x,y],'extensions_x':keycount(s),'extensions_y':keycount(t)};break
    if found:break
   if found:break
  if found: obstructions.append(found)
  # A stronger hostile: forget the order of the endpoint classes and merge
  # each class with its converse, both before and after extension.
  def merged(c,order):
   cc=banks[order][3]
   return min(c,cc[c^((1<<len(pairs(order)))-1)])
  folded=defaultdict(list)
  for key,rows in grouped.items():
   fkey=tuple(sorted(merged(c,n) for c in key))
   for x,sig in rows:
    fs=Counter()
    for (a,b),v in sig.items():fs[tuple(sorted((merged(a,n+1),merged(b,n+1))))]+=v
    folded[fkey].append((x,fs))
  sf=None
  for key,rows in sorted(folded.items()):
   for i,(x,s) in enumerate(rows):
    for y,t in rows[i+1:]:
     if s!=t:
      sf={'n':n,'folded_pair_class':list(key),'word_masks':[x,y],'extensions_x':keycount(s),'extensions_y':keycount(t)};break
    if sf:break
   if sf:break
  if sf:stronger.append(sf)
 require([result[str(n)]['classes'] for n in range(2,7)]==[1,2,4,12,56],'complete small class bank')
 require(obstructions[0]['n']==4 and obstructions[0]['word_masks']==[4,18],'minimal ordinary obstruction')
 require(stronger[0]['n']==4 and stronger[0]['word_masks']==[2,6],'minimal folded obstruction')
 result['pair_class_extension_obstructions']=obstructions
 result['converse_unordered_extension_obstructions']=stronger
 result['checks']=checks
 print(json.dumps(result,sort_keys=True,indent=2))
if __name__=='__main__':main()

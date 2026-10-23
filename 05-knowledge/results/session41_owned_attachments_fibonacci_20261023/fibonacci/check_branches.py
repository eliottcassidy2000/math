#!/usr/bin/env python3
"""Independent branch expansion, fixed blocks versus occurrence-dependent blocks."""
import json
from functools import lru_cache
from itertools import product
F=[0,1,1]
for _ in range(3,24):F.append(F[-1]+F[-2])
def need(x,m):
 if not x:raise RuntimeError(m)
def Z(n):
 s=[]
 for i in range(23,1,-1):
  if F[i]<=n:s.append(i);n-=F[i]
 return tuple(sorted(s))
def r(n):return 2 in Z(n)
def fibre(n,R):
 out=[Z(n)]
 if n>=1:
  out.append((-1,)+Z(n-1))
  if not r(n-1):out.append((1,)+Z(n-1))
 if R and n>=2 and not r(n-2):out.append((-1,1)+Z(n-2))
 return out
def weight(i):return 1 if i==-1 else F[i]
W={-1:'A',1:'B',2:'C'}
for i in range(3,22):
 p=tuple(sorted(list(range(i-1,1,-2))+[-1 if i%2 else 1]))
 W[i]=''.join(W[j] for j in p)
def fixed(n,R):return sorted(set(''.join(W[i] for i in s) for s in fibre(n,R)))
@lru_cache(None)
def language(i,R):
 if i<3:return frozenset([W[i][::-1]])
 words=set()
 for s in fibre(F[i],R):
  if s==(i,):continue
  words.update(''.join(parts) for parts in product(*(language(j,R) for j in reversed(s))))
 return frozenset(words)
def allowed(w,R):return set(''.join(x) for x in product(*[('CB' if R else 'C') if c=='C' else ('BA' if c=='B' else 'A') for c in w]))
for R in [False,True]:
 for i in range(3,8):need(language(i,R)==allowed(W[i][::-1],R),f'full language {i} {R}')
for j in range(6,21):
 base=W[j];chosen='B'+base
 need(chosen in fixed(F[j]+1,False),'B boundary D')
 need(chosen in fixed(F[j]+1,True),'B boundary R')
 need(not any(w[1:]==chosen for w in fixed(F[j]+2,False)),'D dies')
 rescued='A'+chosen
 need(rescued in fixed(F[j]+2,True),'R rescue')
 need(not any(w[1:]==rescued for w in fixed(F[j]+3,True)),'R fixed dies')
@lru_cache(None)
def parses(i,word,R):
 if len(word)!=weight(i):return False
 if i<3:return word==W[i]
 for support in fibre(F[i],R):
  if support==(i,):continue
  cursor=0;ok=True
  for j in reversed(support):
   l=weight(j)
   if not parses(j,word[cursor:cursor+l],R):ok=False;break
   cursor+=l
  if ok:return True
 return False
def root_parses(word,R):
 for support in fibre(len(word),R):
  cursor=0;ok=True
  for j in reversed(support):
   l=weight(j)
   if not parses(j,word[cursor:cursor+l],R):ok=False;break
   cursor+=l
  if ok:return True
 return False
v=W[20][::-1]
patched=v[:34]+'B'+v[35:]
for n in range(1,301):need(root_parses(patched[:n],True),f'old B all-horizon finite {n}')
need(not root_parses(patched[:36],False),'D full obstruction')
rows={str(n):{name:[{'support':list(s),'ascending_word':''.join(W[i] for i in s)} for s in fibre(n,R)] for name,R in [('D',False),('R',True)]} for n in range(34,41)}
print(json.dumps({'status':'PASS','recursive_language_Fn_through_index':7,'fixed_branch_j_range':[6,20],'patched35_recursive_R_prefixes_through':300,'old_supplied36':'not found; prior source supplies35 symbols only','rows34to40':rows},indent=2))

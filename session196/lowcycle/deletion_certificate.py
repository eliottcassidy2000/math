"""Separate exploratory all-order test of cycle-preserving base deletion."""
from apex_transfer import build_cycle_dfa,minimize_dfa
from collections import deque
import time,json
raw,t,o=build_cycle_dfa();p,t,o=minimize_dfa(t,o)
initial=(p[0],frozenset(),0);states=[initial];index={initial:0};parents=[None];bad=[];start=time.time()
for i,(s,deleted,length) in enumerate(states):
 c=o[s]
 if c in (1,2) and not any(o[d]==c for d in deleted):bad.append(i)
 for a in range(4):
  nd=frozenset([s]+[t[d][a] for d in deleted]);nxt=(t[s][a],nd,min(4,length+1))
  if nxt not in index:index[nxt]=len(states);states.append(nxt);parents.append((i,a))
 if i and i%100000==0:print(i,len(states),time.time()-start,flush=True)
def word(i):
 out=[]
 while parents[i] is not None:i,a=parents[i];out.append(a)
 return out[::-1]
print('states',len(states),'elapsed',time.time()-start)
for i in bad:
 s,d,m=states[i];print('bad',o[s],word(i),'lengthcap',m)
 assert (o[s]==1 and m<=2) or (o[s]==2 and m<=3)
from itertools import product
def cv(w):
 s=p[0]
 for a in w:s=t[s][a]
 return o[s]
small_bad=[]
for n in range(4):
 for w in product(range(4),repeat=n):
  c=cv(w)
  if c in (1,2) and all(cv(w[:j]+w[j+1:])!=c for j in range(n)):
   small_bad.append((c,w))
assert small_bad==[(1,(2,)),(1,(3,0)),(2,(1,2,0)),(2,(2,1,2)),(2,(2,3,0)),(2,(3,0,2)),(2,(3,1,0)),(2,(3,2,0)),(2,(3,2,1))]
print('verified: C=1 and base>=3 or C=2 and base>=4 admits one base deletion preserving C')
json.dump({'states':len(states),'bad_witnesses':[(o[states[i][0]],word(i)) for i in bad],'lemma_verified':True},open(__file__.replace('.py','.json'),'w'),indent=2)

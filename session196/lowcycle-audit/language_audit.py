"""Brzozowski-derivative equivalence, independent of Thompson NFA verifier."""
import json,sys,collections,itertools
from pathlib import Path
from independent_audit import automaton,graph_counts,slot_counts,BASE
ZERO=('zero',);EPS=('eps',)
def letter(a):return ('a',a)
def either(*xs):
    flat=set()
    for x in xs:
        if x==ZERO:continue
        if x[0]=='u':flat.update(x[1:])
        else:flat.add(x)
    if not flat:return ZERO
    if len(flat)==1:return next(iter(flat))
    return ('u',)+tuple(sorted(flat))
def chain(*xs):
    flat=[]
    for x in xs:
        if x==ZERO:return ZERO
        if x==EPS:continue
        if x[0]=='c':flat.extend(x[1:])
        else:flat.append(x)
    if not flat:return EPS
    if len(flat)==1:return flat[0]
    return ('c',)+tuple(flat)
def star(x):return EPS if x in (ZERO,EPS) else ('s',x)
def plus(x):return chain(x,star(x))
def word(s):return chain(*(letter(int(a)) for a in s))
def nullable(x):
    k=x[0]
    if k=='eps' or k=='s':return True
    if k=='zero' or k=='a':return False
    if k=='u':return any(nullable(y) for y in x[1:])
    return all(nullable(y) for y in x[1:])
from functools import lru_cache
@lru_cache(None)
def deriv(x,a):
    k=x[0]
    if k in ('eps','zero'):return ZERO
    if k=='a':return EPS if x[1]==a else ZERO
    if k=='u':return either(*(deriv(y,a) for y in x[1:]))
    if k=='s':return chain(deriv(x[1],a),x)
    terms=[]
    for i in range(1,len(x)):
        terms.append(chain(deriv(x[i],a),*x[i+1:]))
        if not nullable(x[i]):break
    return either(*terms)
P=either(EPS,chain(star(either(letter(1),letter(3))),letter(1)))
Q=either(EPS,chain(star(either(letter(2),letter(3))),letter(2)))
R=either(EPS,chain(star(letter(3)),letter(1),P),chain(star(letter(3)),letter(2),Q))
L1=chain(letter(1),star(either(letter(0),letter(1))),either(chain(either(word('20'),word('310'),chain(word('32'),star(letter(3)),letter(1))),P),chain(letter(2),plus(letter(3)),letter(0))))
L2=chain(letter(2),star(either(letter(0),letter(2))),either(chain(either(chain(letter(1),star(letter(3)),letter(2)),word('30')),Q),chain(word('33'),star(letter(3)),letter(0))))
L3=chain(letter(3),either(chain(plus(letter(0)),either(chain(letter(2),R),chain(plus(letter(3)),letter(2),Q))),chain(word('10'),P),word('20'),chain(letter(2),star(letter(3)),letter(1),P)))
expr=either(L1,L2,L3)
vs,T,O,initial,t,o=automaton()
queue=[(initial,expr)];seen=set(queue)
for s,r in queue:
    assert (o[s]==2)==nullable(r),(s,r)
    for a in range(4):
        new=(t[s][a],deriv(r,a))
        if new not in seen:seen.add(new);queue.append(new)
print('Independent derivative product:',len(seen),'states; all-order equivalence verified')
# Count machine words. Cayley-Hamilton on its 32-state adjacency matrix makes
# 32 consecutive zeros sufficient for the residual recurrence for all orders.
v=[0]*len(t);v[initial]=1;seq=[]
for n in range(100):
    seq.append(sum(v[i] for i in range(len(t)) if o[i]==2))
    new=[0]*len(t)
    for i,row in enumerate(t):
        for j in row:new[j]+=v[i]
    v=new
assert seq[:7]==[0,0,0,7,21,57,140]
residual=[seq[n]-5*seq[n-1]+8*seq[n-2]-4*seq[n-3] for n in range(7,100)]
assert residual==[0]*93
assert all(seq[n]==((11*n+6)*2**n)//32-4 for n in range(4,100))
print('Independent inventory: recurrence residual zero for 93 terms, >=32 sufficient by Cayley-Hamilton')
controls={'210':(15,3),'12301':(67,4),'':(1,0),'2':(3,1)}
for w,expected in controls.items():assert graph_counts(tuple(map(int,w)))==slot_counts(tuple(map(int,w)))==expected
print('Independent controls:',controls)
result={'derivative_product_states':len(seen),'inventory':seq[:16],'recurrence_residual_zero_terms':93,'machine_dimension':len(t),'controls':controls}
(BASE/'language_results.json').write_text(json.dumps(result,indent=2)+'\n')

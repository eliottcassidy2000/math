"""Test whether the remaining character (doubled) lifts to a modular-4 guard."""
from pathlib import Path
from runpy import run_path
M=run_path(str(Path(__file__).with_name('gewirtz-model.py')))
G=M['G']
I=tuple(range(56)); tags={I:0};queue=[I]
for p in queue:
    for t,g in zip([0]*6+[1,2],G):
        z=tuple(p[i] for i in g);q=tags[p]^t
        if z not in tags: tags[z]=q;queue.append(z)
        else:assert tags[z]==q
K=[g for g in queue if tags[g] in (0,3)]
assert len(K)==40320
def base(g):return sum(1<<(56*i+j+1) for i,j in enumerate(g))
# Reduce against unit pivots using the two-bit representation x=lo+2*hi.
# Both 2*the unit basis and all reductions with zero lo generate the high-only part.
basis={}; full={}; obstruction=None
def half(v):
    while v>=2 and v.bit_length() in full:v^=full[v.bit_length()]
    if v>=2:full[v.bit_length()]=v
    return v
def reduce(lo,hi,insert=False):
    while lo:
        k=lo.bit_length()
        if k not in basis:
            if insert:
                basis[k]=(lo,hi)
                half(lo)
            return lo,hi
        a,b=basis[k]
        hi^=b^ (a&~lo)
        lo^=a
    return lo,hi
for t,g in enumerate(K,1):
    lo,hi=reduce(base(g),int(tags[g]==3),True)
    if not lo and half(hi)==1:obstruction=t;break
print('mod4 tagged span:',len(basis),len(full),'obstruction at',obstruction,flush=True)
assert obstruction==6837, (len(basis), len(full), obstruction)
if obstruction is None:
    lo,hi=reduce(((1<<(56*56))-1)<<1,0)
    if lo:print('J not in base span');exit()
    r=half(hi)
    print('J reduction',r)

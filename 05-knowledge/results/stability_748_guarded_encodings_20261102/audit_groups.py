#!/usr/bin/env python3
"""Separate combinatorial audit of the additive greedy-group energy account."""
from itertools import combinations,product
from fractions import Fraction as F
import random,json

def need(b,msg):
    if not b:raise RuntimeError(msg)
def audit(n,weights,eta):
    E=list(combinations(range(n),2));q=dict(zip(E,weights))
    W=sum(q[e]*q[f] for e,f in combinations(E,2) if set(e)&set(f))
    remaining={e for e in E if q[e]>0};groups=[];seen_pairs=set();heavy=F(0);light=[];D=F(0)
    while remaining:
        e=max(remaining,key=lambda e:(q[e],e));remaining.remove(e)
        removed={f for f in remaining if set(e)&set(f)};remaining-=removed
        members=removed|{e};pairs={tuple(sorted((a,b))) for a,b in combinations(sorted(members),2) if set(a)&set(b)}
        need(not seen_pairs&pairs,'no wedge budget spent twice');seen_pairs|=pairs
        energy=sum(q[a]*q[b] for a,b in pairs);S=sum(q[f] for f in removed);w=q[e]
        need(energy>=w*S and energy>=S*S/4,'both valid individual group floors')
        if w>eta:heavy+=S
        else:light.append(S)
        D+=S;groups.append((w,S))
    M=sum(w for w,S in groups)
    need(M+D==sum(weights),'literal mass conservation')
    need(W>=D*D/(2*n),'inherited extraction bound')
    joint=eta*heavy+sum(S*S/4 for S in light)
    need(W>=joint,'new additive heavy/light bound')
    if light:need(joint>=eta*heavy+sum(light)**2/(4*len(light)),'light Cauchy step')
    return len(groups)

def main():
    complete=0
    for w in product([F(0),F(1,8),F(1,4)],repeat=6):
        for eta in [F(1,55),F(1,8)]:audit(4,w,eta);complete+=1
    rng=random.Random(748);samples=0
    for n in range(5,13):
        for _ in range(100):
            w=[F(rng.randrange(101),400) for e in combinations(range(n),2)]
            audit(n,w,F(1,55));samples+=1
    print(json.dumps({'status':'passed','complete_four_vertex_controls':complete,
        'seeded_rational_controls':samples,'seed':748,
        'scope':'Independent edge-pair accounting checks; the exact all-size group proof is inherited and the disjoint summation is explicit. Samples are not a stability proof.'},indent=2,sort_keys=True))
if __name__=='__main__':main()

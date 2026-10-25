#!/usr/bin/env python3
"""Finite exact controls; stdlib only. No RH or F1 descent test is claimed."""
import itertools, json, math
from pathlib import Path

def tournament_types(n):
    edges=list(itertools.combinations(range(n),2)); idx={e:i for i,e in enumerate(edges)}
    allmask=(1<<len(edges))-1
    perms=[]
    for p in itertools.permutations(range(n)):
        perms.append([(idx[tuple(sorted((p[a],p[b])))],p[a]>p[b]) for a,b in edges])
    def canonical(mask):
        return min(sum((((mask>>i)&1)^flip)<<j for i,(j,flip) in enumerate(rule)) for rule in perms)
    reps={canonical(m) for m in range(allmask+1)}
    f=sum(canonical(r^allmask)==r for r in reps)
    N=len(reps)
    assert (N-f)%2==0
    # Independent trace: explicit action on reps versus fixed/2-cycle decomposition.
    action={r:canonical(r^allmask) for r in reps}
    for m in range(1,9):
        actual=0
        for r in reps:
            x=r
            for _ in range(m): x=action[x]
            actual+=x==r
        assert actual==(f if m%2 else N)
    return {'n':n,'N':N,'self_converse':f,'two_cycles':(N-f)//2,'quotient':(N+f)//2}

rows=[]
for p in (3,5,7):
    C=(p+1)*p**(p-1)-p*(p+1)//2
    # F_{p^2} punctures: Frobenius has p fixed points and (p^2-p)/2 pairs.
    orbit_lengths=[1]*p+[2]*((p*p-p)//2)
    actual=[]
    for n in range(1,9):
        deleted=sum(d for d in orbit_lengths if n%d==0)
        assert deleted==p**math.gcd(n,2)
        count=p**n-deleted
        assert count>=0
        actual.append(count)
        if n%2==0: assert count==p**n-p*p
        else: assert count!=p**n-p*p
    # Fixed F_p-vector ambient dimension d: p^2 rational marked points.
    for d in (2,3,4):
        assert all(p**(d*m)-p*p>=0 for m in range(1,9))
    rows.append({'p':p,'C_p':C,'forbidden_first_count':C*(p-p*p),'F_p2_deletion_counts_n1_to8':actual,'N_at_1':C*(1-p*p)})
# Coordinatewise power maps fail to preserve the obvious complement after integer lift.
# p=3 marks {0,1,2}^2; (-1,0) is outside, but squaring lands at (1,0).
assert (-1,0) not in set(itertools.product(range(3),repeat=2))
assert (1,0) in set(itertools.product(range(3),repeat=2))
result={'prime_controls':rows,'tournament_types':[tournament_types(n) for n in range(1,6)],'lambda_power_map_hostile':{'deleted_grid':'{0,1,2}^2','source':[-1,0],'square_image':[1,0],'outcome':'naive psi_2 does not restrict to complement'},'scope':'Finite exact checks support separately supplied proofs; no novelty or RH claim.'}
Path(__file__).with_name('verification.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))

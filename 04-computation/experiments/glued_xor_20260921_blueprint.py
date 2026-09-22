"""Exact controls for the glued blueprint audit; no infinite-orbit inference."""
from fractions import Fraction as F
from hashlib import sha256
from itertools import permutations
from math import gcd
from pathlib import Path
import json
import re

ROOT = Path(__file__).resolve().parents[2]
CHECKS = 0


def check(condition, message):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(message)


def odd_step(n, b):
    value = 3*n+b
    if not value:
        raise ValueError('Zero has no finite 2-adic valuation')
    k = 0
    while value % 2 == 0:
        value //= 2
        k += 1
    return value, k


def mm(a, b):
    return tuple(tuple(sum(a[i][k]*b[k][j] for k in range(2))
                       for j in range(2)) for i in range(2))


def mp(a, n):
    r = ((F(1),F(0)),(F(0),F(1)))
    for _ in range(n):
        r = mm(r,a)
    return r


def projective(a):
    flat = sum(a, ())
    first = next(v for v in flat if v)
    return tuple(v/first for v in flat)


def run():
    tri = lambda t: t*(t+1)//2
    for a in range(-8,9):
        for b in range(-8,9):
            N = (a+1)*(b+1)
            check(tri(N-2) == tri(a)+tri(b)+tri(a*b)+a*(b*b-1)+b*(a*a-1),
                  'triangular identity')
            check(a+b+2 == (a+1)+(b+1), 'transported addition')
            check(a*b+a+b+1 == (a+1)*(b+1), 'transported multiplication')

    fibre = [1]
    for _ in range(10):
        fibre.append(4*fibre[-1]+1)
    check(all(odd_step(x,1)[0] == 1 for x in fibre), 'same target, not a path')
    check(all(gcd(a,b)==1 for a,b in zip(fibre,fibre[1:])), 'adjacent gcd')
    check(gcd(fibre[1],fibre[3])==5, 'nonadjacent shared prime hostile')
    n = F(5,7)
    check((3*n+1)/2 == F(11,7) and (3*F(11,7)+1)/8==n,
          'uncancelled rational cycle')

    R = lambda k: (10**(k+1)-7)//3
    check(R(1)==31 and R(2)==331, 'digit indices')
    check(gcd(R(1),R(16))==31, 'sharp horizon witness')
    check(R(8)==17*19607843, 'literal factorization')
    hcounts = {}
    for size in range(2,9):
        def edge(a,b):
            return (a==size-1 and b==0) or (a<b and not(a==0 and b==size-1))
        count = sum(all(edge(a,b) for a,b in zip(p,p[1:]))
                    for p in permutations(range(size)))
        hcounts[str(size)] = count
        check(count == (1 if size==2 else 1+2**(size-2)), 'endpoint-reversal path count')

    B = ((F(1,4),F(-13,4)),(F(1,4),F(3,4)))
    reflection = ((F(-1),F(-2)),(F(0),F(1)))
    identity = mp(B,0)
    minus_identity = ((F(-1),F(0)),(F(0),F(-1)))
    check(mp(B,3)==minus_identity and mp(B,6)==identity, 'vector lift order')
    check(mp(reflection,2)==identity, 'reflection order')
    check(mm(mm(reflection,B),reflection)==mp(B,5), 'reversal relation')
    group = {mm(mp(B,k),mp(reflection,e)) for k in range(6) for e in range(2)}
    check(len(group)==12 and len({projective(g) for g in group})==6,
          'projectivization loses central sign')

    prefixes = 0
    distinct_prefixes = 0
    longest = 0
    sample = []
    for b in (-1,1):
        for start in range(1,1000,2):
            n = start
            q,product,total = F(1),F(1),F(0)
            nodes=[]
            seen=set()
            for L in range(1,301):
                if n in seen:
                    break
                seen.add(n)
                nodes.append(n)
                nxt,k=odd_step(n,b)
                check(k>=1 and nxt>0 and nxt%2==1, 'positive odd exact step')
                total+=q
                product*=1+F(b,3*n)
                q*=F(2**k,3)
                check(q*nxt == start*product, 'multiplicative identity')
                check(q*nxt == start+F(b,3)*total, 'additive identity')
                if b==-1:
                    check(total<3*start and product<1, 'strict finite negative budget')
                else:
                    check(product>1, 'positive correction product')
                    images=nodes[1:]
                    check(all(gcd(x,6)==1 for x in images), 'mod-six image restriction')
                    sorted_images=sorted(images)
                    check(all(x>=3*j-2 for j,x in enumerate(sorted_images,1)),
                          'distinct-state spacing')
                    actual=sum((F(1,x) for x in images),F(0))
                    bound=sum((F(1,3*j-2) for j in range(1,len(images)+1)),F(0))
                    check(actual<=bound, 'reciprocal rearrangement')
                    # Written proof gives C=(4/3)e^(4/9)<3. This rational consequence
                    # avoids making a floating-point exponential a proof gate.
                    check(product**9<=3**9*L, 'rational envelope control')
                    distinct_prefixes+=1
                prefixes+=1
                longest=max(longest,L)
                if start in (1,5,17,27,97,999) and L in (1,2,7,20):
                    sample.append({'b':b,'start':start,'L':L,'endpoint':nxt,
                                   'q':str(q),'product':str(product),'partial_budget':str(total)})
                n=nxt

    budgets=[]
    for start,period in ((1,1),(5,2),(17,7)):
        n=start
        q,total=F(1),F(0)
        word=[]
        for _ in range(period):
            total+=q
            n,k=odd_step(n,-1)
            word.append(k)
            q*=F(2**k,3)
        check(n==start and 0<q<1, 'known minus cycle closes')
        check(total/(1-q)==3*start, 'geometric cycle exhausts budget')
        budgets.append({'start':start,'period':period,'word':word,
                        'block_ratio':str(q),'block_sum':str(total),
                        'infinite_budget':str(total/(1-q))})

    audit=ROOT/'04-computation/lean/CollatzBlueprintAudit/AxiomAudit.lean'
    declarations=re.findall(r'^#print axioms (\S+)',audit.read_text(encoding='utf-8'),re.M)
    check(len(declarations)==35, 'inherited theorem inventory changed; re-audit attribution')
    check(declarations[31]=='CollatzBlueprintAudit.twentyThree_tally_control',
          '32nd listed certificate is not a discrepancy theorem')
    source=Path(__file__)
    return {'status':'PASS','source_sha256':sha256(source.read_bytes()).hexdigest(),
            'universe':{'transport_grid':[-8,8],'inverse_fibre_seed':1,'fibre_last_index':10,
                        'endpoint_reversal_sizes':[2,8],'parameters':[-1,1],
                        'positive_odd_starts':[1,999],'prefix_cap':300,
                        'stop_at_first_repeat':True},
            'checks_passed':CHECKS,'orbit_prefixes_checked':prefixes,
            'positive_distinct_prefixes_checked':distinct_prefixes,'longest_prefix':longest,
            'endpoint_reversal_paths':hcounts,'fibre':fibre,'known_minus_cycle_budgets':budgets,
            'sample_prefixes':sample,'vector_group_order':12,'projective_group_order':6,
            'formal_attribution':{'public_theorems':35,'entry_32':declarations[31],
                                  'general_discrepancy_formalized':False},
            'infinite_results_require_written_proofs':True,'collatz_convergence_proved':False}


if __name__=='__main__':
    result=run()
    Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n',
                                                  encoding='utf-8',newline='\n')
    print('PASS:',result['checks_passed'],'exact controls;',result['orbit_prefixes_checked'],
          'actual prefixes; infinite assertions remain written/cited proofs')

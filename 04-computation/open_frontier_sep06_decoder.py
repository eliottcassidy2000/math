#!/usr/bin/env python3
"""Exact controls for the nonunit cross-divisor closure criterion.

The all-height argument is in the matching result note. No producer imports,
floating point, numerical phase grids, or claimed exhaustion of actual rows.
"""
from fractions import Fraction as F
from itertools import combinations
from math import gcd, lcm
from pathlib import Path
import hashlib
import json
import subprocess

Q = 91**6
CAP = {12: 1, 11: 2, 10: 4, 9: 9, 8: 30, 7: 90}
checks = 0


def require(test, message):
    global checks
    checks += 1
    if not test:
        raise RuntimeError(message)


def ceildiv(x, y):
    return -((-x) // y)


def inside(a, b, target, budget):
    r = (pow(a, -1, b) * target) % b
    s = (target - a*r) // b
    lower = max(ceildiv(-budget-r, b), ceildiv(s-budget, a))
    upper = min((budget-r)//b, (s+budget)//a)
    return lower <= upper


def crossing(A, B, Y, budget):
    A, B = sorted((A, B))
    D = gcd(A, B)
    a, b = A//D, B//D
    delta = gcd(D, Y)
    c, target = D//delta, Y//delta
    require(b <= budget, 'native internal height')
    return c <= budget and inside(a, b, target, budget)


def atlas_edge(A, B):
    D = gcd(A, B)
    p, q = sorted((A//D, B//D))
    if p+q > 356:
        return False
    n = p+q
    prime = 2
    while prime*prime <= n:
        power = 0
        while n % prime == 0:
            n //= prime
            power += 1
        if power and (prime % 3 != 2 or power > 2):
            return False
        prime += 1
    return n == 1 or n % 3 == 2


def components(row):
    unseen = set(range(len(row)))
    output = []
    while unseen:
        start = min(unseen)
        unseen.remove(start)
        block, queue = {start}, [start]
        while queue:
            i = queue.pop()
            for j in list(unseen):
                if atlas_edge(row[i], row[j]):
                    unseen.remove(j)
                    block.add(j)
                    queue.append(j)
        output.append(sorted(block))
    return output


def clearance(row, phase):
    return min(min((phase*x) % 1, 1-(phase*x) % 1) for x in row)


def actual_control(name, V, U, v, u, eta, zeta, scale=None):
    a, b = len(V), len(U)
    K = max(U)
    t = scale if scale is not None else 2*Q*K+1
    row = tuple(t*x for x in V) + tuple(U)
    require(a+b == 13 and a <= b, name+' split')
    require(gcd(*V) == gcd(*U) == gcd(*row) == 1, name+' primitive')
    require(len(set(row)) == 13 and min(row) > 0, name+' distinct')
    require(sum(row) <= Q*Q, name+' physical box')
    require(1 not in U, name+' genuinely nonunit larger component')
    require(components(row) == [list(range(a)), list(range(a,13))],
            name+' actual graph')
    for part in (V, U):
        for x, y in combinations(part, 2):
            require(max(x,y)//gcd(x,y) <= Q, name+' internal height')
    count = 0
    for same, other in ((tuple(t*x for x in V), U),
                        (U, tuple(t*x for x in V))):
        for A, B in combinations(same, 2):
            for Y in other:
                require(not crossing(A, B, Y, Q), name+' complete crossing')
                count += 1
    require(count == 11*a*b//2, name+' support universe')
    # Independent non-entry path: a mixed two-small/one-large relation has
    # distinguished coefficient >=t if every U coordinate is coprime to t.
    # A one-small/two-large row cannot balance because t min(V)>2QK.
    require(all(gcd(t,x) == 1 for x in U), name+' coefficient dominance')
    require(t*min(V) > 2*Q*K, name+' magnitude dominance')
    D = gcd(u, K)
    d = gcd(D, v)
    delta = gcd(D, t*v)  # g=1 in these controls
    aa, bb = u//D, K//D
    R = Q*(aa+bb)-(aa-1)*(bb-1)
    require(v in V and u in U and u < K, name+' selected labels')
    require(v % D == 0, name+' cross divisibility')
    require(CAP[b]*(D//d) <= Q, name+' uniform coefficient gate')
    require(7*(b+1)*lcm(D,v) <= a*Q, name+' uniform width gate')
    require(7*(b+1)*K*v <= a*delta*R, name+' exact native gate')
    require(clearance(V, eta) > F(1,14), name+' smaller safe phase')
    require(clearance(U, zeta) > F(1,14), name+' larger safe phase')
    desired = t*zeta-eta
    j = (desired+F(1,2)).numerator // (desired+F(1,2)).denominator
    phase = (eta+j)/t
    margin = clearance(row, phase)
    require(margin > F(1,14), name+' exact physical gluing witness')
    return dict(name=name, split=[a,b], t=t, g=1, V=V, U=U,
                D=D, v=v, K=K, radius=R, physical_sum=sum(row),
                supports=count, phase=str(phase), clearance=str(margin))


def main():
    table = []
    for a in range(1,7):
        b = 13-a
        tree = 355**(a-1)
        cutoff = min(Q//CAP[b], a*Q//(7*(b+1)*tree))
        table.append([a,b,tree,cutoff])
        require((tree <= F(a*Q,7*(b+1))) == (a <= 5), 'tree comparison')
    require([r[3] for r in table] ==
            [6240321451,38086468,175558,725,2,0], 'exact cutoff table')

    # Complete modest scalar-box universe, not a physical-row census.
    boxes = 0
    for B in range(2,9):
        for a in range(1,B):
            for b in range(a+1,B+1):
                if gcd(a,b) != 1:
                    continue
                S = {a*r+b*s for r in range(-B,B+1)
                     for s in range(-B,B+1)}
                R = B*(a+b)-(a-1)*(b-1)
                require(R >= B*b, 'endpoint lower radius')
                require(2*R > B*(a+b), 'half support strict')
                for x in range(-B*(a+b)-1,B*(a+b)+2):
                    require(inside(a,b,x,B) == (x in S), 'box exact membership')
                    if abs(x) <= R:
                        require(x in S, 'complete central interval')
                require(R+1 not in S, 'first central gap')
                boxes += 1

    arithmetic = 0
    for D in range(1,19):
        for v in range(1,19):
            for t in range(1,13):
                for g in range(1,13):
                    if gcd(t,g) != 1:
                        continue
                    delta = gcd(g*D,t*v)
                    d = gcd(D,v)
                    require(delta % d == 0, 'retained cross divisor')
                    require(g*D//delta <= g*(D//d), 'clearing coefficient bound')
                    require(F(delta, v*D) >= F(1,lcm(D,v)), 'width divisor gain')
                    if v % D == 0:
                        require(g % (g*D//delta) == 0, 'cross divisibility gate divides g')
                    arithmetic += 1

    # Boundary controls: coefficient gate; central-radius subtraction;
    # physical maximum and coprime clocks; collective-versus-pair gcd.
    require(not crossing(8,12,1,3), 'coefficient gate hostile c=4>3')
    require(not inside(2,3,14,3) and 14 <= 3*(2+3), 'unsubtracted radius hostile')
    bad_time = F(7,85)
    require(abs(F(1,12)-bad_time) < F(1,840), 'selected-max hostile proximity')
    require(clearance(tuple(range(1,11))+(85,),bad_time) == 0,
            'whole maximum remains necessary')
    require({F(2*j,4)%1 for j in range(4)} == {F(0),F(1,2)},
            'noncoprime grid collapse')
    require(gcd(6,2,3) == 1 and gcd(6,3) == 3, 'walk collective gcd hostile')
    require(lcm(10,10) <= F(5*200,7*9) < 10*10, 'strict cross-divisor improvement')

    controls = []
    for V,U in [
        ((2,3,6),(3,5,6,9,10,12,15,18,20,30)),
        ((2,3,4,6),(3,5,6,9,10,12,15,18,30)),
        ((2,3,4,6,9),(3,5,6,9,10,12,15,30)),
        ((2,3,4,6,9,12),(3,5,6,9,10,15,30)),
    ]:
        # 60Q+1 is coprime to each displayed larger coordinate.
        controls.append(actual_control('nonunit-%s+%s'%(len(V),len(U)),
                        V,U,3,3,F(1,5),F(1,7)))
    V = (1013861907,1036929995,1060507875,1084612635,1096868145)
    v = V[0]
    U = tuple(sorted(V+(3*v,4*v,9*v)))
    stress = actual_control('billion-divisor-5+8',V,U,v,v,F(1,2),
                            F(1,2)+F(1,24*v),22*Q+1)
    require(stress['D'] > 76388115, 'stress divisor beyond old 11+2 uniform constant')
    require(stress['physical_sum'] == 66123361394683029422040, 'stress physical sum')
    controls.append(stress)

    # Retain the complete inherited necessary profiles, not merely their
    # maximum clocks. Sparse worktrees may read the frozen data from HEAD.
    relative = '05-knowledge/results/lrc14_joint_shadow_empty_core_next_sep06.json'
    base = Path(__file__).resolve().parent.parent
    path = base/relative
    raw = path.read_bytes() if path.exists() else subprocess.check_output(
        ['git','show','HEAD:'+relative], cwd=base)
    require(hashlib.sha256(raw).hexdigest() ==
            '935f3f687b6d7c89cc099e536f536238fd753bcc4c1747906d213cef387ca93f',
            'inherited profile data pin')
    bank = {int(d): {(c,tuple(gs)) for c,gs in level['profiles']}
            for d,level in json.loads(raw)['levels'].items()}
    for entry in controls:
        row = tuple(entry['t']*x for x in entry['V']) + tuple(entry['U'])
        maxima, rejects = [], 0
        for tails in range(1,7):
            greatest = 1
            for positions in combinations(range(13),13-tails):
                occupied = set(positions)
                c = gcd(*(row[i] for i in occupied))
                gs = tuple(sorted(gcd(c,row[i]) for i in range(13)
                                  if i not in occupied))
                greatest = max(greatest,c)
                accepted = (c,gs) in bank[tails]
                if not accepted:
                    rejects += 1
                if entry['name'] != 'billion-divisor-5+8':
                    require(accepted, 'small control full hereditary profile')
            maxima.append(greatest)
        entry['subset_gcd_maxima'] = maxima
        entry['inherited_profile_rejections'] = rejects
        if entry['name'] == 'billion-divisor-5+8':
            require(maxima == [1,183,183,33123,33123,5929017], 'preclosed stress scope')
            require(rejects > 0, 'stress is not an unpaid residual example')
        else:
            require(maxima == [1,1,1,3,3,3], 'small control gcd maxima')

    manifest = dict(Q=Q, table=table, boxes=boxes,
                    divisor_controls=arithmetic, actual_controls=controls)
    digest = hashlib.sha256(json.dumps(manifest,sort_keys=True,
                            separators=(',',':')).encode()).hexdigest()
    print('PROVED analytic criterion; FINITE-EXACT controls only')
    print('UNIFORM_TABLE', json.dumps(table))
    print('SCALAR_BOXES', boxes, 'DIVISOR_CONTROLS', arithmetic)
    for entry in controls:
        print('ACTUAL', json.dumps(entry,sort_keys=True))
    print('EXPLICIT_GATES', checks)
    print('SEMANTIC_SHA256', digest)


if __name__ == '__main__':
    main()

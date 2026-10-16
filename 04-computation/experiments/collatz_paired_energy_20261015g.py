#!/usr/bin/env python3
"""Frozen-source audit of paired energy, joint debt, and global core allocation.

All acceptance uses integers/Fractions. Decimal displays are not evidence.
This reuses the authenticated native rows; it does not extend their universe.
"""
import hashlib
import json
import math
from fractions import Fraction as F
from pathlib import Path

from collatz_reset_debt_analysis_20261015f import initial_norms, ceilroot

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / '05-knowledge/results'


def need(ok, message):
    if not ok:
        raise AssertionError(message)


def frac(x):
    return [x.numerator, x.denominator]


def global_budget(old, new, b, budget, rho=F(7, 10)):
    """Exact ball-intersection test for minimum conformal global-core residual."""
    eta = 1-rho*rho
    U = F(sum(u*u for u in old), b)
    V = F(sum(v*v for v in new), b)
    excess = V-rho*rho*U
    beta = [max(F(0), v-rho*rho*u) for u, v in zip(old, new)]
    B2 = sum(z*z for z in beta)/b
    disc = B2-eta*excess
    need(disc >= 0, 'feasible conformal global core')
    if excess <= 0:
        return True, F(0), B2, excess
    # The unconstrained quadratic minimum lies at beta/eta, within [0,v].
    if budget*budget >= B2/(eta*eta):
        return True, None, B2, excess
    left = excess+eta*budget*budget
    return left*left <= 4*budget*budget*B2, None, B2, excess


def analyse(rec):
    k, b = rec['k'], rec['bin_width']
    X, rho = 1 << k, F(7, 10)
    history = [initial_norms(k)] + [s['row_norms'] for s in rec['stages']]
    energies = [sum(x*x for x in u) for u in history]
    need(energies[0] == rec['initial_envelope'][0], 'same initial source')
    floor = 49*b*X*X
    need(energies[-1] < floor, 'frozen stopping rule')
    pairs = []
    for j in range(len(history)-1):
        end = min(j+2, len(history)-1)
        old, new = energies[j], energies[end]
        # A final singleton is followed by the declared identity padding.
        below = new <= floor
        half = 2*new <= old
        need(below or half, 'paired half-energy or floor in declared census')
        need(3*new <= old+2*floor, 'stronger empirical affine 1/3 envelope')
        u, v = history[j], history[end]
        # Match loss M_i=(u_i-v_i)/2 gives exactly the energy loss.
        matched = sum((a-c)*(a+c) for a,c in zip(u,v))
        need(matched == old-new, 'weighted native matched-mass identity')
        pairs.append(dict(start=j, end=end, numerator_before=old,
                          numerator_after=new, common_denominator=b,
                          below_floor=below, halves=half,
                          above_floor_ratio=frac(F(new,old)) if not below else None))
    debt = [F(0)]*len(history[0])
    scalar_debt = F(0)
    global_failures = []
    for j, (u,v) in enumerate(zip(history, history[1:]), 1):
        residual = [F(max(0,10*c-7*a),3) for a,c in zip(u,v)]
        debt = [rho*d+(1-rho)*r for d,r in zip(debt,residual)]
        R2 = sum(r*r for r in residual)/b
        scalar_debt = rho*scalar_debt+(1-rho)*ceilroot(R2.numerator,R2.denominator)
        ok, _, B2, excess = global_budget(u,v,b,7*X)
        if not ok:
            global_failures.append(dict(stage=j, beta_square=frac(B2),
                                        excess=frac(excess)))
    for _ in range(k+8-(len(history)-1)):
        debt = [rho*d+(1-rho)*r for d,r in zip(debt,history[-1])]
        scalar_debt = rho*scalar_debt+(1-rho)*ceilroot(energies[-1],b)
    debt_square = sum(d*d for d in debt)/b
    bound_vector = [rho**(k+8)*u+d for u,d in zip(history[0],debt)]
    bound_square = sum(z*z for z in bound_vector)/b
    need(debt_square <= scalar_debt*scalar_debt, 'joint-bin debt improves Minkowski bound')
    need(debt_square <= (7*X*(1-rho**(k+8)))**2, 'finite joint debt certificate')
    need(bound_square < F(225,4)*X*X, 'finite final energy consequence')
    return dict(k=k,pairs=pairs,global_seven_failures=global_failures,
                joint_debt_square=frac(debt_square),
                joint_final_bound_square=frac(bound_square),
                joint_final_bound_over_X=math.sqrt(float(bound_square))/X)


def small_controls():
    """Independent rational KKT and brute-box lower-bound checks."""
    checked = 0
    rho = F(3,5)
    for u1 in range(1,5):
        for u2 in range(1,5):
            for v1 in range(u1+1):
                for v2 in range(u2+1):
                    u,v = [F(u1),F(u2)],[F(v1),F(v2)]
                    for budget in [F(0),F(1,2),F(1),F(2),F(4)]:
                        verdict,_,B2,e = global_budget(u,v,1,budget,rho)
                        # Every feasible quarter-grid split must pass the exact test.
                        for a in range(4*v1+1):
                            for c in range(4*v2+1):
                                r = [F(a,4),F(c,4)]
                                feasible = sum((y-z)**2 for y,z in zip(v,r)) <= rho*rho*sum((x-z)**2 for x,z in zip(u,r))
                                if feasible and sum(z*z for z in r) <= budget*budget:
                                    need(verdict,'grid witness cannot contradict optimum')
                        checked += 1
    # A global core does NOT inherit the 0.3 discounted residual recurrence.
    # u=(5,3), v=(5,0), r=(3,0), core norms (2,3), image norms (2,0).
    need(4 <= F(49,100)*13, 'actual global core contracts')
    need(49*34 < 41**2, '5 > .7 sqrt(34)+.3*3')
    # A literal signed identification realizes these two rows.
    rows = [[F(5),F(0)],[F(3,2),F(-3,2)]]
    need((3*3+1)//2 == 5 and (3*5+1)//2 == 8,
         'columns 3 and 8 map to 8 with literal clocks 2 and 0')
    need([sum(abs(x) for x in row) for row in rows] == [5,3], 'native abstract old rows')
    need([abs(sum(row)) for row in rows] == [5,0], 'native abstract new rows')
    return dict(rational_budget_controls=checked, angular_debt_hostile=True)


if __name__ == '__main__':
    paths = [DATA/'collatz_survivor_packets_20261015e.json'] + [
        DATA/f'collatz_compact_resets_20261015f_k{k}_r1.json' for k in range(19,24)]
    inherited = json.loads(paths[0].read_text())
    records = [r for r in inherited['first_reset'] if 8 <= r['k'] <= 18]
    records += [json.loads(p.read_text()) for p in paths[1:]]
    audits = [analyse(r) for r in records]
    above = [(F(*p['above_floor_ratio']),r['k'],p['start'],p['end'])
             for r in audits for p in r['pairs'] if p['above_floor_ratio']]
    worst = max(above)
    result = dict(status='PROVED identities; FINITE-EXACT k8..23; uniform target OPEN',
                  inputs={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
                  pair_transitions=sum(len(r['pairs']) for r in audits),
                  worst_above_floor=dict(ratio=frac(worst[0]),k=worst[1],start=worst[2],end=worst[3]),
                  controls=small_controls(),scales=audits)
    out=DATA/'collatz_paired_energy_20261015g.json'
    out.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps({k:result[k] for k in ['status','pair_transitions','worst_above_floor','controls']},indent=2))
    print('sha256',hashlib.sha256(out.read_bytes()).hexdigest())

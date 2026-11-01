#!/usr/bin/env python3
"""Exact LRC owner certificates, literal FRACTRAN and literal unary TM.
No third-party packages; all acceptance tests survive python -O.
"""
from fractions import Fraction as F
from itertools import combinations
from math import gcd
from pathlib import Path
import hashlib
import json

ROOT = Path(__file__).resolve().parent


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def valid_speeds(speeds):
    values = tuple(speeds)
    if not values or any(type(v) is not int or v <= 0 for v in values):
        raise ValueError('nonempty positive integer speed tuple required')
    return values


def valid_delta(delta):
    if not isinstance(delta, F) or not 0 < delta <= F(1, 2):
        raise ValueError('Fraction threshold in (0,1/2] required')
    return delta


def is_safe(speeds, t, delta):
    values = valid_speeds(speeds)
    valid_delta(delta)
    if not isinstance(t, F):
        raise ValueError('exact Fraction time required')
    a, b = t.numerator, t.denominator
    c, d = delta.numerator, delta.denominator
    return all(d * min(a*v % b, (-a*v) % b) >= c*b for v in values)


def intersect(left, right):
    """Sorted, disjoint, closed rational intervals; retain zero length."""
    out = []
    i = j = 0
    while i < len(left) and j < len(right):
        lo, hi = max(left[i][0], right[j][0]), min(left[i][1], right[j][1])
        if lo <= hi:
            out.append((lo, hi))
        if left[i][1] < right[j][1]:
            i += 1
        elif right[j][1] < left[i][1]:
            j += 1
        else:
            i += 1
            j += 1
    return out


def windows(speeds, delta=F(1,14)):
    values = valid_speeds(speeds)
    valid_delta(delta)
    result = [(F(0), F(1))]
    for v in sorted(set(values)):
        bands = [(F(k, v) + delta/v, F(k+1, v) - delta/v) for k in range(v)]
        result = intersect(result, bands)
        if not result:
            break
    return result


def owners(speeds, t, delta):
    """All labelled tight constraints and signs, including tied owners."""
    out = []
    for index, v in enumerate(speeds):
        phase = t*v
        k = phase.numerator // phase.denominator
        r = phase-k
        if r == delta:
            out.append([index, v, k, '+'])
        if r == 1-delta:
            out.append([index, v, k, '-'])
    return out


def certificate(speeds, delta):
    result = windows(speeds, delta)
    return [{'left': str(l), 'right': str(r),
             'left_owners': owners(speeds, l, delta),
             'right_owners': owners(speeds, r, delta)} for l,r in result]


def wall_cell_oracle(speeds, delta):
    """Independent cell decomposition: evaluate all walls and all midpoints."""
    points = {F(0), F(1)}
    for v in speeds:
        for k in range(v):
            points.add(F(k,v)+delta/v)
            points.add(F(k+1,v)-delta/v)
    points = sorted(points)
    atoms = [(x,x) for x in points if is_safe(speeds,x,delta)]
    for l,r in zip(points,points[1:]):
        if is_safe(speeds,(l+r)/2,delta):
            require(is_safe(speeds,l,delta) and is_safe(speeds,r,delta), 'closedness')
            atoms.append((l,r))
    merged = []
    for l,r in sorted(atoms):
        if merged and l <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1],r))
        else:
            merged.append((l,r))
    return merged


def unit_harmonic_windows(k):
    if type(k) is not int or k < 1:
        raise ValueError('positive integer k required')
    q = k+1
    eps = F(1,q*(q+1))
    intervals = []
    for a in range(1,q):
        if gcd(a,q) == 1:
            u = pow(a,-1,q)
            intervals.append((F(a,q)-eps/u, F(a,q)+eps/(q-u)))
    isolated = [F(b,q+1) for b in range(2,q) if gcd(b,q+1)==1]
    return sorted(intervals+[(t,t) for t in isolated])


def farey_windows(k, delta):
    """Independent THM-826 gap realization, with equality kept as singleton."""
    farey = sorted({F(a,b) for b in range(1,k+1) for a in range(b+1)})
    result = []
    for a,b in zip(farey,farey[1:]):
        lo = a+delta/a.denominator
        hi = b-delta/b.denominator
        if lo <= hi:
            result.append((lo,hi))
    return result


def bounded_search(n,height,delta=F(1,14),primitive_only=True):
    """Complete increasing-row box; no residue/cover/future-multiple pruning."""
    if type(n) is not int or type(height) is not int or not 1<=n<=height:
        raise ValueError('integers 1<=n<=height required')
    valid_delta(delta)
    for row in combinations(range(1,height+1),n):
        g = 0
        for v in row:
            g = gcd(g,v)
        if primitive_only and g!=1:
            continue
        result = windows(row,delta)
        yield row,result[0][0] if result else None


def first_good_period(speeds, cap, delta=F(1,14)):
    for b in range(2,cap+1):
        for a in range(1,b):
            if gcd(a,b)==1 and is_safe(speeds,F(a,b),delta):
                return F(a,b)
    return None


# Genuine ordered fractions; state primes 13 reject, 17 accept.
FRACTRAN = (F(1,6),F(1,35),F(13,33),F(13,77),F(17,11))


def fractran_run(n, program=FRACTRAN, budget=100000):
    if type(n) is not int or n < 1:
        raise ValueError('positive integer FRACTRAN state required')
    for step in range(budget+1):
        for i,f in enumerate(program):
            if n % f.denominator == 0:
                if step == budget:
                    raise RuntimeError('FRACTRAN step budget exhausted')
                n = n//f.denominator*f.numerator
                break
        else:
            return n, step
    raise RuntimeError('unreachable')


def comparison_state(A,B,C,D):
    values = (A,B,C,D)
    if any(type(v) is not int or v<0 for v in values):
        raise ValueError('nonnegative integer registers required')
    return 11*2**A*3**B*5**C*7**D


def fractran_compare(A,B,C,D):
    end,steps = fractran_run(comparison_state(A,B,C,D), budget=min(A,B)+min(C,D)+1)
    require((end%13==0) != (end%17==0),'exclusive terminal output')
    return end%17==0,steps,end


def point_registers(v,t,delta,rounded=False):
    r = v*t.numerator % t.denominator
    c,d = delta.numerator,delta.denominator
    b = t.denominator
    if rounded:
        K = (c*b+d-1)//d
        return r,K,b-r,K
    return d*r,c*b,d*(b-r),c*b


def first_primes(n):
    values = []
    x = 2
    while len(values) < n:
        if all(x % p for p in values if p*p <= x):
            values.append(x)
        x += 1
    return values


def compile_conjunction(inequalities):
    """Actual FRACTRAN family; each pair is (left >= right)."""
    pairs = tuple(inequalities)
    if any(len(pair) != 2 for pair in pairs) or any(type(x) is not int or x < 0 for pair in pairs for x in pair):
        raise ValueError('nonnegative integer comparison registers required')
    primes = first_primes(2*len(pairs)+3)
    registers = list(zip(primes[::2][:len(pairs)], primes[1::2][:len(pairs)]))
    control,reject,accept = primes[-3:]
    program = tuple([F(1,p*q) for p,q in registers] +
                    [F(reject,control*q) for p,q in registers] + [F(accept,control)])
    state = control
    for (a,b),(p,q) in zip(pairs,registers):
        state *= p**a*q**b
    return program,state,{'register_primes':registers,'control':control,'reject':reject,'accept':accept}


def whole_point_machine(speeds,t,delta,rounded=False):
    values = valid_speeds(speeds)
    valid_delta(delta)
    if not isinstance(t,F):
        raise ValueError('exact Fraction time required')
    inequalities = []
    for v in values:
        A,B,C,D = point_registers(v,t,delta,rounded=rounded)
        inequalities.extend(((A,B),(C,D)))
    program,state,layout = compile_conjunction(inequalities)
    bound = sum(min(a,b) for a,b in inequalities)+1
    terminal,steps = fractran_run(state,program,budget=bound)
    require(steps==bound,'whole row exact FRACTRAN steps')
    require((terminal%layout['accept']==0)!=(terminal%layout['reject']==0),'whole row exclusive halt')
    return terminal%layout['accept']==0,steps,program,state,layout,inequalities


# Single-tape TM, blank '_', head initially at first input symbol.
# (state,read) -> (state,write,direction); A and R halt.
TM = {
 ('H','1'):('H','1',1), ('H','X'):('H','X',1), ('H','#'):('B','#',1),
 ('B','X'):('B','X',1), ('B','1'):('C','X',-1), ('B','_'):('A','_',0),
 ('C','1'):('C','1',-1), ('C','X'):('C','X',-1), ('C','#'):('D','#',-1),
 ('D','X'):('D','X',-1), ('D','1'):('H','X',1), ('D','_'):('R','_',0),
}


def tm_compare(A,B):
    if any(type(v) is not int or v<0 for v in (A,B)):
        raise ValueError('nonnegative integer registers required')
    tape = dict(enumerate('1'*A+'#'+'1'*B))
    state,head,steps = 'H',0,0
    bound = 4*(A+B+2)*(min(A,B)+2)
    while state not in ('A','R'):
        require(steps<bound,'TM bound exceeded')
        read = tape.get(head,'_')
        rule = TM.get((state,read))
        if rule is None:
            raise RuntimeError('undefined TM transition outside input promise')
        state,write,move = rule
        tape[head] = write
        head += move
        steps += 1
    return state=='A',steps


def run_tests():
    cases = 0
    max_fr_steps = 0
    for A in range(9):
        for B in range(9):
            for C in range(9):
                for D in range(9):
                    yes,steps,end = fractran_compare(A,B,C,D)
                    require(yes == (A>=B and C>=D), ('FRACTRAN',A,B,C,D,end))
                    require(steps == min(A,B)+min(C,D)+1,'exact fraction count')
                    max_fr_steps = max(max_fr_steps,steps)
                    cases += 1
    tm_cases = 0
    max_tm_steps = 0
    for A in range(65):
        for B in range(65):
            yes,steps = tm_compare(A,B)
            require(yes == (A>=B),('TM',A,B))
            max_tm_steps = max(max_tm_steps,steps)
            tm_cases += 1
    rational_cases = 0
    for b in range(2,17):
        for a in range(b):
            t = F(a,b)
            for v in range(1,13):
                for delta in (F(1,14),F(1,5),F(1,2)):
                    regs = point_registers(v,t,delta)
                    require(fractran_compare(*regs)[0] == is_safe((v,),t,delta),'rational FRACTRAN')
                    rounded_regs = point_registers(v,t,delta,rounded=True)
                    require(fractran_compare(*rounded_regs)[0] == is_safe((v,),t,delta),'rounded rational FRACTRAN')
                    rational_cases += 1
    harmonic = []
    for k in range(1,81):
        delta = F(1,k+2)
        expected = unit_harmonic_windows(k)
        actual = windows(tuple(range(1,k+1)),delta)
        require(actual == expected,('harmonic',k))
        require(farey_windows(k,delta) == expected,('Farey',k))
        count_positive = sum(l<r for l,r in actual)
        count_isolated = sum(l==r for l,r in actual)
        phi_q = sum(gcd(a,k+1)==1 for a in range(1,k+1))
        phi_next = sum(gcd(a,k+2)==1 for a in range(1,k+2))
        require((count_positive,count_isolated)==(phi_q,phi_next-2),'component count')
        measure = sum((r-l for l,r in actual),F(0))
        hprim = sum((F(1,u) for u in range(1,k+1) if gcd(u,k+1)==1),F(0))
        require(measure == 2*hprim/((k+1)*(k+2)),'primitive harmonic measure')
        if k in (1,3,5,12,30,80):
            harmonic.append({'k':k,'positive':count_positive,'isolated':count_isolated,'measure':str(measure)})
    oracle_cases = 0
    for size in range(1,5):
        for row in combinations(range(1,10),size):
            for delta in (F(1,14),F(1,5),F(1,3),F(1,2)):
                require(windows(row,delta)==wall_cell_oracle(row,delta),('oracle',row,delta))
                oracle_cases += 1
    bounded_rows = 0
    for row,witness in bounded_search(13,16):
        require(witness is not None and is_safe(row,witness,F(1,14)),'bounded LRC14 smoke census')
        bounded_rows += 1
    require(bounded_rows==560,'bounded exact universe')
    V26 = tuple(26*i for i in range(1,13))+(339,)
    incoherent = (81,91,131,151,157,196,258,274,313,328,330,339,348)
    zero_free = (43,55,61,70,73,79,83,99,103,104,109,113,156)
    hostiles = []
    for name,row,first in [('V26',V26,F(2,27)),('incoherent',incoherent,F(3,26)),('zero_free',zero_free,F(2,27))]:
        require(first_good_period(row,25) is None,('q25 hostile',name))
        require(first_good_period(row,27)==first,('first good period',name))
        cells = windows(row)
        require(cells,'hostile has windows')
        l,r = cells[0]
        require(is_safe(row,l,F(1,14)),'owner witness')
        require(l.denominator<=14*max(row),'owner denominator cap')
        for v in row:
            require(fractran_compare(*point_registers(v,first,F(1,14)))[0],'hostile machine witness')
        hostiles.append({'name':name,'first_good':str(first),'first_window':[str(l),str(r)],
                         'left_owners':owners(row,l,F(1,14)),'components':len(cells)})
    whole_row_runs = []
    whole_row_artifacts = []
    for name,t,rounded in [('V26-safe',F(2,27),False),('V26-zero-blocked',F(1,26),False),
                           ('V26-safe-rounded',F(2,27),True),('V26-zero-blocked-rounded',F(1,26),True)]:
        yes,steps,program,state,layout,ineq = whole_point_machine(V26,t,F(1,14),rounded=rounded)
        require(yes == is_safe(V26,t,F(1,14)),'whole-point FRACTRAN')
        whole_row_runs.append({'name':name,'accepted':yes,'steps':steps,'fractions':len(program),
                               'input_bit_length':state.bit_length()})
        whole_row_artifacts.append({'name':name,'speeds':V26,'t':str(t),'delta':'1/14',
             'fractions':[str(f) for f in program],'layout':layout,'inequalities':ineq,
             'register_encoding':'rounded' if rounded else 'cross-multiplied',
             'input_bit_length':state.bit_length(),
             'input_sha256_big_endian':hashlib.sha256(state.to_bytes((state.bit_length()+7)//8,'big')).hexdigest()})
    # Freeze machine descriptions, not preprocessor-dependent opaque integers.
    generated={'five_fraction_primitive':[str(f) for f in FRACTRAN],
               'tm_transitions':[{'state':a,'read':b,'next':c,'write':d,'move':e}
                                 for (a,b),(c,d,e) in sorted(TM.items())],
               'whole_row_programs':whole_row_artifacts}
    (ROOT/'machines.json').write_text(json.dumps(generated,indent=2,sort_keys=True)+'\n')
    deep = windows(tuple(range(1,13)))
    require([l for l,r in deep if l==r]==[F(3,14),F(5,14),F(9,14),F(11,14)],'deep-well isolated sidecar')
    tight = windows(tuple(range(1,14)))
    require(tight == [(F(a,14),F(a,14)) for a in (1,3,5,9,11,13)],'tight AP singletons')
    require(sum(r-l for l,r in tight)==0 and len(tight)==6,'zero measure but safe')
    require(windows((1,2),F(1,2))==[],'negative control')
    require(windows((1,),F(1,2))==[(F(1,2),F(1,2))],'half singleton')
    require(owners((1,3),F(1,4),F(1,4)) == [[0,1,0,'+'],[1,3,0,'-']],'tie owners')
    # Explicit demonstration that first-eligible order is essential.
    altered=(F(17,11),)+FRACTRAN[:-1]
    end,steps=fractran_run(comparison_state(0,1,0,0),altered)
    require(end%17==0,'reorder hostile must falsely accept')
    bad_inputs = 0
    for fn,args in [(windows,((0,1),F(1,14))),(windows,((1,),F(0))),
                    (windows,((1,),0.1)),(comparison_state,(-1,0,0,0)),
                    (tm_compare,(True,1)),(unit_harmonic_windows,(0,))]:
        try:
            fn(*args)
        except ValueError:
            bad_inputs += 1
        else:
            raise RuntimeError('malformed input accepted')
    out = {'status':'FINITE-EXACT / PASS; LRC(14) OPEN',
           'fractran_register_cases':cases,'fractran_max_steps':max_fr_steps,
           'tm_register_cases':tm_cases,'tm_max_steps':max_tm_steps,
           'whole_row_fractran_runs':whole_row_runs,
           'rational_fractran_cases_per_encoding':rational_cases,'rational_register_encodings':2,'independent_window_cases':oracle_cases,
           'bounded_box':{'speeds':13,'height':16,'primitive':True,'distinct':True,'rows':bounded_rows,'failures':0},
           'harmonic_k_range':[1,80],'harmonic_selected':harmonic,'canonical_hostiles':hostiles,
           'deep_well_isolated':['3/14','5/14','9/14','11/14'],
           'tight_AP13_measure':'0','tight_AP13_singletons':6,
           'reordering_hostile':'acceptance fraction first falsely accepts (0,1,0,0)',
           'rejected_malformed_inputs':bad_inputs}
    print(json.dumps(out,indent=2,sort_keys=True))


if __name__ == '__main__':
    run_tests()

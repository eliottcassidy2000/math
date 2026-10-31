#!/usr/bin/env python3
"""Exact one-step PRIMEGAME classification on all F=U+S shapes.

Stdlib only. All checks remain active under python -O. The finite representatives
classify arbitrary outside primes by the spectator-equivariance proof in README.
The balance classification itself is inherited, not a new result of this file.
"""
from collections import Counter
from hashlib import sha256
from itertools import combinations
import json
from math import gcd, isqrt, prod
from pathlib import Path

ROOT = Path(__file__).resolve().parent
ACTIVE = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29)
DATA = ACTIVE[:4]
CONTROL = ACTIVE[4:]
SPECTATORS = (31, 37, 41)
BASIS = ACTIVE + SPECTATORS
RULES = ((17,91), (78,85), (19,51), (23,38), (29,33),
         (77,29), (95,23), (77,19), (1,17), (11,13),
         (13,11), (15,2), (1,7), (55,1))
# Independently hand-written multiset rewrite rules on prime labels.
HAND = (((7,13),(17,)), ((5,17),(2,3,13)), ((3,17),(19,)),
        ((2,19),(23,)), ((3,11),(29,)), ((29,),(7,11)),
        ((23,),(5,19)), ((19,),(7,11)), ((17,),()),
        ((13,),(11,)), ((11,),(13,)), ((2,),(3,5)),
        ((7,),()), ((),(5,11)))
BALANCED_PROFILES = {(1,), (3,), (2,1,1)}
CHECKS = Counter()


def require(condition, detail, category='other'):
    CHECKS[category] += 1
    if not condition:
        raise RuntimeError(detail)


def factor(n):
    if n < 1:
        raise ValueError('positive integer required')
    result = {}
    p = 2
    while p*p <= n:
        while n % p == 0:
            result[p] = result.get(p, 0) + 1
            n //= p
        p += 1
    if n > 1:
        result[n] = result.get(n, 0) + 1
    return result


def integer(factors):
    return prod(p**e for p,e in factors.items())


def profile(factors):
    return tuple(sorted(factors.values(), reverse=True))


def fus(factors):
    """None, not a spurious tau-2 extension, at the excluded unit."""
    if not factors:
        return None
    r = len(factors)
    prime = profile(factors) == (1,)
    squarefree = all(e == 1 for e in factors.values())
    return (prod(e+1 for e in factors.values())-2,
            r-int(prime), 2**r-1-int(squarefree))


def direct_fus(n):
    if n == 1:
        return None
    divs = set()
    for d in range(1, isqrt(n)+1):
        if n % d == 0:
            divs.update((d, n//d))
    divs.difference_update((1,n))
    # Direct divisor carrier, not evaluation of the inherited tau/omega formula.
    return (len(divs),
            sum(all(d % j for j in range(2,isqrt(d)+1)) for d in divs),
            sum(all(d % (j*j) for j in range(2,isqrt(d)+1)) for d in divs))


def balanced(factors):
    counts = fus(factors)
    return counts is not None and counts[0] == counts[1]+counts[2]


def canonical_key(factors):
    return (tuple(factors.get(p,0) for p in ACTIVE),
            tuple(sorted((e for p,e in factors.items() if p not in ACTIVE),
                         reverse=True)))


def valid(factors):
    # Equivalent to: zero or one occupied control prime, each exponent <=1.
    return sum(factors.get(p,0) for p in CONTROL) <= 1


def step_integer(n, rules=RULES):
    for index,(a,b) in enumerate(rules,1):
        if n % b == 0:
            return n//b*a, index
    raise RuntimeError('last rule must be unconditional')


def step_hand(factors):
    for index,(take,give) in enumerate(HAND,1):
        if all(factors.get(p,0) >= take.count(p) for p in take):
            output = Counter(factors)
            output.subtract(take)
            output.update(give)
            return {p:e for p,e in output.items() if e}, index
    raise RuntimeError('unconditional hand rule missing')


def representatives():
    values = list(ACTIVE+(31,)) + [p**3 for p in ACTIVE+(31,)]
    for p in ACTIVE:
        others = tuple(q for q in ACTIVE if q != p)
        values.extend(p*p*q*r for q,r in combinations(others,2))
        values.extend(p*p*q*31 for q in others)
        values.append(p*p*31*37)
    values.extend(31*31*q*r for q,r in combinations(ACTIVE,2))
    values.extend(31*31*q*37 for q in ACTIVE)
    values.append(31*31*37*41)
    return sorted(values)


def preserving_criterion(factors, rule):
    """Universal theorem, applied only to balanced first-enabled inputs."""
    exceptional = {1:7*13*17**2, 3:3*17*19**2, 4:2*19*23**2,
                   5:3*11*29**2, 6:29**3, 7:23**3, 8:19**3, 12:8}
    if rule in exceptional:
        return integer(factors) == exceptional[rule]
    if rule == 10:
        return ((factors.get(13,0) == 1 and factors.get(11,0) == 0) or
                (factors.get(13,0) == 2 and factors.get(11,0) == 1))
    if rule == 11:
        return factors.get(11,0) == 1
    return False


def cycle_member(factors):
    """The entire recurrent set inside the balanced-input carrier."""
    for control in (11,13):
        if factors.get(control,0) == 1:
            rest={p:e for p,e in factors.items() if p != control}
            if not rest:
                return True
            if (profile(rest)==(2,1) and
                    all(p in (2,5) or p not in ACTIVE for p in rest)):
                return True
    return False


def rename(factors, mapping):
    return {mapping.get(p,p):e for p,e in factors.items()}


def summarize(rows):
    result = {'total':len(rows),
              'input_profiles':dict(sorted(Counter(r['input_profile'] for r in rows).items())),
              'outcomes':dict(sorted(Counter(r['outcome'] for r in rows).items())),
              'restricted_orbits':dict(sorted(Counter(r['restricted_orbit'] for r in rows).items())),
              'by_rule':[]}
    for rule in range(1,15):
        selected = [r for r in rows if r['rule']==rule]
        result['by_rule'].append({'rule':rule, 'fraction':'%s/%s'%RULES[rule-1],
            'total':len(selected),
            'outcomes':dict(sorted(Counter(r['outcome'] for r in selected).items())),
            'output_profiles':dict(sorted(Counter(r['output_profile'] for r in selected).items()))})
    return result


def encode_profile(factors):
    return ','.join(map(str,profile(factors))) or 'unit'


def encode_factors(factors):
    return '*'.join(str(p) if e==1 else f'{p}^{e}' for p,e in sorted(factors.items())) or '1'


def main():
    values = representatives()
    require(len(values)==len(set(values))==538, 'representative count', 'census')
    require(len({canonical_key(factor(n)) for n in values})==538,
            'duplicate spectator orbit', 'census')
    for (a,b),(take,give) in zip(RULES,HAND):
        require(gcd(a,b)==1, 'fraction must be reduced', 'rules')
        require(factor(a)==dict(Counter(give)) and factor(b)==dict(Counter(take)),
                'hand and rational rule disagree', 'rules')

    # Independent construction: all labelled choices on 13 primes, then quotient.
    alternate = set(BASIS) | {p**3 for p in BASIS}
    for p,q,r in combinations(BASIS,3):
        alternate.update((p*p*q*r,p*q*q*r,p*q*r*r))
    require(len(alternate)==884, 'independent labelled universe', 'census')
    require({canonical_key(factor(n)) for n in alternate} ==
            {canonical_key(factor(n)) for n in values}, 'independent orbit census', 'census')

    # This renames the WHOLE program, not the fixed named-prime machine.
    permutation = dict(zip(ACTIVE, reversed(ACTIVE)))
    permutation.update({31:43,43:31,37:47,47:37,41:53,53:41})
    renamed_rules = tuple((integer(rename(factor(a),permutation)),
                           integer(rename(factor(b),permutation))) for a,b in RULES)
    rows=[]
    for n in values:
        incoming=factor(n)
        require(balanced(incoming) and profile(incoming) in BALANCED_PROFILES,
                f'unbalanced representative {n}', 'inputs')
        require(direct_fus(n)==fus(incoming), f'direct input counts {n}', 'direct_divisors')
        m,rule=step_integer(n)
        outgoing=factor(m)
        hand,hand_rule=step_hand(incoming)
        require((integer(hand),hand_rule)==(m,rule), f'exponent step {n}', 'step_equivalence')
        require(direct_fus(m)==fus(outgoing), f'direct output counts {n}', 'direct_divisors')
        require(balanced(outgoing)==(profile(outgoing) in BALANCED_PROFILES),
                f'output balance profile {n}', 'outputs')
        require(balanced(outgoing)==preserving_criterion(incoming,rule),
                f'preserving theorem {n}', 'preserving_theorem')
        active_part={p:e for p,e in incoming.items() if p in ACTIVE}
        spectator_part={p:e for p,e in incoming.items() if p not in ACTIVE}
        active_output,active_rule=step_integer(integer(active_part))
        require((m,rule)==(active_output*integer(spectator_part),active_rule),
                f'spectator semiconjugacy {n}', 'spectators')
        require({p:e for p,e in outgoing.items() if p not in ACTIVE}==spectator_part,
                f'spectator changed {n}', 'spectators')
        renamed_n=integer(rename(incoming,permutation))
        require(step_integer(renamed_n,renamed_rules)==
                (integer(rename(outgoing,permutation)),rule),
                f'whole-program prime-renaming conjugacy {n}', 'renaming')
        if valid(incoming):
            require(valid(outgoing), f'valid control carrier not closed at {n}', 'valid_carrier')
            require(balanced(outgoing)==(rule in (10,11) or n==8),
                    f'valid preserving criterion {n}', 'valid_carrier')
        if cycle_member(incoming):
            require(m != n and step_integer(m)[0]==n and cycle_member(outgoing),
                    f'balanced 2-cycle {n}', 'restricted_orbits')
            restricted='cycle_2'
        else:
            cursor=n
            for k in range(1,4):
                cursor=step_integer(cursor)[0]
                if not balanced(factor(cursor)):
                    restricted=f'exit_at_{k}'
                    break
            else:
                raise RuntimeError(f'claimed balanced exit bound failed {n}')
            require(not valid(incoming) or k<=2,
                    f'valid balanced exit bound {n}', 'restricted_orbits')
        rows.append({'n':n, 'input_factors':encode_factors(incoming),
            'input_profile':encode_profile(incoming), 'rule':rule,
            'fraction':'%s/%s'%RULES[rule-1], 'm':m,
            'output_factors':encode_factors(outgoing),
            'output_profile':encode_profile(outgoing),
            'outcome':'unit' if m==1 else 'balanced' if balanced(outgoing) else 'unbalanced',
            'active_only':not spectator_part, 'valid_control':valid(incoming),
            'output_valid_control':valid(outgoing),
            'restricted_orbit':restricted,
            'input_fus':fus(incoming), 'output_fus':fus(outgoing)})

    universes={
        'active': [r for r in rows if r['active_only']],
        'all_spectator_types': rows,
        'active_valid': [r for r in rows if r['active_only'] and r['valid_control']],
        'all_valid_spectator_types': [r for r in rows if r['valid_control']]}
    expected={'active':(380,24,354,2), 'all_spectator_types':(538,39,497,2),
              'active_valid':(98,15,81,2), 'all_valid_spectator_types':(181,29,150,2)}
    summary={name:summarize(selected) for name,selected in universes.items()}
    for name,(total,yes,no,unit) in expected.items():
        s=summary[name]
        require((s['total'],s['outcomes'].get('balanced',0),s['outcomes'].get('unbalanced',0),
                 s['outcomes'].get('unit',0))==(total,yes,no,unit),
                f'final census {name}', 'census')

    # Hostiles / boundaries do not infer whole-orbit behavior from one step.
    require(fus(factor(55))==(2,2,2) and not balanced(factor(55)), 'literal 55', 'hostiles')
    require(step_integer(7)==(1,13) and step_integer(17)==(1,9) and fus({}) is None,
            'unit boundary', 'hostiles')
    require(fus(factor(2))==fus(factor(11)) and
            not balanced(factor(step_integer(2)[0])) and balanced(factor(step_integer(11)[0])),
            'FUS profile quotient is not closed', 'hostiles')
    require(all(balanced(factor(2**e))==(e in (1,3)) for e in range(1,31)),
            'pure powers of 2 balance iff exponent 1 or3', 'hostiles')
    require(step_integer(11)[1]!=step_integer(13)[1],
            'fixed-machine prime-label permutation is not a labelled symmetry', 'hostiles')
    require(not balanced(factor(8*31)) and balanced(factor(8)),
            'projection to active factor need not preserve balance', 'hostiles')

    columns=list(rows[0])
    def cell(v):
        if v is None: return 'undefined'
        if isinstance(v,tuple): return ','.join(map(str,v))
        if isinstance(v,bool): return str(int(v))
        return str(v)
    table='\t'.join(columns)+'\n'+''.join('\t'.join(cell(r[c]) for c in columns)+'\n' for r in rows)
    (ROOT/'transitions.tsv').write_text(table)
    result={'status':'FINITE-EXACT representatives; universal coverage proved in README',
            'basis':ACTIVE, 'spectator_representatives':SPECTATORS,
            'independent_labelled_universe':884,
            'universes':summary,
            'balanced_preserving_representatives':[{'n':r['n'],'rule':r['rule'],'m':r['m']}
                                                  for r in rows if r['outcome']=='balanced'],
            'unit_boundary':[{'n':r['n'],'rule':r['rule']} for r in rows if r['outcome']=='unit'],
            'checks':dict(sorted(CHECKS.items())), 'checks_total':sum(CHECKS.values()),
            'table_sha256':sha256(table.encode()).hexdigest()}
    text=json.dumps(result,indent=2,sort_keys=True)+'\n'
    (ROOT/'classification.json').write_text(text)
    print(text,end='')


if __name__=='__main__':
    main()

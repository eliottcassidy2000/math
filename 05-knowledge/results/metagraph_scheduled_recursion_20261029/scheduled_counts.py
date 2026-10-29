"""Exact A000568 prefix recurrence with forced-zero stages contracted.

The inherited gcd-response quotient and odd-factorial weights are unchanged.
An agenda evaluates a state only at lengths that fit its remaining capacity.
No tournament isomorphism classes or actual complement transitions are merged.
"""
from collections import defaultdict
from math import factorial, gcd


def count_prefix(bound, colors=1):
    if type(bound) is not int or bound < 0:
        raise ValueError('bound must be a nonnegative integer')
    if type(colors) is not int or colors < 1:
        raise ValueError('colors must be a positive integer')
    if bound < 3:
        # Continue through the same unit-tail formula below.
        first = 1
    else:
        first = bound if bound % 2 else bound-1
    width = max(1,bound.bit_length())
    digit = (1 << width)-1
    masks = [(1 << (width*c))-1 for c in range(bound//2+2)]
    facts = [factorial(i) for i in range(bound+1)]
    twos = [i-i.bit_count() for i in range(bound+1)]
    odds = [f >> v for f,v in zip(facts,twos)]
    powers = [colors**k for k in range(bound+1)]
    agenda = [{} for _ in range(first+1)]
    agenda[first][0] = {0:1}
    stages = []
    transitions = 0
    skipped = 0
    active_visits = 0
    peak_active = 1
    for p in range(first,2,-2):
        current = agenda[p]
        agenda[p] = None
        visits = sum(map(len,current.values()))
        active_visits += visits
        peak_active = max(peak_active,visits)
        shift = width*((p-1)//2)
        row = sum(gcd(p,j) << (width*((j-1)//2)) for j in range(1,p,2))
        stage_transitions = 0
        for d,bucket in current.items():
            limit = (bound-d)//p
            if limit < 1:
                raise ArithmeticError('agenda scheduled an impossible positive insertion')
            for k in range(limit+1):
                nd = d+k*p
                remaining = bound-nd
                q = min(p-2, remaining if remaining % 2 else remaining-1)
                q = max(1,q)
                mask = masks[(min(q,remaining)+1)//2]
                target = agenda[q].setdefault(nd,defaultdict(int))
                stage_transitions += len(bucket)
                skipped += len(bucket)*((p-2-q)//2)
                if k == 0:
                    for response,value in bucket.items():
                        target[response&mask] += value
                    continue
                factor,remainder = divmod(odds[nd],odds[d]*p**k*odds[k])
                if remainder:
                    raise ArithmeticError('nonintegral odd-factorial transfer')
                factor *= powers[k]
                base = k*(p-1)//2+p*k*(k-1)//2-twos[k]
                offset = k*row
                for response,value in bucket.items():
                    e = base+k*((response>>shift)&digit)
                    target[(response+offset)&mask] += (value*factor) << e
        transitions += stage_transitions
        stages.append(dict(length=p,active_states=visits,transitions=stage_transitions))
    totals = [0]*(bound+1)
    terminal_states = sum(map(len,agenda[1].values()))
    for d,bucket in agenda[1].items():
        for k in range(bound-d+1):
            nd = d+k
            factor,remainder = divmod(odds[nd],odds[d]*odds[k])
            if remainder:
                raise ArithmeticError('nonintegral terminal transfer')
            evaluation = sum(value << (k*response) for response,value in bucket.items())
            totals[nd] += (evaluation*factor*powers[k]) << (k*(k-1)//2-twos[k])
            transitions += len(bucket)
    values=[]
    for total,denominator in zip(totals,odds):
        value,remainder=divmod(total,denominator)
        if remainder:
            raise ArithmeticError('nonintegral Burnside coefficient')
        values.append(value)
    return values,dict(stages=stages,active_visits=active_visits,
                       terminal_states=terminal_states,peak_active=peak_active,
                       transitions=transitions,skipped_branch_stage_visits=skipped)


if __name__ == '__main__':
    import argparse,json
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('bound',type=int)
    parser.add_argument('--colors',type=int,default=1)
    args=parser.parse_args()
    values,stats=count_prefix(args.bound,args.colors)
    print(json.dumps(dict(bound=args.bound,values=values,stats=stats)))

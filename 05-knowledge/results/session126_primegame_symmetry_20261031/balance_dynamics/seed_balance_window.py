#!/usr/bin/env python3
"""Finite seed-2 balance window and exact macro token-potential controls.

Macro words are inherited from the sibling fractran lane's reconstructed
divisor-loop proof. The token minimum and complete balance-window argument
are proved in README section 6a. Stdlib; explicit checks survive -O.
"""
from collections import Counter
import json
from pathlib import Path
from verify_balance_dynamics import (balanced, direct_fus, encode_factors,
    factor, fus, integer, step_hand, step_integer)

ROOT=Path(__file__).resolve().parent
CHECKS=Counter()


def check(condition, message, group):
    CHECKS[group]+=1
    if not condition:
        raise RuntimeError(message)


def macro_word(n,d):
    q,r=divmod(n,d)
    prefix=((1,2)*d+(10,)+(5,6)*d+(11,))*q+(1,2)*r+(1,)
    if r:
        suffix=(3,)+(4,7)*n+(8,)+(5,6)*(r-1)+(11,)
    else:
        suffix=(9,)+(12,)*n+(13,)*(d-1)+(14,)+(5,6)*n+(11,)
    return prefix+suffix


def main():
    macro_trials=0
    macro_steps=0
    for n in range(2,51):
        for d in range(1,n):
            state={5:n,7:d,13:1}
            minimum=sum(state.values())
            for instruction in macro_word(n,d):
                state,chosen=step_hand(state)
                check(chosen==instruction, f'first-enabled macro rule {n},{d}', 'macro_steps')
                minimum=min(minimum,sum(state.values()))
                macro_steps+=1
            expected={5:n,7:d-1,13:1} if n%d else {5:n+1,7:n,13:1}
            expected={p:e for p,e in expected.items() if e}
            check(state==expected, f'macro endpoint {n},{d}', 'macro_endpoints')
            check(minimum==n+d-1-int(bool(n%d)), f'exact token minimum {n},{d}', 'macro_minima')
            check(minimum>=n, f'token lower bound {n},{d}', 'macro_minima')
            macro_trials+=1
    check((macro_trials,macro_steps)==(1225,267171),'macro census','census')

    boundary=5**5*7**4*13
    n=2
    trace=[]
    hits=[]
    for t in range(146):
        state=factor(n)
        counts=fus(state)
        check(counts==direct_fus(n), f'direct seed divisor counts {t}', 'seed_trace')
        row={'step':t,'n':n,'factors':encode_factors(state),'omega_total':sum(state.values()),
             'fus':counts,'balanced':balanced(state)}
        trace.append(row)
        if row['balanced']:
            hits.append(row)
        if n==boundary:
            check(t==145,'boundary update index','seed_trace')
            break
        m,j=step_integer(n)
        out,k=step_hand(state)
        check((m,j)==(integer(out),k), f'independent seed transition {t}', 'seed_trace')
        n=m
    else:
        raise RuntimeError('did not reach C(5,4) within verified finite window')
    expected=[(0,2),(2,825),(4,1925),(5,2275),(13,156),(14,132),(16,308),(17,364),
              (43,380),(45,950),(53,1450),(56,850),(59,580),(62,340),(69,8),(70,60)]
    check([(r['step'],r['n']) for r in hits]==expected,'complete balanced hit list','seed_trace')
    result={'status':'FINITE-EXACT window; universal tail exclusion proved in README',
            'macro_trials':macro_trials,'macro_steps':macro_steps,
            'macro_token_minimum':'n+d-1 if d divides n; n+d-2 otherwise',
            'boundary':{'step':145,'n':boundary,'canonical_state':'C(5,4)'},
            'balanced_visit_count':len(hits),'balanced_visits':hits,
            'last_balanced_visit':hits[-1],
            'checks':dict(sorted(CHECKS.items())),'checks_total':sum(CHECKS.values())}
    out=json.dumps(result,indent=2,sort_keys=True)+'\n'
    (ROOT/'seed_balance_window.json').write_text(out)
    (ROOT/'seed_balance_trace.json').write_text(json.dumps(trace,indent=2,sort_keys=True)+'\n')
    print(out,end='')


if __name__=='__main__':
    main()

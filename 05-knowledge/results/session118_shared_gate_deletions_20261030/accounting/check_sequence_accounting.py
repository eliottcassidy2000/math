"""Integer inventory checks; scalar schedules are NOT asserted circuit-realizable."""
from itertools import product
import json

def require(ok, context):
    if not ok:
        raise RuntimeError(context)

def aggregate_grid():
    cases = equality_cases = 0
    # Bounded nonnegative aggregate counts. Prefix realizability is not assumed.
    for G,B,A,P,Z,H,D in product(range(3), repeat=7):
        q = G-2*B-A-P-Z-H
        if q < 0:
            continue
        for U,I,J,sT,E in product(range(3), repeat=5):
            for K in range(3 if G else 1):
                h = B+A+q+Z+H+U+I+J
                n = G+B+A+h
                I0 = n-J
                require(I == I0-2*G-A+P-U, 'usage inventory')
                require(0 <= I <= h-q, 'terminal free-coordinate budget')
                require(G >= n-2*h and P >= n-4*h, 'frequency bounds')
                for r in [3,4]:
                    s0 = sT+2*G+r*B+3*A+2*P+H+D+K+E
                    if r == 3:
                        extra = sT+5*(I+J+U)+3*q+2*A+3*Z+4*H+D+K+E
                    else:
                        extra = sT+4*(I+J+U)+2*q+A+2*Z+3*H+D+K+E
                    require(s0 == 4*n-(12-r)*h+extra, 'gate-slack identity')
                    zero = all(x == 0 for x in [sT,I,J,U,q,A,Z,H,D,K,E])
                    require((extra == 0) == zero, 'equality iff')
                    if zero:
                        require((B,G,P) == (h,n-2*h,n-4*h), 'tight count classification')
                        equality_cases += 1
                    cases += 1
    return {'exact_identity_cases':cases,'aggregate_equality_cases':equality_cases,
            'count_range':'each row count 0..2; slack/end-state parameters 0..2'}

def prefix_schedules():
    schedules = transitions = 0
    for r in [3,4]:
        for h in range(1,9):
            for n in range(4*h,4*h+21):
                m,I,q,s = n,n,0,4*n-(12-r)*h
                word = ['GP','P2']*(n-4*h)+['GP']*(2*h)+['AP2']*h
                counts = {'GP':0,'P2':0,'AP2':0}
                for row in word:
                    require(m-q >= h+1, 'numerical loop guard')
                    if row == 'GP':
                        require(I >= 2 and m-2*q >= 2, 'GP scalar precondition')
                        q += 1; I -= 2; s -= 2
                    elif row == 'P2':
                        require(q >= 1, 'P2 pair inventory')
                        m -= 1; q -= 1; I += 1; s -= 2
                    else:
                        require(q >= 2, 'AP2 pair inventory')
                        m -= 3; q -= 2; s -= r
                    counts[row] += 1
                    require(q >= 0 and m >= 2*q and 0 <= I <= m-2*q and s >= 0,
                            'prefix nonnegative state')
                    transitions += 1
                require((m,I,q,s) == (h,0,0,0), 'terminal state')
                require(counts == {'GP':n-2*h,'P2':n-4*h,'AP2':h}, 'tight frequencies')
                schedules += 1
    return {'schedules':schedules,'transitions':transitions,'h_range':[1,8],
            'n_range':'4h..4h+20','scope':'prefix-feasible scalar inventories, NOT DAG realizations'}

if __name__ == '__main__':
    print(json.dumps({'status':'PASS','aggregate':aggregate_grid(),
                      'prefix':prefix_schedules()},indent=2,sort_keys=True))

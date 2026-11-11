#!/usr/bin/env python3
"""Small exact controls; no large graph or orbit census."""
import argparse
import hashlib
import itertools
import json
import math
from pathlib import Path

INF = math.inf


def require(p, why):
    if not p:
        raise RuntimeError(why)


def numeric_closure(seed, cap):
    present = [False] * (cap + 1)
    for n in range(1, cap + 1):
        present[n] = n in seed or any(present[a] and present[n-a]
                                      for a in range(1, (n+1)//2))
    return present


def monoid(seed, cap):
    h = [False] * (cap + 1)
    h[0] = True
    for n in range(1, cap + 1):
        h[n] = any(n >= s and h[n-s] for s in seed)
    return h


def direct_state(seed, m, cap=512):
    hvals = monoid(seed, cap)
    sums = {a+b for a,b in itertools.combinations(seed, 2)}
    ivals = [any(n >= q and hvals[n-q] for q in sums) for n in range(cap+1)]
    h = [min((n for n in range(r, cap+1, m) if hvals[n]), default=INF) for r in range(m)]
    i = [min((n for n in range(r, cap+1, m) if ivals[n]), default=INF) for r in range(m)]
    return h, i, hvals, ivals


def add_seed(h, i, m, c):
    ell = m // math.gcd(m,c)
    p = h[:]
    p[0] = m
    hp = [min(j*c+h[(r-j*c)%m] for j in range(ell)) for r in range(m)]
    ip = [min(min(j*c+i[(r-j*c)%m], (j+1)*c+p[(r-(j+1)*c)%m])
              for j in range(ell)) for r in range(m)]
    return hp, ip


def hole_summary(seed, m, i):
    if any(x == INF for x in i):
        return None
    e = sum(s < i[s%m] for s in seed)
    twice = 2*sum(i)//m - (m-1) - 2 - 2*e
    require(twice % 2 == 0, 'hole formula integral')
    candidates = []
    for r in range(m):
        x = i[r]-m
        while x > 0 and x in seed:
            x -= m
        if x > 0:
            candidates.append(x)
    last = max(candidates, default=0)
    return {'holes': twice//2, 'last_hole': last, 'conductor': last+1}


def kernel(a,n):
    ans = 0
    power = 1
    while power*a <= n:
        r = n-power*a
        if 0 <= r < 2*power:
            ans += (-1)**r.bit_count()
        power *= 2
    return ans


def row_budget(n):
    if n == 1:
        return 1
    k = n.bit_length()-1
    leading = n >> (k-1)
    return 2*k-1 if leading == 2 else 2*k+1


def synchronous_step(current):
    return current | {a+b for a,b in itertools.combinations(sorted(current),2)}


def chase_step(M,W,D):
    require(W > 2,'compiler only before stop')
    g=min(D-{0})
    require(1 <= g <= W-2,'second-maximum gap bound')
    return (2*M-g, W-g,
            {u+v-g for u,v in itertools.combinations(sorted(D),2) if u+v <= W})


def run():
    seed_sets = [set(s) for size in range(1,9)
                 for s in itertools.combinations(range(1,9),size)]
    closure_checks = 0
    updates = 0
    for seed in seed_sets:
        literal = numeric_closure(seed,256)
        m = min(seed)
        h,i,hvals,ivals = direct_state(seed,m,256)
        for n in range(1,257):
            require(literal[n] == (n in seed or ivals[n]), 'mixed-seed closure')
            require(literal[n] == (n in seed or n >= i[n%m]), 'Apéry reconstruction')
            closure_checks += 1
        if len(seed) >= 2 and math.gcd(*seed) == 1:
            s = hole_summary(seed,m,i)
            holes = [n for n in range(1,257) if not literal[n]]
            require(s['holes'] == len(holes) and s['last_hole'] == max(holes,default=0), 'hole extraction')
        # Ascending and descending insertion: the latter retains a nonminimal modulus.
        orders = [sorted(seed),sorted(seed,reverse=True)]
        if len(seed) <= 4:
            orders = itertools.permutations(sorted(seed))
        for order in orders:
            m = order[0]
            h = [0]+[INF]*(m-1)
            i = [INF]*m
            current = {m}
            for c in order[1:]:
                h,i = add_seed(h,i,m,c)
                current.add(c)
                hh,ii,_,_ = direct_state(current,m,256)
                require(h == hh and i == ii, 'one-seed finite update')
                updates += 1
    selected = []
    for ss in [{2,3},{2,3,4},{4,7,11},{6,10,15},{1,7,11},{6,9,20}]:
        m = min(ss)
        h,i,_,_ = direct_state(ss,m)
        summary = hole_summary(ss,m,i)
        selected.append({'seeds':sorted(ss),'modulus':m,'h':h,'i':i,**summary})
    interval_checks = 0
    for c in range(1,13):
        for M in range(2*c,2*c+15):
            for z in range(M+1,2*M):
                a=max(c,z-M)
                b=z-a
                require(c <= a < b <= M and a+b == z, 'interval parent witness')
                interval_checks += 1
    sync = []
    for seed in [{1,2},{2,3},{2,5},{3,4},{4,7,11}]:
        current=set(seed)
        depth={s:0 for s in seed}
        _,i,_,_=direct_state(seed,min(seed))
        c=hole_summary(seed,min(seed),i)['conductor']
        t0=max(0,2*c//min(seed)-1)
        records=[]
        for t in range(7):
            if t:
                current |= {a+b for a,b in itertools.combinations(sorted(current),2)}
                for n in current:
                    depth.setdefault(n,t)
            if t >= t0:
                bound=(2*c-1)*(1<<(t-t0))+1
                require(all(n in current for n in range(c,bound+1)), 'expanding synchronous interior')
            records.append({'generation':t,'size':len(current),'maximum':max(current)})
        sync.append({'seeds':sorted(seed),'conductor':c,'t0_bound':t0,'generations':records})

    frontier_steps = 0
    top_block_controls = 0
    for mask in range(1,1<<12):
        current={n+1 for n in range(12) if mask>>n & 1}
        M=max(current)
        blocks=[a for a in current if {a+1,a+2} <= current]
        if not blocks:
            continue
        nxt=synchronous_step(current)
        for a in blocks:
            W=M-a
            D={M-x for x in current if M-x <= W}
            if W == 2:
                require(max(nxt) == 2*M-1 and {2*M-3,2*M-2,2*M-1} <= nxt,
                        'persistent top-three block')
                top_block_controls += 1
            else:
                MM,WW,DD=chase_step(M,W,D)
                require(max(nxt) == MM and DD == {MM-x for x in nxt if 0 <= MM-x <= WW},
                        'exact shrinking-window compiler')
                frontier_steps += 1
    chases=[]
    for initial in [{1,2,3,7},{2,3,4,9},{1,4,5,6,10},{1,2,3,12,13}]:
        current=set(initial)
        M=max(current)
        a=min(a for a in current if {a+1,a+2} <= current)
        W=M-a
        D={M-x for x in current if M-x <= W}
        bound=W-2
        history=[{'maximum':M,'width':W,'deficits':sorted(D)}]
        while W > 2:
            M,W,D=chase_step(M,W,D)
            current=synchronous_step(current)
            require(M == max(current) and D == {M-x for x in current if 0 <= M-x <= W},
                    'complete chase against full states')
            history.append({'maximum':M,'width':W,'deficits':sorted(D)})
            require(len(history)-1 <= bound,'finite chase bound')
        for j in range(2):
            nxt=synchronous_step(current)
            require(max(nxt) == 2*max(current)-1 and
                    {max(nxt)-2,max(nxt)-1,max(nxt)} <= nxt,'post-chase closed regime')
            current=nxt
        chases.append({'initial':sorted(initial),'tracked_block_start':a,'history':history})

    product_checks=0
    for seed in [{2,3},{2,4},{2,3,6},{3,5,15}]:
        cap=1000
        generated=[False]*(cap+1)
        ordinary=[False]*(cap+1)
        ordinary[1]=True
        pair_products={a*b for a,b in itertools.combinations(seed,2)}
        for n in range(2,cap+1):
            ordinary[n]=any(n%s == 0 and ordinary[n//s] for s in seed)
            generated[n]=n in seed or any(n%a == 0 and generated[a] and generated[n//a]
                                         for a in range(2,math.isqrt(n)+1) if a*a<n)
            ideal=any(n%p == 0 and ordinary[n//p] for p in pair_products)
            require(generated[n] == (n in seed or ideal),'joint multiplicative ideal theorem')
            product_checks += 1

    N=128
    A=[[int(a <= n and not(n%2 == 0 and a == n//2)) for a in range(1,N+1)]
       for n in range(1,N+1)]
    inverse=[[0]*N for _ in range(N)]
    for n in range(N):
        inverse[n][n]=1
        for a in range(n):
            inverse[n][a]=-sum(A[n][b]*inverse[b][a] for b in range(a,n))
    for n in range(1,N+1):
        for a in range(1,N+1):
            expected=kernel(a,n) if a <= n else 0
            require(inverse[n-1][a-1] == expected and abs(expected) <= 1,'inverse coefficient')
        require(sum(abs(x) for x in inverse[n-1]) == row_budget(n),'exact row budget')
    f=[0]+[(n*n+3*n)%11-5 for n in range(1,N+1)]
    g=[0]+[sum(f[1:n+1])-(f[n//2] if n%2 == 0 else 0) for n in range(1,N+1)]
    recovered=[0]
    for n in range(1,N+1):
        recovered.append(g[n]-g[n-1]+(-1)**n*recovered[n//2])
    require(recovered == f,'signed reconstruction')
    return {'status':'FINITE-EXACT controls for proved all-size formulas',
            'universes':{'seed_sets':'all255 nonempty subsets of1..8; numeric targets1..256',
                         'insertion_orders':'every order for≤4 seeds; increasing/decreasing otherwise',
                         'inverse_order':N,'synchronous_generation_cap':6,
                         'frontier_states':'all nonempty subsets of1..12 containing a three-block',
                         'joint_product_controls':'four named seed sets; targets2..1000'},
            'closure_membership_checks':closure_checks,'finite_seed_updates':updates,
            'interval_parent_checks':interval_checks,'selected_apery_states':selected,
            'synchronous_controls':sync,'inverse_first_rows':inverse[:12],
            'frontier_one_step_checks':frontier_steps,'persistent_top_block_checks':top_block_controls,
            'complete_frontier_chases':chases,'joint_product_checks':product_checks,
            'inverse_entry_values':sorted(set(x for row in inverse for x in row)),
            'inverse_row_budgets':{str(n):row_budget(n) for n in [1,2,3,4,6,8,12,16,24,32,48,64,96,128]},
            'source_one_cancellation':{str(n):kernel(1,n) for n in range(1,17)},
            'operator_pair_count_claim':False,'shortest_circuit_claim':False}


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--write',action='store_true')
    args=parser.parse_args()
    data=json.dumps(run(),indent=2,sort_keys=True)+'\n'
    target=Path(__file__).with_name('receipt.json')
    if args.write:
        target.write_text(data)
    else:
        require(target.read_text()==data,'receipt drift')
    print(json.dumps({'passed':True,'receipt_sha256':hashlib.sha256(data.encode()).hexdigest()},sort_keys=True))

#!/usr/bin/env python3
"""Native clipping with a fixed source-address tree and adaptive block invoices."""
from collections import defaultdict
from pathlib import Path
import json
import check as c


def run():
    N = 4096
    V = {int(l.split()[0]) for l in Path(__file__).with_name('frontier.txt').read_text().splitlines()
         if l and not l.startswith('#')}
    rows = []
    for s in range(1, 4*N, 2):
        if s % 8 == 5:
            continue
        y, t, d, port = s, 0, 0, 0
        seen = set()
        while y not in V:
            c.check(y not in seen and d < 10000, 'unresolved native route')
            seen.add(y)
            port = y
            y, dt = c.edge(y)
            t += dt
            d += 1
        a, b = c.charge(s, N)
        rows.append((s, y, t, d, port, a, b))
    K = max((s-1)//2 for s, *_ in rows).bit_length()
    def profile(rr):
        E = []
        for k in range(K+1):
            bins = defaultdict(lambda: [0,0])
            for s,u,t,d,p,a,b in rr:
                z = bins[u, ((s-1)//2) % (1 << k)]
                sign = (-1)**t
                z[0] += sign*a
                z[1] += sign*b
            E.append(sum(map(c.Q, bins.values())))
        return E
    F = [r for r in rows if r[3] <= 67 and r[2] <= 64]
    retained = [r for r in F if r[3] <= 36]
    ids = {r[0] for r in retained}
    removed = [r for r in F if r[0] not in ids]
    ef, ek = profile(F), profile(retained)
    delta = [(ef[k]-ef[k+1])-(ek[k]-ek[k+1]) for k in range(K)]
    pair_delta = [0]*K
    pairs = 0
    removed_ids = {r[0] for r in removed}
    for r in removed:
        for z in F:
            if r[0] == z[0] or (z[0] in removed_ids and z[0] < r[0]) or r[1] != z[1]:
                continue
            a = tuple((-1)**r[2] * x for x in r[5:])
            b = tuple((-1)**z[2] * x for x in z[5:])
            beta = c.Q(c.add(a,b))-c.Q(a)-c.Q(b)
            xor = ((r[0]-1)//2) ^ ((z[0]-1)//2)
            k = (xor & -xor).bit_length()-1
            pair_delta[k] += beta
            pairs += 1
    c.check(pair_delta == delta == [-34,-19,-1,-11,13,-6,-8,2,2,0,-3,0,0], 'native boundary delta')
    bf,bk = c.budgets(ef), c.budgets(ek)
    c.check(bf['binary'] == 439 and bk['binary'] == 495 and pairs == 951, 'incoming exact hostile')
    # For every permitted partition, carry its own boundary invoice.
    counts = 0
    least_upper = None
    least_optimal_invoice = None
    def visit(i, cf, ck, bill, word):
        nonlocal counts, least_upper, least_optimal_invoice
        if i == K:
            c.check(ck <= cf+bill, 'partitionwise clipping bound')
            counts += 1
            candidate = (cf+bill, list(word))
            if least_upper is None or candidate[0] < least_upper[0]:
                least_upper = candidate
            if cf == bf['adaptive_1to3']:
                if least_optimal_invoice is None or bill < least_optimal_invoice[0]:
                    least_optimal_invoice = (bill, list(word))
            return
        for j in range(i+1, min(i+3,K)+1):
            df, dk = ef[i]-ef[j], ek[i]-ek[j]
            visit(j, cf+max(df,0), ck+max(dk,0), bill+max(dk-df,0), word+[j-i])
    visit(0,0,0,0,[])
    c.check(bk['adaptive_1to3'] <= least_upper[0] <= bf['adaptive_1to3']+least_optimal_invoice[0], 'optimized invoice')
    return dict(status='PASS', N=N, original_source_norm=sum(c.Q(r[5:]) for r in rows),
                full_energy=ef, retained_energy=ek, full_budgets=bf, retained_budgets=bk,
                removed_labels=[r[0] for r in removed], removed_norm=sum(c.Q(r[5:]) for r in removed),
                boundary_pairs=pairs, delta=delta, binary_invoice=sum(max(-d,0) for d in delta),
                block_partitions=counts, minimum_upper_and_word=least_upper,
                minimum_invoice_on_full_optimum_and_word=least_optimal_invoice,
                scope='one native joint-cutoff hostile, fixed original labels/tree/denominator; not all-scale control')


if __name__ == '__main__':
    print(json.dumps(run(), indent=2, sort_keys=True))

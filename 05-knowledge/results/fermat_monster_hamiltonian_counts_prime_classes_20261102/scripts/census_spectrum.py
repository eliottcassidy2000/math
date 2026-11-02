#!/usr/bin/env python3
"""census_spectrum.py (klein-2026-11-02) -- Hamiltonian path counts H(T) over all tournaments of orders 3..9 (nauty's
gentourng + the C counter scripts/hampath.c), summarising per order: number of tournaments, maximum, number of distinct
values, the odd values below the maximum that are NOT attained, and for the Fermat primes 3, 5, 17, 257, 65537 the
least order at which they occur (65537: order 11 by the hill-climb witnesses in data/climb_65537_11.out and the bound
max H(10) = 15,745 (OEIS A038375, CITED; also max H(n) <= n * max H(n-1) = 33,570 from the census value 3,357 at n = 9)).
usage: python3 census_spectrum.py [gentourng=/opt/homebrew/bin/gentourng]  (expects scripts/hampath compiled next to it)"""
import subprocess, sys, os
HERE = os.path.dirname(os.path.abspath(__file__)); DATA = os.path.join(HERE, '..', 'data')
GEN = sys.argv[1] if len(sys.argv) > 1 else '/opt/homebrew/bin/gentourng'
HP = os.path.join(HERE, 'hampath')
if not os.path.exists(HP): subprocess.run(['cc', '-O2', '-o', HP, os.path.join(HERE, 'hampath.c')], check=True)
first = {}; never = True
for n in range(3, 10):
    gen = subprocess.run([GEN, '-q', str(n)], capture_output=True, text=True).stdout
    out = subprocess.run([HP, 'census', str(n)], input=gen, capture_output=True, text=True).stdout
    open(os.path.join(DATA, 'census_%d.out' % n), 'w').write(out)
    lines = out.split('\n'); head = lines[0]
    vals = {}
    for l in lines[2:]:
        if l.strip() and ':' in l:
            h, cnt = l.split(':')[0].strip(), l.split()[1]; vals[int(h)] = int(cnt)
    mx = max(vals); missing = [h for h in range(1, mx + 1, 2) if h not in vals]
    print('%s; %d distinct values; odd values below the maximum not attained: %s%s' % (head, len(vals), missing[:12], ' ...' if len(missing) > 12 else ''))
    for F in (3, 5, 17, 257, 65537):
        if F in vals and F not in first: first[F] = (n, vals[F])
    never = never and 7 not in vals and 21 not in vals
print('least order attaining each Fermat prime (order, number of tournaments of that order):', first, '; 65537: order 11 (witnesses in data/climb_65537_11.out; max H(10) = 15745 < 65537 < 95095 = max H(11), OEIS A038375)')
print('7 and 21 never attained at orders 3..9 (consistent with THM-1370):', never)

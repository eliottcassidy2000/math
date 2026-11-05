#!/usr/bin/env python3
"""tcgen_driver.py (klein-2026-11-04): run the provably complete extension generator tcgen order by order.
cls_{n}.txt holds "hc digraph6" (canonical, labelg) for every tournament of order n with hc in {1,2}.
Reports per order: one-cycle count vs Garey F(2n-6), two-cycle count vs (n-2)F(2n-9) - L(2n-10), and the
minimum out-/in-degree census of the two-cycle classes (HYP-9459 (1): every two-cycle tournament of order >= 6 has an
out-degree-1 vertex, equivalently by converse an in-degree-1 vertex)."""
import subprocess, sys, time
LABELG = '/opt/homebrew/bin/labelg'
def fib(k):
    a, b = 0, 1
    for _ in range(k): a, b = b, a + b
    return a
def luc(k):
    a, b = 2, 1
    for _ in range(k): a, b = b, a + b
    return a
def rows_from_d6(s):
    n = ord(s[1]) - 63; bits = []
    for ch in s[2:]:
        v = ord(ch) - 63
        for k in range(5, -1, -1): bits.append((v >> k) & 1)
    rows = [0] * n
    for i in range(n):
        for j in range(n):
            if bits[i * n + j]: rows[i] |= 1 << j
    return rows
start, stop = int(sys.argv[1]), int(sys.argv[2])
for n in range(start, stop + 1):
    t0 = time.time(); m = n - 1
    with open('parents_%d.txt' % m, 'w') as fh:
        for l in open('cls_%d.txt' % m): fh.write(l.split()[1] + '\n')
    nparents = sum(1 for _ in open('parents_%d.txt' % m))
    with open('parents_%d.txt' % m) as fin, open('raw_%d.txt' % n, 'w') as fout:
        subprocess.run(['./tcgen', str(m)], stdin=fin, stdout=fout, stderr=open('tcgen_%d.err' % n, 'w'))
    subprocess.run("cut -d' ' -f2 raw_%d.txt | %s -q > can_%d.txt" % (n, LABELG, n), shell=True, check=True)
    cls = {}; nrec = 0
    with open('raw_%d.txt' % n) as fr, open('can_%d.txt' % n) as fc:
        for lr, lc in zip(fr, fc):
            h = int(lr.split()[0]); c = lc.strip(); nrec += 1
            if c in cls: assert cls[c] == h, "hc conflict"
            else: cls[c] = h
    assert nrec == sum(1 for _ in open('can_%d.txt' % n))
    parents = [None] * nparents; recs = [None] * nrec
    import os; os.remove('raw_%d.txt' % n); os.remove('can_%d.txt' % n); os.remove('parents_%d.txt' % m)
    one = sum(1 for h in cls.values() if h == 1); two = len(cls) - one
    # degree census of the two-cycle (and one-cycle) classes
    no_out1 = no_in1 = 0; no_out1_one = no_in1_one = 0
    for c, h in cls.items():
        rows = rows_from_d6(c); outd = [bin(r).count('1') for r in rows]; ind = [sum((rows[i] >> j) & 1 for i in range(n)) for j in range(n)]
        if 1 not in outd: (no_out1, no_out1_one)[h == 1] if False else None
        if h == 2:
            if 1 not in outd: no_out1 += 1
            if 1 not in ind: no_in1 += 1
        else:
            if 1 not in outd: no_out1_one += 1
            if 1 not in ind: no_in1_one += 1
    with open('cls_%d.txt' % n, 'w') as fh:
        for c, h in sorted(cls.items(), key=lambda x: (x[1], x[0])): fh.write('%d %s\n' % (h, c))
    print('order %2d: %d parents, %d labelled extensions; one-cycle %d (Garey F(2n-6) = %d); two-cycle %d ((n-2)F(2n-9)-L(2n-10) = %d); '
          'two-cycle classes without an out-degree-1 vertex: %d, without an in-degree-1 vertex: %d (one-cycle: %d, %d); %.0f s'
          % (n, len(parents), len(recs), one, fib(2 * n - 6), two, (n - 2) * fib(2 * n - 9) - luc(2 * n - 10), no_out1, no_in1, no_out1_one, no_in1_one, time.time() - t0), flush=True)

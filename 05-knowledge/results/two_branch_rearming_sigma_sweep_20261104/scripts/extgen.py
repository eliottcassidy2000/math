"""extgen.py -- generate sigma-ready strong tournaments (digit 1 or 2) by one-vertex extensions from a seed order, canonical
forms via nauty labelg (digraph6).  Usage: python3 extgen.py seed_order max_order
Seeds: all strong sigma-ready tournaments of the seed order from gentourng (via rs17 'pairs'-like filtering done here)."""
import subprocess, sys, collections
seed_n, max_n = int(sys.argv[1]), int(sys.argv[2])
LABELG = '/opt/homebrew/bin/labelg'; GEN = '/opt/homebrew/bin/gentourng'
def rows_from_code(code, n):
    rows = [0] * n; idx = 0
    for i in range(n):
        for j in range(i + 1, n):
            if code[idx] == '1': rows[i] |= 1 << j
            else: rows[j] |= 1 << i
            idx += 1
    return rows
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
# seeds: order seed_n strong tournaments with a digit-1 or digit-2 sole exit: use 'extend' on order seed_n-1? simpler: run rs17 'joincyc'-style via a tiny loop in C? We use extend from seed_n-1 instead.
def canon(d6lines):
    out = subprocess.run([LABELG, '-q'], input='\n'.join(d6lines) + '\n', capture_output=True, text=True).stdout.split()
    return out
# start: order seed_n-1 all tournaments (any), extend to seed_n keeping ready ones -> that gives ALL ready order-seed_n classes (extension of arbitrary tournaments is exhaustive)
n = seed_n - 1
codes = subprocess.run([GEN, '-q', str(n)], capture_output=True, text=True).stdout.split()
inp = '\n'.join(' '.join(map(str, rows_from_code(c, n))) for c in codes) + '\n'
out = subprocess.run(['./rs17', 'extend', str(n)], input=inp, capture_output=True, text=True).stdout.split('\n')
recs = [l.split() for l in out if l.strip()]
d6 = [r[3] for r in recs]; can = canon(d6)
classes = {}
for r, c in zip(recs, can): classes[c] = (int(r[0]), int(r[1]), int(r[2]))
print('order %d: %d ready classes (digit1 %d, digit2 %d) from exhaustive extension of all order-%d tournaments' % (seed_n, len(classes), sum(v[0] for v in classes.values()), sum(v[1] for v in classes.values()), n), flush=True)
cur = classes
for N in range(seed_n + 1, max_n + 1):
    inp = '\n'.join(' '.join(map(str, rows_from_d6(c))) for c in cur) + '\n'
    out = subprocess.run(['./rs17', 'extend', str(N - 1)], input=inp, capture_output=True, text=True).stdout.split('\n')
    recs = [l.split() for l in out if l.strip()]
    can = canon([r[3] for r in recs]); nxt = {}
    for r, c in zip(recs, can): nxt[c] = (int(r[0]), int(r[1]), int(r[2]))
    d1 = {c: v for c, v in nxt.items() if v[0]}; d2 = {c: v for c, v in nxt.items() if v[1]}
    R1 = sorted({v[2] for v in d1.values()}); R2 = sorted({v[2] for v in d2.values()})
    print('order %d: %d ready classes by extension of order-%d READY classes (digit1 %d, digit2 %d); R1: %d values [%d..%d]; R2: %d values [%d..%d]' % (
        N, len(nxt), N - 1, len(d1), len(d2), len(R1), R1[0] if R1 else 0, R1[-1] if R1 else 0, len(R2), R2[0] if R2 else 0, R2[-1] if R2 else 0), flush=True)
    with open('extgen_order%d.txt' % N, 'w') as fh:
        for c, v in nxt.items(): fh.write('%d %d %d %s\n' % (v[0], v[1], v[2], c))
    cur = nxt

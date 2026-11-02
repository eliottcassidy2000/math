import subprocess, itertools, sys, json
exec(open('nearjoin.py').read().split('gen = lambda')[0])
gen = lambda n: subprocess.run(['/opt/homebrew/bin/gentourng', '-q', str(n)], capture_output=True, text=True).stdout.split()
LIMIT = 100000
pairs = {}
for n in range(2, 9):
    for code in gen(n):
        adj = decode(code, n); h = hpaths(adj, list(range(n)))
        for t0 in range(n):
            b = B(adj, n, t0)
            if (h, b) not in pairs: pairs[(h, b)] = (n, code, t0)
    print('order %d done: %d pairs' % (n, len(pairs)), flush=True)
json.dump({'%d,%d' % k: v for k, v in pairs.items()}, open('nearjoin_pairs_le8.json', 'w'))
covered = {}
for (h, b), wit in pairs.items():
    k = 1
    while h + (b << k) <= LIMIT:
        v = h + (b << k)
        if v % 2 == 1 and v not in covered: covered[v] = (h, b, k)
        k += 1
miss = [v for v in range(1, LIMIT + 1, 2) if v not in covered]
print('base orders <= 8: odd numbers <= %d not of the form h + 2^k B: %d; first 30: %s' % (LIMIT, len(miss), miss[:30]))
for X in (1000, 2000, 5000, 10000, 20000, 50000):
    print('  coverage of odd numbers in (%d/2, %d]: %.3f' % (X, X, 1 - sum(1 for v in miss if X // 2 < v <= X) / (X / 4)))

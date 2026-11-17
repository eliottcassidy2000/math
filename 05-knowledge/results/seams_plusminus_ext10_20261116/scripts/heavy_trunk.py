#!/usr/bin/env python3
"""heavy_trunk.py: which order-8 class is the largest branch of the heaviest-parent tree at 8 -> 9 (and of the
least-string tree), with its count, out-degree sequence and strong flag."""
import subprocess
from collections import defaultdict, Counter
exec(open('nearjoin_check.py').read().split('# (1)')[0])
T = '/private/tmp/claude-501/-Users-e-Documents-GitHub-ephrepos-math/d63fb198-a11f-4364-a0ae-e072089f7b73/scratchpad/work/tooth/'
H = defaultdict(list)
for line in open(T + 'counts_z.txt'):
    n, h = map(int, line.split()); H[n].append(h)
def d6_decode(s):
    s = s.strip(); n = ord(s[1]) - 63; bits = ''.join(format(ord(c) - 63, '06b') for c in s[2:])
    return n, [[int(bits[i*n+j]) for j in range(n)] for i in range(n)]
raw8 = subprocess.run(['gentourng', '-q', '-z', '8'], capture_output=True, text=True).stdout.split()
parents = [list(map(int, l.split())) for l in open(T + 'parents_9.txt')]
heavy = Counter(max(ps, key=lambda p: (H[8][p], -p)) for ps in parents); least = Counter(int(x) for x in open(T + 'cparent_9.txt'))
for name, tree in (('heaviest-parent', heavy), ('least-string', least)):
    p, b = tree.most_common(1)[0]; n, A = d6_decode(raw8[p])
    print(f"{name} tree, largest branch at 8 -> 9: order-8 class {p} with {b} children; H = {H[8][p]}, out-degrees {sorted(sum(r) for r in A)}, strong {strong(A)}")
    top = tree.most_common(5); print("   top five branches (class, children, H):", [(q, c, H[8][q]) for q, c in top])

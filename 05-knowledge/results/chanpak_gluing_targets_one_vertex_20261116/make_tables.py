"""Build the target map (TARGET_MAP.md tables) from rows_census.out and jrow_forms.out.
rows_census.out lines: 'R N p q c n' = some tournament of order n with H = N has induced
subtournaments of counts p and q (odd primes, both != 7) with p + q = N + c, c in {1,-1,-3}.
jrow_forms.out lines: 'F p=.. q=.. N=.. n=..' with F in C (contained), K (cover), D (vertex exchange)."""
import collections, re
def isprime(x): return x > 1 and all(x % d for d in range(2, int(x**.5) + 1))
real = collections.defaultdict(set)
for l in open('rows_census.out'):
    _, N, p, q, c, n = l.split(); real[int(N)].add((int(p), int(q), int(c)))
forms = collections.defaultdict(dict)
for l in open('jrow_forms.out'):
    m = re.match(r'(\w) p=(\d+) q=(\d+) N=(\d+) n=(\d+)', l)
    F, p, q, N, n = m.group(1), *map(int, m.groups()[1:])
    forms[(p, q, N)].setdefault(F, n)
lab = {1: 'p+q-1', -1: 'p+q+1', -3: 'p+q+3'}
print('| N | Goldbach pairs p+q=N+1 (p,q>=11) | induced (J) pairs [forms, least order] | rows -1 / -3 realised (any primes != 7) |')
print('|--:|:--|:--|:--|')
noneJ, noneAll = [], []
for N in range(23, 125, 2):
    gp = [(p, N + 1 - p) for p in range(11, (N + 1)//2 + 1) if isprime(p) and isprime(N + 1 - p)]
    J = sorted(x for x in real.get(N, ()) if x[2] == 1)
    other = sorted((x for x in real.get(N, ()) if x[2] != 1), key=lambda x: (-x[2], x[0]))
    js = ', '.join(f"({p},{q}) [{''.join(sorted(forms[(p, q, N)].keys()))} {min(forms[(p, q, N)].values())}]" for p, q, c in J) or '**none**'
    os_ = ', '.join(f"{p}+{q}{'+1' if c == -1 else '+3'}" for p, q, c in other) or '**none**'
    if not J: noneJ.append(N)
    if not real.get(N): noneAll.append(N)
    print(f"| {N} | {' '.join(f'({p},{q})' for p, q in gp) or '-'} | {js} | {os_} |")
print()
print('targets with no induced (J) pair:', len(noneJ), noneJ)
print('targets with no induced pair in any of the three rows:', len(noneAll), noneAll)

from itertools import product


def need(condition, detail):
    if not condition:
        raise AssertionError(detail)


def word(bits):
    rows = [2, 0]
    for b in map(int, bits):
        n = len(rows)
        rows = [a | ((1 << n) if i or not b else 0)
                for i, a in enumerate(rows)] + [b]
    return rows


def rooted_state(bits):
    h, e, t = 1, 0, 2
    for b in map(int, bits):
        h, e, t = (h+t-e, 1, 2*t-1) if b else (h, 0, 2*t)
    return h, e, t


def flip(rows, edges):
    rows = rows.copy()
    for i, j in edges:
        rows[i] ^= 1 << j
        rows[j] ^= 1 << i
    return rows


def blow(rows, widths):
    labels = [i for i in range(len(rows)) for _ in range(widths.get(i, 1))]
    return [sum(1 << j for j, v in enumerate(labels)
                if (i < j if u == v else rows[u] >> v & 1))
            for i, u in enumerate(labels)]


def ham(rows):
    n = len(rows)
    need(all(not (a >> i & 1) for i, a in enumerate(rows)), rows)
    need(all((rows[i] >> j & 1)+(rows[j] >> i & 1) == 1
             for i in range(n) for j in range(i)), rows)
    dp = [[0]*n for _ in range(1 << n)]
    for i in range(n):
        dp[1 << i][i] = 1
    for s, row in enumerate(dp):
        for i, a in enumerate(row):
            if not a:
                continue
            targets = rows[i] & ~s
            while targets:
                bit = targets & -targets
                targets -= bit
                dp[s|bit][bit.bit_length()-1] += a
    return sum(dp[-1])


BASE_PATCH = [(0, 1), (0, 2), (0, 5), (1, 7), (3, 4), (3, 7)]
CORE = word('010'+'1'*5)
PATCH = [(i+(i>=3), j+(j>=3)) for i, j in BASE_PATCH]
FIXED = flip(CORE, PATCH)
need(FIXED == [80, 765, 1017, 1009, 704, 977, 896, 769, 531, 1], FIXED)
for q in (CORE, FIXED):
    need(q[9] == 1, q)
    need([i for i in range(len(q)) if q[i] >> 3 & 1] == [1,2], q)

for s, m in product(range(4), range(1, 3)):
    ww = '0'+'1'*s+'0'+'1'*(m+4)
    widths = {3: s, 9: m}
    a, b = blow(CORE, widths), blow(FIXED, widths)
    need(a == word(ww), ww)
    N = (12*2**s+2)*2**(m+3)+1
    D = 6 if s == 0 else 2+3*2**s
    need((ham(a), ham(b)) == (N-D, N), (s,m))
    h,e,t = rooted_state(ww)
    need((t,h) == (N,N-D), (s,m,t,h))

seeds = [
    ([46,60,24,48,33,4],5,23,11),
    ([126,252,248,48,96,192,8,89],6,69,34),
    ([8,9,19,52,67,343,399,319,31],7,337,168),
    ([126,124,248,48,224,192,8,75],6,121,60),
]
for rows, mark, h, c2 in seeds:
    need(min(rows[mark].bit_count(),len(rows)-1-rows[mark].bit_count()) == 1, rows)
    need((ham(rows), ham(blow(rows,{mark:2}))) == (h,h+2*c2), rows)

SMALL = [
    ('0101',27,[(0,3),(1,3)]),
    ('01011',53,[(1,3)]),
    ('010111',105,[(0,2),(1,5)]),
    ('0101111',209,[(0,7),(4,6)]),
    ('11101',35,[(1,4)]),
    ('10101',43,[(1,3)]),
    ('101011',85,[(0,6),(1,4)]),
    ('1010111',169,[(4,5)]),
    ('0001',31,[(0,2),(1,4)]),
    ('00011',61,[(0,1),(1,3),(4,6)]),
    ('01101',51,[(0,2),(3,5)]),
    ('011011',101,[(0,1),(1,6),(5,6)]),
    ('0110111',201,[(4,5)]),
    ('01101111',401,[(0,2),(1,6)]),
    ('11001',39,[(0,2),(3,4)]),
    ('110011',77,[(0,4),(5,7)]),
    ('1100111',153,[(0,7),(5,8)]),
    ('11001111',305,[(0,2),(0,8)]),
]
for w, N, P in SMALL:
    h,_,t = rooted_state(w)
    need((ham(word(w)),ham(flip(word(w),P)),t) == (h,N,N), w)

for k in range(1,13):
    for w in product('01',repeat=k):
        w = ''.join(w)
        if not w.endswith('1'):
            continue
        h,_,t = rooted_state(w)
        d8 = any(w.startswith(pre+'1') and set(w[len(pre):]) == {'1'}
                 for pre in ('010','100','1110'))
        d12 = w.startswith('10101') and set(w[4:]) == {'1'}
        d14 = any(w.startswith(pre+'1') and set(w[len(pre):]) == {'1'}
                  for pre in ('000','1100','0110'))
        need(((t-h)==8) == d8 and (t-h)!=10 and (((t-h)==12)==d12)
             and (((t-h)==14)==d14), w)

# Signed inclusion-exclusion of the six possible arc reversals.
diffs = []
for mask in range(64):
    q = flip(CORE,[arc for i,arc in enumerate(PATCH) if mask>>i & 1])
    h = ham(q)
    diffs.append([h,(ham(blow(q,{9:2}))-h)//2])
need([i for i, v in enumerate(diffs) if v[1]==208] == [0,63], diffs)
for i in range(6):
    for mask in range(64):
        if mask>>i & 1:
            diffs[mask] = [a-b for a,b in zip(diffs[mask],diffs[mask^(1<<i)])]
degrees = [tuple(sum(v[i] for mask,v in enumerate(diffs) if mask.bit_count()==k)
                 for i in range(2)) for k in range(1,7)]
need(degrees == [(1512,815),(-1968,-976),(480,153),(-12,8),(-4,0),(0,0)],degrees)
print('PASS: tensor grid, exact defect-eight witnesses, and signed patch coefficients')

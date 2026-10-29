import re, itertools, sys
from runs import tour_from_bits
def ham_cycles(A):
    n = len(A)
    cnt = 0
    for w in itertools.permutations(range(1, n)):
        cyc = (0,) + w
        if all(A[cyc[i]][cyc[(i + 1) % n]] for i in range(n)):
            cnt += 1
    return cnt
for f in sys.argv[1:]:
    n = int(re.findall(r'n(\d+)_', f)[0])
    uni_hc = {}
    nonuni_hc1 = 0
    for l in open(f):
        if not re.match(r'^[01]{3,}\s', l) or ' strong ' not in l:
            continue
        bits = l.split()[0]
        m = re.search(r' M=(.*?)( unglued=| P=|$)', l)
        uni = all(t.startswith('1^') for t in m.group(1).split())
        hc = ham_cycles(tour_from_bits(n, bits))
        if uni:
            uni_hc[hc] = uni_hc.get(hc, 0) + 1
        elif hc == 1:
            nonuni_hc1 += 1
    print(f, 'unimodular strong classes by #Hamiltonian cycles:', uni_hc, '| non-unimodular with a unique HC:', nonuni_hc1)

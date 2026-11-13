"""Test the 'j = 0 continuity' identity  sum_{i>=0} (-1)^i a_i(S,u) = 0  (a_0 = H(S-u)), i.e. the exponential
polynomial H(S[u<-TT_j]) = sum_{l>=1} c_l l^j, extended to j = 0, returns H(S-u).  Also the two-module analogue."""
import random, sys
sys.path.insert(0,'scripts')
from identity_test import words_by_visits
from twomod import visit_table
rng = random.Random(11)
ok = bad = 0
for trial in range(4000):
    n = rng.randint(3, 8)
    adj = [0]*n
    for i in range(n):
        for j in range(i+1, n):
            if rng.random() < 0.5: adj[i] |= 1 << j
            else: adj[j] |= 1 << i
    u = rng.randrange(n)
    t = words_by_visits(adj, n, u, cap=8)
    s = sum((-1)**i * c for i, c in t.items())
    if s == 0: ok += 1
    else:
        bad += 1
        if bad <= 5: print('fail', n, u, t, s)
print('single-module alternating sum = 0:', ok, 'fails:', bad)
# two-module: sum_{i,i'} (-1)^{i+i'} a_{ii'} = 0 ?  and the mixed (j=0 in u only) version
ok2 = bad2 = 0
for trial in range(1500):
    n = rng.randint(4, 8)
    adj = [0]*n
    for i in range(n):
        for j in range(i+1, n):
            if rng.random() < 0.5: adj[i] |= 1 << j
            else: adj[j] |= 1 << i
    u, v = rng.sample(range(n), 2)
    tab = visit_table(adj, n, u, v, CU=8, CV=8)
    for ip in range(0, 5):
        s = sum((-1)**i * tab.get((i, ip), 0) for i in range(0, 9))
        if s == 0: ok2 += 1
        else: bad2 += 1
print('two-module, per v-multiplicity alternating u-sum = 0:', ok2, 'fails:', bad2)

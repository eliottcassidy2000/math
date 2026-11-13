"""Distribution of the P_M(1)P_N(1) coefficient c11 = sum (-1)^(i+i') i i' a_{ii'} over two module positions (u,v),
and of the single-module c1 = sum (-1)^(i-1) i a_i (the constant of the exponential polynomial)."""
import random, sys, collections
sys.path.insert(0,'scripts')
from twomod import visit_table
from identity_test import words_by_visits
rng = random.Random(3)
C11 = collections.Counter(); C1 = collections.Counter()
for trial in range(1200):
    n = rng.randint(3, 8)
    adj = [0]*n
    for i in range(n):
        for j in range(i+1, n):
            if rng.random() < 0.5: adj[i] |= 1 << j
            else: adj[j] |= 1 << i
    u, v = rng.sample(range(n), 2)
    tab = visit_table(adj, n, u, v, CU=9, CV=9)
    c11 = sum((-1)**(i+ip) * i * ip * c for (i, ip), c in tab.items() if i >= 1 and ip >= 1)
    C11[c11] += 1
    t = words_by_visits(adj, n, u, cap=9)
    c1 = sum((-1)**(i-1) * i * c for i, c in t.items() if i >= 1)
    C1[c1] += 1
print('c11 values (count):', sorted(C11.items())[:40])
print('c1 values (count):', sorted(C1.items())[:40])

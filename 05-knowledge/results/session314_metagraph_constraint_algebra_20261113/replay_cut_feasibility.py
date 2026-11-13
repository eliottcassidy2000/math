from itertools import product
from collections import defaultdict

def lawful(scores):
    dp = {(0, 0): 0}
    for s in scores:
        new = {}
        for (k, prior), w in dp.items():
            for now in (0, 1):
                q = (k + now, now)
                val = w + now*s - prior*(1-now)
                new[q] = min(new.get(q, val), val)
        dp = new
    return all(w >= k*(k-1)//2 for (k, _), w in dp.items())

def laurent_support(n):
    m = n//2
    factors = []
    for i in range(m):
        for j in range(i+1, m):
            if j >= i+2:
                factors.append(tuple((k == j) - (k == i) for k in range(m)))
            factors.append(tuple(-int(k in (i, j)) for k in range(m)))
        for _ in range(int(i < m-1 or n%2) + int(n%2 and i < m-1)):
            factors.append(tuple(-int(k == i) for k in range(m)))
    supp = {tuple(n-1-i for i in range(m))}
    for step in factors:
        supp |= {tuple(a+b for a,b in zip(word,step)) for word in supp}
    return supp

for n in range(2,10):
    m=n//2
    cut = set()
    for x in product(range(n),repeat=m):
        scores = x + ((m,) if n%2 else ()) + tuple(n-1-v for v in x[::-1])
        if lawful(scores): cut.add(x)
    assert cut == laurent_support(n)
    print(n,len(cut))

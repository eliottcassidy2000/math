"""R = S (mod 2) for every rooted tournament (THM-4878 (B)): exhaustive labelled check of SEAM triples to order 6
and random orders 7-10. SEAM's response() is used for (H, S, R)."""
import sys, os, itertools, random
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(here, '..', '..', '..', '05-knowledge', 'results', 'seven_return_native_recursions_20261029'))
from seam_transfer import response, tournament
tot = bad = 0
for n in range(1, 7):
    pairs = list(itertools.combinations(range(n), 2))
    for mask in range(1 << len(pairs)):
        adj = tournament(n, mask)
        for a in range(n):
            H, S, R = response(adj, a); tot += 1; bad += (R - S) % 2
print('labelled rooted tournaments of order 1-6:', tot, 'violations:', bad)
random.seed(20261104); bad = 0
for _ in range(300):
    n = random.randint(7, 10); adj = tournament(n, random.getrandbits(n * (n - 1) // 2)); a = random.randrange(n)
    H, S, R = response(adj, a); bad += (R - S) % 2
print('random rooted tournaments of order 7-10 (300):', 'violations:', bad)

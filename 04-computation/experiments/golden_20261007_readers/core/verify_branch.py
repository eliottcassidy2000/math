"""Independent brute-force check of the C sieve (J=0, BRANCH=1, DESC) for small K.
For each class r mod 2^K: three large representatives n = r + 2^K q.  Forward orbit x_s (integers).
For each s <= K and EVERY backward word u (length i < s, b <= a_s odd steps), apply u to the integer x_s
(odd backward step needs v = 2 mod 3: v -> (2v-1)/3; even: v -> 2v).  (s,u) certifies the class iff it is valid and
gives m < n for all three representatives; we also re-check T^i(m) = x_s by forward iteration.
Descent: x_s < n for all three.  Prints #uncertified classes at depth K."""
import random, itertools, sys
def T(x): return x//2 if x%2==0 else (3*x+1)//2
random.seed(5)
def certified(r, K, reps):
    orbs=[]
    for n in reps:
        xs=[n]; a=[0]
        for s in range(K):
            a.append(a[-1]+(xs[-1]&1)); xs.append(T(xs[-1]))
        orbs.append((n,xs,a))
    a=orbs[0][2]
    for s in range(1,K+1):
        if all(xs[s] < n for n,xs,_ in orbs): return 'desc'
        for i in range(0, s):
            for u in itertools.product((0,1), repeat=i):
                if sum(u) > a[s]: continue
                ok=True
                for n,xs,_ in orbs:
                    v=xs[s]
                    for bit in u:
                        if bit:
                            if v%3!=2: ok=False;break
                            v=(2*v-1)//3
                        else: v=2*v
                    if not ok or not (v < n): ok=False;break
                    w=v
                    for _ in range(i): w=T(w)
                    assert w==xs[s]
                if ok: return 'branch'
    return None
for K in range(4, int(sys.argv[1])+1):
    unc=0; nb=0
    for r in range(2**K):
        reps=[r+(random.getrandbits(80)<<K) for _ in range(3)]
        c=certified(r,K,reps)
        if c is None: unc+=1
        elif c=='branch': nb+=1
    print(K, "uncertified", unc, "branch-cert (not desc at that s)", nb, flush=True)

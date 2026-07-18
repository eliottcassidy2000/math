# Independent verifier of THM-1121's r=6 universal weighted atlas.
# Obligations transcribed BY HAND from THM-1121-....md (NOT from codex's script).
import numpy as np

ATLAS = [
 (26,2,81),(27,2,7),(28,2,3),(28,13,7),
 (40,3,10),(41,3,3),(41,19,5),(42,13,1),
 (53,4,16),(55,4,2),(56,13,21),
 (68,5,3),(68,21,21),(69,5,1),(69,16,12),(70,27,20),
 (79,6,3),(81,25,5),(82,19,6),(83,6,1),(83,32,5),
 (84,13,8),(93,7,9),(94,7,13),(97,7,1),(97,30,2),
 (98,15,6),(105,8,4),(106,49,46),(107,33,10),
 (109,42,49),(110,17,25),(111,8,1),(111,17,64),
 (112,43,34),
]

def h(q):            # ceil(q/14)
    return -(-q//14)
def la(x,q):         # lattice distance
    r = x % q
    return min(r, q-r)

print("obligations:", len(ATLAS))
print("total weight:", sum(w for _,_,w in ATLAS))
print("modulus range:", min(q for q,_,_ in ATLAS), "..", max(q for q,_,_ in ATLAS))

# Property 2: every obligation safe for every core speed p in 1..12
min_slack = 10**9
bad=[]
for (q,a,w) in ATLAS:
    for p in range(1,13):
        slack = la(p*a,q) - h(q)
        if slack < 0: bad.append((q,a,p,slack))
        min_slack = min(min_slack, slack)
print("PROP2 core-safety: min slack over 35x12 =", min_slack, "-> SAFE" if min_slack>=0 else "-> VIOLATION")
if bad: print("  violations:", bad[:10])

# Property 3: load over [92,332]
def load(k):
    s=0
    for (q,a,w) in ATLAS:
        if la(k*a,q) < h(q): s+=w
    return s
loads_fin = [load(k) for k in range(92,333)]
print("PROP3 max load over [92,332]:", max(loads_fin), "(codex claims 84)")
print("  6*max =", 6*max(loads_fin), "< 505 ?", 6*max(loads_fin) < 505)

# ---- NEW: does load<=84 extend beyond 332?  Scan a large range with numpy. ----
# load(k) = sum_i w_i * [ la(k*a_i, q_i) < h(q_i) ]  -- periodic in k mod q_i.
# Precompute per-obligation killed-residue boolean array of length q_i.
kill_arrays = []
for (q,a,w) in ATLAS:
    arr = np.zeros(q, dtype=np.int64)
    for s in range(q):
        if la(s*a, q) < h(q):
            arr[s] = w
    kill_arrays.append((q,arr))

def max_load_range(k0,k1,chunk=2_000_000):
    best=0; argbest=-1; firstbreak=None
    k=k0
    while k<k1:
        hi=min(k+chunk,k1)
        ks=np.arange(k,hi)
        tot=np.zeros(hi-k,dtype=np.int64)
        for (q,arr) in kill_arrays:
            tot += arr[ks % q]
        m=int(tot.max()); am=int(ks[int(tot.argmax())])
        if m>best: best=m; argbest=am
        if firstbreak is None:
            over = ks[tot>84]
            if over.size>0: firstbreak=int(over[0])
        k=hi
    return best,argbest,firstbreak

for hi in [1000, 10000, 100000, 1_000_000, 20_000_000]:
    b,ab,fb = max_load_range(92,hi)
    print(f"scan [92,{hi}): max load = {b} at k={ab}; first k with load>84 = {fb}")

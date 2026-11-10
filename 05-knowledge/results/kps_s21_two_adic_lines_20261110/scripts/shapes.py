import sys, collections
def ispow2(v): return v > 0 and v & (v-1) == 0
def popc(v): return bin(v).count("1")
tuples = collections.Counter(); wit = {}
for fn in sys.argv[1:]:
    for line in open(fn):
        p = line.split()
        if p[7] != "1": continue
        t = tuple(int(v) for v in p[3:7])
        tuples[t] += 1
        wit.setdefault(t, (fn, p[0], p[1], p[2]))
print("distinct exact tuples:", len(tuples))
# weight-3 test on a grid: fraction of (a,b) in 1..8 x 1..8 with popcount(h)=3
res = []
for t, c in tuples.items():
    al, be, ga, de = t
    good = 0; tot = 0; vals = set()
    for a in range(1, 9):
        for b in range(1, 9):
            h = al + be*2**a + ga*2**b + de*2**(a+b)
            tot += 1
            if h > 0 and popc(h) == 3: good += 1; vals.add(h)
    res.append((good, len(vals), t, c))
res.sort(reverse=True)
print("top tuples by number of grid points with binary weight 3 (good, distinct values, tuple, multiplicity, witness):")
for good, nv, t, c in res[:40]:
    print(good, nv, t, c, wit[t])

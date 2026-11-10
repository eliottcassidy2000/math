import sys, collections
tuples = set()
for fn in sys.argv[1:]:
    for line in open(fn):
        p = line.split()
        if p[7] == "1": tuples.add(tuple(int(v) for v in p[3:7]))
AM = 26
hit = collections.defaultdict(set)
for t in tuples:
    al, be, ga, de = t
    for a in range(1, AM+1):
        for b in range(1, AM+1):
            h = al + be*2**a + ga*2**b + de*2**(a+b)
            if h > 0 and bin(h).count("1") == 3 and h & 1: hit[h].add(t)
A = 20
print("weight-3 targets 1+2^a+2^b, a<=%d: rows d=a-b, columns b=1.." % A)
for d in range(1, 13):
    row = ""
    for b in range(1, A - d + 1):
        a = b + d; N = 1 + 2**a + 2**b
        row += "#" if N in hit else "."
    print("d=%2d %s" % (d, row))
print("positive tuples (1,be,ga,de>=0):", sorted(t for t in tuples if t[0] == 1 and min(t[1:]) >= 0))

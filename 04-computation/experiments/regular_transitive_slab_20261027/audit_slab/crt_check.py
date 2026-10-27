rows = {}
for line in open("/Users/e/Documents/GitHub/math-wt-pyth-cuboids/04-computation/experiments/regular_transitive_slab_20261027/special15_rows.txt"):
    f = line.split(); rows[f[0]] = [int(x, 16) for x in f[1:]]
n = 15
def check(name, blocksize, nblocks):
    o = rows[name]; A = [[(o[i] >> j) & 1 for j in range(n)] for i in range(n)]
    assert all(A[i][j] + A[j][i] == 1 for i in range(n) for j in range(n) if i != j)
    # vertex v = blocksize*b + a  ->  k in Z_15 with k = b mod nblocks, k = a mod blocksize (CRT)
    lab = {}
    for v in range(n):
        b, a = divmod(v, blocksize)
        k = next(k for k in range(n) if k % nblocks == b and k % blocksize == a)
        lab[v] = k
    S = set(); ok = True
    for i in range(n):
        for j in range(n):
            if i != j and A[i][j]: S.add((lab[j] - lab[i]) % n)
    # circulant iff the arc set is exactly {k -> k+d : d in S} and S, -S partition Z_15 \ {0}
    ok = all(A[i][j] == (1 if (lab[j] - lab[i]) % n in S else 0) for i in range(n) for j in range(n) if i != j)
    print(name, "is the Z_15 circulant with connection set", sorted(S), ":", ok and len(S) == 7 and not (S & {(-d) % n for d in S}))
check("R5[C3]", 3, 5)
check("C3[R5]", 5, 3)

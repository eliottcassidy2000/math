"""D4 check: number of distinct residual maps of the parity-address map q_+ after
prefixes of length k, truncated to the next r output bits (functions of the next r
input bits).  A synchronous deterministic transducer computing q_+ needs at least this
many distinct states at time k.  Compare plus, minus, and 5n+1."""
def T(n, q, s): return n // 2 if n % 2 == 0 else (q * n + s) // 2
def qmap(n, q, s, d):
    out = 0; x = n
    for j in range(d):
        out |= (x & 1) << j; x = T(x, q, s)
    return out
for (q, s, name) in ((3, 1, "plus"), (3, -1, "minus"), (5, 1, "5n+1")):
    row = []
    for k in range(1, 9):
        r = k  # look ahead r bits
        d = k + r
        funcs = set()
        for p in range(1 << k):
            f = tuple((qmap(p + (m << k), q, s, d) >> k) for m in range(1 << r))
            funcs.add(f)
        row.append((k, len(funcs)))
    print(name, "(k, #distinct residual maps on next k bits):", row)

# v_invol.py - involution pairs {h, iota h} per cell: claim = 2^(n-2-t') unordered pairs in cell (t', t'-1)
def rev(h, n):
    r = 0
    for i in range(n):
        if (h >> i) & 1: r |= 1 << (n - 1 - i)
    return r
for n in (12, 13, 16):
    m = n - 1
    cells = {}
    fixed = 0
    for r in range(1 << m):
        h = (1 << m) | r
        hp = (1 << m) | rev((~r) & ((1 << m) - 1), m)
        if h == hp: fixed += 1; continue
        if h < hp:
            d = h ^ hp
            tp = n - d.bit_length()
            s = (d & -d).bit_length() - 1
            cells[(tp, s)] = cells.get((tp, s), 0) + 1
    ok = all(s == tp - 1 for (tp, s) in cells)
    print("n=%d fixed points=%d ; all involution pairs in cells (t', t'-1): %s ; counts/2^(n-2-t'): %s" % (
        n, fixed, ok, " ".join("%d:%.3f" % (tp, c / 2.0 ** (n - 2 - tp)) for (tp, s), c in sorted(cells.items()))))

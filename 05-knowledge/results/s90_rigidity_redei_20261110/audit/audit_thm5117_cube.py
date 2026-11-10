#!/usr/bin/env python3
"""Is THM-5117's "exact affine scalar cube" an exchange cube in the sense of THM-5122 (disjoint pair exchanges in one
deletion core, protected closure)?  Recount of the case k = 2 (closure of order 15) by a subset dynamic programme.
(THM-5117 proves the cube for every k; only k = 2 is recounted here, to confirm the reading of the construction.)

Construction read from THM-5117: module 1 = core [18, 8, 3, 5, 14] (H = 13) plus an unselected sink, cuts {0, 1} or
{1, 2} (the exchange 0 <-> 2); module 2 = core [32, 17, 3, 7, 13, 30] (H = 37) plus an unselected sink, cuts {0, 4} or
{0, 5} (the exchange 4 <-> 5); ordinal composition (earlier module -> later module), then the protected closure
(global sink v, root r beating exactly the selected vertices, v -> r).
Expected from THM-5117 for k = 2: coefficients -2 * 37 = -74 and 124, an affine cube with both increments nonzero.
"""


def H(out):
    n = len(out)
    cnt = [dict() for _ in range(1 << n)]
    for v in range(n):
        cnt[1 << v][v] = 1
    total = 0
    full = (1 << n) - 1
    for mask in range(1, 1 << n):
        row = cnt[mask]
        if not row:
            continue
        if mask == full:
            total = sum(row.values())
            break
        for v, c in row.items():
            nxt = out[v] & ~mask
            while nxt:
                low = nxt & -nxt
                nxt ^= low
                d = cnt[mask | low]
                w = low.bit_length() - 1
                d[w] = d.get(w, 0) + c
        cnt[mask] = None
    return total


def module(core):
    n = len(core)
    return [core[i] | (1 << n) for i in range(n)] + [0]       # unselected sink of the module


def compose(mods):
    out = []
    off = 0
    sizes = [len(m) for m in mods]
    total = sum(sizes)
    for m in mods:
        later = ((1 << total) - 1) & ~((1 << (off + len(m))) - 1)
        for row in m:
            out.append((row << off) | later)
        off += len(m)
    return out


def closure(C, sel):
    n = len(C)
    adj = [C[i] | (1 << n) | (0 if sel >> i & 1 else 1 << (n + 1)) for i in range(n)]
    adj.append(1 << (n + 1))
    adj.append(sel)
    return adj


def main():
    c13 = [18, 8, 3, 5, 14]
    c37 = [32, 17, 3, 7, 13, 30]
    m13 = module(c13)
    m37 = module(c37)
    print("H(core13) =", H(c13), " H(core37) =", H(c37))
    for name, mods, cuts in (("h13 then h37", [m13, m37], [((0, 1), (1, 2)), ((0, 4), (0, 5))]),
                             ("h37 then h13", [m37, m13], [((0, 4), (0, 5)), ((0, 1), (1, 2))])):
        C = compose(mods)
        offs = [0, len(mods[0])]
        vals = {}
        for b0 in (0, 1):
            for b1 in (0, 1):
                sel = 0
                for off, cut, b in zip(offs, cuts, (b0, b1)):
                    for v in cut[b]:
                        sel |= 1 << (off + v)
                vals[(b0, b1)] = H(closure(C, sel))
        g0 = vals[(1, 0)] - vals[(0, 0)]
        g1 = vals[(0, 1)] - vals[(0, 0)]
        mixed = vals[(1, 1)] - vals[(1, 0)] - vals[(0, 1)] + vals[(0, 0)]
        print(f"{name}: closure order {len(C) + 2}; counts {vals}; increments {g0}, {g1}; mixed difference {mixed}; "
              f"affine with both increments nonzero: {mixed == 0 and g0 != 0 and g1 != 0}")


if __name__ == "__main__":
    main()

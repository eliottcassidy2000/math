# Compensated banks: keep a whole line of exterior (tape) choices alive.
# With the lowest upper vertex a non-contact and the highest lower vertex a contact (u1 = z1 = 0) the count is
#     H = h (1 + 2 N_up + 2 N_low) + 2 [ c + a alpha + b beta + 2 G alpha beta ],   a = a_o + a_e, b = b_o + b_e, G = G_e + G_o.
# Seeds with the same (h, G) whose (a, b) lie on a line of direction (p, -q), i.e. a_i = a_0 + p t_i, b_i = b_0 - q t_i,
# have counts differing by 2 (c_i - c_0 + t_i rho) for EVERY exterior on the line  p alpha - q beta = rho.
# A line bank of length L: L such seeds whose values c_i + t_i rho are consecutive.   (p, q) = (0, 0) is the synchronised case.
import sys, collections
seeds = []
for line in open(sys.argv[1]):
    t = line.split(); h, Ge, Go, ao, ae, bo, be, cut = map(int, t[:8]); m = int(t[10])
    seeds.append((h, Ge + Go, ao + ae, bo + be, cut, int(t[9]), tuple(int(v) for v in t[11:11 + m])))
print("seeds:", len(seeds))
dirs = [(0, 0), (1, 0), (0, 1), (1, 1), (1, 2), (2, 1), (1, 3), (3, 1), (2, 3), (3, 2), (1, 4), (4, 1), (3, 4), (4, 3)]
RHO = range(-60, 61)
best = {}
for (p, q) in dirs:
    groups = collections.defaultdict(list)
    for sd in seeds:
        h, G, a, b, c = sd[:5]
        if (p, q) == (0, 0): key = (h, G, a, b); t = 0
        else:
            key = (h, G, q * a + p * b, a % p if p else 0, b % q if q else 0); t = a // p if p else -(b // q)
        groups[key].append((t, c, sd))
    top = (1, None)
    for key, lst in groups.items():
        if len(lst) < 2: continue
        rhos = [0] if (p, q) == (0, 0) else ([r for r in RHO if r >= 0] if q == 0 else ([r for r in RHO if r <= 0] if p == 0 else RHO))
        for rho in rhos:
            vals = {}
            for t, c, sd in lst: vals.setdefault(c + t * rho, sd)
            vs = sorted(vals); run = 1; start = vs[0]; cur = 1; cs = vs[0]
            for i in range(1, len(vs)):
                if vs[i] == vs[i - 1] + 1: cur += 1
                else: cur = 1; cs = vs[i]
                if cur > run: run = cur; start = cs
            if run > top[0]: top = (run, (key, rho, start, [vals[start + i] for i in range(run)]))
    best[(p, q)] = top
    run, info = top
    if info:
        key, rho, start, members = info
        print("direction (p, q) = (%d, %d): longest line bank %d  [h = %d, G = %d, exterior line %d*alpha - %d*beta = %d]" % (p, q, run, key[0], key[1], p, q, rho))
        for sd in members[:8]: print("      a = %d b = %d cut = %d  block %s pattern %d" % (sd[2], sd[3], sd[4], ' '.join(map(str, sd[6])), sd[5]))

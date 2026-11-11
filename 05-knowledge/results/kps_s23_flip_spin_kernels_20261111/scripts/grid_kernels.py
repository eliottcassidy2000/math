"""grid_kernels.py - which two-port kernels cover a whole grid uniformly?

mf_m.txt (S21/S22 scratch): code x y al be ga de exact   with  H(K[TT_a at x, TT_b at y]) = al + be 2^a + ga 2^b + de 2^(a+b)
(exact = 1: the four-term formula holds for all a, b >= 1, i.e. the ports are visited at most twice).
In the coordinates (h, pi) = (H, Pi_2) of the two modules the kernel form is  al h h' + be pi h' + ga h pi' + de pi pi'.

A kernel with (al, be, ga, de) = (1, 2^u, 0, d 2^w) (or be <-> ga) gives  1 + 2^(a+u) + d 2^(a+b+w): the d-grid
N(d, r, s) = 1 + 2^s + d 2^(r+s) for all s >= u + 1, r >= w - u + 1."""
import sys, collections


def odd(x):
    v = 0
    while x % 2 == 0:
        x //= 2
        v += 1
    return x, v


def ispow2(x):
    return x > 0 and x & (x - 1) == 0


files = sys.argv[1:]
forms = collections.Counter()
first = {}
sig = collections.Counter()
for fn in files:
    m = int(fn.split("mf_")[1].split(".")[0])
    for l in open(fn):
        t = l.split()
        code, x, y, al, be, ga, de, ex = map(int, t[:8])
        if not ex:
            continue
        f = (al, be, ga, de)
        forms[f] += 1
        first.setdefault(f, (m, code, x, y))
        det = al * de - be * ga
        sig[(m, "elliptic" if det > 0 else ("parabolic" if det == 0 else "hyperbolic"))] += 1
print("distinct exact forms (al, be, ga, de):", len(forms))
print("signature by order (sign of al de - be ga):", sorted(sig.items()))
grid = {}
for f, n in forms.items():
    al, be, ga, de = f
    if al != 1 or de <= 0:
        continue
    for (p, q) in ((be, ga), (ga, be)):
        if q == 0 and ispow2(p):
            d, w = odd(de)
            u = p.bit_length() - 1
            key = d
            val = (u, w, first[f])
            if key not in grid or (u, w) < grid[key][:2]:
                grid[key] = val
print("grids d with a uniform kernel 1 + 2^(a+u) + d 2^(a+b+w) (d: (u, w, (order, code, x, y))):")
for d in sorted(grid):
    print("   d = %d: u = %d, w = %d, kernel %s  -> all s >= %d, r >= %d" % (d, grid[d][0], grid[d][1], grid[d][2], grid[d][0] + 1, max(1, grid[d][1] - grid[d][0] + 1)))
# diagonal forms (be = ga = 0): one-port pairs with joined module; offsets al and slopes de
diag = sorted({(al, de) for (al, be, ga, de) in forms if be == 0 and ga == 0})
print("diagonal forms (al, de), al = offset, de = slope:", len(diag), "; negative offsets:", [x for x in diag if x[0] < 0][:40])
# forms with al = -1: anti-grids  -1 + ...
anti = {}
for f in forms:
    al, be, ga, de = f
    if al != -1 or de <= 0:
        continue
    for (p, q) in ((be, ga), (ga, be)):
        if q == 0 and ispow2(-p) if p < 0 else False:
            d, w = odd(de)
            anti.setdefault(d, (f, first[f]))
print("anti-grids -1 - 2^(a+u) + d 2^(a+b+w): d ->", {d: anti[d] for d in sorted(anti)})

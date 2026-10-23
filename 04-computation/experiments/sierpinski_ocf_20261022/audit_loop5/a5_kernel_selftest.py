# Self-test of a5_kernel.c against the THM-4757 section 4 two-loop kernel table (eps_a = (-1)^(a/2)):
#   figure-eight (a,b) even: 8(eps_a + eps_b) - 12 eps_a eps_b; theta (a,b,c) all even: 4(eps_a+eps_b+eps_c) - 8 eps_a eps_b eps_c,
#   a,b odd and c even: -4 eps_c, otherwise 0; dumbbell (a,b;c): -32 eps_c if a = b = 2 mod 4 and c even, otherwise 0.
# Usage: python3 a5_kernel_selftest.py KERNEL_BINARY.  Output: a5_kernel_selftest.out.
import subprocess, itertools, sys
binary = sys.argv[1]
def K(V, edges, lams):
    inp = f"{V} {len(edges)}\n" + ''.join(f"{a} {b}\n" for a, b in edges) + f"{len(lams)}\n" + ''.join(' '.join(map(str, l)) + '\n' for l in lams)
    return [int(x) for x in subprocess.run([binary], input=inp, capture_output=True, text=True).stdout.split()]
eps = lambda a: (-1) ** (a // 2)
def f8(a, b): return 8 * (eps(a) + eps(b)) - 12 * eps(a) * eps(b) if a % 2 == 0 and b % 2 == 0 else 0
def th(a, b, c):
    if a % 2 == 0 and b % 2 == 0 and c % 2 == 0: return 4 * (eps(a) + eps(b) + eps(c)) - 8 * eps(a) * eps(b) * eps(c)
    for x, y, z in ((a, b, c), (a, c, b), (b, c, a)):
        if x % 2 and y % 2 and z % 2 == 0: return -4 * eps(z)
    return 0
def db(a, b, c): return -32 * eps(c) if a % 4 == 2 and b % 4 == 2 and c % 2 == 0 else 0
R = range(1, 9)
tests = [('figure-eight', 1, [(0, 0), (0, 0)], list(itertools.product(R, R)), lambda l: f8(*l)),
         ('theta', 2, [(0, 1), (0, 1), (0, 1)], list(itertools.product(range(1, 7), repeat=3)), lambda l: th(*l)),
         ('dumbbell', 2, [(0, 0), (1, 1), (0, 1)], list(itertools.product(range(1, 8), range(1, 8), range(1, 7))), lambda l: db(*l))]
for name, V, edges, lams, ref in tests:
    got = K(V, edges, lams)
    bad = [(l, g, ref(l)) for l, g in zip(lams, got) if g != ref(l)]
    print(f"{name}: {len(lams)} length vectors (lengths up to 8/6/7), {sum(1 for g in got if g)} nonzero, mismatches with the THM-4757 table: {len(bad)} {bad[:3]}")

#!/usr/bin/env python3
"""Audit A, misc: THM-4601 (i) lag-one form; 'every K-uniform collision rule of the tiling compiler is a D-chain
absorption at equal Terras time with zero debt' on the compiler's own exponent-1459 rule (collatz_reset2_rules_20261007.md
section 3); class-decided collisions must be equal-time zero-debt (slope argument, checked on the explicit heads)."""
import random
from a2_barrier import make_source, T

def v2(n):
    return (n & -n).bit_length() - 1

def uword(z, L):
    out = []
    for _ in range(L):
        y = 3*z + 1; a = v2(y); out.append(a); z = y >> a
    return out, z

def part_i():
    rnd = random.Random(1)
    n_res = n_reset = 0
    for _ in range(3000):
        K = rnd.randint(3, 60)
        # random source, any first reset letter
        t = rnd.getrandbits(80) | 1
        x = 2*3**(K - 1)*t - 1
        z = 2*3**(K - 2)*t - 1                      # last ones-run point of h_1 = 2^(K-1) t - 1
        assert x - 1 == 3*z + 1
        h1 = (t << (K - 1)) - 1
        w, end = uword(h1, K - 2)
        assert end == z and w == [1]*(K - 2)
        # D = 1 chain (x, z), state (1, 2), equals the lag-one chain (x, x-1) from Terras time 1 on
        a, b, k = x, z, 1
        c, d, k2 = x, x - 1, 0
        for s in range(40):
            if s >= 1:
                assert a == c and b == d and k == k2
            k += (a & 1) - (b & 1); a, b = T(a), T(b)
            k2 += (c & 1) - (d & 1); c, d = T(c), T(d)
        if x % 8 == 5:
            # immediate merge of (x, x-1): absorbed at Terras time 3 with equal odd steps
            c, d, kk = x, x - 1, 0
            for s in range(3):
                kk += (c & 1) - (d & 1); c, d = T(c), T(d)
            assert c == d and kk == 0
            n_reset += 1
        elif x % 8 == 1:
            n_res += 1
        else:
            raise AssertionError("x = 3 or 7 mod 8 impossible")
    print(f"(i) lag-one form: x-1 = 3z+1; D=1 chain = (x, x-1) chain from time 1 on (3000 sources); "
          f"x=5 mod 8 merges at Terras time 3 ({n_reset} cases); residual x=1 mod 8 ({n_res} cases)")

def part_compiler_1459():
    # section 3 of collatz_reset2_rules_20261007.md: u, v1 (D = 1), v2 (D = 2) with F_u(-1) etc.; source 2^1459 - 1
    u = (2, 3, 2, 2, 2, 2, 1, 5, 1, 1, 2)
    v2w = (2, 2, 2, 1, 2, 2, 2, 1, 1, 3, 1, 1, 1)
    K = 1459
    n = (1 << K) - 1; D = 2; h = (1 << (K - D)) - 1
    ws, zs = uword(n, K - 1 + len(u) + 1)
    wc, zc = uword(h, K - 1 - D + len(v2w) + 1)
    assert ws[:K - 1] == [1]*(K - 1) and tuple(ws[K - 1:K - 1 + len(u)]) == u
    assert wc[:K - 1 - D] == [1]*(K - 1 - D) and tuple(wc[K - 1 - D:K - 1 - D + len(v2w)]) == v2w
    assert zs == zc and wc[-1] == ws[-1] + 2
    # D-chain from the run ends at equal Terras time with zero debt
    x = 2*3**(K - 1) - 1; y = (x + 1)//3**D - 1
    a, b, k = x, y, D
    for s in range(sum(u) + 30):
        if a == b and k == 0:
            break
        k += (a & 1) - (b & 1); a, b = T(a), T(b)
    assert a == b and k == 0
    print(f"tiling-compiler rule 2^1459-1 ~> 2^1457-1: actual words match the compiler; the D=2 chain from the run ends "
          f"absorbs at equal Terras time {s} with zero debt (sum u + 1 = {sum(u) + 1})")

if __name__ == "__main__":
    part_i()
    part_compiler_1459()
    print("DONE")

"""Pure-Python check of a rotation-order witness: order x in [0, p) by (x * num) mod den; excess by Legendre symbols (pow).
Usage: python verify_rotation_alpha.py p num den"""
import sys
p, num, den = map(int, sys.argv[1:4])
order = [x for _, x in sorted(((x * num) % den, x) for x in range(p))]
assert sorted(order) == list(range(p)) and len(set((x * num) % den for x in range(p))) == p   # no ties
leg = [0] + [1 if pow(a, (p - 1) // 2, p) == 1 else -1 for a in range(1, p)]
ex = sum(leg[(order[j] - order[i]) % p] for i in range(p) for j in range(i + 1, p))
nat = sum(leg[(j - i) % p] for i in range(p) for j in range(i + 1, p))
print('p=%d alpha=%d/%d: rotation-order excess %d, natural excess %d, forward arcs %d vs %d, beaten: %s'
      % (p, num, den, ex, nat, (p * (p - 1) // 2 + ex) // 2, (p * (p - 1) // 2 + nat) // 2, ex > nat))

# Audit of HYP-9391 (b) / check D: Legendre moments at every odd prime in the auditor's ysums output (both classes
# mod 4), a_l = -L(l) = -sum_w phi(w(w-1)(w-l)), sums over l != 0, 1, PARI coefficients a = 8.4.a.a, c = 8.6.a.a,
# d = 4.6.a.a:
#   (T2) sum phi(l) a_l^2 = -q - 1 - a          (T4) sum phi(l) a_l^4 = -2q^2 - 3q - 1 - 3q a - c - d
#   (U2) sum a_l^2 = q^2 - 2q - 3                (U4) sum a_l^4 = 2q^3 - 4q^2 - 9q - 3 - d
# Also the Sym-power bookkeeping (the H^1_! traces implied by T4, U4):
#   sum phi(l) Sym^4 = T4 - 3q (T2 + q) + 2q^2  should be -(1 + c + d)  (one local invariant at l = 1, H^1_! = c + d)
#   sum Sym^4 = U4 - 3q (U2 - q(q-2)) - 2q^2 (q-2)  should be -(3 + d)   (three local invariants, H^1_! = d)
# Usage: moments_check.py ysums_file forms_file ; writes moments_check.out
import sys, os
HERE = os.path.dirname(os.path.abspath(__file__))
forms = {}
for line in open(sys.argv[2]):
    p, x84, x86, x46, nb = map(int, line.split()); forms[p] = (x84, x86, x46)
res = {'T2': [], 'T4': [], 'U2': [], 'U4': [], 'Sym4t': [], 'Sym4u': []}
qs = []
for line in open(sys.argv[1]):
    r = list(map(int, line.split()))
    q, M2t, M4t, M2u, M4u = r[0], r[7], r[8], r[9], r[10]
    if q < 5: continue
    a, c, d = forms[q]
    qs.append(q)
    res['T2'].append(M2t == -q - 1 - a)
    res['T4'].append(M4t == -2 * q * q - 3 * q - 1 - 3 * q * a - c - d)
    res['U2'].append(M2u == q * q - 2 * q - 3)
    res['U4'].append(M4u == 2 * q ** 3 - 4 * q * q - 9 * q - 3 - d)
    sym4t = M4t - 3 * q * (M2t + q) + 2 * q * q          # sum phi(l) = -1 for l != 0, 1
    sym4u = M4u - 3 * q * (M2u - q * (q - 2)) - 2 * q * q * (q - 2)
    res['Sym4t'].append(sym4t == -(1 + c + d))
    res['Sym4u'].append(sym4u == -(3 + d))
n1 = sum(1 for q in qs if q % 4 == 1); n3 = sum(1 for q in qs if q % 4 == 3)
out = [f"{len(qs)} odd primes 5..{max(qs)} ({n1} = 1 mod 4, {n3} = 3 mod 4)"]
for k, v in res.items():
    bad = [q for q, ok in zip(qs, v) if not ok]
    out.append(f"  {k}: holds at {sum(v)}/{len(v)}{'; fails at ' + str(bad[:10]) if bad else ''}")
print('\n'.join(out))
open(os.path.join(HERE, 'moments_check.out'), 'w').write('\n'.join(out) + '\n')

# taxicab_splitaxis_kpc1_verify.py — ADVERSARIAL VERIFIER (kind-pasteur-2026-06-10-S1)
# Independent check of worker claim A3 (THM-461 split lemma):
#   For gcd(x,y)=1, q = x^2 - xy + y^2 satisfies:
#     (i)   q is odd,
#     (ii)  every prime p != 3 dividing q is p == 1 (mod 3),
#     (iii) v3(q) <= 1, with v3(q) = 1  iff  3 | x+y,
#   and x^3+y^3 = (x+y)*q, q = N(x + y*omega) (Eisenstein norm),
#   x^2-xy+y^2 ~GL2(Z)~ x^2+xy+y^2 (same representation counts),
#   r_E(m) = 6*B(m) with B(m) = sum_{d|m} chi_{-3}(d), so cofactors m on the
#   split axis have r_E(m) = 6*B(m) >= 6 > 0.
#
# FRESH METHOD (different from worker): factor every q via a numpy
# smallest-prime-factor sieve, range 1 <= y <= x <= 2000;
# plus brute-force cross-checks of the norm-form identities and r = 6B.
import math
import numpy as np

X = 2000
LIM = X * X          # q = x^2 - xy + y^2 <= x^2 for 0 <= y <= x

# SPF sieve (numpy)
spf = np.arange(LIM + 1, dtype=np.int32)
for i in range(2, math.isqrt(LIM) + 1):
    if spf[i] == i:
        sl = spf[i * i::i]
        sl[sl == np.arange(i * i, LIM + 1, i, dtype=np.int32)] = i

cnt = 0
viol_odd = viol_split = viol_v3max = viol_v3iff = viol_identity = 0
for x in range(1, X + 1):
    xx = x * x
    x3 = xx * x
    for y in range(1, x + 1):
        if math.gcd(x, y) != 1:
            continue
        cnt += 1
        q = xx - x * y + y * y
        if x3 + y * y * y != (x + y) * q:
            viol_identity += 1
        if q % 2 == 0:
            viol_odd += 1
        v3 = 0
        m = q
        while m % 3 == 0:
            m //= 3
            v3 += 1
        if v3 > 1:
            viol_v3max += 1
        if (v3 == 1) != ((x + y) % 3 == 0):
            viol_v3iff += 1
        # factor 3-free part via SPF; every prime must be == 1 mod 3
        while m > 1:
            p = int(spf[m])
            if p % 3 != 1:
                viol_split += 1
                break
            while m % p == 0:
                m //= p

print(f"[split] primitive pairs 1<=y<=x<={X}: {cnt}")
print(f"[split] identity x^3+y^3=(x+y)q violations: {viol_identity}")
print(f"[split] q even (should be 0): {viol_odd}")
print(f"[split] prime p!=3, p==2 mod 3 dividing q (should be 0): {viol_split}")
print(f"[split] v3(q) >= 2 (should be 0): {viol_v3max}")
print(f"[split] (v3==1) <-> (3|x+y) violations (should be 0): {viol_v3iff}")

# Norm form / GL2-equivalence / r = 6B cross-checks (brute force, small range)
def chi(d):
    r = d % 3
    return 1 if r == 1 else (-1 if r == 2 else 0)

def B(m):
    s = 0
    for d in range(1, m + 1):
        if m % d == 0:
            s += chi(d)
    return s

T = 400
R = math.isqrt(4 * T) + 2
rplus = [0] * (T + 1)   # x^2+xy+y^2
rminus = [0] * (T + 1)  # x^2-xy+y^2
for a in range(-R, R + 1):
    for b in range(-R, R + 1):
        vp = a * a + a * b + b * b
        if 0 < vp <= T:
            rplus[vp] += 1
        vm = a * a - a * b + b * b
        if 0 < vm <= T:
            rminus[vm] += 1
bad_equiv = [m for m in range(1, T + 1) if rplus[m] != rminus[m]]
bad_6B = [m for m in range(1, T + 1) if rplus[m] != 6 * B(m)]
print(f"[forms] m<={T}: #(r(x^2+xy+y^2) != r(x^2-xy+y^2)): {len(bad_equiv)}")
print(f"[forms] m<={T}: #(r != 6*B(m)): {len(bad_6B)}")
# split-axis m (3^{0,1} * split primes) must have B(m) >= 1, i.e. r_E >= 6:
bad_pos = 0
for m in range(1, T + 1):
    mm = m
    v3 = 0
    while mm % 3 == 0:
        mm //= 3
        v3 += 1
    onaxis = (v3 <= 1)
    t = mm
    while t > 1 and onaxis:
        p = t
        for pp in range(2, math.isqrt(t) + 1):
            if t % pp == 0:
                p = pp
                break
        if p % 3 != 1:
            onaxis = False
        while t % p == 0:
            t //= p
    if onaxis and rplus[m] < 6:
        bad_pos += 1
print(f"[forms] split-axis m <= {T} with r_E(m) < 6: {bad_pos}")

ok = (viol_identity == viol_odd == viol_split == viol_v3max == viol_v3iff == 0
      and not bad_equiv and not bad_6B and bad_pos == 0)
print("VERDICT(A3):", "PASS" if ok else "FAIL")

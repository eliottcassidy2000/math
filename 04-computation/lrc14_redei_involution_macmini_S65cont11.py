from fractions import Fraction as F
import random
random.seed(94)
# CORRECTED: the involution is tau <-> 1 - tau (exact for integer speeds).
viol = 0
for _ in range(100000):
    V = random.randrange(20, 300)
    tau = F(random.randrange(1, 4000), random.randrange(2, 4000))
    v = random.randrange(1, V)
    a, b = (v * tau) % 1, (v * (1 - tau)) % 1
    if min(a, 1-a) != min(b, 1-b): viol += 1
print(f"tau <-> 1-tau clearance equality: violations {viol} / 100k")
# REDEI PARITY LAW on symmetric grids: #lonely grid points == [1/2 lonely] (mod 2)
def lonely(S, tau):
    return all(14 * ((v * tau) % 1) >= 1 and 14 * (1 - (v * tau) % 1) >= 1 for v in S)
par_viol = 0; odd_counts = 0; covering_odd = 0
for trial in range(400):
    N = random.randrange(20, 80)                     # symmetric grid k/N, k = 1..N-1
    S = sorted(random.sample(range(1, 100), 6))
    wit = [k for k in range(1, N) if lonely(S, F(k, N))]
    fix = 1 if (N % 2 == 0 and lonely(S, F(1, 2))) else 0
    if (len(wit) - fix) % 2 != 0: par_viol += 1
    if len(wit) % 2 == 1:
        odd_counts += 1
        if all(any(v % q == 0 for v in S) for q in (2,)): covering_odd += 1
print(f"parity law [#witnesses == fixed-layer (mod 2)]: violations {par_viol} / 400")
print(f"odd witness counts: {odd_counts}; of those with an even speed present: {covering_odd}")
# the depth-1 blindness: lonely at 1/2 <=> all speeds odd
blind_viol = 0
for _ in range(20000):
    S = random.sample(range(1, 200), 6)
    if lonely(S, F(1,2)) != all(v % 2 == 1 for v in S): blind_viol += 1
print(f"[lonely at 1/2 <=> all-odd]: violations {blind_viol} / 20k")

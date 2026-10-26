# Natural numbers with their powers below them (n, n^2, n^3, ...) and doublings beside them (n, 2n, 4n, ...).  Built in
# order n = 1, 2, 3, ...: which numbers and chain values already occur in the structure built from smaller n, and which
# chain values are hit later by a larger number's chain.  Writing n = 2^a b (b odd):
#   n is new  <=>  n is odd and not a perfect power (n > 1);
#   n^2 = 2^(2a) b^2 lies in the doubling chain of b^2, which is LATER than n iff b > 2^a; first instance n = 6
#   (6^2 = 36 = 9 * 2^2, in 9's doubling chain).
# Also: preimages of perfect numbers under s, and repunit coincidences s(p^a) = s(q^b).  mac-mini-2026-10-25-pyth.
from sympy import isprime, perfect_power
LIM = 10 ** 6
owner = {}            # value -> (n, kind, exponent) of the first chain containing it
events = []
for n in range(1, 61):
    chains = [(n ** j, 'pow', j) for j in range(1, 30) if n ** j <= LIM and (n > 1 or j == 1)]
    chains += [(n * 2 ** i, 'dbl', i) for i in range(1, 40) if n * 2 ** i <= LIM]
    pre = [(v, kind, e, owner[v]) for v, kind, e in chains if v in owner]
    new_n = n not in owner
    for v, kind, e in chains:
        owner.setdefault(v, (n, kind, e))
    if n <= 40:
        events.append((n, new_n, [(v, f"{kind}^{e}", f"from {o[0]} ({o[1]}^{o[2]})") for v, kind, e, o in pre if v <= 2000][:4]))
print("n | new? | chain values already present (value, how n reaches it, who placed it first)")
for n, new_n, pre in events:
    print(f"  {n:3d} {'NEW ' if new_n else 'old '} {pre}")
# squares landing in a LATER number's doubling chain
later = []
for n in range(2, 200):
    a = (n & -n).bit_length() - 1; b = n >> a
    if a >= 1 and b * b > n: later.append(n)
print("n whose square lies in the doubling chain of a later odd number b^2 (n = 2^a b, b > 2^a):", later[:20])
# new numbers = odd non-perfect-powers
news = [n for n in range(3, 120) if n % 2 and not perfect_power(n)]
print("new numbers (odd, not perfect powers) <= 120:", news)

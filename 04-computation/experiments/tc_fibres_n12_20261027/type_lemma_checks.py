"""type_lemma_checks.py -- finite facts used in the order-12 type lemma of THM-4812 (a) and in (d).
(1) every hyperplane of F_3^4 (and F_3^3 itself) contains a vector of Hamming weight >= 3;
(2) the window computation: f(n) <= tc * max |Aut| bound only for n <= 22 (f = strong Hamiltonian-path minimum)."""
import itertools
vecs4 = list(itertools.product(range(3), repeat=4))
hyper = set()
for c in vecs4:
    if any(c):
        # normalise the functional up to scalar
        lead = next(x for x in c if x)
        cn = tuple((x * pow(lead, -1, 3)) % 3 for x in c)
        hyper.add(cn)
bad = 0
for c in hyper:
    ker = [x for x in vecs4 if sum(a * b for a, b in zip(c, x)) % 3 == 0]
    assert len(ker) == 27
    if not any(sum(1 for t in x if t) >= 3 for x in ker):
        bad += 1
print(f"(1) hyperplanes of F_3^4: {len(hyper)}, without a weight->=3 vector: {bad}; F_3^3 contains (1,1,1)")
def f(n):
    k, r = divmod(n, 3)
    return 3 * 5 ** (k - 1) if r == 0 else (5 ** k if r == 1 else 9 * 5 ** (k - 1))
def B(n, tcmax):
    best = 0
    for p7 in (0, 1):
        for r5 in (0, 1, 2):
            if (9 if p7 else 1) * 3 ** r5 > tcmax:
                continue
            rest = n - 5 * r5 - 7 * p7
            if rest < 0:
                continue
            best = max(best, 3 ** (rest // 3) * 5 ** r5 * 21 ** p7)
    return best
possible = [n for n in range(3, 200) if f(n) <= 7 * B(n, 7) or f(n) <= 21 * B(n, 21)]
print(f"(2) orders n in [3, 200] where f(n) <= tc * B(n) for tc in {{7, 21}}: {possible}")
print(f"    tc = 7 alone: {[n for n in range(3, 200) if f(n) <= 7 * B(n, 7)]}; tc = 21 alone: {[n for n in range(3, 200) if f(n) <= 21 * B(n, 21)]}")
print(f"    n = 22: f = {f(22)}, f/21 = {f(22)/21:.1f}, B(22, 21) = {B(22, 21)} (one P7 and five C3 orbits), B(22, 7) = {B(22, 7)}")

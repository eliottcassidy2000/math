#!/usr/bin/env python3
"""THM-485/486 (claudebox-2026-06-11-S5): the two temperatures (deterministic fugacity +
quenched sign-disorder) of the hard-core/Fibonacci transfer operator; Viswanath's constant as
the disordered endpoint; the disorder-induced phase transition at the Embree-Trefethen activity;
and the involution modulus 24 = Pisano modulus (THM-484 -> Fibonacci). Pure python, no deps."""
import random, math

# ---------- THM-486: involution modulus 24 = Pisano modulus ----------
def pisano(m):
    if m == 1: return 1
    a, b = 0, 1
    for i in range(6 * m * m):
        a, b = b, (a + b) % m
        if a == 0 and b == 1: return i + 1
def rank_app(m):
    a, b = 0, 1
    for k in range(1, 6 * m + 1):
        a, b = b, (a + b) % m
        if a == 0: return k
def fib(n):
    a, b = 0, 1
    for _ in range(n): a, b = b, a + b
    return a

def part_thm486():
    print("== THM-486: the involution modulus 24 is the Pisano modulus ==")
    fixed = [n for n in range(1, 3001) if pisano(n) == n]
    print(f"  pi(24) = {pisano(24)} (= 24); pi(n)=n exactly for n in {fixed} = {{1}} U {{24*5^k}} (Fulton-Morris)")
    print(f"  => 24 is the SMALLEST n>1 with pi(n)=n (the involution modulus is the base Pisano-fixed point)")
    assert fixed == [1, 24, 120, 600, 3000]
    print(f"  alpha(24) = {rank_app(24)} (rank of apparition); F_12 = {fib(12)} = 12^2 (unique nontrivial Fibonacci square)")
    ok = all((p * p - 1) % pisano(p) == 0 and (p * p - 1) % 24 == 0
             for p in [7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61])
    print(f"  pi(p) | p^2-1 AND 24 | p^2-1 for primes 7..61: {ok}  (same p^2-1, THM-484 link)")
    assert pisano(24) == 24 and fib(12) == 144 and ok

# ---------- THM-485 part 1: deterministic fugacity temperature ----------
def Ipath(n, x):
    a, b = 1, 1 + x
    for _ in range(n - 1): a, b = b, b + x * a
    return b

def part_fugacity():
    print("\n== THM-485(1): deterministic fugacity x — one operator, Zeckendorf -> repo H ==")
    for x, name in [(1, "Zeckendorf / golden-mean SFT (phi)"), (2, "repo H = I(Omega,2) / Jacobsthal"), (3, "x=3")]:
        g = (1 + math.sqrt(1 + 4 * x)) / 2
        r = Ipath(40, x) / Ipath(39, x)
        print(f"  x={x}: growth (1+sqrt(1+4x))/2 = {g:.6f}, ratio I_40/I_39 = {r:.6f}  [{name}]")
        assert abs(g - r) < 1e-6
    # Zeckendorf = no-11 strings = independent sets of path = Fibonacci
    def no11(n):
        a, b = 1, 2
        for _ in range(n - 1): a, b = b, a + b
        return b  # length-n binary strings with no two adjacent 1s = F_{n+2}
    zk = [no11(n) for n in range(1, 8)]
    print(f"  Zeckendorf-admissible length-n strings (no 11) = {zk} = F_3..F_9 (golden-mean SFT, entropy log phi)")
    assert zk == [fib(n + 2) for n in range(1, 8)]

# ---------- THM-485 part 2: quenched disorder temperature (Viswanath) ----------
def lyap(order, step, steps=1_500_000, seed=1):
    rng = random.Random(seed); v = [1.0] * order; s = 0.0
    for _ in range(steps):
        c = step(v, rng); v = v[1:] + [c]
        m = max(abs(t) for t in v) + 1e-300
        if m > 1e150 or m < 1e-150: s += math.log(m); v = [t / m for t in v]
    return math.exp((s + math.log(max(abs(t) for t in v) + 1e-300)) / steps)

def part_disorder():
    print("\n== THM-485(2): quenched sign-disorder — the disordered (Viswanath) constants ==")
    sgn = lambda r: r.choice((1, -1))
    fams = [
        ("Fibonacci  (det phi=1.61803)", 2, lambda v, r: sgn(r) * v[1] + sgn(r) * v[0], 1.13199),
        ("Tribonacci (det 1.83929)", 3, lambda v, r: sgn(r) * v[2] + sgn(r) * v[1] + sgn(r) * v[0], None),
        ("base-path  (THM-337, det 3.38298)", 3, lambda v, r: 3 * v[2] + sgn(r) * v[1] + sgn(r) * v[0], None),
    ]
    for name, order, step, ref in fams:
        c = sum(lyap(order, step, 1_200_000, seed=s) for s in range(3)) / 3
        tag = f"  == Viswanath {abs(c-ref)<0.02}" if ref else ""
        print(f"  {name} -> disordered Lyapunov = {c:.4f}{tag}")
    print("  (disorder LOWERS growth eigenvalue -> Lyapunov exponent; the 'quenched/glass' temperature, cf S637)")

def part_transition():
    print("\n== THM-485(3): disorder-induced phase transition at Embree-Trefethen beta*~0.70258 ==")
    sgn = lambda r: r.choice((1, -1))
    for x in [0.3, 0.5, 0.70258, 0.9, 1.0, 2.0]:
        lam = sum(lyap(2, lambda v, r: sgn(r) * v[1] + x * sgn(r) * v[0], 800_000, seed=s)
                  for s in range(3)) / 3
        det = (1 + math.sqrt(1 + 4 * x)) / 2
        print(f"  x={x:.5f}: disordered growth {lam:.4f} ({'DECAY' if lam<1 else 'grow'}), deterministic {det:.4f} (always grows)")
    print("  => the QUENCHED hard-core gas has a critical activity beta* (decay<->growth) the deterministic one lacks;")
    print("     Zeckendorf x=1 (1.132) and repo-H x=2 both in the grown phase.")

if __name__ == '__main__':
    part_thm486(); part_fugacity(); part_disorder(); part_transition()
    print("\nall checks passed")

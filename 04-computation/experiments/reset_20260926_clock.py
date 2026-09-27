"""Exact, clock-correct reset bounds, precision ledger, and ternary roles.

No convergence assumption; all orbit loops have explicit caps. No dependencies.
"""
from collections import Counter
from fractions import Fraction


def need(condition, context):
    if not condition:
        raise RuntimeError(context)


def v2(n):
    need(n > 0, ("v2 domain", n))
    return (n & -n).bit_length() - 1


def v3(n):
    count = 0
    while n % 3 == 0:
        n //= 3
        count += 1
    return count


def U(n):
    return (3*n+1) >> v2(3*n+1)


def value(state):
    q, r, u = state
    return (q << r) + u


def product(state):
    q, r, u = state
    return (q << r)*u


def split(n):
    if n == 1:
        return None
    r = n.bit_length()-1
    return 1, r, n-(1 << r)


def step(state):
    q, r, u = state
    a = v2(3*u+1)
    v = (3*u+1) >> a
    if a < r:
        return (3*q, r-a, v), a, a, 0, "consume"
    if a > r:
        return (v, a-r, 3*q), a, r, 0, "swap"
    b = v2(3*q+v)
    return None, a, r+b, b, "collision"


def main():
    counts = Counter()
    maximum = Fraction(0)
    for q in range(1, 64, 2):
        for u in range(1, 256, 2):
            for r in range(1, 13):
                state = q, r, u
                y = value(state)
                nxt, a, m, b, tag = step(state)
                need(m == v2(3*y+1), ("actual valuation", state))
                counts[tag] += 1
                if nxt:
                    need(value(nxt) == U(y), ("actual value", state))
                    need(2*m == r+a-nxt[1], ("precision ledger", state))
                    q1, _, u1 = nxt
                    need((q1 % 3 == 0) != (u1 % 3 == 0), ("role xor", state))
                    if tag == "consume":
                        need(v3(q1) == v3(q)+1 and v3(u1) == 0, ("consume clock", state))
                    else:
                        need(v3(u1) == v3(q)+1 and v3(q1) == 0, ("swap clock", state))
                    continue
                need(2*m == r+a+2*b, ("collision ledger", state))
                if not y <= 4*u <= 3*y:
                    continue
                counts["balanced"] += 1
                x = U(y)
                # P5 stops when the collision endpoint is1.
                z = U(x) if x > 1 else 1
                fresh = split(z)
                new_k = product(fresh) if fresh else 0
                old_k = product(state)
                need(16*new_k < 27*old_k, ("27/16 reset", state))
                maximum = max(maximum, Fraction(new_k, old_k))
                if x > 1 and v2(3*x+1) >= 2:
                    counts["balanced_extra_ge2"] += 1
                    need(64*new_k < 27*old_k, ("27/64 reset", state))

    # Positive and hostile controls for both sharp constants.
    sharp_one = []
    sharp_two = []
    for k in range(65):
        y = ((1 << (6*k+8))-31)//9
        state = 3*(y-1)//8, 1, (y+3)//4
        need(value(state) == y and all(t & 1 for t in (state[0], state[2])), ("sharp one domain", k))
        need(step(state)[4] == "collision" and v2(3*y+1) == 2, ("sharp one collision", k))
        x = U(y)
        need(v2(3*x+1) == 1 and U(x) == (1 << (6*k+5))-3, ("sharp one clock", k))
        rat = Fraction(product(split(U(x))), product(state))
        need(rat == Fraction(81*y*y+126*y-527, 48*(y-1)*(y+3)), ("sharp one identity", k))
        sharp_one.append(rat)
    for k in range(1, 66):
        z = (1 << (6*k+1))-1
        y = (16*z-7)//9
        state = (3*y-11)//8, 1, (y+11)//4
        need(value(state) == y and all(t & 1 for t in (state[0], state[2])), ("sharp two domain", k))
        need(step(state)[4] == "collision" and v2(3*y+1) == 2, ("sharp two collision", k))
        x = U(y)
        need(v2(3*x+1) == 2 and U(x) == z, ("sharp two clock", k))
        need(y <= 4*state[2] <= 3*y, ("sharp two balance", k))
        sharp_two.append(Fraction(product(split(z)), product(state)))
    need(all(a < b < Fraction(27,16) for a,b in zip(sharp_one, sharp_one[1:])), "sharp one monotone controls")
    need(all(a < b < Fraction(27,64) for a,b in zip(sharp_two, sharp_two[1:])), "sharp two monotone controls")

    # Preserve one actual source through every consume/swap of a phase.
    phases = Counter()
    for source in range(3, 4096, 2):
        state = split(source)
        initial_r = state[1]
        sum_a = sum_m = 0
        run = 0
        for j in range(512):
            q, r, u = state
            need(v3(q) == run, ("actual consume run", source, j))
            nxt, a, m, b, tag = step(state)
            sum_a += a
            sum_m += m
            if nxt is None:
                need(2*sum_m == initial_r+sum_a+2*b, ("phase ledger", source))
                phases["collision"] += 1
                break
            run = run+1 if tag == "consume" else 0
            need(2*sum_m == initial_r+sum_a-nxt[1], ("open ledger", source))
            state = nxt
        else:
            phases["capped"] += 1

    # Fix arbitrarily many ternary digits and still regenerate unbounded precision.
    for depth in range(1, 7):
        period = 3**depth
        reference = ((4**4-1)//3) % period
        for j in range(9):
            k = 4+j*period
            u = (4**k-1)//3
            nxt, _, _, _, tag = step((1,1,u))
            need(u % period == reference and v3(u) == 0, ("fixed ternary guard", depth, j))
            need(tag == "swap" and nxt == (1,2*k-1,3), ("unbounded precision", depth, j))

    print("exact universe q odd1..63,u odd1..255,r1..12:", dict(sorted(counts.items())))
    print("largest scanned P5 balanced product ratio:", maximum)
    print("sharp27/16 family k0..64; first ratio:", sharp_one[0])
    print("sharp27/64 family k1..65; first ratio:", sharp_two[0])
    print("actual canonical-source phases odd3..4095,cap512:", dict(sorted(phases.items())))
    print("fixed ternary depths1..6,nine members each: same residues and v3(core)=0, unbounded-symbolic precision")
    print("PASS: all arithmetic uses integers/Fraction; no convergence extrapolation")


if __name__ == "__main__":
    main()

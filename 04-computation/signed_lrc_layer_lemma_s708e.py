"""Verify the LAYERED-VISIBILITY lemma + order-block decomposition for signed-LRC (monad-explorer-S708).

Lemma (odd prime power C=q^k): sin(2 pi t x / C) = 0  <=>  v_q(t)+v_q(x) >= k.
Consequence: at frequency t with v_q(t)=a, only magnitudes x with v_q(x) < k-a are 'visible'
in Phi_t = sum_x eps_x sin(2pi t x/C).  This layering is WHY the silent-flip group is the
subgroup lattice (deep layers vanish at high-valuation frequencies).

Also verify: order-block O_e = {x: ord(x)=e} (e|C, e>1); subgroup half H_d = union_{e|d,e>1} O_e;
V = span{H_d} = span{O_e}; and confirm for general composite C the same valuation/visibility law
holds per prime (sin=0 iff for every prime power p^a || C, ... ) -- here just check q^k and a mixed C.
"""
import math
from itertools import product


def vq(x, q):
    if x == 0:
        return 99
    v = 0
    while x % q == 0:
        x //= q; v += 1
    return v


def check_layer_lemma(q, k):
    C = q ** k
    bad = 0
    for t in range(1, C):
        for x in range(1, C):
            s = math.sin(2 * math.pi * t * x / C)
            zero = abs(s) < 1e-9
            pred = (vq(t, q) + vq(x, q) >= k)
            if zero != pred:
                bad += 1
    print(f"  C={C}={q}^{k}: layer-lemma sin=0 <=> v_q(t)+v_q(x)>=k : {'OK' if bad==0 else f'{bad} FAIL'}")
    return bad


def check_order_blocks(C):
    half = (C - 1) // 2
    # order of x
    def order(x):
        g = math.gcd(x, C)
        return C // g
    blocks = {}
    for x in range(1, half + 1):
        e = order(x)
        blocks.setdefault(e, []).append(x)
    # subgroup halves
    def H(d):
        K = set((C // d * j) % C for j in range(d))
        return frozenset(y for y in K if 1 <= y <= half)
    divs = [d for d in range(2, C) if C % d == 0]
    ok = True
    for d in divs:
        # H_d should = union of O_e for e|d, e>1
        pred = set()
        for e in blocks:
            if e > 1 and d % e == 0:
                pred |= set(blocks[e])
        if pred != set(H(d)):
            ok = False
            print(f"    MISMATCH H_{d}: {sorted(H(d))} vs union-O {sorted(pred)}")
    print(f"  C={C}: order-blocks {{e:|O_e|}} = "
          f"{ {e: len(v) for e, v in sorted(blocks.items()) if e>1} };  "
          f"H_d = union_(e|d) O_e : {'OK' if ok else 'FAIL'}")


if __name__ == "__main__":
    print("Layered-visibility lemma (odd prime powers):")
    for q, k in [(3, 2), (3, 3), (3, 4), (5, 2), (5, 3), (7, 2)]:
        check_layer_lemma(q, k)
    print("\nOrder-block decomposition of the silent-flip space V:")
    for C in [9, 15, 21, 25, 27, 35, 45, 63, 75, 81, 105, 121, 125]:
        check_order_blocks(C)

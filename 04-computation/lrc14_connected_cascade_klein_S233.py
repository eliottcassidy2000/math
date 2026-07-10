#!/usr/bin/env python3
"""
THE CONNECTED CASCADE (klein-2026-07-09-S233, HYP-5880) -- corrected objects.

CORRECTION TO THM-684(I) ESTABLISHED HERE: the orthogonality box object of the
character-tuple layer sum Sum_{prod chi = chi_0} prod c-hat(chi_l) chi_l(u_l)
is NOT the product count M_t = #{y in B^t : prod y = prod u}; it is the
COMMON-MULTIPLIER (partial live) count

    A_t(U) = #{c in (Z/q)^x : c*u mod q in B for every u in U},

with the exact identity  Full(U) := Sum_{prod chi = chi_0} prod c-hat chi(u)
= A_t(U)/(q-1)  (c-hat(chi) = (1/(q-1)) Sum_{y in B} conj(chi)(y)).
Consequences: A_2({u1,u2}) = N_{u1/u2} (THM-683's ratio object, verbatim), and
A_13(full support) = LM(q) -- the box counts ARE the partial live counts.

Connected (pure/cumulant) counts by Mobius over sub-supports
(Pure(empty)=1, A_0 := q-1, A_1 = b exactly):

    PureN(U) = Sum_{V subset U} (-b/Q)^{|U|-|V|} A_{|V|}(V)/Q ,   Q = q-1.

Integer forms:  Pure2 = Q*A2 - b^2          (dev2 := Pure2/Q  = centered N_w)
                Pure3 = Q^2*A3 - b*Q*Sum_{3 pairs}A2 + 2 b^3   (dev3 := Pure3/Q^2)

Layer ledger:  LM/Q = Sum_t (b/Q)^{13-t} Sum_{|U|=t} PureN(U);
layer_t = (b/Q)^{13-t} * (1/Q) * Sum_{supports} dev_t;  budget = (b/Q)^13.

Parts: (1) small-prime direct character-sum convention lock (A_t identity
PASS, product-count form FAIL, Mobius = direct pure sums); (2) exact full
sweep (78 pairs + 286 triples, bitmask popcounts) at primes to 5003 for the
generic and near-dilation instances: dev scales, CS/energy domination,
absolute + signed layer totals, the measured q0.
"""
import cmath
from itertools import combinations
from math import comb

GEN = [12, 33, 46, 47, 68, 73, 79, 81, 85, 87, 91, 112, 120]
DIL = [20, 41, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260]


def band(q):
    return (q + 13) // 14, (13 * q) // 14  # [lo, hi] closed


# ---------------------------------------------------------------- Part 1
def convention_lock(q=61, g=2):
    lo, hi = band(q)
    B = list(range(lo, hi + 1))
    b = len(B)
    Q = q - 1
    dlog = {}
    x = 1
    for k in range(Q):
        dlog[x] = k
        x = x * g % q
    assert len(dlog) == Q, "g not primitive"
    ch = lambda a, y: cmath.exp(-2j * cmath.pi * a * dlog[y] / Q)  # conj(chi_a)(y)
    chat = [sum(ch(a, y) for y in B) / Q for a in range(Q)]
    chi = lambda a, u: cmath.exp(2j * cmath.pi * a * dlog[u % q] / Q)

    def A_t(us):
        return sum(1 for c in range(1, q) if all(lo <= c * u % q <= hi for u in us))

    def M_t(us):  # the (incorrect-for-the-identity) product box count
        T = 1
        for u in us:
            T = T * u % q
        if len(us) == 2:
            return sum(1 for y in B if T * pow(y, -1, q) % q in set(B))
        return sum(1 for y1 in B for y2 in B
                   if T * pow(y1 * y2 % q, -1, q) % q in set(B))

    print(f"PART 1 -- convention lock at q={q} (b={b}); direct character sums:")
    for us in [(12, 33, 47), (46, 7, 26)]:
        # t=2 on the first two coords
        u1, u2 = us[0], us[1]
        full2 = sum(chat[a] * chat[(-a) % Q] * chi(a, u1) * chi((-a) % Q, u2)
                    for a in range(Q))
        a2 = A_t([u1, u2])
        m2 = M_t([u1, u2])
        ok = abs(full2 - a2 / Q) < 1e-9
        bad = abs(full2 - m2 / Q) > 1e-6
        print(f"  t=2 {u1, u2}: Full={full2.real:+.6f}  A2/Q={a2 / Q:+.6f} "
              f"[{'PASS' if ok else 'FAIL'}]   M2/Q={m2 / Q:+.6f} "
              f"[product form {'differs, as claimed' if bad else 'UNEXPECTEDLY equal'}]")
        # t=3
        full3 = 0
        pure3d = 0
        for a1 in range(Q):
            for a2_ in range(Q):
                a3 = (-a1 - a2_) % Q
                term = (chat[a1] * chat[a2_] * chat[a3]
                        * chi(a1, us[0]) * chi(a2_, us[1]) * chi(a3, us[2]))
                full3 += term
                if a1 and a2_ and a3:
                    pure3d += term
        a3c = A_t(list(us))
        pairs = [A_t([us[i], us[j]]) for i, j in combinations(range(3), 2)]
        pure3m = (Q * Q * a3c - b * Q * sum(pairs) + 2 * b ** 3) / Q ** 3
        ok3 = abs(full3 - a3c / Q) < 1e-9
        okm = abs(pure3d - pure3m) < 1e-9
        print(f"  t=3 {us}: Full={full3.real:+.6f}  A3/Q={a3c / Q:+.6f} "
              f"[{'PASS' if ok3 else 'FAIL'}]   pure direct={pure3d.real:+.6f} "
              f"Mobius={pure3m:+.6f} [{'PASS' if okm else 'FAIL'}]")
    # the A_13 = LM closure, exact
    vs = [v % q for v in GEN]
    lm = sum(1 for c in range(1, q) if all(lo <= c * v % q <= hi for v in vs))
    print(f"  A_13(full GEN support) = {lm} = LM({q}) by definition "
          f"(the box counts ARE the partial live counts)")


# ---------------------------------------------------------------- Part 2
def sweep(S, name, q):
    lo, hi = band(q)
    b = hi - lo + 1
    Q = q - 1
    r = b / Q
    v = [x % q for x in S]
    assert all(v), (name, q)
    bm = []
    for vi in v:  # bitmask over c in [1,q): bit c set iff c*vi in band
        m = 0
        for c in range(1, q):
            if lo <= c * vi % q <= hi:
                m |= 1 << c
        bm.append(m)
    A2 = {}
    dev2 = {}
    for i, j in combinations(range(13), 2):
        a = (bm[i] & bm[j]).bit_count()
        A2[i, j] = a
        dev2[i, j] = (Q * a - b * b) / Q
    dev3 = {}
    for i, j, k in combinations(range(13), 3):
        a3 = (bm[i] & bm[j] & bm[k]).bit_count()
        p3 = Q * Q * a3 - b * Q * (A2[i, j] + A2[i, k] + A2[j, k]) + 2 * b ** 3
        dev3[i, j, k] = p3 / Q / Q
    d2 = sorted(abs(x) for x in dev2.values())
    d3 = sorted(abs(x) for x in dev3.values())
    budget = r ** 13
    lay2a = r ** 11 / Q * sum(abs(x) for x in dev2.values())
    lay3a = r ** 10 / Q * sum(abs(x) for x in dev3.values())
    lay2s = r ** 11 / Q * sum(dev2.values())
    lay3s = r ** 10 / Q * sum(dev3.values())
    e2 = sum(x * x for x in dev2.values())
    csC = d3[-1] / e2 ** 0.5
    print(f"{name} q={q:5d}: |dev2| max={d2[-1]:8.2f} med={d2[39]:7.2f} | "
          f"|dev3| max={d3[-1]:9.2f} med={d3[143]:8.2f} | "
          f"L2abs={lay2a:.4f} L3abs={lay3a:.4f} (bud {budget:.4f}) | "
          f"L2sig={lay2s:+.4f} L3sig={lay3s:+.4f} | max|dev3|/sqrtE2={csC:.2f}")
    return dict(q=q, name=name, d2max=d2[-1], d2med=d2[39], d3max=d3[-1],
                d3med=d3[143], lay2a=lay2a, lay3a=lay3a, lay2s=lay2s,
                lay3s=lay3s, budget=budget, csC=csC)


def main():
    convention_lock()
    print("\nPART 2 -- exact sweep: all 78 pairs + 286 triples "
          "(dev2 = centered ratio pair; dev3 = connected triple, layer scale)")
    primes = [139, 239, 353, 383, 547, 1009, 2003, 3001, 5003]
    rows = []
    for name, S in (("GEN", GEN), ("DIL", DIL)):
        for q in primes:
            rows.append(sweep(S, name, q))
        print()
    print("PART 3 -- scaling + the measured q0 (t<=3 truncation):")
    for name in ("GEN", "DIL"):
        rs = [r for r in rows if r["name"] == name]
        import math
        g2 = math.log(rs[-1]["d2max"] / rs[0]["d2max"]) / math.log(rs[-1]["q"] / rs[0]["q"])
        g3 = math.log(rs[-1]["d3max"] / rs[0]["d3max"]) / math.log(rs[-1]["q"] / rs[0]["q"])
        print(f"  {name}: max|dev2| ~ q^{g2:.2f}, max|dev3| ~ q^{g3:.2f}")
        for frac, tag in ((1.0, "q0(abs<=budget)"), (0.5, "q0(abs<=budget/2)")):
            qs = [r["q"] for r in rs if r["lay2a"] + r["lay3a"] <= frac * r["budget"]]
            print(f"    {tag}: {'q=' + str(min(qs)) if qs else 'NOT REACHED by q=5003'}"
                  f"   [L2abs+L3abs vs {frac}*(b/Q)^13]")
        w = rs[-1]
        print(f"    at q={w['q']}: abs total={w['lay2a'] + w['lay3a']:.4f} "
              f"signed total={w['lay2s'] + w['lay3s']:+.4f} "
              f"(abs/|signed| = {(w['lay2a'] + w['lay3a']) / abs(w['lay2s'] + w['lay3s']):.1f}x)")


if __name__ == "__main__":
    main()

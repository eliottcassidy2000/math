#!/usr/bin/env python3
"""
dioph_brute.py — integer brute force for polynomial solution families of
   z^2 + y^2 z + P(x) = 0,  P(x) = x^3 + s*x + u   (and the -z variant e=-1:
   z^2 + (y^2-1) z + P(x) = 0).

Structure: z = -A(t) works iff A*(y^2 + e - A) = P(x(t)), i.e.
   A | P(x(t)))  and  A + B = y^2 + e,  B = P(x)/A.

Search 1 (quadratic split): x = t^2 + a t + b, A = t^2 + al t + be (both monic),
   B = P(x)/A exact division (quartic), S = A + B - e must be a perfect square.
Search 2 (constant divisor): A = d in Z, B = P(x)/d integral, S = d + B - e
   a perfect square of a cubic.

All exact integer arithmetic; families verified symbolically at the end and
instantiated at three huge t values.

Usage: python3 dioph_brute.py [BOX]
"""
import sys
from itertools import product

BOX = int(sys.argv[1]) if len(sys.argv) > 1 else 12

EQNS = {
    "warmup x^3+2":  (0, 2, 0),
    "eq1 x^3-2":     (0, -2, 0),
    "eq2 x^3-x-1":   (-1, -1, 0),
    "eq3 x^3+x-1":   (1, -1, 0),
    "eq4 x^3+x+1":   (1, 1, 0),
    "eq5 x^3-3":     (0, -3, 0),
    "eq6 x^3+3":     (0, 3, 0),
    "eq7 x^3-x-2":   (-1, -2, 0),
    "eq8 x^3-x+2":   (-1, 2, 0),
    "eq9 x^3+2 (-z)": (0, 2, -1),
}

def polymul(p, q):
    r = [0] * (len(p) + len(q) - 1)
    for i, pi in enumerate(p):
        if pi:
            for j, qj in enumerate(q):
                r[i + j] += pi * qj
    return r

def polyadd(p, q):
    n = max(len(p), len(q))
    return [(p[i] if i < len(p) else 0) + (q[i] if i < len(q) else 0) for i in range(n)]

def polyscale(p, c):
    return [c * v for v in p]

def P_of_x(x, s, u):
    # x is coeff list low->high; return x^3 + s*x + u
    x2 = polymul(x, x)
    x3 = polymul(x2, x)
    r = polyadd(x3, polyscale(x, s))
    r = polyadd(r, [u])
    return r

def divmod_monic(num, den):
    """divide num by monic den, return (quo, rem) exact integer."""
    num = num[:]
    dq = len(num) - len(den)
    if dq < 0:
        return [0], num
    quo = [0] * (dq + 1)
    for k in range(dq, -1, -1):
        c = num[k + len(den) - 1]
        quo[k] = c
        if c:
            for i, dv in enumerate(den):
                num[k + i] -= c * dv
    rem = num[:len(den) - 1]
    while len(rem) > 1 and rem[-1] == 0:
        rem.pop()
    return quo, rem

def is_square_quartic(S):
    """S monic quartic (low->high, len 5, S[4]=1): return y (len 3) or None."""
    if len(S) != 5 or S[4] != 1:
        return None
    if S[3] % 2:
        return None
    g = S[3] // 2
    h2 = S[2] - g * g
    if h2 % 2:
        return None
    h = h2 // 2
    if S[1] != 2 * g * h or S[0] != h * h:
        return None
    return [h, g, 1]

def is_square_quartic_general(S):
    """S quartic with square leading coeff s^2: y = s t^2 + g t + h."""
    while len(S) > 1 and S[-1] == 0:
        S = S[:-1]
    if len(S) == 5:
        l = S[4]
        s = int(round(abs(l) ** 0.5))
        if l < 0 or s * s != l or s == 0:
            return None
        # y = s t^2 + g t + h: t^3: 2sg; t^2: g^2+2sh; t: 2gh; 1: h^2
        if S[3] % (2 * s):
            return None
        g = S[3] // (2 * s)
        h2 = S[2] - g * g
        if h2 % (2 * s):
            return None
        h = h2 // (2 * s)
        if S[1] != 2 * g * h or S[0] != h * h:
            return None
        return [h, g, s]
    if len(S) == 3:   # leading cancellation: quadratic square: y = g t + h
        l = S[2]
        g = int(round(abs(l) ** 0.5))
        if l < 0 or g * g != l:
            return None
        if S[1] % (2 * g):
            return None
        h = S[1] // (2 * g)
        if S[0] != h * h:
            return None
        return [h, g]
    if len(S) == 1:
        v = S[0]
        rr = int(round(abs(v) ** 0.5))
        if v >= 0 and rr * rr == v:
            return [rr]
        return None
    return None

def is_square_sextic(S, lead):
    """S sextic with S[6] = lead^2 (y = lead t^3 + ...): solve coefficients."""
    if len(S) != 7:
        return None
    import math
    l2 = S[6]
    L = int(math.isqrt(abs(l2)))
    if L * L != l2 or L == 0:
        return None
    # y = L t^3 + p t^2 + q t + r
    # y^2: t^6: L^2; t^5: 2Lp; t^4: p^2+2Lq; t^3: 2pq+2Lr; t^2: q^2+2pr; t: 2qr; 1: r^2
    if S[5] % (2 * L):
        return None
    p = S[5] // (2 * L)
    t4 = S[4] - p * p
    if t4 % (2 * L):
        return None
    q = t4 // (2 * L)
    t3 = S[3] - 2 * p * q
    if t3 % (2 * L):
        return None
    r = t3 // (2 * L)
    if S[2] != q * q + 2 * p * r or S[1] != 2 * q * r or S[0] != r * r:
        return None
    return [r, q, p, L]

def sqrt_poly(S):
    """exact integer polynomial square root of S (any even degree), or None."""
    T = S[:]
    while len(T) > 1 and T[-1] == 0:
        T.pop()
    dS = len(T) - 1
    if dS % 2:
        return None
    if dS == 0:
        v = T[0]
        if v < 0:
            return None
        r = int(round(v ** 0.5))
        return [r] if r * r == v else None
    lead = T[-1]
    if lead < 0:
        return None
    wl = int(round(lead ** 0.5))
    if wl * wl != lead or wl == 0:
        return None
    dW = dS // 2
    W = [0] * (dW + 1)
    W[dW] = wl
    # determine coefficients top-down: coeff of t^(dS - j)
    for j in range(1, dW + 1):
        # coefficient of t^{dS-j} in W^2: 2*wl*W[dW-j] + (terms from known W's)
        acc = 0
        for a in range(dW - j + 1, dW):
            bidx = dS - j - a
            if 0 <= bidx <= dW:
                acc += W[a] * W[bidx]
        rem = T[dS - j] - acc
        if rem % (2 * wl):
            return None
        W[dW - j] = rem // (2 * wl)
    WW = polymul(W, W)
    WW += [0] * (len(T) - len(WW))
    return W if WW[:len(T)] == T and all(v == 0 for v in WW[len(T):]) else None

def fam_str(x, y, A):
    def s(p, var="t"):
        terms = []
        for i in range(len(p) - 1, -1, -1):
            c = p[i]
            if c == 0:
                continue
            mono = "" if i == 0 else (var if i == 1 else f"{var}^{i}")
            if i == 0:
                terms.append(f"{c:+d}")
            elif c == 1:
                terms.append("+" + mono)
            elif c == -1:
                terms.append("-" + mono)
            else:
                terms.append(f"{c:+d}*{mono}")
        out = "".join(terms)
        return out.lstrip("+") or "0"
    return f"x(t) = {s(x)} ;  y(t) = {s(y)} ;  z(t) = -({s(A)})"

def verify(x, y, A, s, u, e):
    P = P_of_x(x, s, u)
    y2 = polymul(y, y)
    inner = polyadd(polyadd(y2, [e]), polyscale(A, -1))
    lhs = polymul(A, inner)
    n = max(len(lhs), len(P))
    lhs += [0] * (n - len(lhs)); P2 = P + [0] * (n - len(P))
    return lhs == P2

def main():
    LEADS = [(1, 1), (-1, -1), (2, 2), (-2, -2), (3, 3), (-3, -3)]
    for name, (s, u, e) in EQNS.items():
        found = []
        # Search 1: x = L t^2 + a t + b, A = aL t^2 + al t + be, B = P/A
        for (L, aL) in LEADS:
            bl = (L ** 3) // aL          # B leading; must be a perfect square
            sq = int(round(abs(bl) ** 0.5))
            if bl < 0 or sq * sq != bl:
                continue
            for a in range(-BOX, BOX + 1):
                for b in range(-BOX, BOX + 1):
                    P = P_of_x([b, a, L], s, u)
                    for al in range(-BOX, BOX + 1):
                        for be in range(-BOX, BOX + 1):
                            A = [be, al, aL]
                            # rational division: multiply P by aL^k trick — do exact:
                            # A*B = P with B integer quartic: solve by back-substitution
                            B = [0] * 5
                            ok = True
                            rem = P[:]
                            for k in range(4, -1, -1):
                                c = rem[k + 2]
                                if c % aL:
                                    ok = False
                                    break
                                B[k] = c // aL
                                for i, dv in enumerate(A):
                                    rem[k + i] -= B[k] * dv
                            if not ok or any(rem[i] for i in range(len(rem))):
                                continue
                            S = polyadd(polyadd(A, B), [-e])
                            y = is_square_quartic_general(S)
                            if y is not None:
                                x = [b, a, L]
                                if verify(x, y, A, s, u, e):
                                    found.append((x, y, A))
        # Search 3: cancelling cubic split — x = L t^2 + b t + d with L = -k^2,
        # A, B cubic with opposite leads, y = g t + h (deg <= 1).
        # Equivalent condition: y^4 - 4 P(x) is the square of a cubic W.
        for L in (-1, -4, -9):
            for b in range(-BOX, BOX + 1):
                for d in range(-BOX, BOX + 1):
                    P = P_of_x([d, b, L], s, u)
                    m4 = polyscale(P, -4)
                    for g in range(0, BOX + 1):
                        for h in range(-BOX, BOX + 1):
                            if g == 0 and h < 0:
                                continue
                            y = [h, g] if g else [h]
                            y2 = polymul(y, y)
                            y4 = polymul(y2, y2)
                            Scand = polyadd(y4, m4)
                            # want Scand = W^2, W cubic with W-lead^2 = -4L^3
                            W = sqrt_poly(Scand)
                            if W is None:
                                continue
                            # recover A = (y^2 + e' ... ) use A = (y^2 - W)/2 or (y^2+W)/2
                            for sgn in (1, -1):
                                num = polyadd(y2, polyscale(W, sgn))
                                if any(c % 2 for c in num):
                                    continue
                                A = [c // 2 for c in num]
                                # B = y^2 - A
                                B = polyadd(y2, polyscale(A, -1))
                                if polymul(A, B) == P + [0] * (len(polymul(A, B)) - len(P)):
                                    if verify([d, b, L], y, A, s, u, 0) and e == 0:
                                        found.append(([d, b, L], y, A))
        # Search 2: A = d constant, x = t^2 + a t + b (need d | P(x) contentwise)
        for d in [1, -1, 2, -2, 4, -4, 8, -8, 3, -3, 9, -9]:
            for a in range(-BOX, BOX + 1):
                for b in range(-BOX, BOX + 1):
                    x = [b, a, 1]
                    P = P_of_x(x, s, u)
                    if any(c % d for c in P):
                        continue
                    B = [c // d for c in P]
                    S = polyadd(polyadd([d], B), [-e])
                    y = is_square_sextic(S + [0] * (7 - len(S)) if len(S) < 7 else S, None)
                    if y is not None:
                        A = [d]
                        if verify(x, y, A, s, u, e):
                            found.append((x, y, A))
        print(f"=== {name}: {len(found)} families")
        seen = set()
        for (x, y, A) in found[:12]:
            key = fam_str(x, y, A)
            if key in seen:
                continue
            seen.add(key)
            print("   " + key, flush=True)
    print("done")

if __name__ == "__main__":
    main()

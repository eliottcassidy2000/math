#!/usr/bin/env python3
"""Independent audit of THM-4849 (D) and (E) and of HYP-9446 (1), (2).

Everything is computed from the definitions:
  - reduced denominators d(wv) = D(wv)/gcd(C(wv), D(wv)) from exact integers;
  - "cancels p^e" = p^e | gcd(C(w'v), D(wv));
  - "locked at p under v" = every marking w' of the cycle cancels p (HYP-9446's definition), with p | D(wv)
    and p not dividing 6(2^v - 3).
No use is made of the theorem being audited when computing the left-hand sides.
"""
import random
from fractions import Fraction as Fr
from math import gcd
import sympy
from sympy import factorint, isprime, primerange

def C_of(w, a=3):
    k = len(w); s = 0; t = 0
    for j in range(k):
        t += a ** (k - 1 - j) * 2 ** s
        s += w[j]
    return t

def D_of(w, a=3):
    return 2 ** sum(w) - a ** len(w)

def vp(n, p):
    if n == 0:
        return 10 ** 9
    c = 0
    while n % p == 0:
        n //= p
        c += 1
    return c

def vp_frac(x, p):
    if x == 0:
        return 10 ** 9
    return vp(x.numerator, p) - vp(x.denominator, p)

def ordmod(a, n):
    if gcd(a, n) != 1:
        return None
    from sympy.ntheory import n_order
    return n_order(a, n)

def rots(w):
    return [tuple(w[i:] + w[:i]) for i in range(len(w))]

out = []
def say(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    out.append(s)

def check(c, m):
    if not c:
        say("FAIL:", m)
        raise SystemExit(1)

# ---------------------------------------------------------------------------------------------------------
def check_identity_symbolic():
    H, v, k, x, D, C, Dv = sympy.symbols('H v k x D C Dv')
    # with C = x*D(w), D(w) = 2^H - 3^k, D(wv) = 2^(H+v) - 3^(k+1)
    twoH, threek, twov = sympy.symbols('twoH threek twov')
    Dw = twoH - threek
    Dwv = twoH * twov - 3 * threek
    Cw = x * Dw
    lhs = 3 * Cw + twoH
    rhs = twoH * (1 - (twov - 3) * x) + x * Dwv
    check(sympy.expand(lhs - rhs) == 0, "identity (D)")
    say("(D) identity C(wv) = 2^H(1 - (2^v-3)x) + x D(wv) with C(wv) = 3C(w) + 2^H: verified symbolically")

def check_attachment_law(ntests=3000, seed=7):
    rnd = random.Random(seed)
    tested = 0
    vacuous_case = 0
    for _ in range(ntests):
        k = rnd.randint(1, 6)
        w = tuple(rnd.randint(1, 7) for _ in range(k))
        v = rnd.randint(1, 9)
        Dw, Cw = D_of(w), C_of(w)
        if Dw == 0:
            continue
        wv = w + (v,)
        Dwv, Cwv = D_of(wv), C_of(wv)
        check(Cwv == 3 * Cw + 2 ** sum(w), "C(wv) = 3C(w) + 2^H")
        x = Fr(Cw, Dw)
        if Dwv == 0:
            continue
        for p, e in factorint(abs(Dwv)).items():
            if p == 2 or Dw % p == 0:
                # hypothesis p does not divide 2 D(w)
                continue
            red = Fr(Cwv, Dwv).denominator
            lhs = vp(red, p)
            if (2 ** v - 3) % p == 0:
                vacuous_case += 1  # predicted impossible
                continue
            xi = Fr(1, 2 ** v - 3)
            rhs = e - min(e, vp_frac(x - xi, p))
            check(lhs == rhs, f"attachment law w={w} v={v} p={p}")
            tested += 1
    say(f"(D) attachment law v_p(d(wv)) = e - min(e, v_p(x - xi_v)) on {tested} (word, letter, prime) tests: OK; "
        f"cases with p | D(wv), p not | D(w), p | 2^v-3 encountered: {vacuous_case} (the case is vacuous)")
    check(vacuous_case == 0, "vacuous case occurred")

def check_vacuity_argument():
    # p | 2^v - 3 and p | D(wv) imply p | D(w): exhaustive small check
    cnt = 0
    for k in range(1, 6):
        for H in range(k, 14):
            # words of length k, sum H
            def comps(n, parts):
                if parts == 1:
                    yield (n,)
                    return
                for a in range(1, n - parts + 2):
                    for rest in comps(n - a, parts - 1):
                        yield (a,) + rest
            for w in comps(H, k):
                for v in range(1, 12):
                    m = 2 ** v - 3
                    Dwv = D_of(w + (v,))
                    g = gcd(abs(m), abs(Dwv))
                    if g > 1:
                        check(D_of(w) % g == 0 or all(D_of(w) % p == 0 for p in factorint(g)), "vacuity")
                        cnt += 1
    say(f"(D) vacuity: whenever a prime p divides both 2^v - 3 and D(wv), it divides D(w) ({cnt} instances checked)")

# ---------------------------------------------------------------------------------------------------------
def locking_definition(w, v, p, e):
    """(1): every point of the cycle of w is p-integral and = xi_v mod p^e"""
    Dw = D_of(w)
    pts = [Fr(C_of(r), Dw) for r in rots(w)]
    if any(vp_frac(x, p) < 0 for x in pts):
        return False
    xi = Fr(1, 2 ** v - 3)
    return all(vp_frac(x - xi, p) >= e for x in pts)

def locking_letters(w, v, p, e):
    L = ordmod(2, p ** e)
    return all((vi - v) % L == 0 for vi in w)

def cancels_all(w, v, p, e):
    Dwv = D_of(w + (v,))
    return all(gcd(C_of(r + (v,)), Dwv) % p ** e == 0 for r in rots(w))

def check_locking(ntests=12000, seed=11):
    rnd = random.Random(seed)
    stats = {"tests": 0, "locked": 0, "higher_e": 0, "locked_higher_e": 0}
    edge = {"p|2^v-3": 0, "p|2^v-3, (2) true, (1) false": 0}
    while stats["tests"] < ntests:
        k = rnd.randint(1, 7)
        # bias toward congruent letters so that locked cases occur
        L = rnd.choice([3, 4, 5, 6, 8, 10, 11, 12, 18, 20])
        c = rnd.randint(1, L)
        if rnd.random() < 0.5:
            w = tuple(c + L * rnd.randint(0, 2) for _ in range(k))
            v = c + L * rnd.randint(0, 2)
        else:
            w = tuple(rnd.randint(1, 9) for _ in range(k))
            v = rnd.randint(1, 12)
        if sum(w) + v > 62:
            continue
        Dwv = D_of(w + (v,))
        if Dwv == 0 or D_of(w) == 0:
            continue
        if abs(Dwv) > 10 ** 19:
            continue
        for p, emax in factorint(abs(Dwv)).items():
            if p in (2, 3):
                continue
            for e in range(1, emax + 1):
                if (2 ** v - 3) % p == 0:
                    edge["p|2^v-3"] += 1
                    if locking_letters(w, v, p, e) and not locking_definition(w, v, p, e):
                        edge["p|2^v-3, (2) true, (1) false"] += 1
                    continue
                a = locking_definition(w, v, p, e)
                b = locking_letters(w, v, p, e)
                cc = cancels_all(w, v, p, e)
                check(a == b, f"locking iff fails: w={w} v={v} p^e={p}^{e}")
                check(a == cc, f"(1) vs all markings cancel p^e: w={w} v={v} p^e={p}^{e}")
                stats["tests"] += 1
                if e >= 2:
                    stats["higher_e"] += 1
                if a:
                    stats["locked"] += 1
                    if e >= 2:
                        stats["locked_higher_e"] += 1
    say(f"(E) locking theorem (1) <=> (2) <=> 'every marking cancels p^e' on {stats['tests']} (word, letter, p^e) tests "
        f"({stats['higher_e']} with e >= 2): OK; locked: {stats['locked']} ({stats['locked_higher_e']} with e >= 2)")
    say(f"(E) edge cases excluded by the hypothesis p not | 2^v-3: {edge}")

def edge_example():
    w, v = (3,), 3
    say(f"(E) the hypothesis p not | 2^v - 3 is needed: w = {w}, v = {v}, p = 5: D(wv) = {D_of(w + (v,))}, "
        f"letters = v mod ord_5(2) = 4 holds, but the cycle {{1/5}} is not 5-integral (locking (1) false)")

# ---------------------------------------------------------------------------------------------------------
def interval_table(qmax=100):
    rows = []
    for q in primerange(3, qmax):
        k = (q - 1) // 2
        m = k + 1
        w = (1,) * (k - 1) + (k + 2,)
        Dw1 = D_of(w + (1,))
        # locked primes by definition: p | D(w.1), p not | 6, every marking cancels p
        g = 0
        for r in rots(w):
            g = gcd(g, C_of(r + (1,)))
        g = gcd(g, Dw1)
        locked = sorted(p for p in factorint(g) if p not in (2, 3)) if g > 1 else []
        G = gcd(2 ** m - 1, 3 ** m - 1)
        gp = sorted(factorint(G)) if G > 1 else []
        check(locked == gp, f"interval corollary at q={q}")
        rows.append((q, m, locked, G))
    say("(E) interval corollary: locked primes of the R_q interval cycle under (1) (by definition) = primes of "
        "gcd(2^m-1, 3^m-1), m=(q+1)/2, for every odd prime q < 100:")
    for q, m, locked, G in rows:
        say(f"   q={q:3d} m={m:3d} locked {locked} gcd {G}")

# ---------------------------------------------------------------------------------------------------------
def connection_sets(q):
    k = (q - 1) // 2
    for mask in range(2 ** k):
        yield frozenset((x if (mask >> (x - 1)) & 1 == 0 else q - x) for x in range(1, k + 1))

def canon(S, q):
    return min(tuple(sorted((x + t) % q for x in S)) for t in range(q))

def gaps(S, q):
    s = sorted(S)
    return tuple(((s[(i + 1) % len(s)] - s[i]) % q) or q for i in range(len(s)))

def hyp1_table(qs=(5, 7, 11, 13, 17, 19, 23)):
    say("HYP-9446 (1): locked tournament cycles (definition: some p, some v with p | D(wv), p not | 6(2^v-3), "
        "every marking cancels p):")
    for q in qs:
        k = (q - 1) // 2
        I = frozenset(range(1, k + 1))
        rot = {frozenset((m * x) % q for x in I) for m in range(1, q)}
        cyc = {}
        for S in connection_sets(q):
            cyc.setdefault(canon(S, q), []).append(S)
        nlocked = 0
        nonrot = 0
        examples = []
        for neck, Ss in cyc.items():
            w = gaps(set(Ss[0]), q)
            Cs = [C_of(r) for r in rots(w)]
            # candidate primes: p | C(w'v) for all w' implies p | 3(C(w') - C(w'')); p != 3
            dg = 0
            for c in Cs[1:]:
                dg = gcd(dg, c - Cs[0])
            cands = [p for p in factorint(dg) if p not in (2, 3)] if dg > 1 else []
            hits = []
            for p in cands:
                L = ordmod(2, p)
                for v in range(1, L + 1):
                    if (2 ** v - 3) % p == 0:
                        continue
                    if D_of(w + (v,)) % p != 0:
                        continue
                    if all((3 * c + 2 ** q) % p == 0 for c in Cs):
                        hits.append((p, v))
            if hits:
                nlocked += 1
                isrot = any(S in rot for S in Ss)
                if not isrot:
                    nonrot += 1
                if len(examples) < 8:
                    examples.append((w, hits, "rot" if isrot else "NONROT"))
        say(f"   q={q}: {len(cyc)} tournament cycles, {nlocked} locked, of which not rotational: {nonrot}; examples {examples}")

def hyp1_example():
    w = (1, 1, 1, 1, 1, 1, 1, 4, 1, 1, 10)
    q = 23
    # recover a connection set from the word: positions s_0..s_(k-1)
    pos = [0]
    for v in w[:-1]:
        pos.append(pos[-1] + v)
    P = set(pos)
    trans = [t for t in range(q) if (lambda S: 0 not in S and all(((x in S) != ((-x) % q in S)) for x in range(1, q)))({(x + t) % q for x in P})]
    k = 11
    I = frozenset(range(1, k + 1))
    rot = {frozenset((m * x) % q for x in I) for m in range(1, q)}
    Ss = [frozenset((x + t) % q for x in P) for t in trans]
    Dw1 = D_of(w + (1,))
    allc = all(gcd(C_of(r + (1,)), Dw1) % 7 == 0 for r in rots(w))
    say(f"HYP-9446 (1) example {w}: base points {len(trans)} (connection set {sorted(Ss[0]) if Ss else None}), "
        f"rotational {any(S in rot for S in Ss)}, letters mod 3 {sorted(set(x % 3 for x in w))}, "
        f"7 | D(w.1): {Dw1 % 7 == 0}, every marking cancels 7 under (1): {allc}")

def hyp2():
    qs = list(primerange(3, 6000))
    good = [q for q in qs if gcd(2 ** ((q + 1) // 2) - 1, 3 ** ((q + 1) // 2) - 1) == 1]
    say(f"HYP-9446 (2): odd primes q < 6000: {len(qs)}; with gcd(2^m-1,3^m-1) = 1: {len(good)}; first 15: {good[:15]}")
    # m odd implies gcd = 1 ?
    bad = [q for q in qs if ((q + 1) // 2) % 2 == 1 and gcd(2 ** ((q + 1) // 2) - 1, 3 ** ((q + 1) // 2) - 1) > 1]
    say(f"   first q with m odd and gcd > 1: {bad[:3]}; gcd there = {gcd(2 ** ((bad[0] + 1) // 2) - 1, 3 ** ((bad[0] + 1) // 2) - 1)}")
    # congruence conditions: p locks R_q (interval cycle under (1)) iff p | gcd
    for p, mod, res in ((5, 8, 7), (7, 12, 11), (13, 24, 23), (11, 20, 19)):
        L = ordmod(2, p) * ordmod(3, p) // gcd(ordmod(2, p), ordmod(3, p))
        lhs = [q for q in qs if gcd(2 ** ((q + 1) // 2) - 1, 3 ** ((q + 1) // 2) - 1) % p == 0]
        rhs = [q for q in qs if q % mod == res]
        say(f"   p={p}: lcm(ord_p 2, ord_p 3) = {L}; p | gcd  <=>  q = {res} mod {mod}: {lhs == rhs} "
            f"(first: {lhs[:6]})")
        check(lhs == rhs, "congruence")

def main():
    check_identity_symbolic()
    check_vacuity_argument()
    check_attachment_law()
    check_locking()
    edge_example()
    interval_table()
    hyp1_table()
    hyp1_example()
    hyp2()

if __name__ == "__main__":
    main()

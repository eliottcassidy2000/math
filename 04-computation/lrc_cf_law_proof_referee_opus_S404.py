"""
opus-2026-07-19-S404: referee for THM-1291 (the CF active-leg law, proved).

Checks, all exact:
(R1) LEMMA/(A) referee, 200k random (a,q,mtilde): u* := min{u>=1 : delta(u) <= mtilde}
     is ALWAYS a convergent denominator of a/q (standard CF).  delta(u) = q||u a/q||.
(R2) Prop 0 referee on random families: at any maximizer, ALL active speeds lie in a
     single +- residue class mod q.
(R3) Corpus verdicts: for the S403 corpus + controls: u*, whether u* is in V
     (hypothesis H), smallest active speed v0, and the law chain
     [H => v0 = u* = convergent, all straddling sums = kq, D = k*delta(u*)].
     Expected: H holds on all near-floor families; fails exactly on
     mixed-control (u*=2 not in V) and primes13 (u*=1 not in V).
(R4) (B)-violation hunt: random families, collect cases u* not in V; confirm (A)
     still holds in every case (u* convergent) -- H is necessary, (A) never fails.
"""
import random
from math import gcd
from fractions import Fraction
random.seed(1291)

def convergent_denoms(a, q):
    num, den = q, a
    quots = []
    while den:
        quots.append(num // den); num, den = den, num % den
    h0, h1 = 1, quots[0]; k0, k1 = 0, 1
    dens = [1, h1] if False else None
    # convergents of q/a are h/k -> of a/q are k/h; denominators of a/q convergents = h's
    Hs = [1]
    hprev, h = 1, quots[0]
    kprev, k = 0, 1
    Hs = [h] if False else None
    # redo cleanly: convergents p/qd of a/q: p=k, qd=h
    hs = [hprev, h]  # h_{-1}=1? standard: h_{-1}=1,h_0=quots[0]
    for c in quots[1:]:
        hprev, h = h, c*h + hprev
        hs.append(h)
    return set([1] + hs)  # include q0=1; hs are denominators q1,q2,...

def delta(u, a, q):
    r = (u*a) % q
    return min(r, q-r)

def ustar(a, q, mt):
    for u in range(1, q+1):
        if delta(u, a, q) <= mt:
            return u
    return q

def exact_max(V, qmax):
    bg, bq, wit = 0, 1, None
    for q in range(2, qmax+1):
        for aa in range(1, q):
            if gcd(aa, q) != 1: continue
            m = q
            for v in V:
                r = (v*aa) % q
                r = min(r, q-r)
                if r < m:
                    m = r
                    if m*bq < bg*q: break
            if m*bq > bg*q:
                bg, bq, wit = m, q, (aa, q)
    return bg, bq, wit

# (R1)
bad = 0
for _ in range(200000):
    q = random.randint(3, 3000)
    a = random.randint(1, q-1)
    if gcd(a, q) != 1: continue
    mt = random.randint(1, q//2)
    u = ustar(a, q, mt)
    if u not in convergent_denoms(a, q):
        bad += 1
        if bad < 5: print("R1 FAIL:", a, q, mt, u, sorted(convergent_denoms(a,q)))
print(f"(R1) u* is a convergent denominator: failures {bad}/~200k")

# (R2)+(R4)
r2bad, r4cases, r4fail = 0, 0, 0
for _ in range(400):
    n = random.randint(4, 13)
    V = random.sample(range(1, 120), n)
    bg, bq, (aa, q) = *exact_max(V, 120)[:2], exact_max(V, 120)[2]
    if bg == 0: continue
    act = [v for v in V if delta(v % q if v % q else q, aa, q) == bg and v % q]
    act = [v for v in V if (v % q) and delta(v % q, aa, q) == bg]
    if not act: continue
    r0 = act[0] % q
    if not all((v % q) in (r0, (q - r0) % q) for v in act):
        r2bad += 1
        if r2bad < 4: print("R2 FAIL:", V, aa, q, act)
    u = ustar(aa, q, bg)
    if u not in V:
        r4cases += 1
        if u not in convergent_denoms(aa, q):
            r4fail += 1
print(f"(R2) active class is single +- residue class: failures {r2bad}/400")
print(f"(R4) H-violations found: {r4cases}; (A)-failures among them: {r4fail}")

# (R3) corpus
FAMS = [
 ("AP13", list(range(1,14)), 300), ("GW", list(range(1,12))+[13,24], 300),
 ("{1..12,14}", list(range(1,13))+[14], 300),
 ("3/41", list(range(1,12))+[13,36], 300),
 ("ladder4", list(range(1,12))+[13,48], 300),
 ("ladder5", list(range(1,12))+[13,60], 300),
 ("loose {1..12,15}", list(range(1,13))+[15], 300),
 ("primes13", [2,3,5,7,11,13,17,19,23,29,31,37,41], 300),
 ("mixed", [3,7,10,14,22,25,31,38,44,52,57,63,70], 300),
 ("N19", list(range(1,18))+[19,54], 300),
 ("N31", list(range(1,30))+[31,120], 300),
]
print("\n(R3) corpus:")
for name, V, qm in FAMS:
    bg, q, (aa, _) = *exact_max(V, qm)[:2], exact_max(V, qm)[2]
    u = ustar(aa, q, bg)
    conv = u in convergent_denoms(aa, q)
    H = u in V
    act = [v for v in V if (v % q) and delta(v % q, aa, q) == bg]
    v0 = min(act) if act else None
    lawstr = ""
    if H:
        above = [v for v in act if (v*aa) % q == bg]
        below = [v for v in act if (v*aa) % q == q - bg]
        pairs = [(vi, vj) for vi in above for vj in below]
        sums_ok = all((vi+vj) % q == 0 for vi, vj in pairs)
        Ds = sorted(set(Fraction(bg, q)*(vi+vj) for vi, vj in pairs))
        lawstr = f"; straddle sums all = kq: {sums_ok}; D set {Ds} (delta(u*)={delta(u,aa,q)})"
    print(f"  {name}: t*={aa}/{q}, mt={bg}, u*={u} conv:{conv} H(u* in V):{H} v0={v0}{lawstr}")

# Verifier (lane architecture): independent exact check of P1 (multiplicative reduction with atoms 49, 63, 343).
# Algorithm differs from the lane's forward sieve: backward recursion on divisors.  member[n] = 1 iff n = 1 or some
# generator g divides n with member[n/g] = 1.  Generators that can divide n are read off the factorisation of n
# (smallest-prime-factor sieve).  Also tests the session-context list (49, 63, every 7p with p != 3, p = 7 allowed)
# and the necessity of each of 49, 63, 343 separately.
LIM = 2_000_000
spf = list(range(LIM + 1))
i = 2
while i * i <= LIM:
    if spf[i] == i:
        for j in range(i * i, LIM + 1, i):
            if spf[j] == j:
                spf[j] = i
    i += 1

def prime_factors(n):
    ps = []
    while n > 1:
        p = spf[n]
        ps.append(p)
        while n % p == 0:
            n //= p
    return ps

def misses(atoms, allow_7p_with_p7=False):
    member = bytearray(LIM + 1)
    member[1] = 1
    out = []
    for n in range(3, LIM + 1, 2):
        ps = prime_factors(n)
        ok = False
        for p in ps:                      # (a) odd primes != 7
            if p != 7 and member[n // p]:
                ok = True; break
        if not ok:
            for a in atoms:               # (b) atoms
                if n % a == 0 and member[n // a]:
                    ok = True; break
        if not ok and n % 7 == 0:         # (c) 7p, p prime, p != 3 (and p != 7 unless allowed)
            for p in ps:
                if p == 3: continue
                if p == 7 and not allow_7p_with_p7: continue
                g = 7 * p
                if n % g == 0 and member[n // g]:
                    ok = True; break
        if ok:
            member[n] = 1
        else:
            out.append(n)
    return out

m_full = misses([49, 63, 343])
print("atoms 49, 63, 343 : misses up to", LIM, "=", m_full)
m_no343 = misses([49, 63])
print("atoms 49, 63      : misses =", m_no343)
m_ctx = misses([49, 63], allow_7p_with_p7=True)
print("session-context list (49, 63, 7p with p != 3, p = 7 allowed): misses =", m_ctx)
m_no49 = misses([63, 343])
print("atoms 63, 343     : number of misses =", len(m_no49), "first", m_no49[:10])
m_no63 = misses([49, 343])
print("atoms 49, 343     : number of misses =", len(m_no63), "first", m_no63[:10])
pred = {7, 21}
a = 3
while 7 ** a <= LIM:
    pred.add(7 ** a)
    if 3 * 7 ** a <= LIM:
        pred.add(3 * 7 ** a)
    a += 2
print("prediction {7, 21} + {7^a, 3*7^a : a odd >= 3} for the no-343 misses:", set(m_no343) == pred)
print("CHECK", "PASS" if m_full == [7, 21] and set(m_no343) == pred and m_no343 == m_ctx else "FAIL")

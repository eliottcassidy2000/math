# S94 density lane. Exact checks of L1 (readout form), L2 (affine fibre lemma), shifted collision equation.
import itertools, sys, random
def q(h):
    r = 0
    # q(2h)=q(h), q(2h+1)=q(h)+ceil(h/2)
    bits = []
    while h:
        bits.append(h & 1); h >>= 1
    # rebuild from top
    v = 0; r = 0
    for b in reversed(bits):
        if b: r = r + (v + 1)//2
        v = 2*v + b
    return r
def eta(h): return 2*h + 1 - 2*q(h)
def rev(x, w):
    r = 0
    for i in range(w):
        if (x >> i) & 1: r |= 1 << (w-1-i)
    return r
def wgt(d): return 1 if d == 1 else 2**(d-2)
ok = True
# (0) quadratic form q(h) = sum_{i<j} h_i h_j wgt(j-i)
bad = 0
for h in range(1 << 12):
    b = [(h >> i) & 1 for i in range(12)]
    s = sum(b[i]*b[j]*wgt(j-i) for i in range(12) for j in range(i+1, 12))
    if s != q(h): bad += 1
print("Q0 quadratic form of q on h < 2^12: violations", bad); ok &= (bad == 0)
# L1 readout form
bad = 0; cnt = 0
for h1 in range(0, 64):
    x1, y1 = eta(h1), 2*h1+1
    for s in range(0, 6):
        A, B = x1 + (1 << s)*y1, x1
        for h2 in range(0, 200):
            N = x1*eta(h2) + (1 << s)*y1*(2*h2+1)
            if N != A*(2*h2+1) - 2*B*q(h2): bad += 1
            cnt += 1
print("L1 readout identity: instances", cnt, "violations", bad); ok &= (bad == 0)
# L2 affine fibre lemma, all pairs of w-bit words, w <= 9
for w in range(1, 10):
    bad = 0; cnt = 0
    for x in range(1 << w):
        for xp in range(1 << w):
            eps = [((x >> i) & 1) - ((xp >> i) & 1) for i in range(w)]
            xplus = sum(1 << i for i in range(w) if eps[i] == 1)
            xminus = sum(1 << i for i in range(w) if eps[i] == -1)
            q0 = q(xplus) - q(xminus)
            lin = 0
            for i in range(w):
                if eps[i] == 0 and ((x >> i) & 1):
                    lin += sum(eps[j]*wgt(abs(i-j)) for j in range(w) if j != i)
            if q(x) - q(xp) != q0 + lin: bad += 1
            cnt += 1
    print("L2 affine fibre lemma w=%d: pairs %d violations %d" % (w, cnt, bad)); ok &= (bad == 0)
# shifted collision equation for the readout N_sigma(h) = eta(h) + 2^sigma (2h+1):
# (c + 2^(s+sigma+1)) X - P Xi = gamma, gamma = 2(q(x)-q(x')) + (P mod 2) eps_top + (c mod 2) eps_bottom
for n in (12, 14):
    for sigma in (1, 2, 3, 5):
        d = {}
        for h in range(1 << (n-1), 1 << n):
            d.setdefault(eta(h) + (1 << sigma)*(2*h+1), []).append(h)
        npairs = 0; bad = 0; maxdiff = 0
        for N, L in d.items():
            for a in range(len(L)):
                for b in range(len(L)):
                    if a == b: continue
                    h, hp = L[a], L[b]
                    dd = h ^ hp
                    top = dd.bit_length() - 1
                    s = (dd & -dd).bit_length() - 1
                    w = top - s + 1
                    tp = n - 1 - top
                    P = h >> (top + 1)
                    Z = h & ((1 << s) - 1)
                    x = (h >> s) & ((1 << w) - 1); xp = (hp >> s) & ((1 << w) - 1)
                    c = (1 << (s+1)) - rev(Z, s)
                    X = x - xp; Xi = rev(x, w) - rev(xp, w)
                    gam = 2*(q(x) - q(xp)) + (P & 1)*(((x >> (w-1)) & 1) - ((xp >> (w-1)) & 1)) + (c & 1)*((x & 1) - (xp & 1))
                    if (c + (1 << (s + sigma + 1)))*X - P*Xi != gam: bad += 1
                    if abs(gam) > (1 << w): bad += 1
                    npairs += 1
                    maxdiff = max(maxdiff, abs(h - hp))
        print("shifted equation n=%d sigma=%d: ordered colliding pairs %d violations %d max|h-h'| %d (bound 2^(n-1-sigma)=%d)"
              % (n, sigma, npairs, bad, maxdiff, 1 << (n-1-sigma)))
        ok &= (bad == 0) and (maxdiff < (1 << (n-1-sigma)) or npairs == 0)
print("ALL CHECKS PASS" if ok else "SOME CHECK FAILED")

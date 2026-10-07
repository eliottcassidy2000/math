import sys, re
sys.path.insert(0, '.')
from fractions import Fraction as Fr
exec(open('t4593.py').read().split('# ---------- sanity')[0])   # import Tmap, step
base = '/Users/e/Documents/GitHub/math-wt-chessboard-20261006/04-computation/experiments/golden_20261007_readers/core/A/'
for fn, Dmax in [('mers_D7_witness.txt', 7), ('mers_D61_witness.txt', 61)]:
    txt = open(base + fn).read()
    m = re.search(r'by Terras time (\d+).*bits: ([01]+)', txt, re.S)
    t0 = int(m.group(1)); bits = [int(c) for c in m.group(2)]
    sts = [(-D, Fr(1, 3**D) - 1) for D in range(1, Dmax + 1, 2)]
    absorbed = False; first_all = None
    for t in range(t0):
        sts = [step(s, bits[t]) for s in sts]
        if any(s == (0, Fr(0)) for s in sts): absorbed = True
        if first_all is None and len(set(sts)) == 1: first_all = t + 1
    print(fn, 'partners', len(sts), 'bits', len(bits), 'forced prefix', bits[:4], 'none absorbed:', not absorbed,
          'all equal at t0:', len(set(sts)) == 1, 'first all-equal time:', first_all, 'state k =', sts[0][0])
    # direct check: build a 2-adic p (mod 2^(t0+200)) with these first t0 bits by inverse Terras, then compare actual orbits
    # p with given parity vector: solve bit by bit
    P = 0; mod = 1
    # construct residue rho mod 2^t0 having parity vector bits[:t0]
    rho = 0
    for j in range(t0):
        # choose bit j of rho so that parity at step j matches
        for cand in (rho, rho + (1 << j)):
            x = cand; ok = True
            for t in range(j + 1):
                if x % 2 != bits[t]: ok = False; break
                x = (x // 2) if x % 2 == 0 else (3 * x + 1) // 2
            if ok: rho = cand; break
        else:
            raise SystemExit('no residue?')
    import random
    p = rho + (random.getrandbits(80) << t0)
    qs = []
    for D in range(1, Dmax + 1, 2):
        # q_D = 3^{-D} p + 3^{-D} - 1 must be a 2-adic integer: work mod 2^(t0+80) via inverse of 3^D
        M = 1 << (t0 + 80)
        inv = pow(3**D, -1, M)
        qs.append(((p + 1) * inv - 1) % M)
    # orbits mod 2^(t0+80) are exact for t0 steps if we keep track (each step loses at most 1 bit)
    vals = [p % (1 << (t0 + 80))] + qs
    prec = t0 + 80
    for t in range(t0):
        vals = [(v // 2) if v % 2 == 0 else (3 * v + 1) // 2 for v in vals]
        prec -= 1
        vals = [v % (1 << prec) for v in vals]
    print('   direct 2-adic orbits (mod 2^80 after t0 steps): all partners equal:', len(set(vals[1:])) == 1, ' driver differs:', vals[0] != vals[1])

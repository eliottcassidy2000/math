# -*- coding: utf-8 -*-
# kps-2026-07-11-S127: WHICH INVARIANT governs the seven-sector moment extremality -- E2 or E3?
#
# The wide-direction residue (THM-701/703/opus-S220) = consec MINIMIZES the empty-count functional
# F2 := 4 m1 - m2 = E[N(5-N)]  (degree-2; k=9,10) and a degree-3 analog (k=8).  Phi <= 1 - (1/6)F2.
# mac-mini's stated tool = ADDITIVE-ENERGY (E2 = #{a+b=c+d}) bricks. But opus HYP-5683 proved that for
# the LRC DENSITY FLOOR the governing invariant is SCHUR TRIPLES (E3 = #{a+b=c}), NOT E2 -- because
# loneliness is scale- but NOT translation-invariant, and only E3 shares that symmetry (E2 is translation-
# invariant, blind to the tight-AP vs translated-AP distinction).  meas(S7)/Phi has the SAME symmetry
# (scale-inv, NOT translation-inv; THM-536).  SYMMETRY-MATCH PREDICTION: F2 should track E3, not E2.
#
# TEST: over a battery of 0-anchored 8-cores spanning the (E2,E3) plane -- including opus's separator
# shape (translated AP: max E2, low E3) -- correlate F2 (and cap-Phi) against E2 vs E3.
import random
from itertools import product

def sector(y): return min(6, int((y % 1.0) * 7))
def moments(E, N=96001):
    # N(x) = # empty inner sectors {1..6}; returns (Phi, m1=E[N], m2=E[N(N-1)], F2=4m1-m2)
    s1 = s2 = p0 = p1 = 0
    for i in range(N):
        x = (i + 0.5) / N
        cov = set(sector(e * x) for e in E) & {1,2,3,4,5,6}
        n = 6 - len(cov)
        s1 += n; s2 += n*(n-1)
        if n == 0: p0 += 1
        elif n == 1: p1 += 1
    m1 = s1/N; m2 = s2/N
    return (p0/N + (p1/N)/3, m1, m2, 4*m1 - m2)

def E2(F):
    S = list(F); from_ = {}
    for a in S:
        for b in S: from_[a+b] = from_.get(a+b,0)+1
    return sum(v*v for v in from_.values())
def E3(F):
    Sset = set(F); c = 0
    for a in F:
        for b in F:
            if a+b in Sset: c += 1
    return c   # includes free triples a+0=a when 0 in F

def pear(xs, ys):
    n=len(xs); mx=sum(xs)/n; my=sum(ys)/n
    sx=sum((x-mx)**2 for x in xs)**.5; sy=sum((y-my)**2 for y in ys)**.5
    return sum((x-mx)*(y-my) for x,y in zip(xs,ys))/(sx*sy) if sx*sy>0 else 0

def main():
    cap9 = 0.49426
    consec = list(range(8))
    battery = {'consec {0..7}': consec}
    # opus's separator family: translated APs (max E2, LOW E3) -- scale-inv kills raw translate, so use
    # 0-anchored near-translates and sum-free-flavored 0-cores
    battery['dilate {0,2,..,14}'] = [2*i for i in range(8)]          # = consec by scale-inv (control)
    battery['transl-AP {0,3,5,7,9,11,13,15}'] = [0,3,5,7,9,11,13,15] # 0 + odd block: few Schur triples
    battery['odd+0 {0,1,3,5,7,9,11,13}'] = [0,1,3,5,7,9,11,13]
    battery['sidon {0,1,3,7,12,20,30,44}'] = [0,1,3,7,12,20,30,44]
    battery['geom {0,1,2,4,8,16,32,64}'] = [0,1,2,4,8,16,32,64]
    random.seed(23)
    rows=[]
    for name,E in battery.items():
        phi,m1,m2,f2 = moments(E)
        rows.append((name,E2(E),E3(E),phi,f2))
    for _ in range(30):
        E=[0]+sorted(random.sample(range(1,33),7))
        phi,m1,m2,f2=moments(E); rows.append(('rand',E2(E),E3(E),phi,f2))
    print(' set                                E2     E3   Phi     F2=4m1-m2   cap9-Phi')
    for name,e2,e3,phi,f2 in rows:
        if name!='rand':
            print(f'  {name:34s} {e2:5d} {e3:5d}  {phi:.4f}   {f2:.4f}    {cap9-phi:+.4f}')
    e2s=[r[1] for r in rows]; e3s=[r[2] for r in rows]; f2s=[r[4] for r in rows]; phis=[r[3] for r in rows]
    print()
    print(f'  corr(F2, E2) = {pear(f2s,e2s):+.3f}     corr(F2, E3) = {pear(f2s,e3s):+.3f}')
    print(f'  corr(Phi,E2) = {pear(phis,e2s):+.3f}     corr(Phi,E3) = {pear(phis,e3s):+.3f}')
    print('  (consec MINIMIZES F2 / MAXIMIZES Phi; the governing invariant is the one F2 anti-correlates with strongest)')
    # the separator verdict: does max-E2-low-E3 shape look tight (like consec) or loose?
    cons=[r for r in rows if r[0].startswith('consec')][0]
    sep =[r for r in rows if 'transl-AP' in r[0]][0]
    print(f'\n  SEPARATOR: consec E2={cons[1]} E3={cons[2]} F2={cons[4]:.3f} ; transl-AP E2={sep[1]} E3={sep[2]} F2={sep[4]:.3f}')
    print(f'  transl-AP has {"HIGHER" if sep[4]>cons[4] else "LOWER"} F2 than consec => it is {"LOOSER (E3 governs)" if sep[4]>cons[4] else "as tight (E2 governs)"}')

if __name__ == '__main__':
    main()

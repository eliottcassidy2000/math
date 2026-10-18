"""Audit: evaluate m_N = 1/t at the CM points printed in Beukers arXiv:2403.03301v3, App. C 9.11 (s7) and 9.14
(Apery2 (x) central), with the referee's eta code from a3_cm_classification.py (t = h/(1 + a h + b h^2), 1/t = m_N)."""
import sys, os
import mpmath as mp
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from a3_cm_classification import m_value
mp.mp.dps = 60
s = mp.sqrt
I = mp.mpc(0, 1)
pts7 = {  # tau : Beukers' printed t
    '(1+sqrt(-1/7))/2': ((1 + I * s(mp.mpf(1) / 7)) / 2, '-1'),
    '(2+sqrt(-3))/7': ((2 + I * s(3)) / 7, '1/2'),
    '(3+sqrt(-19))/14': ((3 + I * s(19)) / 14, '1/8'),
    '(1+sqrt(-5/7))/2': ((1 + I * s(mp.mpf(5) / 7)) / 2, '-1/8'),
    '(1+3sqrt(-3))/14': ((1 + 3 * I * s(3)) / 14, '1/24'),
    'sqrt(-1/7)': (I * s(mp.mpf(1) / 7), '1/27'),
    '(1+sqrt(-13/7))/2': ((1 + I * s(mp.mpf(13) / 7)) / 2, '-1/64'),
    'sqrt(-4/7)': (I * s(mp.mpf(4) / 7), '1/125'),
    '(1+sqrt(-61/7))/2': ((1 + I * s(mp.mpf(61) / 7)) / 2, '-1/10648'),
}
pts5 = {
    '(3+sqrt(-11))/10': ((3 + I * s(11)) / 10, '1/4'),
    '(1+sqrt(-3/5))/2': ((1 + I * s(mp.mpf(3) / 5)) / 2, '-1/3'),
    '(1+sqrt(-7/5))/2': ((1 + I * s(mp.mpf(7) / 5)) / 2, '-1/28'),
    '(1+sqrt(-19))/10': ((1 + I * s(19)) / 10, '1/36'),
    '(-1+2i)/5': ((-1 + 2 * I) / 5, '1/18'),
    'sqrt(-2/5)': (I * s(mp.mpf(2) / 5), '1/72'),
    'sqrt(-3/5)': (I * s(mp.mpf(3) / 5), '1/588'),
    '(1+sqrt(-23/5))/2': ((1 + I * s(mp.mpf(23) / 5)) / 2, '-1/828'),
    '(1+sqrt(-47/5))/2': ((1 + I * s(mp.mpf(47) / 5)) / 2, '-1/15228'),
}
for N, pts in ((7, pts7), (5, pts5)):
    for lab, (tau, printed) in pts.items():
        W = -1 / (N * tau)
        z = tau if tau.imag >= W.imag else W
        mv = m_value(z, N)
        tv = 1 / mv
        num, _, den = printed.partition('/')
        pv = mp.mpf(int(num)) / mp.mpf(int(den or 1))
        flag = 'OK' if abs(tv - pv) < 1e-30 else f'MISMATCH (printed {printed}, computed t = 1/{mp.nstr(mv.real, 12)})'
        print(f'N={N} tau={lab:20s} m_N = {mp.nstr(mv, 15):>28}  -> {flag}')

#!/usr/bin/env python3
"""The universal 1/7 buffering law referee (HYP-6840, opus-S301).
buffered(a<-b) = #{exits u of a in one window : b's OPEN bad set contains u}
             = w_a/7 + O(gcd(w_a,w_b))   -- fraction 1/7 for EVERY partner.
(Corrects the session's first w_a/(7d) guess: the step-14*w_b AP takes w_a/d
values in the window, each with multiplicity d; the d cancels.)"""
from fractions import Fraction as F
from math import gcd
import random
random.seed(7)
def circ(x):
    fx = x - (x.numerator // x.denominator); return min(fx, 1-fx)
c = 7; worst = 0
for trial in range(40):
    wa = random.randint(20, 400); wb = random.randint(20, 400)
    if wa % 7 == 0 or wb % 7 == 0 or wa == wb: continue
    d = gcd(wa, wb)
    buf = sum(1 for j in range(wa) if circ(F(wb)*((F(c,14)+c*j)/wa)/c) < F(1,14))
    dev = abs(buf - wa/7)
    worst = max(worst, dev - d)
    print(f"wa={wa:4} wb={wb:4} d={d:3} buffered={buf:4}  wa/7={wa/7:7.2f}  dev={dev:5.2f}")
print(f"max (deviation - gcd) over battery: {worst:.2f}  (law: dev = O(gcd))")

"""
CUSP FORM of M on speed-space. Integer points are downward SPIKES (resonance alignment with the units).
(A) the AP spike: M({1..12,x}) vs x -- the cusp shape.
(B) cusp DEPTHS at integer sets: is the AP the deepest? which sets spike?
(C) decoupling cusp: multiple speeds -> inf, the (6/7)^k factorization of L.
"""
import numpy as np
def M_grid(S,Q=300000,n=14):
    t=np.arange(1,Q)/Q; f=np.ones(Q-1)
    for v in S:
        f=np.minimum(f,np.abs(((v*t+0.5)%1)-0.5))
    return f.max()
def L_grid(S,Q=300000,n=14):
    t=np.arange(Q)/Q; safe=np.ones(Q,bool)
    for v in S:
        safe &= (np.abs(((v*t+0.5)%1)-0.5)>1/n)
    return safe.mean()
print("(A) the AP cusp: M({1..12,x}) vs x near 13 (downward SPIKE at the integer):")
base=list(range(1,13))
for x in [12.0,12.5,12.9,12.99,13.0,13.01,13.1,13.5,14.0,15.0,26.0,39.0]:
    M=M_grid(base+[x]); spike = "<-- SPIKE" if M<0.072 else ""
    print(f"   x={x:>6}: M={M:.6f}  (1/14={1/14:.5f},1/13={1/13:.5f}) {spike}")
print("   => integer 13 (and multiples 26,39: 13|x => resonance) spike to 1/14; non-integer/coprime stay 1/13")
print("\n(B) cusp depths: M at integer 13-sets (AP deepest = 1/14?):")
sets={"AP {1..13}":list(range(1,14)),"{1..12,14}":list(range(1,13))+[14],
      "{1..12,26}":list(range(1,13))+[26],"{2..14}":list(range(2,15)),
      "{1..11,13,14}":[1,2,3,4,5,6,7,8,9,10,11,13,14],"{1,3,5,..,25} odd":list(range(1,26,2)),
      "{1..13} doubled? {2,4,..,26}":list(range(2,27,2))}
for nm,S in sets.items():
    print(f"   {nm:>26}: M={M_grid(S):.6f}  depth below 1/13 = {1/13-M_grid(S):+.6f}")
print("\n(C) decoupling cusp: send k speeds to infinity, L -> (6/7)^k * L(remaining)?")
core=[1,2,3,4,5]
for k in range(0,5):
    big=[1000+331*j for j in range(k)]   # k 'generic' large speeds
    S=core+big; L=L_grid(S, 400000); Lc=L_grid(core,400000)
    print(f"   core{{1..5}} + {k} huge: L={L:.6f}   (6/7)^{k}*L(core)={ (6/7)**k * Lc:.6f}")

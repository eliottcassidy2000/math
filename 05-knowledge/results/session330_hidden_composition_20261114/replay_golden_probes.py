"""Exact same-order probe census and literal Hamiltonian path controls."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent.parent /
                      "session328_fermat_gluing_and_phase_repair_20261114"))
from collections import Counter
from itertools import combinations
from profile_tools import canon, extend, hp_all, profile, strong

W1 = (2,12,9,49,39,71,31)
W2 = (32,33,11,19,7,92,31)
V1 = (10,36,33,22,7,89,31)
V2 = (16,17,11,35,108,71,15)
U1 = (2,4,81,7,43,15,59)
U2 = (2,4,65,23,39,15,59)

def glue(U,W):
    a,b=len(U),len(W)
    return tuple(t|(((1<<b)-1)<<a) for t in U)+tuple((t<<a)|(1<<(a+b)) for t in W)+((1<<a)-1,)

def depth(h,z):
    d=0
    while z>h:
        h,z=z-h,h
        d+=1
    return d,h,z

classes={(0,)}
for n in range(1,8):
    P={profile(A) for A in classes if strong(A)}
    if n>=3:
        Z={}
        for h,q in P:
            p=(q-4*(-1)**n*h)
            assert p%12==0
            if n%2==0:
                assert q%8==0
            Z[h,q]=depth(h,p//12)
        D=[abs(h*r-k*q) for (h,q),(k,r) in combinations(P,2) if h*r!=k*q]
        print(n,len(classes),len(P),min(D,default=None),dict(sorted(Counter(t[0] for t in Z.values()).items())))
        if n==5 or n==6:
            assert min(D)==24
        if n==7:
            assert min(D)==12 and len(P)==131
            assert [((h,q),(k,r)) for (h,q),(k,r) in combinations(sorted(P),2)
                    if h%3 and k%3 and abs(h*r-k*q)==12]==[((89,1240),(91,1268))]
            assert (Z[89,1240],Z[91,1268])==((2,45,44),(2,46,45))
            assert {p:t for p,t in Z.items() if t[0]>=7}=={(61,944):(7,6,1),(69,1068):(7,7,1)}
    if n<7:
        classes={canon(extend(A,c)) for A in classes for c in range(1<<n)}
assert all(strong(A) for A in (W1,W2,U1,U2))
assert list(map(profile,(W1,W2,U1,U2)))==[(61,944),(69,1068),(37,740),(45,732)]
for U in (U1,U2):
    H,Q=profile(U)
    Z=[hp_all(glue(U,W))[-1] for W in ((0,),W1,W2)]
    assert Z==[777,883*Q-761*H,999*Q-861*H]
    Y1,Y2=Z[1:]
    assert (999*Y1-883*Y2,861*Y1-761*Y2)==(24*H,24*Q)
    print((H,Q), Z)
assert all(strong(A) for A in (V1,V2)) and list(map(profile,(V1,V2)))==[(89,1240),(91,1268)]
U=(2,4,9,3)
H,Q=profile(U)
Z1,Z2=(hp_all(glue(U,W))[-1] for W in (V1,V2))
assert (Z1,Z2)==(1151*Q-973*H,1177*Q-995*H)
assert (1177*Z1-1151*Z2,995*Z1-973*Z2)==(24*H,24*Q)
assert Z1%3 and Z2%3 and Z2%11==6*H%11!=0
print('unit probe', (H,Q), (Z1,Z2))
print('same-order golden probe readout: PASS')

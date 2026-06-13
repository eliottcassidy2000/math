"""
monad-explorer-2026-06-07-S5   (HYP-2298 Q1: the TOWER of bridge lattices)
=========================================================================
Double-glued bridge lattice  L_{a,b} = Z[zeta6]-module on {1, w_a, w_b, w_a w_b},
  w_t = ((2t-1)+i sqrt(4t-1))/(2t),  |w_t|=1,  |1-w_t|^2 = 1/t   (THM-433 engine).

PREDICTION (from the THM-433 divider mechanism):  a UNIT VECTOR of L_{a,b} is one of
  rosettes:            zeta6^j, zeta6^j w_a, zeta6^j w_b, zeta6^j w_a w_b   (4*6=24)
  single dividers:     alpha (1-w_a)  [N(alpha)=a],   alpha (1-w_b) [N=b]
  product divider:     alpha (1-w_a)(1-w_b)  [N(alpha)=ab]   since |(1-w_a)(1-w_b)|^2=1/(ab)
  => #units(L_{a,b}) = 24 + r_E(a) + r_E(b) + r_E(ab)        (if no OTHER families appear)
We TEST this: float-find all unit vectors in a box, then for EACH check it matches a
predicted exact form alpha*prod(1-w) with the right Eisenstein norm.  Reports whether
the prediction is exact or whether extra (unpredicted) unit vectors exist.
"""
from itertools import product
from math import isqrt, gcd

def is3sq(d):  # d == 3*square  -> w_t degenerates into Q(sqrt-3)
    if d%3: return False
    q=d//3; r=isqrt(q); return r*r==q

def rE(t):
    if t<=0: return 0
    exc=0
    for dv in range(1,t+1):
        if t%dv==0: exc += (1 if dv%3==1 else (-1 if dv%3==2 else 0))
    return 6*exc

def wval(t):
    import cmath
    d=4*t-1
    return complex((2*t-1)/(2*t), (d**0.5)/(2*t))

def find_units(a,b,R):
    import cmath
    z6=complex(0.5,(3**0.5)/2); wa=wval(a); wb=wval(b)
    gens=[1, wa, wb, wa*wb]              # Z[zeta6]-module generators
    units=[]
    rng=range(-R,R+1)
    for coords in product(rng,repeat=8):
        if not any(coords): continue
        # element = sum_k (coords[2k]+coords[2k+1]*z6) * gens[k]
        z=0
        for k in range(4):
            z += (coords[2*k]+coords[2*k+1]*z6)*gens[k]
        if abs(abs(z)-1.0)<1e-9:
            units.append(z)
    return units, z6, wa, wb

def eis_norm_of(z, base):
    """if z/base is an Eisenstein integer m+n z6, return its norm, else None."""
    import cmath
    z6=complex(0.5,(3**0.5)/2)
    w=z/base
    # solve w = m + n z6 ; z6=(1+ i sqrt3)/2 -> Im(w)=n sqrt3/2 -> n=2 Im/sqrt3 ; m=Re-n/2
    n=2*w.imag/(3**0.5); m=w.real-n/2
    mr,nr=round(m),round(n)
    if abs(m-mr)<1e-6 and abs(n-nr)<1e-6:
        return mr*mr+mr*nr+nr*nr
    return None

print("="*78)
print("DOUBLE-GLUED BRIDGE LATTICE  L_{a,b}  unit-vector counts (float-find R, then")
print("structural exact-match each vector to alpha*prod(1-w) ).")
print("  predicted: #units = 24 + r_E(a) + r_E(b) + r_E(ab)")
print("="*78)
pairs=[(2,3),(3,4),(2,5),(3,5),(2,3),(3,13),(4,9),(2,4)]
seen=set()
for (a,b) in pairs:
    if (a,b) in seen: continue
    seen.add((a,b))
    da,db=4*a-1,4*b-1
    if is3sq(da) or is3sq(db) or a==b:
        print(f"  L_{{{a},{b}}}: degenerate/invalid, skip"); continue
    for R in (2,3):
        units,z6,wa,wb=find_units(a,b,R)
        # classify
        bases={'1':1,'wa':wa,'wb':wb,'wawb':wa*wb,
               '(1-wa)':(1-wa),'(1-wb)':(1-wb),'(1-wa)(1-wb)':(1-wa)*(1-wb)}
        classmap={}
        extra=0
        for z in units:
            tag=None
            for name,base in bases.items():
                nn=eis_norm_of(z,base)
                if nn is not None:
                    # expected norm
                    exp={'1':1,'wa':1,'wb':1,'wawb':1,'(1-wa)':a,'(1-wb)':b,'(1-wa)(1-wb)':a*b}[name]
                    if nn==exp:
                        tag=name; break
            if tag is None: extra+=1
            else: classmap[tag]=classmap.get(tag,0)+1
        pred=24+rE(a)+rE(b)+rE(a*b)
        if R==3:
            print(f"  L_{{{a},{b}}} (d_a={da},d_b={db}):  #units={len(units):>3}  predicted={pred:>3}  "
                  f"r_E({a})={rE(a)},r_E({b})={rE(b)},r_E({a*b})={rE(a*b)}  extra={extra}  "
                  f"{'OK' if len(units)==pred and extra==0 else 'MISMATCH'}")
            # show class breakdown
            bd=", ".join(f"{k}:{v}" for k,v in classmap.items())
            print(f"        classes: {bd}")

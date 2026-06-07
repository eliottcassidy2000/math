from itertools import product
from math import isqrt
import cmath
def wval(t):
    d=4*t-1; return complex((2*t-1)/(2*t),(d**0.5)/(2*t))
def find(a,b,R):
    z6=complex(0.5,(3**0.5)/2); wa=wval(a); wb=wval(b)
    gens=[1,wa,wb,wa*wb]; out=[]
    for c in product(range(-R,R+1),repeat=8):
        if not any(c): continue
        z=sum((c[2*k]+c[2*k+1]*z6)*gens[k] for k in range(4))
        if abs(abs(z)-1.0)<1e-9: out.append((c,z))
    return out,z6,wa,wb
# stability check L_{2,3}
for R in (2,3,4):
    out,_,_,_=find(2,3,R)
    print(f"L_2,3 R={R}: #units={len(out)}")
# identify extra family at R=3
out,z6,wa,wb=find(2,3,3)
known={'1':1,'wa':wa,'wb':wb,'wawb':wa*wb,'(1-wa)':1-wa,'(1-wb)':1-wb,'(1-wa)(1-wb)':(1-wa)*(1-wb)}
exp={'1':1,'wa':1,'wb':1,'wawb':1,'(1-wa)':2,'(1-wb)':3,'(1-wa)(1-wb)':6}
def eisnorm(w):
    n=2*w.imag/(3**0.5); m=w.real-n/2; mr,nr=round(m),round(n)
    if abs(m-mr)<1e-6 and abs(n-nr)<1e-6: return mr*mr+mr*nr+nr*nr,(mr,nr)
    return None,None
extras=[]
for c,z in out:
    tagged=False
    for nm,base in known.items():
        nn,_=eisnorm(z/base)
        if nn==exp[nm]: tagged=True; break
    if not tagged: extras.append((c,z))
print(f"\n{len(extras)} extra vectors. Trying cross-term bases:")
# candidate extra bases
cross={'wa-wb':wa-wb,'wb-wa':wb-wa,'1-wawb':1-wa*wb,'wa-wb*...':wa-wb,
       '(wa-wb)':wa-wb,'(1-wa wb)':1-wa*wb,'wa(1-wb)':wa*(1-wb),'wb(1-wa)':wb*(1-wa),
       'wa-wawb':wa-wa*wb,'wb-wawb':wb-wa*wb,'1-wa-wb+...':None}
for c,z in extras:
    print(f"  coords={c}")
    for nm,base in cross.items():
        if base is None: continue
        nn,al=eisnorm(z/base)
        if nn is not None and nn<=20:
            print(f"     = (norm {nn} eis {al}) * {nm}   |{nm}|^2={abs(base)**2:.4f}")

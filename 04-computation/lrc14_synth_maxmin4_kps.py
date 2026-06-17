from fractions import Fraction as Fr
def maxmin(S):
    cset=set()
    Sl=sorted(S)
    for v in Sl:
        for k in range(v):
            cset.add(Fr(2*k+1,2*v))
    for i in range(len(Sl)):
        a=Sl[i]
        for j in range(i+1,len(Sl)):
            b=Sl[j]; s=a+b; d=b-a
            for ka in range(a):
                for kb in range(b):
                    t=Fr(ka+kb,s)
                    if 0<t<1: cset.add(t)
                    if d!=0:
                        t2=Fr(kb-ka,d)
                        if 0<t2<1: cset.add(t2)
    best=Fr(0)
    arg=None
    for t in cset:
        m=Fr(1)
        for v in Sl:
            x=v*t; lo=x.numerator//x.denominator
            d1=x-lo; r=(lo+1)-x
            dv=d1 if d1<r else r
            if dv<m: m=dv
            if m<=best: break
        if m>best: best=m; arg=t
    return best,arg
import sys
for name,S in [("AP",list(range(1,14))),("sporadic",[1,2,3,4,5,6,7,8,9,10,11,13,24]),("12->36",[1,2,3,4,5,6,7,8,9,10,11,13,36])]:
    b,t=maxmin(S)
    print(f"{name:10s} max-min={b} = {float(b):.6f}  at tau={t}")
    sys.stdout.flush()

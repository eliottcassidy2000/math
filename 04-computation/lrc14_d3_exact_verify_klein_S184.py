#!/usr/bin/env python3
"""klein-S184 -- exact D3 (Farey-cell piecewise integration) for a given shape, to verify the
exhaustive global minimizer clears bar in exact arithmetic."""
from fractions import Fraction as Fr
import sys
TH=Fr(1,7); M=Fr(6,7); BAR=Fr(83549,252252)
FC={}
def farey(n):
    if n in FC: return FC[n]
    fs=set()
    for d in range(1,n+1):
        for m in range(0,d+1): fs.add(Fr(m,d))
    FC[n]=sorted(fs); return FC[n]
def moments(E):
    E=sorted(E);k=len(E); breaks=farey(E[-1]-E[0]); m1=Fr(0);m2=Fr(0);m3=Fr(0)
    for c in range(len(breaks)-1):
        a,b=breaks[c],breaks[c+1]; mid=(a+b)/2
        cj=[(e*mid).__floor__() for e in E]
        order=sorted(range(k),key=lambda i: E[i]*mid-cj[i])
        ph=[(E[order[r]],Fr(-cj[order[r]])) for r in range(k)]
        gaps=[(Fr(ph[r+1][0]-ph[r][0]),ph[r+1][1]-ph[r][1]) for r in range(k-1)]
        gaps.append((Fr(ph[0][0]-ph[k-1][0]),Fr(1)+ph[0][1]-ph[k-1][1]))
        pts={a,b}
        for (s,ic) in gaps:
            if s!=0:
                xc=(TH-ic)/s
                if a<xc<b: pts.add(xc)
        pts=sorted(pts)
        for t in range(len(pts)-1):
            lo,hi=pts[t],pts[t+1]; m2m=(lo+hi)/2; A=Fr(0);Bc=Fr(0)
            for (s,ic) in gaps:
                if s*m2m+ic>TH: A+=s; Bc+=ic-TH
            # W(x)=A x + Bc on [lo,hi]; integrate W, W^2, W^3
            def ipow(p):  # int_lo^hi (A x + Bc)^p dx
                if p==1: return A/2*(hi*hi-lo*lo)+Bc*(hi-lo)
                if p==2: return A*A/3*(hi**3-lo**3)+A*Bc*(hi*hi-lo*lo)+Bc*Bc*(hi-lo)
                if p==3: return A**3/4*(hi**4-lo**4)+3*A*A*Bc/3*(hi**3-lo**3)+3*A*Bc*Bc/2*(hi*hi-lo*lo)+Bc**3*(hi-lo)
            m1+=ipow(1); m2+=ipow(2); m3+=ipow(3)
    return m1,m2,m3
def D3(E):
    m1,m2,m3=moments(E); den=m2-m3/M
    return m1/M+(m1-m2/M)**2/den if den>0 else m1/M
if __name__=="__main__":
    E=eval(sys.argv[1]) if len(sys.argv)>1 else (0,2,4,6,8,9,10,12,14,16,18)
    d3=D3(E)
    print(f"shape {tuple(E)}: D3 = {d3} = {float(d3):.6f}; bar={float(BAR):.6f}; margin {float(d3-BAR):+.6f} -> {'CLEARS (EXACT)' if d3>=BAR else 'BELOW'}")

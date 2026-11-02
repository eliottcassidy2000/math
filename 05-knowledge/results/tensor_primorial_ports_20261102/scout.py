"""Floating exploration only; never an exact stability certificate."""
import json
from math import exp
from scipy.optimize import differential_evolution

def scout(R):
    n=2*R+1
    k=4/R*(1-.75**R)
    best=None
    for h in range(R):
        L=R-h
        def cost(z):
            x,s=z
            M=h/4+L*x; D=L*s; V=M+D
            if V<=1: return 10+(1-V)
            if D>2*L*L*x: return 10+D-2*L*L*x
            qmax=.25 if h else x
            B=max(L*max(x*s,x*s/2+s*s/4),2*D*D/n,2*V*V/n-qmax*V)
            F=1-.75**h*(1-x)**L
            return max(B+F-k*V,B+k*M-k*V,B+(1-exp(-4*V))/4-k*V)
        out=differential_evolution(cost,[(0,.25),(0,1)],seed=123+R+h,tol=1e-9,popsize=10,maxiter=150)
        row={'R':R,'h':h,'x':float(out.x[0]),'s':float(out.x[1]),'defect':float(out.fun)}
        if best is None or row['defect']<best['defect']: best=row
    return best

if __name__=='__main__':
    for R in [5,8,12,16,20,21,22,23]: print(json.dumps(scout(R)),flush=True)

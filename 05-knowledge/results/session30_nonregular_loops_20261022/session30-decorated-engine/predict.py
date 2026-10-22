"""Parameter-free formal first/second-order predictions. No uniform error claim."""
import math
import numpy as np
from decorated import evaluate

def predict(S,derivative_radius=.025,derivative_points=16):
 n=len(S);M=np.asarray(S,dtype=float)/n
 xs=np.linalg.eigvalsh(1j*M)**2
 rho=float(np.exp(.5*np.sum(np.log1p(xs)-np.log1p(-xs))))
 Ds=[sum(float(np.sum(xs**m))*(2*m)**k/m for m in range(1,151,2)) for k in range(1,5)]
 d1,d2,d3,d4=Ds
 b1=d1;b2=d2+d1*d1;b3=d3+3*d1*d2+d1**3;b4=d4+4*d1*d3+3*d2*d2+6*d1*d1*d2+d1**4
 p1=(b2-b1)/2;p2=(3*b4-2*b3-3*b2+2*b1)/24
 a1=float(evaluate(M,1));a2=float(evaluate(M,2))
 # Cauchy coefficients of A1(exp(t)); derivatives are Euler z-derivatives.
 points=np.exp(2j*np.pi*np.arange(derivative_points)/derivative_points)
 vals=np.array([evaluate(M,1,np.exp(derivative_radius*u)) for u in points])
 ad1=float((np.mean(vals/points)/derivative_radius).real)
 ad2=float((2*np.mean(vals/points**2)/derivative_radius**2).real)
 shift=(ad2+(2*d1-1)*ad1)/2
 first=a1+p1/n
 second=a2+.5*a1*a1+(a1*p1+shift)/n+p2/n**2
 return {'n':n,'rho':rho,'A1':a1,'A2':a2,'D_A1':ad1,'D2_A1':ad2,'P1':p1,'P2':p2,'Y1':first,'Y2':second,'y1':-n*first,'y2':-n*(first+second),'kappa1':rho*(1+first),'kappa2':rho*(1+first+second)}

if __name__=='__main__':
 import json,time
 out=[]
 for q in (7,11,19,23):
  r={i*i%q for i in range(1,q)}
  S=np.array([[0 if i==j else 1 if (j-i)%q in r else -1 for j in range(q-1)] for i in range(q-1)])
  t=time.time();rec=predict(S);rec['family']='deleted Paley';rec['seconds']=time.time()-t;out.append(rec)
 print(json.dumps(out,indent=2))

import importlib.util,json
from pathlib import Path
sp=importlib.util.spec_from_file_location('selector',Path(__file__).with_name('select.py'));sel=importlib.util.module_from_spec(sp);sp.loader.exec_module(sel)
def need(c,m):
 if not c:raise RuntimeError(m)
def compose(x,y):
 h,q,A,G,T,S,E=x;j,k,B,J,U,V,W=y
 return h*j,q*k,A*k+h*B,J*q+j*G,A*J+T*j+h*U+E*V,S*j,h*W
s=[(1,2,b,b,0,b,b) for b in [0,1]]
checks=0
for k in range(4):
 q,A,G,T,D=sel.data(k)
 for a in range(5):
  for b in range(5-a):
   for P in (range(1<<(a-1),1<<a) if a else [0]):
    for R in range(1<<b):
     for j in set([0,k]):
      z=(45,q,A,G,T+2*j,45,45)
      for bit in reversed(sel.digits(P,a)):z=compose(s[bit],z)
      for bit in sel.digits(R,b):z=compose(z,s[bit])
      N=2*z[2]+z[0]-z[4]
      need(N==sel.baseline(k,P,R,b)-2*j,'selector mismatch')
      need(5760*4096**k*2**(a+b)<=N-45,'resource bound')
      need(N>=D-2*k,'monotonic bound')
      checks+=1
need(sel.select(9907) is None,'minimum hostile')
for N in [9909,25031,55275]:need(sel.select(N) is not None,'positive selector')
print(json.dumps({'status':'PASS','context_checks':checks,'scope':'k=0..3,total padding<=4,j=0 or k','positive_targets':[9909,25031,55275],'negative_target':9907},indent=2))

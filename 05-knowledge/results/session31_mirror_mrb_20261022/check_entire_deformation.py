import json
from pathlib import Path
import mpmath as mp
mp.mp.dps=65
h=lambda z:mp.log(z)/z
A=lambda t:2/mp.pi*mp.exp(-t*mp.log(8)/8)*mp.sin(mp.pi*t/8)
def Csum(t):
    return mp.nsum(lambda k:(-1)**k*mp.expm1(t*h(k)),[1,mp.inf])
def Cint(t,j=0):
    def fun(x):
        if not x: return (t if j==0 else 1 if j==1 else 0)/mp.pi
        w=h(1+1j*x)
        return mp.im(w**j*mp.exp(t*w))/mp.sinh(mp.pi*x)
    return mp.quad(fun,[0,1,3,8,20,mp.inf])
rows=[]
for s in ['0','.5','.999999','1','1.000001','2','-1']:
    t=mp.mpf(s); a=Csum(t); b=Cint(t)
    assert abs(a-b)<mp.mpf('1e-55'), (s,a,b)
    rows.append({'t':s,'C':mp.nstr(a,58),'D8':mp.nstr(A(t)-a,58),'method_error':mp.nstr(abs(a-b),8)})
dp=mp.diff(A,1)-Cint(mp.mpf(1),1)
assert abs(dp-mp.diff(lambda t:A(t)-Csum(t),1))<mp.mpf('1e-55')
out={'status':'two-method highprecision verification, not interval arithmetic','rows':rows,'D8_prime_1':mp.nstr(dp,58),'eta_prime_1':mp.nstr(mp.euler*mp.log(2)-mp.log(2)**2/2,58),'nearby_numerical_zero':mp.nstr(mp.findroot(lambda t:A(t)-Cint(t),(.99999,1.00001)),58)}
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))

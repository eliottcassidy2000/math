"""Exact quadratic pullback and independent numerical normalization check."""
from context import save,need
from math import comb,isqrt
from fractions import Fraction
import sympy as s
import mpmath as mp

def run():
    root=s.sqrt(645);u=(127+5*root)/2
    trace=s.simplify(u+1/u);cube_trace=s.simplify(u**3+u**-3)
    X=-432*u**3;J=s.simplify(X+864+432**2/X)
    factor=s.simplify(320*root*(u**3+1)/(u**3-1))
    need(trace==127 and cube_trace==2048002 and J==-960**3 and factor==8127,'exact norm, trace, pullback and derivative')
    # Norm -1 impossible modulo3. Check all smaller positive unit y values.
    small=[(y,isqrt(645*y*y+4)) for y in range(1,5)]
    need(all(x*x!=645*y*y+4 for y,x in small),'no smaller norm-one unit')
    N=100
    A=[comb(2*k,k)*comb(3*k,k)*comb(6*k,3*k) for k in range(N)]
    B=[sum(A[j]*comb(k+j,k-j)*(-432)**(k-j) for j in range(k+1)) for k in range(N)]
    # Companion recurrence (a,b,c)=(432,60,0), independent of convolution.
    recur=[1]
    for k in range(N-1):
        num=-(2*k+1)*(432*k*k+432*k+312)*recur[k]-(432**2)*k**3*(recur[k-1] if k else 0)
        need(num%(k+1)**3==0,'integral recurrence');recur.append(num//(k+1)**3)
    need(B==recur,'independent complete companion coefficients')
    mp.mp.dps=120;uu=(127+5*mp.sqrt(645))/2;values=[]
    for exponent in [1,3]:
        base=432*uu**exponent;z=1/base
        val=24*sum(mp.mpf(B[k])*(-1)**k*(-3669+320*mp.sqrt(645)*(mp.mpf(k)+mp.mpf('.5')))/base**(mp.mpf(k)+mp.mpf('.5')) for k in range(N))
        q=6912*z;K=320*mp.sqrt(645)+7338
        tail=32*mp.sqrt(z)*K*q**N*((N+1)/(1-q)+q/(1-q)**2)
        values.append({'U_exponent':exponent,'value':mp.nstr(val,90),'absolute_error_against_1_over_pi':mp.nstr(abs(val-1/mp.pi),10),'analytic_absolute_tail_bound':mp.nstr(tail,10)})
    classical=mp.sqrt(15)/3200*sum(mp.mpf(A[k])*(5418*k+263)/(-960**3)**k for k in range(N))
    need(abs(classical-1/mp.pi)<mp.mpf('1e-110'),'classical numerical positive control')
    need(abs(mp.mpf(values[1]['value'])-1/mp.pi)<mp.mpf('1e-85'),'corrected branch numerical control')
    out={'scope':'PROVED exact algebra and coefficient transform; cited classical discriminant43 identity; VERIFIED high-precision evaluations with a proved analytic tail estimate (not an interval-arithmetic roundoff certificate). U means (127+5sqrt645)/2.','trace':str(trace),'cube_trace':str(cube_trace),'J':str(J),'derivative_coefficient':str(factor),'first_A':A[:8],'first_B':B[:8],'terms':N,'decimal_precision':mp.mp.dps,'evaluations':values,'classical_value':mp.nstr(classical,90)}
    save('modular_bridge.json.gz',out)
    print(out)
if __name__=='__main__':run()

"""Exact even Collatz conjugates and context separation by halving guards."""
from support import need,save
def v2(n):return (n&-n).bit_length()-1
def C(n):return 3*n+1 if n%2 else n//2
def scaled(x,a):
    need(x>0 and x%(2**a)==0,'scaled source belongs to declared even realm')
    return x//2 if x%(2**(a+1))==0 else 3*x+2**a
def context(n,m):
    r=v2(abs(m-n));k=1
    while 5**k<2**(r+1):k+=1
    t=(-5**k*n)%(2**(r+1));suffix=2**k*t
    a=10**k*n+suffix;b=10**k*m+suffix
    need(0<=suffix<10**k and v2(a)>=k+r+1 and v2(b)==k+r,'same appended suffix separates exact halving clocks')
    return {'n':n,'m':m,'digits':k,'suffix':str(suffix).zfill(k),'sources':[a,b],'clocks':[v2(a),v2(b)]}
def main():
    for a in range(1,9):
        for n in range(1,10001):need(scaled(2**a*n,a)==2**a*C(n),'exact integer conjugacy with native step clock')
    for n in range(1,10001,2):
        k=v2(3*n+1);u=(3*n+1)//2**k
        need(u%3==pow(2,-k,3) and u%2==1,'ternary output phase retains halving exponent')
    tests=[context(n,m) for n in range(1,65) for m in range(n+1,65)]
    selected=[context(10,100),context(1,257),context(701,2107)]
    orbit=[6]
    while orbit[-1]!=2:orbit.append(scaled(orbit[-1],1))
    need(orbit==[6,20,10,32,16,8,4,2],'actual even-realm6/10/4 trajectory')
    save('even_towers.json',{'status':'FINITE_EXACT_WITH_ALL_INTEGER_PROOFS','conjugacy_cases':80000,
      'separating_context_pairs':len(tests),'selected_contexts':selected,'even_orbit':orbit,
      'scope':'an exact conjugate of Collatz; no new convergence claim; all positive source integers are distinguishable by append/halving guards'})
    print('PASS even realms1..8; 2016 distinct-source clock selectors; literal6/10/4 orbit',flush=True)
if __name__=='__main__':main()

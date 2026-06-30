"""
Concentration rate of H. CV^2(H) = G(n)/n! - 1, G(n)=Sum_k k! [x^n] W^k, W=x(1+x^2)/(1-x^2).
[x^n] W^k = Sum_{a+b=(n-k)/2} C(k,a) C(k+b-1,b)  (n-k even).
Dominant: one part=3 => CV^2 ~ 2(n-2)/(n(n-1)) ~ 2/n. Test n*CV^2 -> 2.
"""
from math import comb, factorial
from fractions import Fraction as F
def G(n):
    tot=0
    for k in range(1,n+1):
        if (n-k)%2: continue
        s=(n-k)//2
        inner=sum(comb(k,a)*comb(k+(s-a)-1, s-a) for a in range(s+1))
        tot+=factorial(k)*inner
    return tot
print(f"{'n':>3} {'G(n)/n!':>14} {'CV^2=G/n!-1':>14} {'2(n-2)/(n(n-1))':>16} {'n*CV^2':>9}")
for n in list(range(3,21))+[25,30,40,50]:
    g=G(n); r=F(g,factorial(n)); cv=r-1
    lead=F(2*(n-2),n*(n-1))
    print(f"{n:>3} {float(r):>14.6f} {float(cv):>14.6f} {float(lead):>16.6f} {float(n*cv):>9.5f}")
print("\nG(n) sequence:", [G(n) for n in range(1,13)])
print("CV^2 ~ 2/n: H concentrates; relative spread sqrt(Var)/E ~ sqrt(2/n).")

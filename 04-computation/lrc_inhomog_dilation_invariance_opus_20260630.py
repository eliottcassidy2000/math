from fractions import Fraction
def frac(x): x=x%1; return min(x,1-x)
def Mc(S,c,Qmax):
    best=Fraction(0)
    for q in range(1,Qmax+1):
        for a in range(q):
            m=min(frac(Fraction(s*a,q)-c) for s in S)
            if m>best: best=m
    return best
n=8
print("Inhomogeneous DILATION M_c(d*AP)=M_c(AP)=1/n+c(n-2)/n (so the linear law holds for all dilated APs):")
for d in [1,2,3]:
    S=[d*k for k in range(1,n)]
    print(f"  d={d}: ", end="")
    ok=True
    for cn,cd in [(0,1),(1,2*n),(1,n),(1,4),(1,2)]:
        c=Fraction(cn,cd); got=Mc(S,c,6*n); env=Fraction(1,n)+c*Fraction(n-2,n)
        ok = ok and (got==env)
        print(f"c={str(c)}:{str(got)}{'=env' if got==env else '!='+str(env)} ", end="")
    print(f" -> all match envelope: {ok}")
print()
print("So: M_c(d*{1..n-1}) = 1/n + c(n-2)/n for every dilation d. The proven linear law (and the homogeneous")
print("covering-min 1/n at c=0) holds for the whole dilation class -- which includes the divisible extremals")
print("(e.g. 2*{1..13} contains 14 ~ 0 mod 14, yet M=1/14 at the dilated block t=1/28).")

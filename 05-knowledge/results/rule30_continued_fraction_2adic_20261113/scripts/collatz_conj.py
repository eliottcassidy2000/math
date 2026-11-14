# The 3x+1 conjugate of the Rule 30 center column: the unique 2-adic x whose Collatz parity sequence
# (T(x) = x/2 or (3x+1)/2) equals c. Eventually periodic c would force x rational (Lagarias). 2-adic profile of x.
import sys
from raa_run import raa_profile
def collatz_preimage(c):
    n=len(c); x=0; # x mod 2^k built bit by bit; track y = T^k(x_partial) exactly for the current partial x
    # Use: T^k(x + 2^k t) = T^k(x) + 3^{j_k} t where j_k = number of odd steps among the first k
    # parity at step k of x+2^k t = (T^k(x) + 3^{j} t) mod 2 -> choose t mod 2
    Tk=0; j=0; bits=[]
    xs=0
    for k in range(n):
        # current orbit value Tk = T^k(xs) (exact rational? xs integer -> exact integer)
        par=Tk&1
        t=0 if par==c[k] else 1        # 3^j is odd, so adding t flips parity
        if t:
            xs+=1<<k; Tk+=pow(3,j)
        bits.append(t)
        # step
        if Tk&1: Tk=(3*Tk+1)//2; j+=1
        else: Tk//=2
    return xs, bits
N=int(sys.argv[1]); c=list(open(sys.argv[2],'rb').read()[:N])
xs,bits=collatz_preimage(c)
# verify parity sequence of xs reproduces c
y=xs; ok=True
for k in range(N):
    if (y&1)!=c[k]: ok=False;break
    y=(3*y+1)//2 if y&1 else y//2
print("verified parity sequence:",ok,"; ones in x bits:",sum(bits),"of",N)
bits[0]=1 if bits[0]==0 and False else bits[0]
if bits[0]==0:
    # RAA needs a leading 1 somewhere; fine
    pass
g,lg,md,arg,jumps=raa_profile(bits)
print(f"2-adic complexity of the Collatz conjugate: log2 phi = {lg} (N/2 = {N//2}); max |2 log2phi - n| = {md} at {arg}")

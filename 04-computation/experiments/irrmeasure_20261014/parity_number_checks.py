"""Checks for THM-4634 (parity-number dictionary).  T(x) = x/2 (even), (3x+1)/2 (odd) on Z_2; Q(x) = parity vector;
r(x) = sum v_i 2^{-i-1} in [0,1] (the parity vector read as a binary expansion)."""
from fractions import Fraction as Fr
import itertools
def T(x): return x//2 if x%2==0 else (3*x+1)//2
# (a) Terras: x mod 2^k -> first k parities is a bijection on Z/2^k
for k in range(1,15):
    seen=set()
    for x in range(2**k):
        v=[]; y=x
        for _ in range(k): v.append(y%2); y=T(y)
        seen.add(tuple(v))
    assert len(seen)==2**k, k
print("(a) Terras bijection x mod 2^k <-> first k parities: verified for k <= 14")
# (b) parity numbers of n <= 2000 that reach the cycle {1,2}: r(n) = prefix + 2^-L r(c), r(1)=2/3, r(2)=1/3; denominator 3*2^L
def parity_number(n):
    v=[]; x=n; L=0
    while x not in (1,2):
        v.append(x%2); x=T(x); L+=1
    r=sum(Fr(b,2**(i+1)) for i,b in enumerate(v)) + Fr(1,2**L)*(Fr(2,3) if x==1 else Fr(1,3))
    return r,L
ok=True
for n in range(1,2001):
    r,L=parity_number(n)
    d=r.denominator; 
    # lowest-terms denominator must be exactly 3*2^L with L the hitting time of {1,2} (the last prefix digit is odd: T^(L-1) n in {4, 1/3})
    j=0; dd=d
    while dd%2==0: dd//=2; j+=1
    if dd!=3 or j!=L: ok=False; print("denominator anomaly", n, r, j, L)
print("(b) parity numbers of n <= 2000: lowest-terms denominators are exactly 3*2^L (L = hitting time of {1,2}):", ok, "; r(1) =", parity_number(1)[0], " r(2) =", parity_number(2)[0], " r(27) =", parity_number(27)[0])
# (c) rational cycles from binary words: for a word w of length P with k ones at positions d_0 < ... < d_{k-1} (odd steps),
# x_0 = sum_i 3^{k-1-i} 2^{d_i} / (2^P - 3^k) is a T-periodic point with parity word w (as a 2-adic rational, odd denominator);
# check for all words of length P <= 7 using exact rationals (T extends to rationals with odd denominators)
def Tq(x): return x/2 if (x.numerator%2==0) else (3*x+1)/2   # parity of a rational with odd denominator = parity of numerator
cnt=0; integer_cycles=set()
for P in range(1,8):
    for w in itertools.product((0,1),repeat=P):
        k=sum(w)
        if k==0: continue
        D=2**P-3**k
        if D==0: continue
        pos=[i for i,b in enumerate(w) if b==1]
        S=sum(3**(k-1-i)*2**pos[i] for i in range(k))
        x0=Fr(S,D)
        assert x0.denominator%2==1
        x=x0; word=[]
        for _ in range(P): word.append(x.numerator%2); x=Tq(x)
        assert x==x0 and tuple(word)==w, (w,x0,x,word)
        cnt+=1
        if x0.denominator==1:
            cyc=[]; y=x0
            for _ in range(P): cyc.append(y); y=Tq(y)
            integer_cycles.add(tuple(sorted(set(int(c) for c in cyc))))
print("(c) %d binary words of length <= 7 realised as rational T-cycles (word recovered exactly); integer cycles found:"%cnt, sorted(integer_cycles, key=lambda c:(len(c),c)))
# (d) doubling map: r(T x) = frac(2 r(x)) and the trivial cycle <-> {1/3, 2/3}
print("(d) doubling orbit of 2/3:", [Fr(2,3)*2**i % 1 for i in range(4)], "; of 4/7 = rho(1/5) (period 3, word 100 -> x0 = 1/5; rho(4/5) = 1/7, rho(2/5) = 2/7):", [Fr(4,7)*2**i % 1 for i in range(4)])

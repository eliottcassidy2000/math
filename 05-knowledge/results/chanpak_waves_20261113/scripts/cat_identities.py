# Verify: (1) one-hub closed form H(U); (2) split law; (3) two-segment formula; (4) factor migration;
# all against the 4D recurrence (itself checked against brute-force Hamiltonian path counts in rootchain_check.py).
import random
def step(s,m):
    H,E,L,h=s
    if m=='a': return (H,0,2*L,h)
    if m=='b': return (H+L-E,h,2*L-h,h)
    if m=='c': return (H+L-E,H-E,H+L,H)
def run(word,s=(1,1,1,1)):
    for m in word: s=step(s,m)
    return s
def word_of(U,k):  # MSB first, a=1 b=0
    return ''.join('a' if (U>>(k-1-i))&1 else 'b' for i in range(k))
def Hclosed(U):
    if U==0: return 1
    k=U.bit_length(); s=1
    for t in range(k):
        if not (U>>t)&1: s+=U>>(t+1)
        if t+1<k and ((U>>(t+1))&1) and not ((U>>t)&1): s+=1
    return s
bad=0
for U in range(1<<13):
    k=max(1,U.bit_length())
    for pad in (0,2):   # leading b's are irrelevant
        w='b'*pad+word_of(U,k)
        if run(w)[0]!=Hclosed(U): bad+=1
print("(1) closed form H(U): mismatches",bad)
def rev_comp(W,m): return sum((1-((W>>t)&1))<<(m-1-t) for t in range(m))
def Hm(W,m):  # count of the m-bit word W (leading zeros = leading b's are fine)
    return Hclosed(W)
bad=0
for _ in range(3000):
    V=random.randint(1,1<<12); m=random.randint(1,10); W=random.randrange(1<<m)
    U=(V<<m)|W
    rhs=Hclosed(V)+Hm(W,m)-1+V*rev_comp(W,m)+((V&1) and not (W>>(m-1))&1)
    if Hclosed(U)!=rhs: bad+=1
print("(2) split law H(2^m V+W)=H(V)+H(W)-1+V*revcomp_m(W)+[V odd][top bit of W is 0]: mismatches",bad)
def Q(y): return 1+sum(1<<j for j,ch in enumerate(y) if ch=='b')
def P(x): return 1+sum(1<<(len(x)-1-i) for i,ch in enumerate(x) if ch=='a')
bad=0
for _ in range(3000):
    x=''.join(random.choice('ab') for _ in range(random.randint(0,9)))
    y=''.join(random.choice('ab') for _ in range(random.randint(0,9)))
    Hx=run(x)[0]; Hy=run(y)[0]
    eps=(1 if (len(x)==0 or x[-1]=='b') else 0)*(1 if (len(y)==0 or y[0]=='a') else 0)
    if run(x+'c'+y)[0]!=Hx*Hy+P(x)*Q(y)-eps: bad+=1
print("(3) two-segment formula N(x c y)=H(x)H(y)+P(x)Q(y)-eps: mismatches",bad)
bad=0
for _ in range(3000):
    W=''.join(random.choice('abc') for _ in range(random.randint(0,7)))
    W=W.replace('ac','bc')  # keep c after b/c (mode 1) not required, any word ok
    y=''.join(random.choice('abc') for _ in range(random.randint(0,6)))
    j=random.randint(1,4)
    s1=run(W+'a'*j+'c'+y); s2=run(W+'c'+'b'*j+y)
    if s1!=s2: bad+=1
print("(4) factor migration W a^j c y == W c b^j y (full root state): mismatches",bad)

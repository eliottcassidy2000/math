"""Forward binary construction of protected address pairs; no external packages."""
from math import gcd

def need(x):
    if not x:
        raise RuntimeError("certificate failed")

def order(m):
    if m == 1:
        return 1
    t, v = 1, 2 % m
    while v != 1:
        t, v = t + 1, v * 2 % m
    return t

def apply(word, x=0, z=1, p=0, M=None):
    for b in word:
        b = int(b)
        if b == 0:
            z += x + p
        x, p = 2*x+b, b
        if M is not None:
            x, z = x % M, z % M
    return (x, z, p)

def invert_suffix(s, x, z, m):
    i2 = pow(2, -1, m) if m > 1 else 0
    before = '1' + s
    for i in range(len(s)-1, -1, -1):
        b = int(s[i])
        x = (x-b)*i2 % m
        if not b:
            z = (z-x-int(before[i])) % m
    return x,z

def steer(M, prefix, suffix, x, z, free=''):
    need(prefix.startswith('1') and suffix and z % 2)
    a = (M & -M).bit_length()-1
    t2, m = 2**a, M // 2**a
    d = order(m)
    n = min(a,len(suffix))
    need((x-int(suffix,2)) % 2**n == 0)
    if len(suffix) < a:
        suffix = format(x % t2, '0'+str(a)+'b')
    pre = prefix if prefix.endswith('1') else prefix+'1'
    x0,z0,_ = apply(pre, M=m)
    x1,z1 = invert_suffix(suffix,x,z,m)
    y0,y1 = (x0+1)%m,(x1+1)%m
    C = '01' + '1'*((-2)%d)
    Q = '1' + C + '1'*((-1)%d)
    E = '001' + '1'*((-3)%d)
    V = E + Q*((-3)%m)
    zero = C*((2*y0)%m)
    onward = C*((-2*y1)%m)
    _,z0p,_ = apply(zero, x0,z0,1,m)
    _,drift,_ = apply(onward,-1,0,1,m)
    count = (-2*(z1-z0p-drift))%m
    word = pre+zero+V*count+onward
    reset = '1'*(max(1,(a+d-1)//d)*d)
    r = a+((-3-a)%d)
    W = '1'*r+'001'
    word += reset
    _,z2,_=apply(word, M=t2)
    _,suffix_drift,_=apply(suffix,-1,0,1,t2)
    v = (z-z2-suffix_drift)%t2
    b = t2//2
    if b>1:
        word += W*(m*((v//2)*pow(-m,-1,b)%b))
    word += reset
    pad = W*(m*max(1,b))+reset
    for bit in free:
        word += pad if bit=='0' else '1'*len(pad)
    word += suffix
    U,Z,_ = apply(word,M=M)
    need((U,Z)==(x%M,z%M))
    return word

if __name__ == '__main__':
    import random
    rng = random.Random(2536)
    trials = 0
    for m in (1,3,5,7,9,15,25,27,35,45,49,77,81,125,143,225,255):
        for a in range(7):
            M=m*2**a
            for _ in range(12):
                p='1'+format(rng.randrange(32),'05b')
                s=format(rng.randrange(8),'03b')
                x=rng.randrange(M)
                k=2**min(a,len(s))
                x=(x-x%k+int(s,2)%k)%M
                z=2*rng.randrange(M)+1
                steer(M,p,s,x,z,'01' if _==0 else '')
                trials += 1
    # With the literal three-vertex core (rows 4,1,2), verify the selector
    # at the two exact forbidden values, as congruence targets only.
    for M in (25,36,72,360,1260):
        for N in (7,21):
            e,f={0:(0,1),1:(1,0),2:(0,0)}[N%3]
            # u=e,w=f mod M adjusted through F values only when necessary.
            a=(M & -M).bit_length()-1
            dy=2**a
            odd=M//dy
            c=5+6*e+3*f+12*e*f+2*e+f
            # 3* (F(u)-1) = N-c mod M. Pick F(u) odd using true CRT.
            v=next(j for j in range(1,2*M+1,2) if (3*(j-1)+c-N)%M==0)
            U=steer(M,'10',str(e),e,v)
            V=steer(M,'11',str(f),f,1)
            u,z,_=apply(U,M=M)
            w,t,_=apply(V,M=M)
            need((5+3*(z+t-2)+6*u+3*w+12*u*w+2*e+f-N)%M==0)
    print(f'{trials} forward address certificates; forbidden 7/21 modular controls pass.')

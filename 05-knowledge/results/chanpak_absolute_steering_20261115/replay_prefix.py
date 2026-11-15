from math import gcd, isqrt
from random import Random


def need(t):
    if not t:
        raise AssertionError("prefix programming check failed")


def p(n):
    ans = 0
    while n:
        n, b = divmod(n, 2)
        if not b:
            ans += n + n % 2
    return ans


def advance(x, q, w, M=None):
    for c in w:
        if c == "0":
            q += (x + x % 2)//2
        x = 2*x + int(c)
        if M is not None:
            x %= 2*M
            q %= M
    return x, q


def period(m):
    if m == 1:
        return 1
    k, z = 1, 2 % m
    while z != 1:
        k, z = k+1, 2*z % m
    return k


def compiler(M):
    a = (M & -M).bit_length()-1
    m = M >> a
    o = period(m)
    D = "1"*((-2) % o) + "01"
    R = "1"*((-1) % o) + D + "1" + D*((-2) % m) if m > 1 else ""
    K = o*((a+1+o)//o)
    L0 = o*((a+o)//o)
    resets = []
    for j in range(a):
        d = j+2
        top = L0 + m*K
        zeros = [top,top+d]+[L0+i*K for i in range((-1-pow(2,d,m)) % m)]
        word = "".join("0" if i in zeros else "1" for i in range(top+d,-1,-1))
        resets.append(word*m)
        x, q = advance(2*M-1, 0, resets[-1], M)
        need(x == 2*M-1 and q % m == 0 and q % (1 << a) == (-m << j) % (1 << a))
    rr = len(R)*(m-1)
    need(m == 1 or advance(2*m-1,0,R,m) == (2*m-1,pow(4,-1,m)))
    l = (2*M-1).bit_length()
    length = 2*l + rr + a + 1 + sum(map(len,resets))
    return M,a,m,R,resets,l,length


def append(comp, x, q, xgoal, qgoal):
    M,a,m,R,resets,l,length = comp
    def bits(x):
        return format(x % (2*M), f"0{l}b")
    final = bits(xgoal+(1 << l))
    expected = (qgoal-advance(2*M-1,0,final,M)[1]) % M
    init = bits(-1-(1 << l)*x)
    xr, qr = advance(x,q,init,M)
    chosen = R*((4*(expected-qr)) % m) + "1"*(a+1)
    xr, qr = advance(xr,qr,chosen,M)
    need(xr == 2*M-1 and qr % m == expected % m)
    code = ((expected-qr)*pow(-m,-1,1 << a)) % (1 << a)
    r = "".join(w for i,w in enumerate(resets) if code >> i & 1)
    out = init + chosen + r
    out += "1"*(length-l-len(out)) + final
    need(len(out) == length)
    return out


def program(comp, prefix, x, q):
    M = comp[0]
    refer = p(prefix)//2
    suffix = append(comp,prefix,refer,x,q)
    if prefix == 0:
        result = int(suffix,2)
    else:
        result = (prefix << len(suffix)) + int(suffix,2)
    need(advance(prefix,refer,suffix,M) == (x % (2*M),q % M))
    return result


def bracket(N,M):
    comp = compiler(M)
    L = comp[-1]
    r = max(1,isqrt(N//(1 << L)).bit_length()-1)
    A = 4034 << r
    C = 2855-4034+A
    target = (N-C)//2
    r1 = (A//2) % M
    g = gcd(45,M)
    ui = target*pow(r1,-1,g) % g if g > 1 else 0
    qi = ((target-r1*ui)//g)*pow(45//g,-1,M//g) % (M//g) if M//g > 1 else 0
    D = 1 << L
    lowprefix = max(0,((N-C)//(A+45)+1)//D-1)
    highprefix = max(0,((N-C + A*D-1)//(A*D)))
    out = []
    for prefix in (lowprefix,highprefix):
        u = program(comp,prefix,ui,qi)
        H = C + A*u + 45*p(u)
        need((H-N) % (2*M) == 0)
        out.append(H)
    need(out[0] <= N <= out[1])
    need(max(abs(h-N) for h in out) <= 120000*isqrt(N*D))
    return L, out


def main():
    rng = Random(12911)
    total=0
    for M in [1,2,3,4,5,6,7,8,9,10,11,12,13,16,18,20,24,27,28,32,45,50,64]:
        comp=compiler(M)
        for i in range(40):
            v=rng.randrange(1 << i)
            x=rng.randrange(2*M)
            q=rng.randrange(M)
            r=program(comp,v,x,q)
            need(p(r)//2 % M == q)
            total+=1
    for M in [1,2,3,5,6,8,9,16,45]:
        L=compiler(M)[-1]
        for i in range(10):
            N=(1 << (2*L+20))+2*rng.randrange(1 << (2*L+19))+1
            bracket(N,M)
            total+=1
    print("PASS: independent full-integer recurrence, fixed prefix and two-sided absolute brackets;", total,"controls")
    for i in [8,12,16,20]:
        print("dyadic modulus exponent",i,"bounded suffix length",compiler(1 << i)[-1])


if __name__ == "__main__":
    main()

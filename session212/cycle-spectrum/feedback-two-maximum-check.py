from itertools import product

def exhaustive(n):
    full = (1 << n) - 1
    last = 1 << (n-1)
    ends = []
    for a in range(1, full):
        b = full ^ a
        ends.append((a & -a, 1 << (a.bit_length()-1),
                     b & -b, 1 << (b.bit_length()-1)))
    best = -1
    for x,y in product(range(full+1), repeat=2):
        # Choose the favorable direction of the interapex edge.
        extra = (bool(x&1) and not bool(y&last)) or (bool(y&1) and not bool(x&last))
        total = extra + sum(bool(x&i and not y&j and y&k and not x&l)
                            for i,j,k,l in ends)
        best = max(best,total)
    return best

def formula(n):
    return ((1 << (n//2))-2)*((1 << ((n+1)//2))-2)+1+n%2

def compressed(n, a, b, c):
    s = [2]*a + [1]*b + [0]*c
    literal = 1+b+max(b-1,0)+sum((1<<(j-i-1))*s[i]*(2-s[j])
             for j in range(len(s)) for i in range(j))
    A,B,D = 1<<a,1<<b,1<<c
    value = (1<<n)-4*(A+D)+5 if b == 0 else (1<<n)-(2*B+2)*(A+D)+B+b+3
    assert value == literal
    return value

if __name__ == '__main__':
    wanted = [1,3,5,14,37,86]
    for n,expected in enumerate(wanted, 2):
        got = exhaustive(n)
        assert got == expected
        print(n,got)
    for n in range(6,29):
        m = n-2
        assert formula(n) == max(compressed(n,a,b,m-a-b)
                                 for a in range(m+1) for b in range(m-a+1))
    print('compressed checks passed')

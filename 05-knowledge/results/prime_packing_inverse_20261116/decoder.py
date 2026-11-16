"""Exact inverse for fixed-root binary tournament continuations; no external modules."""
from math import isqrt, log1p

CORE = (2,16,3,23,165,79,31,111)
ROOT = 1

def f(z):
    ans = 1
    while z:
        u = z >> 1
        if not z & 1:
            ans += u + (u & 1)
        z = u
    return ans

def address_count(z, A, g):
    return A*(z+1)+g*f(z)

def residues(q):
    rev = [0]*q
    for r in range(1,q):
        rev[r] = (rev[r>>1]>>1) + (q>>1)*(r&1)
    for r in range(q):
        yield r,q-1-rev[r],f(r)-1,int(2*r<q)

def inverse(N, A, g, q=None):
    """Yield every address reaching N, including z=0 if the initial count matches."""
    if N <= 0 or A <= 0 or g <= 0:
        return
    if q is None:
        q = 1 << (max(1,isqrt(N)) - 1).bit_length()
    if q < 2:
        q = 2
    if q & (q-1):
        raise ValueError('q must be a power of two')
    for r,s,t,d in residues(q):
        K,C = A*q+g*s, A*(r+1)+g*t
        lo = max(0,-((C+2*g-N)//(K+g)))
        hi = (N-C-g)//K
        for b in range(lo,hi+1):
            if C+K*b+g*(f(b)+d*(b&1)) == N:
                yield q*b+r

def append_word(rows,root,z):
    word = [1-int(bit) for bit in reversed(bin(z)[2:])] if z else []
    # Address is reverse-complemented. Physical notation: 0:N, 1:S.
    rows = list(rows)
    for bit in word:
        n=len(rows)
        for j in range(n):
            if bit or j!=root:
                rows[j] |= 1<<n
        rows.append(0 if bit else 1<<root)
    return rows

def hamiltonian(rows):
    n=len(rows)
    D=[[0]*n for _ in range(1<<n)]
    for v in range(n):
        D[1<<v][v]=1
    full=(1<<n)-1
    for mask in range(1,full):
        for v,val in enumerate(D[mask]):
            if val:
                rem=rows[v]&~mask
                while rem:
                    bit=rem&-rem
                    D[mask|bit][bit.bit_length()-1]+=val
                    rem-=bit
    return sum(D[full])

def verify():
    for k in range(1,8):
        q=1<<k
        for r,s,t,d in residues(q):
            for b in range(28):
                z=q*b+r
                if f(z)!=f(b)+t+b*s+d*(b&1):
                    raise ArithmeticError((b,r,q))
    for A,g in [(2,1),(28,13),(30,13),(52,11)]:
        bank={}
        for z in range(1400//A):
            bank.setdefault(address_count(z,A,g),set()).add(z)
        for N in range(1,1400):
            if set(inverse(N,A,g)) != bank.get(N,set()):
                raise ArithmeticError((A,g,N))
    for z in range(50):
        h=hamiltonian(append_word(CORE,ROOT,z))
        if h!=address_count(z,52,11):
            raise ArithmeticError((z,h))
    for N in (63,1799,10007,10009,30011,30013,1000000007,1000000009):
        out=list(inverse(N,52,11))
        print(N,out)
    for A,g in [(28,13),(30,13),(52,11)]:
        print('address mass', A,g,log1p(g/A)/g)
    print('PASS')

if __name__=='__main__':
    verify()

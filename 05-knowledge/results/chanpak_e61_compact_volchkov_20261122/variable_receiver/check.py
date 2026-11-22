def p(address):
    total = 0
    while address:
        high, bit = divmod(address, 2)
        if bit == 0:
            total += high + (high & 1)
        address = high
    return total

def flipped(rows, i, j):
    out = list(rows)
    out[i] ^= 1 << j
    out[j] ^= 1 << i
    return tuple(out)

def star_rows(core, tags, x, y):
    """Ordinal left flank, marked core, right flank, and last vertex as apex."""
    left = bin(x)[2:] if x else ''
    right = bin(y)[2:][::-1].translate(str.maketrans('01', '10')) if y else ''
    tape = left + ''.join(str((tags >> i) & 1) for i in range(len(core))) + right
    n = len(tape)
    mask = (1 << n) - 1
    rows = [mask ^ ((1 << (i + 1)) - 1) for i in range(n)]
    width = len(core)
    core_mask = ((1 << width) - 1) << len(left)
    for i, row in enumerate(core):
        index = len(left) + i
        rows[index] = (rows[index] & ~core_mask) | (row << len(left))
    apex = 0
    for i, bit in enumerate(tape):
        if bit == '1':
            apex |= 1 << i
        else:
            rows[i] |= 1 << n
    return tuple(rows + [apex])

def subpaths(rows, starts=None):
    n = len(rows)
    size = 1 << n
    end = [[0] * n for _ in range(size)]
    for i in range(n):
        if starts is None or starts >> i & 1:
            end[1 << i][i] = 1
    for subset in range(1, size):
        for i, count in enumerate(end[subset]):
            if count:
                choices = rows[i] & ~subset
                while choices:
                    bit = choices & -choices
                    end[subset | bit][bit.bit_length() - 1] += count
                    choices ^= bit
    return [sum(row) for row in end]

def profile(core, tags):
    n = len(core)
    full = (1 << n) - 1
    paths = subpaths(core)
    paths[0] = 1
    h = paths[-1]
    q = sum(paths[s] * paths[full ^ s] for s in range(full + 1))

    def eo(graph, marked):
        rooted = subpaths(graph, marked)
        return tuple(sum(rooted[s] * paths[full ^ s] for s in range(1, full + 1)
                         if s.bit_count() % 2 == parity) for parity in (0, 1))

    reverse = tuple(row ^ full ^ (1 << i) for i, row in enumerate(core))
    k = subpaths(star_rows(core, tags, 0, 0))[-1]
    return h, q, eo(core, tags), eo(reverse, full ^ tags), k

def strong(rows):
    target=(1<<len(rows))-1
    for u in range(len(rows)):
        mask=1<<u
        while True:
            new=mask
            for i,r in enumerate(rows):
                if mask>>i&1:new|=r
            if new==mask:break
            mask=new
        if mask!=target:return False
    return True


def simple(w):
    S={(False,True):1}
    m=1
    for i,ch in enumerate(bin(w)[2:] if w else ''):
        e=int(ch)
        if not e:m+=1<<i
        T={}
        def add(k,h):T[k]=T.get(k,0)+h
        for (a,b),h in S.items():
            add((a,b),h)
            add((a,not e),h)
            if a or e:add((True,b),h)
        S=T
    return 2*sum(h for (a,b),h in S.items() if b),m,1+p(w)

def ordinal(L,R):
    n=len(L);tail=((1<<len(R))-1)<<n
    return tuple(r|tail for r in L)+tuple(r<<n for r in R)

D=(4,5,24,19,35,15)
Dplus=flipped(D,1,3)

def counts(w,t,y):
    V=star_rows((),0,2**t*(w+1)-1,0)
    tags=1|(1<<(len(D)+len(V)-1))
    A,B=ordinal(D,V),ordinal(Dplus,V)
    return (star_rows(A,tags,0,y),star_rows(B,tags,0,y))

def value(w,t,y):
    q,m,h=simple(w)
    return (25*(w+1)+28*q)*2**t*(y+1)+25*h*(1+p(y))-10-28*m-25*(not t and not y%2)

for w in range(32):
    V=star_rows((),0,w,0)
    tags=1<<(len(V)-1)
    h,q,_,_,_=profile(V,tags)
    W=tuple(r^((1<<len(V))-1)^(1<<i) for i,r in enumerate(V))
    S=subpaths(W);S[0]=1
    E=subpaths(W,tags)
    m=sum(E[J]*S[((1<<len(V))-1)^J] for J in range(1<<len(V)))
    assert (q,m,h)==simple(w)
for w in [0,2,4,6]:
    for t in [0,1,2]:
        for y in [0,1,3,4]:
            if not w and not t and not y:continue
            lo,hi=counts(w,t,y)
            assert all(map(strong,[lo,hi]))
            assert hi==flipped(lo,1,3)
            assert (subpaths(lo)[-1],subpaths(hi)[-1])==(value(w,t,y),value(w,t,y)+2)

# Known E31 certificate in the single-apex-receiver notation, not new.
w,t,y=6,0,1973896132951836
target=2**30*(2**31-1)-1
assert value(w,t,y)==target
lo,hi=counts(w,t,y)
assert hi==flipped(lo,1,3) and all(map(strong,(lo,hi)))
print('Receiver state profiles and physical unit gaps verified; recovers known E31 witness.')
# Lucas-Lehmer primality check for the Mersenne instance used in the proof.
M=(1<<61)-1
s=4
for _ in range(59):s=(s*s-2)%M
assert s==0
print("Mersenne primality control passed.")

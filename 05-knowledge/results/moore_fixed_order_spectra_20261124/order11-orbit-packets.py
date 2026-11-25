"""Replays the local voltage test, compact coded chart, and exact path controls."""
from itertools import combinations
from math import factorial

P = 11
N = 270
D = (0, 1, 2, 5)
assert {(x-y) % P for x in D for y in D} == set(range(P))

def gap(d):
    return min(d % P, -d % P)

def pair_girth(s, t, shifts):
    adj = [set() for _ in range(2*P)]
    def edge(x,y):
        adj[x].add(y)
        adj[y].add(x)
    for a,d in enumerate((s,t)):
        if d:
            for x in range(P):
                edge(a*P+x, a*P+(x+d)%P)
    for x in range(P):
        for d in shifts:
            edge(x, P+(x+d)%P)
    return all(len(adj[x]&adj[y]) <= int(y not in adj[x])
               for x,y in combinations(range(2*P),2))

def audit():
    expected = {(1,2):5,(2,4):1,(3,4):2,(3,5):4,(1,5):3}
    for s in range(6):
        for t in range(s,6):
            assert pair_girth(s,t,(0,)) == (s != t or not s)
            good = {d for d in range(1,6) if pair_girth(s,t,(0,d))}
            if s:
                assert good == ({expected[s,t]} if (s,t) in expected else set())
            elif t:
                assert good == set(range(1,6))-{t,gap(2*t)}
            else:
                assert good == set(range(1,6))
            triples = [d for d in combinations(range(1,11),2)
                       if pair_girth(s,t,(0,*d))]
            assert len(triples) == (30 if s==t==0 else 0)
    for s in range(1,6):
        cyclic = {d for pair,d in expected.items() if s in pair}
        assert set(range(1,6)) - {s,gap(2*s)} - cyclic == {gap(4*s)}
    assert [5*t for t in range(55) if t%15==5] == [25,100,175,250]
    assert [c for c in (25,100,175,250) if c<=5*(N-c)] == [25,100,175]
    print('voltage pentagon and residual orbit cuts verified')

# Two labelled vertices per remote orbit. For i<j the 4 cross-pair bits in
# lexicographic clone order form one nibble (low bit first). Orientation bit
# zero agrees with the fixed lexicographic order on (orbit label, clone).
DOUBLES, TRIPLES = list(combinations(range(P),2)), list(combinations(range(P),3))

def owner(i,j):
    d = (j-i)%N
    return i if d<=N//2 else j  # i<j; smaller index breaks antipodal tie

def blocks():
    return list(combinations(range(N),2))

def bits(n,k):
    assert 0<=n<1<<k
    return [(n>>i)&1 for i in range(k)]

def take(a,k):
    return sum(next(a)<<i for i in range(k))

def encode(internal,ports,wrap,data):
    q = {}
    extra = [[] for _ in range(N)]
    empty = [[] for _ in range(N)]
    for i,j in blocks():
        p = (i,j)
        a = owner(*p)
        S = data.get(p,())
        if not S:
            q[p] = 0
            empty[a].append(p)
        elif len(S)==1:
            q[p] = 1+S[0]
        else:
            q[p] = 10+len(S)
            seq,k = (DOUBLES,6) if len(S)==2 else (TRIPLES,8)
            extra[a].extend(bits(seq.index(S),k))
    for a in range(N):
        # Five (index,phase) contacts with the port orbits, 3+4 raw bits
        # compressed to 6 by numbering the 55 vertices 0..54 per color.
        register = bits(internal[a],3)
        for s in ports[a]:
            register.extend(bits(s,6))
        assert len(extra[a])<=30
        register += extra[a] + [0]*(30-len(extra[a]))
        register += wrap[3*a:3*a+3]
        assert len(register)<=66 and len(empty[a])>=82
        for p,b in zip(empty[a], register):
            q[p] = 14*b
    return q

def decode(q):
    data, descriptions = {}, []
    multi = [[] for _ in range(N)]
    empty = [[] for _ in range(N)]
    for p in blocks():
        a = owner(*p)
        v = q[p]
        if v in (0,14):
            empty[a].append(v//14)
        elif 1<=v<=11:
            data[p] = (v-1,)
        else:
            assert v in (12,13)
            multi[a].append((p,v))
    wrap = []
    for a in range(N):
        b = iter(empty[a])
        s = take(b,3)
        ports = tuple(take(b,6) for _ in range(5))
        descriptions.append((s,ports))
        used = 0
        for p,v in multi[a]:
            seq,k = (DOUBLES,6) if v==12 else (TRIPLES,8)
            data[p] = seq[take(b,k)]
            used += k
        assert not take(b,30-used)
        if a<45:
            wrap.extend(bits(take(b,3),3))
    return descriptions,wrap,data

def codec_control():
    # A deterministic storage stress test, not a putative Moore graph.
    # Cyclic distance 1 gives one owned triple per row and distances 2..3
    # two owned doubles. Distances 4..25 give the single connections.
    local = {1:(0,1,3),2:(2,8),3:(5,9)}
    data = {}
    for i,j in blocks():
        d = gap_like(i,j)
        if d<=25:
            data[i,j] = local.get(d,((i+j)%11,))
    internal = [i%6 for i in range(N)]
    ports = [tuple((7*i+13*j)%55 for j in range(5)) for i in range(N)]
    wrap = [(i*i+3*i+4)%7<3 for i in range(135)]
    assert decode(encode(internal,ports,wrap,data)) == (list(zip(internal,ports)),wrap,data)
    print('540 vertices: all 36,315 cross-pair nibbles round trip')

def gap_like(i,j):
    return min(j-i,N-(j-i))

def count(n,arc):
    dp = [[0]*n for _ in range(1<<n)]
    out = [[j for j in range(n) if i!=j and arc(i,j)] for i in range(n)]
    for i in range(n):
        dp[1<<i][i] = 1
    for a in range(1<<n):
        for i,m in enumerate(dp[a]):
            if m:
                for j in out[i]:
                    if not a>>j&1:
                        dp[a | 1<<j][j] += m
    return sum(dp[-1])

def h11(S):
    return count(11,lambda i,j:(j-i)%11 in S)

def heights():
    four = set()
    for m in range(32):
        four.add(h11({(11-d if (m>>(d-1))&1 else d) for d in range(1,6)}))
    assert four == {92411,93027,93467,95095}
    pale = {1,3,4,5,9}
    assert h11(pale)==95095
    assert all(h11({(11-d if gap(d)==s else d) for d in pale}) == 93027
               for s in range(1,6))
    assert count(5,lambda i,j:(i<j) != ((j-i)%5 in (2,3)))==11
    assert factorial(540)<1<<4129 and 95095*factorial(540)<1<<4146
    print('25-site identity packet:',11**5)
    print('36-site identity packet + kernel:',[11**5*h for h in sorted(four)])

if __name__ == '__main__':
    audit()
    codec_control()
    heights()

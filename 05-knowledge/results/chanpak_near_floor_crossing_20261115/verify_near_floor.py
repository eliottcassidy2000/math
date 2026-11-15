#!/usr/bin/env python3
"""Exact span-2 near-floor induction certificate. Standard-library Python only."""

from collections import Counter
from itertools import product, permutations

W = 2
ALPH = tuple(range(4)) # mask bits correspond to distances 1,2
NEW = -2


def check(ok, label):
    if not ok:
        raise RuntimeError(label)


def rev(mask, distance):
    return 1 <= distance <= W and bool(mask & (1 << (distance-1)))


def fragments(state, letter, retain=W):
    """Each item is one possible partial-path transition, with repetitions."""
    def finish(seq):
        def age(i):
            return 0 if i == NEW else i+1 if 0 <= i < retain-1 else -1
        v = tuple(sorted((age(h), age(t)) for h,t in seq))
        if sum(h == -1 for h,t in v) <= 1:
            yield v

    yield from finish(state+((NEW,NEW),))
    for i,(h,t) in enumerate(state):
        if t < 0 or not rev(letter,t+1):
            yield from finish(state[:i]+((h,NEW),)+state[i+1:])
            for j,(h2,t2) in enumerate(state):
                if j != i and h2 >= 0 and rev(letter,h2+1):
                    yield from finish(tuple(p for k,p in enumerate(state)
                                            if k not in (i,j))+((h,t2),))
        if h >= 0 and rev(letter,h+1):
            yield from finish(state[:i]+((NEW,t),)+state[i+1:])


def components(state, letter):
    """SCCs in their directed order; (ages, whether an old vertex exists)."""
    S = list(state)
    beats = {s for s in range(W) if rev(letter,s+1)}
    hit = next((i for i,(a,_) in enumerate(S) if beats.intersection(a)),len(S))
    j = len(S)
    while j>hit and not S[j-1][1] and all(a in beats for a in S[j-1][0]):
        j -= 1
    middle = S[hit:j]
    S[hit:j] = [(tuple(sorted({NEW}.union(*(set(a) for a,_ in middle)))),
                 any(old for _,old in middle))]
    dest = tuple((tuple(0 if a==NEW else a+1 for a in ages if a < W-1),
                  old or (W-1 in ages)) for ages,old in S)
    return dest if dest[0][0] else None


def matrix_model():
    ps = [()]
    ix = {():0}
    trans = [dict() for _ in ALPH]
    for v in ps:
        for a in ALPH:
            for z,k in Counter(fragments(v,a)).items():
                if z not in ix:
                    ix[z] = len(ps)
                    ps.append(z)
                trans[a][ix[z],ix[v]] = k
    return ps,trans


P,M = matrix_model()
D = len(P)
INIT = (1,)+(0,)*(D-1)
check(D==27, "profile dimension")
FINAL = tuple(i for i,z in enumerate(P) if len(z)==1)


def dot(x,y):
    return sum(a*b for a,b in zip(x,y))


def hcount(v):
    return sum(v[i] for i in FINAL)


def push(a,v):
    z = [0]*D
    for (i,j),c in M[a].items():
        z[i] += c*v[j]
    return tuple(z)


def pull(a,r):
    z = [0]*D
    for (i,j),c in M[a].items():
        z[j] += c*r[i]
    return tuple(z)


def close(v,a,b):
    h,t = v[0]
    if h not in (0,1):
        return False
    if t<0:
        return True
    if t>h:
        return not rev(b if h==0 else a,t-h)
    return rev(b if t==0 else a,h-t)


# A nonnegative final row for H - 2*hc, keyed by the two final letters.
G = {(a,b):tuple(sum(1 for v in fragments(s,a,W+1) for z in fragments(v,b,W+2)
                    if len(z)==1 and not close(z,a,b)) for s in P)
     for a,b in product(ALPH,repeat=2)}


def dom(v,w):
    return all(a<=b for a,b in zip(v,w))


def pareto(V):
    out = []
    for v in sorted(set(V),key=lambda x:(sum(x),x)):
        if not any(dom(w,v) for w in out):
            out.append(v)
    return set(out)


def step(F):
    out = {}
    for s,V in F.items():
        for a in ALPH:
            t = components(s,a)
            if t is not None:
                out.setdefault(t,set()).update(push(a,v) for v in V)
    return {t:pareto(V) for t,V in out.items()}


def final_edges(s):
    for a in ALPH:
        t = components(s,a)
        if t is None:
            continue
        for b in ALPH:
            u = components(t,b)
            if u is not None and len(u)==1:
                yield a,b


def exact_small():
    ntested = 0
    for n in range(3,7):
        for word in product((0,),(0,1),*([ALPH]*(n-2))):
            adj = tuple(tuple(u!=v and (not rev(word[v],v-u) if u<v else
                                        rev(word[u],u-v)) for v in range(n))
                        for u in range(n))
            H = cyc_0 = 0
            for z in permutations(range(n)):
                if all(adj[u][v] for u,v in zip(z,z[1:])):
                    H += 1
                    cyc_0 += z[0]==0 and adj[z[-1]][0]
            R = {0}
            S = {0}
            for _ in range(n):
                R = R.union(v for u in R for v in range(n) if adj[u][v])
                S = S.union(v for u in S for v in range(n) if adj[v][u])
            v = INIT
            st = ()
            for a in word[:-2]:
                v = push(a,v)
                if st is not None:
                    st = components(st,a)
            a,b = word[-2:]
            check(hcount(push(b,push(a,v))) == H, ('H',word))
            check(dot(G[a,b],v) == H-2*cyc_0, ('G',word))
            if st is not None:
                st = components(st,a)
                if st is not None:
                    st = components(st,b)
            check((st is not None and len(st)==1) == (len(R)==len(S)==n),
                  ('strong',word))
            ntested += 1
    return ntested


def f(n):
    return (3,5,9)[(n-3)%3]*5**((n-3)//3)


def fit(F,mid,last):
    """Find affine templates v + (3^k-1)/2 * (w-v) by three exact samples."""
    check(F.keys()==mid.keys()==last.keys(), 'fit states')
    out = {}
    for s in F:
        rows = {(v,w) for v in F[s] for w in mid[s]
                if tuple(4*b-3*a for a,b in zip(v,w)) in last[s]}
        check({v for v,w in rows}==F[s] and {w for v,w in rows}==mid[s] and
              {tuple(4*b-3*a for a,b in zip(v,w)) for v,w in rows}==last[s] and
              len(rows)==len(F[s])==len(mid[s])==len(last[s]), 'fit bijection')
        out[s] = {(v,delta(w,v)) for v,w in rows}
    return out


def delta(v,w):
    return tuple(a-b for a,b in zip(v,w))


def shift3(v,c):
    return tuple(a+b for a,b in zip(v,c)),tuple(3*b for b in c)


def certify(F,m):
    """Symbolic two-step induction on both initial point and slope."""
    exact = {s:set() for s in F}
    two_edges = 0
    for s,X in F.items():
        for v,c in X:
            check(all(a>=0 for a in c), 'negative template slope')
            for w,d in X- {(v,c)}:
                check(any(a>b and e>=p for a,b,e,p in zip(v,w,c,d)),
                      'persistent antichain separation')
            for a in ALPH:
                t = components(s,a)
                if t is None:
                    continue
                for b in ALPH:
                    u = components(t,b)
                    if u is None:
                        continue
                    p = push(b,push(a,v))
                    q = push(b,push(a,c))
                    found = False
                    for z,k in F[u]:
                        x,y = shift3(z,k)
                        if x==p and y==q:
                            exact[u].add((z,k))
                        if dom(x,p) and dom(y,q):
                            found = True
                    check(found,'two-step symbolic domination')
                    two_edges += 1
    check(all(exact[s]==F[s] for s in F),'physical two-step attainment')
    pairs_H = {(hcount(v),hcount(c)) for s,X in F.items() if len(s)==1 for v,c in X}
    pairs_G = {(dot(G[a,b],v),dot(G[a,b],c)) for s,X in F.items() for a,b in final_edges(s)
               for v,c in X}
    check(all(v>=m and c>=2*m for v,c in pairs_H) and (m,2*m) in pairs_H,
          'exact H minimum')
    check(all(v>=3*m-2 and c>=6*m for v,c in pairs_G) and (3*m-2,6*m) in pairs_G,
          'exact G minimum')
    return sum(map(len,F.values())),two_edges


def equality_fibres(n,levels):
    """Lossless search of all words with G <= f(n)-2 using backward row minima."""
    lo = f(n)-2
    K = n-2
    rows = {}
    for s in levels[K]:
        rows[K,s] = pareto({G[a,b] for a,b in final_edges(s)})
    for k in range(K-1,-1,-1):
        for s in levels[k]:
            rows[k,s] = pareto({pull(a,r) for a in ALPH
                                if (t:=components(s,a)) is not None
                                for r in rows[k+1,t]})
    def go(k,s,v):
        if not any(dot(r,v)<=lo for r in rows[k,s]):
            return 0
        if k==K:
            c=0
            for a,b in final_edges(s):
                if dot(G[a,b],v)<=lo:
                    check(dot(G[a,b],v)==lo, ('lower counterexample',n))
                    check(hcount(push(b,push(a,v)))==f(n), ('equality counterexample',n))
                    c+=1
            return c
        return sum(go(k+1,t,push(a,v)) for a in ALPH if (t:=components(s,a)) is not None)
    return go(0,(),INIT)


def main():
    print('Independent labeled graph checks H, hc, SCC:',exact_small())
    levels = [{():{INIT}}]
    while len(levels)<=21:
        levels.append(step(levels[-1]))
    h = lambda n:min(hcount(v) for s,X in levels[n].items() if len(s)==1 for v in X)
    g = lambda n:min(dot(G[a,b],v) for s,X in levels[n-2].items() for a,b in final_edges(s) for v in X)
    for n in range(3,19):
        check(h(n)==(3375 if n==16 else 10125 if n==18 else f(n)), ('finite H',n))
        check(g(n)==h(n)-2, ('finite G',n))
    for n,m in [(16,3375),(17,5625)]:
        F = fit(levels[n],levels[n+2],levels[n+4])
        print('Base, H, (templates, checked transitions):',n,m,certify(F,m))
    for n,m in [(16,3375),(17,5625)]:
        for k in range(3):
            j = n+2*k
            check(3**k*m>=f(j) and ((3**k*m==f(j))==(j==17)), ('strictness base',j))
    check(27>25, 'ratio improves under +6')
    eq=[(n,equality_fibres(n,levels)) for n in (*range(3,16),17)]
    check(all(c>0 for n,c in eq), 'nonempty equality strata')
    print('Exact-minimum equality ribbon word multiplicities:',eq)
    print('PASS')


if __name__ == '__main__':
    main()

#!/usr/bin/env python3
"""Exact, standalone inverse certificates; no asserts (also works with -O)."""
from fractions import Fraction
from functools import lru_cache
import json


def require(x, message):
    if not x:
        raise RuntimeError(message)


@lru_cache(None)
def g(n):
    if n <= 2:
        return 1
    m = (n + 1)//2
    return g(m) + (m - m % 2 if n % 2 else 0)


def word(n):
    letters = []
    while n > 2:
        letters.append('Z' if n % 2 else 'O')
        n = (n + 1)//2
    return ''.join(reversed(letters))


def step(s, letter):
    D,E,U,d = s
    return (D,0,2*U,d) if letter == 'O' else (D+U-E,d,2*U-d,d)


def graph_word(w):
    a = [2,0]
    for c in w:
        n = len(a)
        a.append(1 if c == 'Z' else 0)
        for i in range(n):
            if c == 'O' or i:
                a[i] |= 1 << n
    return tuple(a)


def counts(a):
    n=len(a); size=1<<n
    dp=[[0]*n for _ in range(size)]
    H=[0]*size; H[0]=1
    for i in range(n): dp[1<<i][i]=1
    for mask in range(1,size):
        H[mask]=sum(dp[mask])
        for i in range(n):
            z=dp[mask][i]
            if not z: continue
            todo=a[i]&~mask
            while todo:
                b=todo&-todo; todo-=b
                dp[mask|b][b.bit_length()-1]+=z
    return H,dp


def reachable(a,r):
    seen=1<<r
    while True:
        more=seen
        for i in range(len(a)):
            if seen>>i&1: more |= a[i]
        if more == seen: return seen == (1<<len(a))-1
        seen=more


def profile(a,r):
    H,dp=counts(a); n=len(a); full=(1<<n)-1
    rev=tuple(sum(1<<j for j in range(n) if a[j]>>i&1) for i in range(n))
    _,back=counts(rev)
    U=sum(back[m][r]*H[full^m] for m in range(1,full+1) if m>>r&1)
    return H[full],dp[full][r],U,H[full^(1<<r)]


def seeds(nmax):
    out=set()
    witness13=None
    for n in range(2,nmax+1):
        pairs=[(i,j) for i in range(n) for j in range(i+1,n)]
        for bits in range(1<<len(pairs)):
            a=[0]*n
            for k,(i,j) in enumerate(pairs):
                if bits>>k&1: a[i]|=1<<j
                else: a[j]|=1<<i
            for r in range(n):
                if reachable(a,r):
                    s=profile(a,r);out.add(s)
                    if s[0]==13 and witness13 is None:
                        witness13={'adjacency_bitmasks':a,'root':r,'profile':s}
    return out,witness13


def bank_image(bank, bound):
    """All final-Z counts <= bound, every retained-root length, not a depth cut."""
    todo=[s for s in bank if s[2]<=bound];seen=set(todo);image=set()
    while todo:
        s=todo.pop()
        for letter in 'OZ':
            t=step(s,letter)
            if letter=='Z' and t[0]<=bound: image.add(t[0])
            if t[2]<=bound and t not in seen:
                seen.add(t);todo.append(t)
    return image,len(seen)


def main():
    limit=100000
    for n in range(2,limit+1):
        s=(1,0,2,1)
        for c in word(n): s=step(s,c)
        require(s==(g(n),n%2,n,1),('word',n,s))
        require(1<=g(n)<=n,('range',n))
        if n%2: require(2*g(n)>=n+1,('odd bound',n))
        require(g(n)+4*n==g(2*n)+4*n,('redundancy',n))
    for n in range(2,257):
        a=graph_word(word(n));s=profile(a,0)
        require(s==(g(n),n%2,n,1),('physical profile',n,s))
    for k in range(2,31):
        a,b=2**k+1,3*2**(k-1)
        require(g(a)+2*a==g(b)+2*b==3*2**k+3,('collision',k))
    bound=2001
    compact={3}|{g(n)+2*n for n in range(2,(bound-1)//2+1) if g(n)+2*n<=bound}
    oneapex={1}|{g(u) for u in range(3,2*bound,2) if g(u)<=bound}
    require(13 not in oneapex,'one-apex hostile')
    # Complete per-target inverse, checked against the one-pass image.
    for h in range(3,bound+1,2):
        choices=[n for n in range(max(2,(h+2)//3),(h-1)//2+1) if g(n)+2*n==h]
        require(bool(choices)==(h in compact and h!=3),('compact inverse',h))
        choices=[u for u in range(h,2*h,2) if g(u)==h]
        require(bool(choices)==(h in oneapex),('apex inverse',h))
    banks={}; witness=None
    for order in [3,4,5]:
        bank,w=seeds(order)
        require(all(U>=2*d and 0<=E<=D<=U for D,E,U,d in bank),'seed inequalities')
        image,states=bank_image(bank,bound)
        # Seeds themselves may already be strong/count-realizing; union them safely.
        image|={D for D,E,U,d in bank if D<=bound}
        banks[str(order)]={'seed_profiles':len(bank),'expanded_profiles':states,
          'covered_odds_through_2001':len(image),'misses_through_201':[h for h in range(1,202,2) if h not in image]}
        if w: witness=w
    J=16
    lower=sum((Fraction(1,2**(j+2)+1) for j in range(J)),Fraction())
    upper=lower+Fraction(1,2**(J+1))
    require(upper<Fraction(432,1000),'density upper')
    print(json.dumps({'status':'PROVED formulas; FINITE-EXACT replay; no universal odd coverage',
      'word_state_checks':limit-1,'literal_graph_profile_checks':255,
      'compact_image_through_2001':len(compact),
      'compact_misses_through_101':[h for h in range(1,102,2) if h not in compact],
      'oneapex_image_through_2001':len(oneapex),
      'oneapex_misses_through_101':[h for h in range(1,102,2) if h not in oneapex],
      'density_constant_bounds':{'lower':str(lower),'upper':str(upper),'decimal_upper':float(upper)},
      'seed_banks':banks,'ordinary_tournament_with_13':witness},indent=2,sort_keys=True))


if __name__=='__main__': main()

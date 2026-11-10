#!/usr/bin/env python3
"""Independent integer walk sieve and endpoint-deletion audit. No producer imports."""
import itertools, json, sys
from collections import Counter

def require(ok, msg):
    if not ok: raise RuntimeError(msg)

def cayley(p, D):
    return [[int(i != j and (j-i)%p in D) for j in range(p)] for i in range(p)]

def walk_sieve(a):
    # Inclusion-exclusion forces an (n-1)-step walk to visit every vertex.
    # With n slots this is exactly a Hamilton path. Integers are arbitrary precision.
    n=len(a); total=0
    for mask in range(1,1<<n):
        vs=[v for v in range(n) if mask>>v&1]
        neighbors=[[j for j,w in enumerate(vs) if a[v][w]] for v in vs]
        counts=[1]*len(vs)
        for _ in range(n-1): counts=[sum(counts[j] for j in ns) for ns in neighbors]
        total+=(-1 if (n-len(vs))%2 else 1)*sum(counts)
    return total

def subset_paths(a):
    # Remove first endpoint from a path rather than extending producer paths.
    n=len(a); out=[sum(a[v][w]<<w for w in range(n)) for v in range(n)]
    dp=[[0]*n for _ in range(1<<n)]; h=[0]*(1<<n)
    for s in range(1,1<<n):
        vs=s
        while vs:
            vb=vs&-vs; v=vb.bit_length()-1; vs-=vb; rest=s^vb
            if not rest: dp[s][v]=1
            else:
                ws=rest&out[v]; value=0
                while ws:
                    wb=ws&-ws; ws-=wb; value+=dp[rest][wb.bit_length()-1]
                dp[s][v]=value
        h[s]=sum(dp[s])
    return h

def inflate(a):
    n=len(a); b=[r+[1,0] for r in a]
    b.extend([[0]*n+[0,1],[1]*n+[0,0]])
    return b

def brute(a):
    return sum(all(a[x][y] for x,y in zip(q,q[1:])) for q in itertools.permutations(range(len(a))))

def main():
    controls=[]
    for p in (3,5):
        for bits in range(1<<((p-1)//2)):
            D=[i if bits>>(i-1)&1 else p-i for i in range(1,(p+1)//2)]
            a=cayley(p,D); b=inflate(a)
            require(walk_sieve(a)==brute(a)==subset_paths(a)[-1], 'small core control')
            require(walk_sieve(b)==brute(b)==subset_paths(b)[-1], 'small inflated control')
            controls.append([p,D,brute(a),brute(b)])
    rows=[]
    for bits in range(32):
        D=[i if bits>>(i-1)&1 else 11-i for i in range(1,6)]
        a=cayley(11,D); h=subset_paths(a); H=walk_sieve(a)
        require(H==h[-1], 'core independent methods disagree')
        full=(1<<11)-1; F2=sum(h[s]*h[full^s] for s in range(1,full))
        inflated=subset_paths(inflate(a))[-1]
        require(inflated==3*H+F2, 'cover identity disagrees with direct inflation')
        rows.append(dict(D=sorted(D),H=H,F2=F2,inflated=inflated))
    minimum=min(r['H'] for r in rows); gamma=min(r['inflated'] for r in rows)
    require(sorted(set(r['H'] for r in rows))==[92411,93027,93467,95095], 'H spectrum')
    require(minimum==92411 and gamma==1848759, 'minimum/gamma')
    witness=next(r for r in rows if r['D']==[2,3,4,5,10])
    require(witness['inflated']==gamma, 'specified sharp witness')
    minimizing={tuple(r['D']) for r in rows if r['inflated']==gamma}
    orbit={tuple(sorted(a*d%11 for d in witness['D'])) for a in range(1,11)}
    require(minimizing==orbit and len(orbit)==10, 'minimum multiplier orbit')
    producer_masks=sorted(sum(1<<(i-1) for i in range(1,6) if 11-i in D) for D in minimizing)
    require(producer_masks==[1,8,9,12,14,17,19,22,23,30], 'producer mask convention')
    require(minimum**2+22==8539792943>gamma, 'multi-orbit bound')
    # Non-strong fixed singleton padding is a hostile to dropping strongness.
    a=cayley(11,witness['D']); padded=[r+[1] for r in a]+[[0]*12]
    require(subset_paths(padded)[-1]==witness['H']<gamma, 'non-strong hostile')
    seven=[]
    for bits in range(8):
        D=[i if bits>>(i-1)&1 else 7-i for i in range(1,4)]
        h=subset_paths(cayley(7,D)); require(walk_sieve(cayley(7,D))==h[-1], 'seven core sieve')
        seven.append((h[-1],sum(h[s]*h[127^s] for s in range(1,127))))
    mixed_values=[3*a*r['H']+x*r['H']+a*r['F2']+x*r['F2'] for a,x in seven for r in rows]
    require(min(mixed_values)==3049061631, 'mixed minimum')
    require(mixed_values.count(min(mixed_values))*2==120, 'mixed fixed orientations')
    mixed=dict(minimum=min(mixed_values),minimizing_fixed_orientations=120,seven_pairs=seven)
    if '--inject-failure' in sys.argv: require(gamma==1848760,'deliberate wrong expected gamma')
    print(json.dumps(dict(status='PASS',method='integer walk inclusion-exclusion; first-endpoint deletion DP',controls=controls,H_histogram=dict(sorted(Counter(r['H'] for r in rows).items())),gamma=gamma,multi_orbit_bound=minimum**2+22,witness=witness,rows=rows,mixed=mixed),indent=2,sort_keys=True))
if __name__=='__main__': main()

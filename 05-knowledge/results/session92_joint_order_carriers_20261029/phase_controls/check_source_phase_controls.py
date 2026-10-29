"""Independent, bounded replay; stdlib only, no source producer imports."""
from fractions import Fraction as Q
from itertools import permutations, combinations
import json

CHECKS = 0

def require(value, message):
    global CHECKS
    CHECKS += 1
    if not value:
        raise RuntimeError(message)

PAIR = ((66,12,41,17,39,75,30),(66,20,9,51,37,71,30))

def matrix(rows):
    n = len(rows)
    return [[bool(rows[i] & (1 << j)) for j in range(n)] for i in range(n)]

def validate(a):
    n = len(a)
    require(all(not a[i][i] for i in range(n)), 'zero diagonal')
    require(all(a[i][j] != a[j][i] for i,j in combinations(range(n),2)), 'tournament')

def complement_off_path(a):
    n = len(a)
    return [[False if i == j else (j == i+1 if abs(i-j)==1 else not a[i][j])
             for j in range(n)] for i in range(n)]

def exceptional_sink(a, p):
    return [row + [i != p] for i,row in enumerate(a)] + [[j == p for j in range(len(a))]+[False]]

def triangle_substitute(a,p):
    owners = [(i,j) for i in range(len(a)) for j in range(3 if i==p else 1)]
    return [[(y == (x+1)%3 if i==j and x!=y else False) if i==j else a[i][j]
             for j,y in owners] for i,x in owners]

def hamilton_dp(a):
    n=len(a); dp=[[0]*n for _ in range(1<<n)]
    for i in range(n): dp[1<<i][i]=1
    for mask,row in enumerate(dp):
        for i,value in enumerate(row):
            if value:
                for j in range(n):
                    if not mask>>j&1 and a[i][j]: dp[mask|1<<j][j]+=value
    return sum(dp[-1])

def hamilton_permutations(a):
    return sum(all(a[p[i]][p[i+1]] for i in range(len(a)-1)) for p in permutations(range(len(a))))

def good(speeds,e):
    d=Q(1,3)-e; pieces=[(Q(0),Q(1))]
    for v in speeds:
        new=[]
        for left,right in pieces:
            for k in range(v):
                lo=max(left,(k+d)/v); hi=min(right,(k+1-d)/v)
                if lo<=hi: new.append((lo,hi))
        pieces=[]
        for lo,hi in sorted(new):
            if pieces and lo<=pieces[-1][1]: pieces[-1]=(pieces[-1][0],max(hi,pieces[-1][1]))
            else: pieces.append((lo,hi))
    return pieces

def expected(which,e):
    if which==1: return [(Q(4,9)-e/3,Q(4,9)+e/6),(Q(5,9)-e/6,Q(5,9)+e/3)]
    return [(Q(2,9)-e/6,Q(2,9)+e/3),(Q(7,9)-e/3,Q(7,9)+e/6)]

def verify_whole_interval(speeds, limit, formula):
    # Every endpoint is affine in e. Its ordering is constant between crossings.
    # On each open cell the intersection endpoints are fixed affine expressions.
    # Two interior evaluations prove equality with the proposed affine formula;
    # all cell boundaries are checked separately (including singleton components).
    lines={(Q(0),Q(0)),(Q(1),Q(0))}
    for v in speeds:
        for k in range(v):
            lines.add(((Q(k)+Q(1,3))/v,-Q(1,v)))
            lines.add(((Q(k+1)-Q(1,3))/v,Q(1,v)))
    cuts={Q(0),limit}
    for (a,b),(c,d) in combinations(lines,2):
        if b!=d:
            e=(c-a)/(b-d)
            if 0<e<limit: cuts.add(e)
    cuts=sorted(cuts); samples=set(cuts)
    for a,b in zip(cuts,cuts[1:]): samples.update(((2*a+b)/3,(a+2*b)/3))
    for e in sorted(samples): require(good(speeds,e)==formula(e),f'exact cell at {speeds}, {e}')
    return len(cuts)-1,len(samples)

def main():
    controls=[]
    for rows in PAIR:
        a=matrix(rows); validate(a)
        require(all(a[i][i+1] for i in range(6)), 'native path')
        b=complement_off_path(a);validate(b)
        require(complement_off_path(b)==a, 'native complement involution')
        counts=[]
        for c in (exceptional_sink(b,1),triangle_substitute(b,1)):
            validate(c); h=hamilton_dp(c)
            require(h==hamilton_permutations(c),'subset DP/permutation cross-check')
            counts.append(h)
        controls.append(counts)
    require(controls==[[307,777],[293,981]],'inherited tournament counts')
    cells=[]
    for which,speeds in ((1,(1,3,6)),(2,(2,3,6))):
        require(good(speeds,Q(0))==expected(which,Q(0)), 'active phase points')
        cells.append(verify_whole_interval(speeds,Q(1,72),lambda e,w=which:expected(w,e)))
        cells.append(verify_whole_interval(speeds+(7,),Q(1,84),lambda e,w=which:[] if w==1 else expected(w,e)))
        for e in (Q(0),Q(1,72),Q(1,1000)):
            require(sum((b-a for a,b in good(speeds,e)),Q(0))==e,'shared length germ')
    print(json.dumps({'checks':CHECKS,'tournament_counts_by_pair_sink_C3':controls,
        'interval_cells_and_samples':cells,'local_length_germs':[0,1],
        'speed7_survival':[False,True],'status':'FINITE-EXACT inherited source replay'},sort_keys=True))

if __name__=='__main__':main()

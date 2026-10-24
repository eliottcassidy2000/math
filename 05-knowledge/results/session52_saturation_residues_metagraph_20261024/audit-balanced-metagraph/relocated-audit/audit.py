#!/usr/bin/env python3
"""Independent reverse coefficient, full reflection orbit, and literal controls.
This module does not import the producer and uses no tournament-class census.
"""
from collections import Counter
from functools import lru_cache
from itertools import combinations, permutations
from math import gcd
from pathlib import Path
import argparse, hashlib, json

def require(condition, message):
    if not condition: raise RuntimeError(message)

def reverse_coefficients(n):
    """Eliminate highest indexed degree variable, using its exact coefficient.
    State contains already allocated outdegrees including all mandatory path arcs.
    Return coefficients separately by degree of vertex zero, then classify n-1.
    """
    p=n//2
    @lru_cache(None)
    def coefficient(degree):
        r=len(degree)
        if r==1: return tuple(int(degree[0]==a) for a in (p-1,p))
        answer=[0,0]
        # Last vertex has no undetermined edge to its predecessor.
        free=r-2
        for final in (p-1,p):
            win_count=final-degree[-1]
            if not 0<=win_count<=free: continue
            for wins in combinations(range(free),win_count):
                winset=set(wins)
                residual=tuple(degree[i]+(i<free and i not in winset) for i in range(r-1))
                if any(d>p or d+(r-2)<p-1 for d in residual): continue
                part=coefficient(residual)
                answer[0]+=part[0];answer[1]+=part[1]
        return tuple(answer)
    result={}
    for end_label,end_degree in [('L',p-1),('H',p)]:
        output=[0,0]
        for wins in combinations(range(n-2),end_degree):
            winset=set(wins)
            residual=tuple(1+int(i<n-2 and i not in winset) for i in range(n-1))
            part=coefficient(residual)
            output[0]+=part[0];output[1]+=part[1]
        for start_label,value in zip('LH',output):result[start_label+end_label]=value
    return dict(sorted(result.items())),coefficient.cache_info().currsize

def reflection_factors(n):
    """Build factors from literal full-edge reflection orbits, then project first half."""
    p=n//2
    unseen={(a,b) for a in range(n) for b in range(a+2,n)}
    factors=[]
    while unseen:
        a,b=min(unseen)
        orbit={(a,b),(n-1-b,n-1-a)}
        unseen-=orbit
        forward=[0]*p;backward=[0]*p
        for x,y in orbit:
            if x<p:forward[x]+=1
            if y<p:backward[y]+=1
        factors.append((tuple(forward),tuple(backward)))
    return factors

def folded_blue(n):
    p=n//2
    factors=reflection_factors(n)
    # Order differs from producer half-polynomial and is determined by edge orbits.
    poly={(1,)*p:1}
    for factor in factors:
        nextpoly=Counter()
        for a,c in poly.items():
            for term in factor:
                d=tuple(x+y for x,y in zip(a,term))
                if max(d)<=p:nextpoly[d]+=c
        poly=nextpoly
    blue={k:0 for k in ('LL','LH','HL','HH')}
    for d,c in poly.items():
        if all(x in (p-1,p) for x in d):blue['HL' if d[0]==p else 'LH']+=c
    return blue, len(factors)

def matrix(n,word):
    A=[[0]*n for _ in range(n)]
    for i in range(n-1):A[i][i+1]=1
    for k,(i,j) in enumerate((i,j) for i in range(n) for j in range(i+2,n)):
        if (word>>k)&1:A[i][j]=1
        else:A[j][i]=1
    return A

def scores(A):return tuple(map(sum,A))
def transpose(A):return list(map(list,zip(*A)))
def C_matrix(A):
    n=len(A);B=transpose(A)
    for i in range(n-1):B[i][i+1]=1;B[i+1][i]=0
    return B

def R_matrix(A):
    n=len(A);return [[A[n-1-j][n-1-i] for j in range(n)] for i in range(n)]

def extend(A):
    n=len(A);p=n//2
    out=[r+[int(sum(r)==p-1)] for r in A]
    out.append([int(sum(r)==p) for r in A]+[0])
    return out

def cycle_flip(A):
    n=len(A);B=transpose(A)
    for i in range(n):B[i][(i+1)%n]=1;B[(i+1)%n][i]=0
    return B

def encode(A):
    n=len(A)
    return sum(A[i][j]<<k for k,(i,j) in enumerate((i,j) for i in range(n) for j in range(i+2,n)))

def literal(n):
    """Full actual tournament matrices in binary word order, not Gray-code scores."""
    p=n//2;m=(n-1)*(n-2)//2;mask=(1<<m)-1
    endpoint=Counter();blue=Counter();internal_words=set();reflection_moved=0
    fixed_C=0;CR_fixed=0;digest=hashlib.sha256();hostile={}
    for w in range(1<<m):
        A=matrix(n,w);d=scores(A)
        if not all(x in (p-1,p) for x in d):continue
        tag=('H' if d[0]==p else 'L')+('H' if d[-1]==p else 'L')
        endpoint[tag]+=1
        B=C_matrix(A);e=scores(B);R=R_matrix(A)
        require(encode(B)==w^mask,'native C must reverse off-path arcs only')
        require(e==tuple(n-1-d[i]+(i==0)-(i==n-1) for i in range(n)),'C degree transport')
        inside=all(x in (p-1,p) for x in e)
        require(inside==(tag=='HL'),'internal iff high-to-low')
        require(scores(R)==tuple(n-1-x for x in d[::-1]),'literal R reverse-converse scores')
        require(C_matrix(R)==R_matrix(B),'C and R commute')
        isblue=R==A
        if isblue:blue[tag]+=1
        if inside:
            internal_words.add(w)
            lifted=extend(A)
            require(scores(lifted)==(p,)*(n+1),'lift is regular')
            require(all(lifted[i][(i+1)%(n+1)] for i in range(n+1)),'lift closes actual fixed Hamilton cycle')
            require(extend(B)==cycle_flip(lifted),'C equivariance on every actual arc')
            require(all(sum(lifted[i][j] for j in range(n+1) if j not in (i,(i-1)%(n+1),(i+1)%(n+1)))==p-1 for i in range(n+1)),'off-cycle Eulerian outdegrees')
            if A==B:fixed_C+=1
            if R==B:CR_fixed+=1
            if not isblue:reflection_moved+=1
        if tag!='HL' and 'wrong_endpoint_internal' not in hostile:
            hostile['wrong_endpoint_internal']={'word':w,'tag':tag,'source_scores':d,'C_scores':e}
        if transpose(A)!=B and 'full_converse_is_not_C' not in hostile:
            hostile['full_converse_is_not_C']={'word':w,'source_scores':d,'C_scores':e,'converse_scores':scores(transpose(A))}
        digest.update(f'{w}|{tag}|{int(isblue)}|{int(inside)}\n'.encode())
    result={'n':n,'endpoint':dict(sorted(endpoint.items())),'blue_endpoint':dict(sorted(blue.items())),
            'internal_ports':len(internal_words),'literal_digest':digest.hexdigest(),'C_fixed_internal':fixed_C,
            'CR_fixed_internal':CR_fixed,'reflection_moved_internal':reflection_moved,'hostile_witnesses':hostile}
    if n>=4:
        clines={min(w,w^mask) for w in internal_words}
        edge_orbits={min(w,w^mask,encode(R_matrix(matrix(n,w))),encode(R_matrix(matrix(n,w)))^mask) for w in internal_words}
        require(len(clines)*2==len(internal_words),'C line free pairing')
        require(CR_fixed==0,'CR has no fixed word because there is a reflection-fixed off-path edge')
        require(4*len(edge_orbits)==len(internal_words)+sum(blue.values())-blue['LH'],'Klein quotient formula')
        result['native_C_lines']=len(clines);result['optional_R_edge_orbits']=len(edge_orbits)
    return result

def cycle_rotation_fixed(n,step):
    """Count invariant orientations directly on actual anticycle edge orbits.
    This does not use a divisor formula, thus checks composite fixed-set formula.
    """
    N=n+1;p=n//2;edges={tuple(sorted((i,j))) for i in range(N) for j in range(i+1,N) if (j-i)%N not in (1,N-1)}
    orbit_factors=[]
    while edges:
        a,b=min(edges);orbit=[];x,y=a,b
        while (x,y) not in orbit:
            orbit.append((x,y));x=(x+step)%N;y=(y+step)%N
        require((x,y)==(a,b),'odd order cannot reverse one edge inside its rotation orbit')
        unoriented={tuple(sorted(pair)) for pair in orbit}
        require(len(unoriented)==len(orbit),'no edge reversal stabilizer')
        edges-=unoriented
        forward=[0]*N;backward=[0]*N
        for x,y in orbit:forward[x]+=1;backward[y]+=1
        orbit_factors.append((tuple(forward),tuple(backward)))
    counts={(0,)*N:1}
    for terms in orbit_factors:
        after=Counter()
        for ex,c in counts.items():
            for term in terms:
                nxt=tuple(a+b for a,b in zip(ex,term))
                if max(nxt)<=p-1:after[nxt]+=c
        counts=after
    return counts.get((p-1,)*N,0),len(orbit_factors)

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',default='normal.json');parser.add_argument('--max-order',type=int,default=12);parser.add_argument('--literal-max',type=int,default=8);args=parser.parse_args()
    result={'status':'PASS','method':'independent descending coefficient elimination; R polynomial from actual full-edge orbits; binary-order actual matrix enumeration','counts':[],'literal':[],'rotation_controls':[]}
    for n in range(2,args.max_order+1,2):
        ep,states=reverse_coefficients(n);blue,orbits=folded_blue(n)
        row={'n':n,'endpoint':ep,'ports':sum(ep.values()),'internal':ep['HL'],'blue_endpoint':blue,'blue_internal':blue['HL'],'blue_boundary':blue['LH'],'reverse_dp_states':states,'reflection_edge_orbits':orbits}
        require(ep['LL']==ep['HH'],'endpoint reversal')
        if n>=4:
            require(ep['HL']%2==0 and blue['HL']%2==0,'C parity')
            row['native_internal_C_lines']=ep['HL']//2
            row['optional_R_internal_edge_orbits']=(ep['HL']+blue['HL'])//4
            require((ep['HL']+blue['HL'])%4==0,'R edge quotient integrality')
        result['counts'].append(row);print(json.dumps(row),flush=True)
        if n<=args.literal_max:
            lit=literal(n)
            require(all(lit['endpoint'].get(k,0)==v for k,v in ep.items()),'literal endpoint coefficients')
            require(all(lit['blue_endpoint'].get(k,0)==v for k,v in blue.items()),'literal blue coefficients')
            result['literal'].append(lit);print('literal',n,'PASS',flush=True)
        N=n+1
        for step in range(1,N):
            if step!=1 and gcd(step,N)==1:continue
            val,orbitcount=cycle_rotation_fixed(n,step)
            d=gcd(step,N)
            if d==1:require(val==1<<((n-2)//2),'circulant fixed count')
            result['rotation_controls'].append({'n':n,'step':step,'gcd':d,'fixed':val,'edge_orbits':orbitcount})
        if all(N%a for a in range(2,int(N**.5)+1)):require(ep['HL']%N==(1<<((n-2)//2))%N,'prime rotation congruence')
    Path(args.output).write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
if __name__=='__main__':main()

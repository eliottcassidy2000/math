#!/usr/bin/env python3
"""Independent automaton-series and literal-query audit of rank ten."""
from pathlib import Path
from itertools import product
from fractions import Fraction
import hashlib,json,subprocess,sys
ROOT=Path(__file__).resolve().parent
SRC=ROOT.parent/'coarse_response_rank'
KEYS=list(product(range(8),repeat=2))
BASIS=[(0,1),(0,2),(1,0),(1,1),(1,2),(1,3),(1,4),(2,0),(2,1),(2,2)]

def need(ok,msg):
    if not ok:raise RuntimeError(msg)

def kernel(a,c,reverse):
    rows=[0]*5
    for i in range(3):
        rows[i]|=1<<((i+(2 if reverse else 1))%3)
        for mask,p in [(a,3),(c,4)]:
            if mask>>i&1:rows[i]|=1<<p
            else:rows[p]|=1<<i
    rows[3]|=1<<4
    # State is (old multiplicities, used-port mask, last label).
    levels=[{} for _ in range(12)]
    for v in range(5):
        e=[0,0,0];ports=0
        if v<3:e[v]=1
        else:ports=1<<(v-3)
        levels[1][tuple(e)+(ports,v)]=1
    series={}
    for length in range(1,12):
        for state,z in levels[length].items():
            e=state[:3];ports,last=state[3:]
            if ports==3:series[e]=series.get(e,0)+z
            if length==11:continue
            for v in range(5):
                if not rows[last]>>v&1:continue
                nxt=list(e);pp=ports
                if v<3:
                    if nxt[v]==3:continue
                    nxt[v]+=1
                else:
                    bit=1<<(v-3)
                    if pp&bit:continue
                    pp|=bit
                key=tuple(nxt)+(pp,v);levels[length+1][key]=levels[length+1].get(key,0)+z
    return series

def numerators():
    out=[]
    for a,c in KEYS:
        f=kernel(a,c,False);g=kernel(a,c,True);D={}
        for e in product(range(4),repeat=3):
            value=0
            for k,coef in enumerate([1,-3,3,-1]):
                if min(e)>=k:
                    x=tuple(v-k for v in e);value+=coef*(f.get(x,0)-g.get(x,0))
            if value:D[e]=value
        out.append(D)
    return out

def hp(rows):
    n=len(rows);N=1<<n;dp=[[0]*n for _ in range(N)]
    for i in range(n):dp[1<<i][i]=1
    for mask in range(1,N):
        for last,value in enumerate(dp[mask]):
            if not value:continue
            nxt=rows[last]&~mask
            while nxt:
                bit=nxt&-nxt;nxt-=bit;dp[mask|bit][bit.bit_length()-1]+=value
    return sum(dp[-1])

def literal(blocks,a,c,reverse):
    sizes=list(map(len,blocks));offsets=[0,sizes[0],sizes[0]+sizes[1]];n=sum(sizes);rows=[0]*(n+2);rows[n]=1<<(n+1)
    for k,B in enumerate(blocks):
        dest=(k+(2 if reverse else 1))%3
        for v,r in enumerate(B):
            u=offsets[k]+v;rows[u]=(r<<offsets[k])|(((1<<sizes[dest])-1)<<offsets[dest])
            for mask,p in [(a,n),(c,n+1)]:
                if mask>>k&1:rows[u]|=1<<p
                else:rows[p]|=1<<u
    return hp(rows)

def det(A):
    # Fraction-free Bareiss elimination.
    A=[row[:] for row in A];previous=1;sign=1;n=len(A)
    for k in range(n-1):
        pivot=next((i for i in range(k,n) if A[i][k]),None)
        need(pivot is not None,'singular determinant control')
        if pivot!=k:A[pivot],A[k]=A[k],A[pivot];sign=-sign
        p=A[k][k]
        for i in range(k+1,n):
            for j in range(k+1,n):
                numerator=A[i][j]*p-A[i][k]*A[k][j]
                need(numerator%previous==0,'Bareiss exact division');A[i][j]=numerator//previous
        for i in range(k+1,n):A[i][k]=0
        previous=p
    return sign*A[-1][-1]

def main():
    names=['PROOF.md','verify.py','normal.json','optimized.json'];pins={n:hashlib.sha256((SRC/n).read_bytes()).hexdigest() for n in names}
    expected={'PROOF.md':'8bb69ba348858d5c0808a9b781e3f27f5cf4baef583e728ce63e26e788cba140','verify.py':'ce99a320bb47a1632aa3e2f40adb5689e01658d0c1e6c8b1253f579c05e97137','normal.json':'309b514d51f3b93ca0f6461263222d2d1e435cb79eb3a169f77d8171f02abcc0','optimized.json':'309b514d51f3b93ca0f6461263222d2d1e435cb79eb3a169f77d8171f02abcc0'}
    need(pins==expected,'reviewed source pins')
    data=json.loads((SRC/'normal.json').read_text());need((SRC/'normal.json').read_bytes()==(SRC/'optimized.json').read_bytes(),'producer paired receipts')
    D=numerators();support=sorted(set().union(*(set(x) for x in D)));need(len(support)==48 and [list(e) for e in support]==data['kernel_monomials'],'independent full numerator support')
    basis=[D[KEYS.index(k)] for k in BASIS];coefficients=data['twice_reconstruction_coefficients_in_lexicographic_query_order'];need(len(coefficients)==64,'complete coefficient matrix')
    for i,row in enumerate(coefficients):
        need(len(row)==10 and all(isinstance(v,int) and abs(v)<=4 for v in row),'coefficient bounds')
        for e in product(range(4),repeat=3):need(2*D[i].get(e,0)==sum(row[j]*basis[j].get(e,0) for j in range(10)),'full automaton numerator identity')
    allvalues=[];rank=[];checks=0
    for record in data['full64_literal_controls']:
        values=[];gaps=[]
        for a,c in KEYS:
            pair=[literal(record['child_rows'],a,c,rev) for rev in (False,True)];values.append(pair);gaps.append(pair[0]-pair[1]);checks+=2
        need(values==record['full64_response_pairs'],'all exact physical counts')
        observed=[gaps[KEYS.index(k)] for k in BASIS];need(observed==record['basis_contrasts'],'physical basis vector')
        for gap,row in zip(gaps,coefficients):need(2*gap==sum(x*y for x,y in zip(observed,row)),'all physical contrast reconstructions')
        if record['kind']=='rank_witness':rank.append(observed)
        allvalues.append(values)
    determinant=det(rank);need(determinant==-3057647616,'literal rank ten determinant')
    for flags in ([],['-O']):
        replay=subprocess.check_output([sys.executable,*flags,str(SRC/'verify.py')]);need(replay==(SRC/'normal.json').read_bytes(),'normal/optimized producer exact replay')
    need(pins=={n:hashlib.sha256((SRC/n).read_bytes()).hexdigest() for n in names},'frozen source')
    return dict(status='PASS',method='Bounded multiplicity port automaton, full numerator recovery by (1-XYZ)^3, independent forward physical DP and fraction-free Bareiss determinant',full_polynomial_identities=64,polynomial_coordinates_per_identity=64,full_kernel_support=48,literal_physical_counts=checks,exact_determinant=determinant,complete_rank_matrix=rank,independent_full_numerators=[[[*e,v] for e,v in sorted(x.items())] for x in D],source_sha256=pins,scope='Exact rank ten for universal linear coarse-contrast recovery, not a separation lower bound or full named reconstruction.')
if __name__=='__main__':print(json.dumps(main(),indent=2,sort_keys=True))

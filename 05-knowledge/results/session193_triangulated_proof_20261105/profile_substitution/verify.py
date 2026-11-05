#!/usr/bin/env python3
"""Exact multi-port profile transfer, geometric fan transport, and infinite-rank minors.
Only Python's standard library is required. No assertion is removed by -O.
"""
from fractions import Fraction
from itertools import combinations, permutations
from math import comb, factorial, gcd
from pathlib import Path
import hashlib
import json
from seed_data import STATES, EXPECTED, EDGES
ROOT = Path(__file__).resolve().parent

def need(ok, message):
    if not ok:
        raise ValueError(message)

def tournament(rows):
    n=len(rows)
    need(all(r >= 0 and r < 1<<n and not r>>i&1 for i,r in enumerate(rows)), 'row bounds')
    need(all(((rows[i]>>j&1)+(rows[j]>>i&1))==1 for i,j in combinations(range(n),2)), 'tournament')

def path_data(rows):
    n=len(rows); N=1<<n
    end=[[0]*n for _ in range(N)]
    for i in range(n): end[1<<i][i]=1
    for mask in range(1,N):
        for i,c in enumerate(end[mask]):
            if not c: continue
            targets=rows[i]&~mask
            while targets:
                b=targets&-targets; targets-=b
                end[mask|b][b.bit_length()-1]+=c
    h=[sum(row) for row in end]; h[0]=1
    return end,h

def all_bins(h,D):
    """p[d][M]: covers of induced mask M by d labelled possibly-empty paths."""
    N=len(h); p=[[0]*N for _ in range(D+1)]; p[0][0]=1
    for d in range(1,D+1):
        for mask in range(N):
            sub=mask
            while True:
                p[d][mask]+=h[sub]*p[d-1][mask^sub]
                if sub==0: break
                sub=(sub-1)&mask
    return p

def endpoint_data(rows,D):
    n=len(rows); full=(1<<n)-1
    end,h=path_data(rows)
    converse=[full^(1<<i)^row for i,row in enumerate(rows)]
    start,_=path_data(converse)
    return end,start,all_bins(h,D),h

def core_operator(rows,B,D):
    """Rows (a_d,b_d,c_d) give P_W(d)=d[a P_S(d)-b P_S(d-1)+c P_S(d+1)]."""
    end,start,p,h=endpoint_data(rows,D+1); full=len(h)-1
    coefficients={}
    for d in range(1,D+1):
        left=right=cross=0
        for u in B:
            for mask in range(1,full+1):
                left+=start[mask][u]*p[d-1][full^mask]
                right+=end[mask][u]*p[d][full^mask]
        for q in B:
            for u in B:
                if q==u: continue
                for first in range(1,full+1):
                    if not end[first][q]: continue
                    rest=full^first; second=rest
                    while second:
                        cross+=end[first][q]*start[second][u]*p[d-1][rest^second]
                        second=(second-1)&rest
        coefficients[d]=(p[d][full]-cross,left,right)
    return coefficients, sum(end[-1][q] for q in B)

def apply_operator(coeff,values,d):
    a,b,c=coeff[d]
    return d*(a*values[d]-b*values[d-1]+c*values[d+1])

def wrap(A,B,S):
    """Literal order: A, S, w. Every S name has the same exterior incidences."""
    a=len(A); s=len(S); w=a+s; rows=list(A)+[0]*(s+1)
    for i,row in enumerate(S): rows[a+i]=((1<<a)-1)|(row<<a)
    for q in B: rows[q]|=1<<w
    rows[w]=sum(1<<i for i in range(w) if i not in B)
    return rows,list(range(a))+['s'+str(i) for i in range(s)]+['w']

def changed(A,Z):
    a,names=A; z,newnames=Z; pos={v:i for i,v in enumerate(newnames)}
    need(set(names)<=set(newnames),'physical old-name inclusion')
    return [(u,v) for i,u in enumerate(names) for j,v in enumerate(names) if i<j and
            ((a[i]>>j)&1)!=((z[pos[u]]>>pos[v])&1)]

def direct_profile(rows):
    """Independent permutation-and-cuts enumeration of every ordered cover count."""
    n=len(rows); F=[0]*(n+1)
    for word in permutations(range(n)):
        back=sum(not(rows[a]>>b&1) for a,b in zip(word,word[1:]))
        good=n-1-back
        for j in range(back+1,n+1): F[j]+=comb(good,j-1-back)
    return F

def bins_from_F(F,D):
    return [sum(comb(d,j)*F[j] for j in range(1,min(d,len(F)-1)+1)) for d in range(D+1)]

def cycle_count(rows):
    n=len(rows); end=[[0]*n for _ in range(1<<n)]; end[1][0]=1
    for mask in range(1,1<<n):
        for i,c in enumerate(end[mask]):
            targets=rows[i]&~mask
            while c and targets:
                b=targets&-targets; targets-=b
                end[mask|b][b.bit_length()-1]+=c
    return sum(end[-1][i] for i in range(n) if rows[i]&1)

def strong(rows):
    for s in range(len(rows)):
        seen=1<<s
        while True:
            nxt=seen
            for i,row in enumerate(rows):
                if seen>>i&1: nxt|=row
            if seen==nxt: break
            seen=nxt
        if seen!=(1<<len(rows))-1: return False
    return True

def omega(u,v): return u[0]*v[1]-u[1]*v[0]+u[2]*v[3]-u[3]*v[2]

def graph(points):
    rows=[0]*len(points)
    for i,j in combinations(range(len(points)),2):
        a=omega(points[i],points[j]); need(a!=0,'geometric tie')
        x,y=(i,j) if a>0 else (j,i); rows[x]|=1<<y
    return rows

def compile_substitution(outer,inners):
    need(all(v[0]==1 for v in outer+sum(inners,[])),'common positive affine chart')
    D=E=0
    for a,b in combinations(range(len(outer)),2):
        va,vb=outer[a],outer[b]
        for u in inners[a]:
            for v in inners[b]:
                B=(vb[3]-va[3])*(u[2]+v[2])+(va[2]-vb[2])*(u[3]+v[3])
                D=max(D,abs(B)); E=max(E,abs(omega(u,v)))
    M=D+E+1; points=[]
    for (_,p,q,r),block in zip(outer,inners):
        points.extend((1,M*M*p+M*(r*Q-q*R)+P,M*q+Q,M*r+R) for _,P,Q,R in block)
    return points,{'D':D,'E':E,'M':M}

def word_points(word):
    K=3*len(word)+1
    return [(1,0,K,0),(1,0,0,K)]+[(1,3*i,-1 if c&2 else 1,1 if c&1 else -1) for i,c in enumerate(word,1)]

def substitute_graph(Q,factors):
    labels=[(a,i) for a,S in enumerate(factors) for i in range(len(S))]
    rows=[0]*len(labels)
    for x,y in combinations(range(len(labels)),2):
        a,i=labels[x]; b,j=labels[y]
        forward=(factors[a][i]>>j&1) if a==b else (Q[a]>>b&1)
        u,v=(x,y) if forward else (y,x); rows[u]|=1<<v
    return rows

def labelled_tournaments(n):
    pairs=list(combinations(range(n),2))
    for code in range(1<<len(pairs)):
        rows=[0]*n
        for k,(a,b) in enumerate(pairs):
            if code>>k&1: a,b=b,a
            rows[a]|=1<<b
        yield rows

def determinant(a):
    a=[[Fraction(x) for x in row] for row in a]; result=Fraction(1)
    for j in range(len(a)):
        r=next((r for r in range(j,len(a)) if a[r][j]),None)
        if r is None: return 0
        if r!=j: a[r],a[j]=a[j],a[r]; result=-result
        v=a[j][j]; result*=v
        for r in range(j+1,len(a)):
            q=a[r][j]/v
            for k in range(j,len(a)): a[r][k]-=q*a[j][k]
    need(result.denominator==1,'integer determinant')
    return int(result)

def nested_response(t):
    """Coefficients of evaluation at 1 after D^t, without sampling source n."""
    coeff={1:1}
    for _ in range(t):
        out={}
        for d,c in coeff.items():
            out[d]=out.get(d,0)+c*d*d
            out[d+1]=out.get(d+1,0)+c*d
            if d>1: out[d-1]=out.get(d-1,0)-c*d
        coeff=out
    return coeff

def transitive_alpha(rows):
    best=0
    for mask in range(1<<len(rows)):
        vertices=[i for i in range(len(rows)) if mask>>i&1]
        if len(vertices)<=best: continue
        triangle=any((rows[a]>>b&1 and rows[b]>>c&1 and rows[c]>>a&1) or
                     (rows[b]>>a&1 and rows[c]>>b&1 and rows[a]>>c&1)
                     for a,b,c in combinations(vertices,3))
        if not triangle: best=len(vertices)
    return best

def main():
    receipt={'scope':'Exact finite controls and instances of separately proved all-order identities.',
             'fan':{},'nested':{},'small_universe':{}}
    # Exhaustive tiny tournaments, every nonempty incoming port set, arbitrary tiny clocks.
    bank=[rows for n in range(1,4) for rows in labelled_tournaments(n)]
    full_checks=0
    for A in bank:
        for Bmask in range(1,1<<len(A)):
            B=tuple(i for i in range(len(A)) if Bmask>>i&1)
            for S in bank:
                rows,_=wrap(A,B,S); tournament(rows); N=len(rows)
                coeff,_=core_operator(A,B,N)
                sf=direct_profile(S); values=bins_from_F(sf,N+1)
                af=direct_profile(rows); actual=bins_from_F(af,N)
                need(gcd(*af)%gcd(*sf)==0,'entire ordered-profile content is preserved as a divisor')
                need(all(actual[d]==apply_operator(coeff,values,d) for d in range(1,N+1)),
                     'complete multi-port polynomial versus literal permutation cuts')
                need(actual[1]==path_data(rows)[1][-1],'independent vertex subset Hamilton count')
                full_checks+=1
    substitution_checks=0
    for Q in bank:
        for address in range(len(Q)):
            for S in bank:
                factors=[[0] for _ in Q]; factors[address]=S
                rows=substitute_graph(Q,factors)
                need(gcd(*direct_profile(rows))%gcd(*direct_profile(S))==0,
                     'arbitrary quotient context preserves ordered-profile content divisor')
                substitution_checks+=1
    receipt['substitution_content_checks']=substitution_checks
    receipt['small_universe']={'labelled_clock_tournaments':len(bank),'full_polynomial_checks':full_checks,
                              'core_orders':[1,2,3],'clock_orders':[1,2,3],'all_nonempty_port_sets':True}
    # All 14 physical states, exact integral rank-four seed geometry.
    seed=json.loads((ROOT/'seed_geometry.json').read_text())
    C3=[(1,0,1,0),(1,0,0,1),(1,0,-1,-1)]; K1=[(1,0,0,0)]
    T2,_=compile_substitution(C3,[C3,C3,C3])
    clocks={'K1':K1,'TT2':[(1,0,0,0),(1,1,0,0)],'C3':C3,'H13':word_points((1,2,0)),'T2':T2}
    contents={}
    for name,pts in clocks.items():
        bins=all_bins(path_data(graph(pts))[1],len(pts))
        contents[name]=gcd(*(bins[d][-1] for d in range(1,len(pts)+1)))
    need(contents=={'K1':1,'TT2':1,'C3':3,'H13':1,'T2':81},'exact source profile contents')
    receipt['profile_content']=contents
    clock_values={name:all_bins(path_data(graph(pts))[1],2)[2][-1] for name,pts in clocks.items()}
    geometric_checks=physical_checks=cycle_checks=move_checks=0
    for name,(A,B) in STATES.items():
        outer=seed[name]['points']; need(graph(outer)==wrap(A,B,[0])[0],'seed exact pairing signs')
        need(strong(A),'all inherited finite cores are strong')
        if name in ('RL_bridge','RL_reroot','RL'):
            need(len(outer)-transitive_alpha(graph(outer))==3,'seed leaves two-apex class exactly')
        coeff,cycles=core_operator(A,B,2)
        a,unused,b=coeff[1]; need((a,b,cycles)==EXPECTED[name],'inherited core/port coefficients')
        margin=min(abs(omega(u,v)) for u,v in combinations(outer,2))
        receipt['fan'][name]={'core':A,'ports':B,'coefficients':list(EXPECTED[name]),
                              'profile_operator_d1_d2':coeff,'seed_minimum_pairing':margin,'clock_tests':{}}
        for clock,points in clocks.items():
            S=graph(points); inners=[K1 for _ in outer]; inners[len(A)]=points
            lifted,scale=compile_substitution(outer,inners)
            rows,names=wrap(A,B,S)
            need(graph(lifted)==rows and strong(rows),'all labelled substituted arcs and strongness')
            geometric_checks+=1
            h=path_data(S)[1][-1]; p=clock_values[clock]; predicted=a*h+b*p
            entry={'H_by_transfer':predicted,'cycles_by_transfer':cycles*h,'scale':scale}
            if clock in ('K1','TT2','C3'):
                need(path_data(rows)[1][-1]==predicted,'independent physical fan Hamilton count')
                need(cycle_count(rows)==cycles*h,'physical cycle DP')
                physical_checks+=1; cycle_checks+=1
                if clock=='C3':
                    actualp=all_bins(path_data(rows)[1],2)[2][-1]
                    sv=all_bins(path_data(S)[1],3)
                    need(actualp==apply_operator(coeff,[sv[d][-1] for d in range(4)],2),'fan second profile coordinate')
                    entry['P2']=actualp
            receipt['fan'][name]['clock_tests'][clock]=entry
    for clock,points in clocks.items():
        S=graph(points); h=path_data(S)[1][-1]; p=clock_values[clock]
        wrapped={name:wrap(A,B,S) for name,(A,B) in STATES.items()}
        for a,b,flips in EDGES:
            need(changed(wrapped[a],wrapped[b])==flips,'literal physical fan move transport'); move_checks+=1
        H={name:a*h+b*p for name,(a,b,_) in EXPECTED.items()}
        for parent,left,right in [('P','L','R'),('L','LL','LR'),('R','RL','RR')]:
            need(H[left]==3*H[parent]+2*h and H[right]==3*H[parent]+4*h,'homogeneous ternary fan law')
        for final,intermediate in [('LL',['L','LL_birth']),('LR',['L','LR_birth','LR_flip']),
                                   ('RL',['R','RL_bridge','RL_reroot']),('RR',['R','RR_flip','RR_birth'])]:
            need(all(H[z]<H[final] for z in intermediate),'transported nonovershooting physical route')
    # C3 nesting, exact context minors, and feedback witnesses.
    for m in range(9):
        matrix=[[sum(c*j**n for j,c in nested_response(t).items()) for n in range(1,m+2)] for t in range(m+1)]
        expected=factorial(m+1)
        for k in range(1,m+1): expected*=factorial(k)**2
        need(determinant(matrix)==expected,'all-order Vandermonde determinant formula control')
        need(nested_response(m)[m+1]==factorial(m),'nonzero top evaluation coefficient')
        receipt['nested'][str(m)]={'coefficients':nested_response(m),'minor_determinant':expected}
    nested_checks=0
    for n in range(1,4):
        points=[(1,i,0,0) for i in range(n)]; rows=graph(points)
        for t in range(5):
            expected=sum(c*j**n for j,c in nested_response(t).items())
            need(path_data(rows)[1][-1]==expected,'nested polynomial versus vertex DP')
            need(transitive_alpha(rows)==n+t,'nested feedback exact finite control')
            nested_checks+=1
            if t<4:
                points,_=compile_substitution(C3,[points,K1,K1]); rows=graph(points)
    # Genuine nesting order changes both the labelled graph and its count.
    seedP=seed['P']['points']; Pgraph=graph(seedP)
    pp=all_bins(path_data(Pgraph)[1],2)[2][-1]
    nested_geom,_=compile_substitution(C3,[seedP,K1,K1])
    need(pp==644 and path_data(graph(nested_geom))[1][-1]==675,'noncommuting nested order control')
    need(receipt['fan']['P']['clock_tests']['C3']['H_by_transfer']==159,'reverse nested order control')
    receipt['noncommuting_orders']={'P_after_D':159,'D_after_P':675,'P_K1_P2':pp}
    # Exact hostile: replacing a transitive clock by C3 multiplies the digits by three.
    pC3=clock_values['C3']; need(pC3==12,'triangle bin value')
    need((9*3+11*pC3,29*3+33*pC3,31*3+33*pC3)==(159,483,489),'scaled-digit hostile')
    receipt['checks']={'geometry':geometric_checks,'vertex_dp':physical_checks,'cycles':cycle_checks,
                       'physical_moves':move_checks,'nested_dp_feedback':nested_checks,'tiny_full_polynomials':full_checks,
                       'general_substitution_content':substitution_checks}
    raw=(json.dumps(receipt,indent=2,sort_keys=True)+'\n').encode()
    (ROOT/'receipt.json').write_bytes(raw)
    print(json.dumps({'status':'PASS',**receipt['checks'],'receipt_sha256':hashlib.sha256(raw).hexdigest()},sort_keys=True))
if __name__=='__main__': main()

#!/usr/bin/env python3
"""Independent exact atlas and support bijections, no nauty or repo scripts.
A bit in a graph is an edge; in a tournament it reverses i->j, i<j.
Pairs are itertools.combinations(range(n),2), low bit first.
Run: python3 -B walsh_support.py --max-n 6 --out walsh_support.json
"""
from __future__ import annotations
import argparse, itertools, json, math
from collections import Counter
from fractions import Fraction
from pathlib import Path
import sympy as sp


def need(condition, message):
    if not condition: raise AssertionError(message)


def vp(a,p):
    e=0
    while a%p==0: e+=1; a//=p
    return e


def permdata(n):
    pairs=list(itertools.combinations(range(n),2)); index={e:k for k,e in enumerate(pairs)}
    out=[]
    for p in itertools.permutations(range(n)):
        images=[]; source_inversions=0; target_inversions=0
        for k,(i,j) in enumerate(pairs):
            a,b=p[i],p[j]; dest=index[tuple(sorted((a,b)))]; images.append(1<<dest)
            if a>b: source_inversions |= 1<<k; target_inversions |= 1<<dest
        out.append((images,source_inversions,target_inversions))
    return pairs,out


def transform(mask,images):
    res=0
    while mask:
        low=mask&-mask; res |= images[low.bit_length()-1]; mask-=low
    return res


def atlas(n):
    pairs,perms=permdata(n); size=1<<len(pairs); fact=math.factorial(n)
    tours=[]; seen=set()
    for t in range(size):
        if t in seen: continue
        orb={transform(t,a)^c for a,b,c in perms}; seen.update(orb)
        tours.append({'mask':t,'aut':fact//len(orb),'orbit':len(orb)})
    graphs=[]; seen=set()
    for x in range(size):
        if x in seen: continue
        signs={}; even=True
        for a,b,c in perms:
            y=transform(x,a); sign=(-1)**((x&b).bit_count()%2)
            if y in signs and signs[y]!=sign: even=False
            signs[y]=sign
        seen.update(signs)
        if even: graphs.append({'mask':x,'aut':fact//len(signs),'orbit':len(signs),'signs':signs})
    M=[[sum(s*(-1)**((t['mask']&x).bit_count()%2) for x,s in g['signs'].items()) for g in graphs] for t in tours]
    return pairs,perms,tours,graphs,M


def invert_mod(A,p):
    n=len(A); B=[[v%p for v in A[i]]+[int(i==j) for j in range(n)] for i in range(n)]
    for j in range(n):
        k=next(k for k in range(j,n) if B[k][j]); B[j],B[k]=B[k],B[j]
        u=pow(B[j][j],-1,p); B[j]=[(x*u)%p for x in B[j]]
        for k in range(n):
            if k!=j:
                u=B[k][j]; B[k]=[(a-u*b)%p for a,b in zip(B[k],B[j])]
    return [r[n:] for r in B]


def cofactor_matching(M,p=None):
    """First-row/first-eligible-column matching with inverse downdates.
    Exact O(N^3) field operations once the inverse is available.
    p=None: Q. Otherwise p must be prime and matrix invertible mod p.
    """
    n=len(M); rows=list(range(n)); cols=list(range(n)); A=[list(r) for r in M]
    if p is None: Q=[[Fraction(v) for v in r] for r in sp.Matrix(A).inv().tolist()]
    else: Q=invert_mod(A,p)
    match=[]; transcript=[]
    while rows:
        def nz(v): return bool(v if p is None else v%p)
        j=next(j for j in range(len(cols)) if nz(A[0][j]) and nz(Q[j][0]))
        match.append((rows[0],cols[j])); transcript.append([rows[0],cols[j],str(Q[j][0])])
        keep=[k for k in range(len(cols)) if k!=j]
        if p is None:
            h=Q[j][0]; Qnew=[[Q[k][l]-Q[k][0]*Q[j][l]/h for l in range(1,len(rows))] for k in keep]
        else:
            h=pow(Q[j][0],-1,p); Qnew=[[(Q[k][l]-Q[k][0]*Q[j][l]*h)%p for l in range(1,len(rows))] for k in keep]
        A=[[A[i][k] for k in keep] for i in range(1,len(rows))]; Q=Qnew
        rows=rows[1:]; cols=[cols[k] for k in keep]
    need(len({j for i,j in match})==n,'matching columns')
    need(all(M[i][j] and (p is None or M[i][j]%p) for i,j in match),'matching support')
    return match,transcript


def matching_support(M,forbidden_rows=(),forbidden_cols=(),predicate=None):
    n=len(M); rs=[i for i in range(n) if i not in forbidden_rows]; cs=set(range(n))-set(forbidden_cols)
    pred=predicate or bool; owner={}
    def augment(i,seen):
        for j in sorted(cs):
            if not pred(M[i][j]) or j in seen: continue
            seen.add(j)
            if j not in owner or augment(owner[j],seen): owner[j]=i; return True
        return False
    for i in rs:
        if not augment(i,set()): return None
    return sorted((i,j) for j,i in owner.items())


def graph_edges(mask,pairs): return [list(e) for k,e in enumerate(pairs) if mask>>k&1]

def tournament_scores(mask,pairs,n):
    s=[0]*n
    for k,(i,j) in enumerate(pairs): s[j if mask>>k&1 else i]+=1
    return sorted(s)


def audit(n):
    pairs,perms,T,G,M=atlas(n); N=len(T); need(len(G)==N,'class counts')
    m=len(pairs); A=sp.Matrix(M); det=int(A.det()); need(det and abs(det)&(abs(det)-1)==0,'power-two determinant')
    Dg=sp.diag(*[x['aut'] for x in G]); Dt=sp.diag(*[x['aut'] for x in T])
    need(A*Dg*A.T==2**m*Dt,'weighted row orthogonality')
    Ot=sp.diag(*[x['orbit'] for x in T]); Og=sp.diag(*[x['orbit'] for x in G])
    need(A.T*Ot*A==2**m*Og,'weighted column orthogonality')
    B=[[Fraction(G[j]['aut']*M[i][j]**2,2**m*T[i]['aut']) for j in range(N)] for i in range(N)]
    need(all(sum(r)==1 for r in B),'stochastic row sums')
    need(all(sum(B[i][j] for i in range(N))==1 for j in range(N)),'stochastic column sums')
    rational,trace=cofactor_matching(M)
    need(matching_support(M) is not None,'independent augmenting matching')
    # Every deleted row/column pair has a support matching (including zero entries).
    for i in range(N):
        for j in range(N): need(matching_support(M,(i,),(j,)) is not None,'strict Hall deletion')
    prime_data={}
    for p in [3,5,7]:
        mp,tp=cofactor_matching(M,p)
        need(all(vp(T[i]['aut'],p)==vp(G[j]['aut'],p) for i,j in mp),'p valuation matching')
        need(Counter(vp(x['aut'],p) for x in T)==Counter(vp(x['aut'],p) for x in G),'valuation distribution')
        prime_data[str(p)]={'valuation_distribution':dict(sorted(Counter(vp(x['aut'],p) for x in T).items())),
            'matching':[[T[i]['mask'],G[j]['mask'],M[i][j]] for i,j in mp], 'cofactor_trace':tp,
            'nonzero_entries':sum(v%p!=0 for r in M for v in r)}
        # Each nonzero entry obeys the valuation inequality.
        need(all(not M[i][j]%p or vp(T[i]['aut'],p)<=vp(G[j]['aut'],p) for i in range(N) for j in range(N)),'valuation filter')
    odd=lambda a:a//(a&-a)
    both35=matching_support(M,predicate=lambda x:math.gcd(x,15)==1)
    summary={'n':n,'N':N,'det':det,'det_exponent':abs(det).bit_length()-1,
      'support_entries':sum(bool(v) for r in M for v in r),'possible_entries':N*N,
      'rank_mod2':int(any(v%2 for r in M for v in r)),
      'tournament_aut_histogram':dict(sorted(Counter(t['aut'] for t in T).items())),
      'graph_aut_histogram':dict(sorted(Counter(g['aut'] for g in G).items())),
      'joint_odd_part_equal':Counter(odd(t['aut']) for t in T)==Counter(odd(g['aut']) for g in G),
      'coprime15_matching_exists':both35 is not None,
      'rational_matching':[[T[i]['mask'],G[j]['mask'],M[i][j]] for i,j in rational],
      'rational_cofactor_trace':trace,
      'prime_data':prime_data,
      'tournaments':[{**t,'score_sequence':tournament_scores(t['mask'],pairs,n)} for t in T],
      'even_graphs':[{k:v for k,v in g.items() if k!='signs'}|{'edges':graph_edges(g['mask'],pairs)} for g in G],
      'matrix':M,'stochastic_matrix':[[str(v) for v in r] for r in B]}
    print(f'n={n} N={N} det={det} support={summary["support_entries"]}/{N*N} odd-part-equal={summary["joint_odd_part_equal"]} mod15-unit-matching={both35 is not None}',flush=True)
    for p,data in prime_data.items(): print(f'  p={p}: valuation distribution={data["valuation_distribution"]}; all matched entries p-units and aut valuations preserved',flush=True)
    if n<=5:
        print('  rational matching (tournament mask, graph mask, Phi):',summary['rational_matching'],flush=True)
        print('  mod3 matching:',prime_data['3']['matching'],flush=True)
        print('  mod5 matching:',prime_data['5']['matching'],flush=True)
    return summary


def main():
    parser=argparse.ArgumentParser(); parser.add_argument('--max-n',type=int,default=6); parser.add_argument('--out',default='walsh_support.json'); args=parser.parse_args()
    need(0<=args.max_n<=6,'bounded exhaustive implementation supports n<=6')
    data=[audit(n) for n in range(args.max_n+1)]
    Path(args.out).write_text(json.dumps(data,indent=2)+'\n')
    print('ALL EXACT CHECKS PASSED',flush=True)

if __name__=='__main__':main()

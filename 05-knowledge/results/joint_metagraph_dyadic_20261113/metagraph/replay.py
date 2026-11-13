#!/usr/bin/env python3
"""Exact native complement-line walk compression; no third-party dependency."""
from collections import Counter, defaultdict
from fractions import Fraction as F
from itertools import combinations, permutations, product
from functools import lru_cache, reduce
from pathlib import Path
import hashlib, json

CHECKS = 0
def require(test, message):
    global CHECKS
    CHECKS += 1
    if not test:
        raise RuntimeError(message)

def mm(A, B):
    if not A: return []
    out = [[F(0) for _ in B[0]] for _ in A]
    for i, row in enumerate(A):
        for k, x in enumerate(row):
            if x:
                for j, y in enumerate(B[k]):
                    if y: out[i][j] += x*y
    return out

def column_rank_sequence(P, L):
    basis={}
    def insert(col):
        v=col[:]
        for pivot,b in sorted(basis.items()):
            if v[pivot]:
                scale=v[pivot];v=[x-scale*y for x,y in zip(v,b)]
        pivot=next((i for i,x in enumerate(v) if x),None)
        if pivot is not None:
            scale=v[pivot];basis[pivot]=[x/scale for x in v]
    V=L; dims=[]
    while True:
        for col in zip(*V):insert(list(col))
        dims.append(len(basis))
        if len(dims)>1 and dims[-1]==dims[-2]:return dims
        V=mm(P,V)

def add(A, B):
    return [[x+y for x,y in zip(a,b)] for a,b in zip(A,B)]

def eye(n): return [[F(i==j) for j in range(n)] for i in range(n)]

def make(n, w):
    T = [[0]*n for _ in range(n)]
    for i in range(n-1): T[i][i+1]=1
    for k,(i,j) in enumerate((e for e in combinations(range(n),2) if e[1]>e[0]+1)):
        T[i][j]=(w>>k)&1
        T[j][i]=1-T[i][j]
    return tuple(tuple(r) for r in T)

def converse(T): return tuple(tuple(T[j][i] for j in range(len(T))) for i in range(len(T)))

@lru_cache(None)
def canon(T):
    # Sorting score buckets is invariant; permutations within them exhaust all
    # bijections respecting that sort. This code is not the unrestricted min.
    d=tuple(map(sum,T)); groups=defaultdict(list)
    for i,s in enumerate(d): groups[s].append(i)
    E=list(combinations(range(len(T)),2))
    return min(sum(T[p[i]][p[j]]<<k for k,(i,j) in enumerate(E))
               for parts in product(*(permutations(groups[s]) for s in sorted(groups)))
               for p in [sum((tuple(x) for x in parts),())])

def merged(T): return min(canon(T),canon(converse(T)))

def profile(T):
    d=tuple(sorted(map(sum,T))); n=len(T)
    return min(d,tuple(sorted(n-1-x for x in d)))

def hp(T):
    n=len(T);D=[[0]*n for _ in range(1<<n)]
    for i in range(n):D[1<<i][i]=1
    for s in range(1,1<<n):
        for j in range(n):
            if s>>j&1:
                D[s][j]+=sum(D[s^(1<<j)][i]*T[i][j] for i in range(n) if i!=j)
    return sum(D[-1])

def autom(T):
    n=len(T);ds=list(map(sum,T));E=list(combinations(range(n),2))
    return sum(all(ds[i]==ds[p[i]] for i in range(n)) and
               all(T[i][j]==T[p[i]][p[j]] for i,j in E)
               for p in permutations(range(n)))

def wire(n):
    m=(n-1)*(n-2)//2; total=1<<m; allmask=total-1
    tables=[make(n,w) for w in range(total)]
    nodes=sorted({merged(T) for T in tables}); ix={c:i for i,c in enumerate(nodes)}
    names=[ix[merged(T)] for T in tables]; reps={names[w]:tables[w] for w in range(total)}
    A=[[0]*len(nodes) for _ in nodes]; C=[[[0]*len(nodes) for _ in nodes] for _ in range(2)]
    for w,T in enumerate(tables):
        i,j=names[w],names[w^allmask]
        blue=all(T[a][b]==T[n-1-b][n-1-a] for a,b in combinations(range(n),2))
        A[i][j]+=1; C[int(blue)][i][j]+=1
    d=list(map(sum,A)); keys=sorted({profile(T) for T in tables}); keyix={s:i for i,s in enumerate(keys)}
    blocks=[keyix[profile(reps[i])] for i in range(len(nodes))]
    require(sum(d)==total, 'all darts')
    require(all(A[i][j]==A[j][i] for i in range(len(nodes)) for j in range(len(nodes))), 'symmetric native adjacency')
    require(all(C[c][i][i]%2==0 for c in range(2) for i in range(len(nodes))), 'loop doubling')
    require(all(C[0][i][j]%2==0 for i in range(len(nodes)) for j in range(len(nodes))), 'black twin multiplicity')
    for i,T in reps.items():
        factor=1+(canon(T)!=canon(converse(T)))
        require(d[i]*autom(T)==factor*hp(T), 'source weight H/Aut')
    # Independent second path: enumerate every physical path of one representative,
    # then divide by its free automorphism action and retain both ordinary sheets.
    A2=[[0]*len(nodes) for _ in nodes]
    for i,T in reps.items():
        counts=Counter(); factor=1+(canon(T)!=canon(converse(T))); aut=autom(T)
        for p in permutations(range(n)):
            if not all(T[a][b] for a,b in zip(p,p[1:])):continue
            U=[list(r) for r in converse(T)]
            for a,b in zip(p,p[1:]):U[a][b],U[b][a]=1,0
            counts[ix[merged(tuple(tuple(r) for r in U))]]+=1
        for j,v in counts.items():
            require((v*factor)%aut==0,'free path orbit divisor')
            A2[i][j]=v*factor//aut
    require(A==A2,'tiling and physical-path native matrices agree')
    return nodes,reps,A,C,d,keys,blocks

def analyze(n):
    nodes,reps,A,C,d,keys,blocks=wire(n);N=len(nodes);b=len(keys)
    mass=[sum(d[i] for i in range(N) if blocks[i]==j) for j in range(b)]
    I=[[F(blocks[i]==j) for j in range(b)] for i in range(N)]
    E=[[F(d[i],mass[j]) if blocks[i]==j else F(0) for i in range(N)] for j in range(b)]
    P=[[F(x,d[i]) for x in row] for i,row in enumerate(A)]
    krylov=column_rank_sequence(P,I)
    require(krylov=={3:[2,2],4:[3,3],5:[7,10,10],6:[14,26,34,34]}[n],
            'exact rational Krylov rank sequence')
    project=lambda X:mm(I,mm(E,X))
    residual=lambda X:[[x-y for x,y in zip(a,z)] for a,z in zip(X,project(X))]
    PI=mm(P,I);Q=mm(E,PI);B=residual(PI)
    powers=[I];Fs=[eye(b)]
    for k in range(1,9):
        powers.append(mm(P,powers[-1]));Fs.append(mm(E,powers[-1]))
    Ks=[];V=B
    for k in range(7):
        Ks.append(mm(E,mm(P,V)));V=residual(mm(P,V))
    rec=[eye(b),Q]
    for k in range(2,9):
        X=mm(Q,rec[k-1])
        for j in range(k-1): X=add(X,mm(Ks[j],rec[k-2-j]))
        rec.append(X)
    require(rec==Fs,'memory recurrence through length eight')
    K0=Ks[0]
    for j in range(b):
        for k in range(b):
            require(mass[j]*K0[j][k]==sum(d[i]*B[i][j]*B[i][k] for i in range(N)), 'weighted Gram identity')
    require(add(mm(Q,Q),K0)==Fs[2], 'two-step consequence')
    lump=all(not x for row in B for x in row)
    require(lump==all(not x for row in K0 for x in row),'zero defect iff lumpable')
    witness=None
    for i in range(N):
        for j in range(i):
            if blocks[i]==blocks[j] and PI[i]!=PI[j]:
                diff=next(k for k in range(b) if PI[i][k]!=PI[j][k])
                witness={'codes':[nodes[j],nodes[i]],'profile':keys[blocks[i]],
                         'representative_upper_masks':[sum(reps[t][u][v]<<k for k,(u,v) in enumerate(combinations(range(n),2))) for t in (j,i)],
                         'path_counts':[hp(reps[t]) for t in (j,i)],
                         'automorphism_orders':[autom(reps[t]) for t in (j,i)],
                         'source_dart_weights':[d[t] for t in (j,i)],
                         'target_profile':keys[diff],
                         'one_step_probabilities':[str(PI[t][diff]) for t in (j,i)],
                         'projected_rows':[[str(x) for x in PI[t]] for t in (j,i)],
                         'literal_row_node_codes':nodes,
                         'literal_rows':[A[t] for t in (j,i)]}
                break
        if witness:break
    equal_weight_hostile=None
    if n==5:
        ids=[i for i in range(N) if profile(reps[i])==(1,1,2,3,3)]
        require(len(ids)==2,'equal-weight hostile profile fibre')
        require(all(hp(reps[i])==9 and autom(reps[i])==1 and d[i]==9 and
                    canon(reps[i])==canon(converse(reps[i])) for i in ids), 'equal-weight hostile common invariants')
        require(sorted(PI[i][0] for i in ids)==[F(0),F(1,9)],'equal-weight hostile actual return')
        equal_weight_hostile={'codes':[nodes[i] for i in ids],
            'representative_upper_masks':[sum(reps[t][u][v]<<k for k,(u,v) in enumerate(combinations(range(n),2))) for t in ids],
            'common_profile':[1,1,2,3,3],'common_H':9,'common_Aut':1,'common_dart_degree':9,
            'both_self_converse':True,'literal_row_node_codes':nodes,'literal_rows':[A[i] for i in ids],
            'probabilities_to_transitive':[str(PI[i][0]) for i in ids]}
    nonzero=next(((j,k) for j in range(b) for k in range(b) if K0[j][k]),None)
    defect=None
    if nonzero:
        j,k=nonzero
        defect={'source_profile':keys[j],'target_profile':keys[k],
                'true_two_step':str(Fs[2][j][k]),'memoryless_two_step':str(mm(Q,Q)[j][k]),
                'difference':str(K0[j][k])}
    tt=blocks.index(0)
    neighbor=next(i for i in range(N) if A[tt][i])
    tt_control={'neighbor_code':nodes[neighbor], 'neighbor_profile':keys[blocks[neighbor]],
                'neighbor_degree':d[neighbor], 'neighbor_block_mass':mass[blocks[neighbor]],
                'true_two_step':str(Fs[2][0][0]), 'memoryless_two_step':str(mm(Q,Q)[0][0])}
    hidden_receipt=None
    if n==5:
        columns=[]; selected=[]
        for k in range(b):
            ids=[i for i in range(N) if blocks[i]==k]
            for i in ids[:-1]:
                col=[F(0)]*N;col[i]=F(1);col[ids[-1]]=-F(d[i],d[ids[-1]])
                columns.append(col);selected.append(i)
        J=list(map(list,zip(*columns)))
        Dh=[row for i,row in enumerate(residual(mm(P,J))) if i in selected]
        # selected is already increasing here; require this before restriction.
        require(selected==sorted(selected),'hidden coordinate order')
        def determinant(X):
            return sum((-1)**sum(p[i]>p[j] for i in range(len(p)) for j in range(i+1,len(p))) *
                       reduce(lambda a,i:a*X[i][p[i]],range(len(p)),F(1))
                       for p in permutations(range(len(X))))
        den=[sum((-1)**k*determinant([[Dh[i][j] for j in S] for i in S])
                 for S in combinations(range(len(selected)),k)) for k in range(len(selected)+1)]
        require(den==[F(1),-F(6590,186615),-F(16714,186615),F(2126,186615)],'order-five cubic memory denominator')
        for k in range(3,len(Ks)):
            Z=[[sum(den[j]*Ks[k-j][u][v] for j in range(4)) for v in range(b)] for u in range(b)]
            require(all(not x for row in Z for x in row),'cubic kernel recurrence')
        numer=[sum(den[j]*Ks[k-j][0][0] for j in range(k+1)) for k in range(3)]
        require(numer==[F(20735,373230),-F(5340,373230),F(300,373230)],'transitive memory entry numerator')
        hidden_receipt={'matrix':[[str(x) for x in row] for row in Dh],
                        'denominator_coefficients':[str(x) for x in den],
                        'transitive_numerator_coefficients':[str(x) for x in numer]}
    return {'n':n,'merged_nodes':N,'score_blocks':b,'darts':sum(d),
            'source_block_masses':mass,'score_profiles':keys,'lumpable':lump,
            'witness':witness,'first_nonzero_two_step_defect':defect,
            'krylov_rank_sequence':krylov,'equal_weight_hostile':equal_weight_hostile,
            'transitive_return_control':tt_control,'hidden_memory':hidden_receipt,
            'raw_matrix_sha256':hashlib.sha256(json.dumps(A,separators=(',',':')).encode()).hexdigest()}

def uniform_hostile():
    for n in range(5,33):
        T=make(n,0);U=[list(row) for row in T]
        for a,b in combinations((1,2,3),2):U[a][b],U[b][a]=U[b][a],U[a][b]
        def cycles(V):
            return {s for s in combinations(range(n),3) if
                    V[s[0]][s[1]]*V[s[1]][s[2]]*V[s[2]][s[0]] or
                    V[s[0]][s[2]]*V[s[2]][s[1]]*V[s[1]][s[0]]}
        old=cycles(T);new=cycles(U)
        require(old=={(i,i+1,i+2) for i in range(n-2)},'staircase triangle list')
        require(new==(old-{(0,1,2),(2,3,4)})|{(0,1,3),(1,3,4)},'triangle reversal update')
        require(profile(T)==profile(U),'all-order hostile scores')
        td=lambda S:sorted(sum(i in s for s in S) for i in range(n))
        require(td(old)!=td(new),'converse-invariant separator')
    return {'checked_orders':'5..32','changed_triangle':'(1,2,3), zero based',
            'removed_triangles':[[0,1,2],[2,3,4]],'added_triangles':[[0,1,3],[1,3,4]]}

if __name__=='__main__':
    rows=[analyze(n) for n in range(3,7)]
    require([r['merged_nodes'] for r in rows]==[2,3,10,34],'known complete-node controls')
    require([r['lumpable'] for r in rows]==[True,True,False,False],'minimal score-compression obstruction')
    universal_control=uniform_hostile()
    print(json.dumps({'universal_hostile_controls':universal_control,'status':'FINITE-EXACT','universe':'all native fixed-path tilings, n=3..6; no filters',
                      'walk_coefficients_checked':'0..8','rows':rows,'always_active_checks':CHECKS,
                      'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},sort_keys=True,indent=2))

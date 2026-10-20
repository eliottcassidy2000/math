#!/usr/bin/env python3
"""Exact fixed-cube Walsh certificate at n=10, K=C7 and quotient prime q=3.
Universe is every K-fixed tournament/graph; quotient Gamma=C6 x S3.
Global type/evenness interpretation uses the proved ambient-Sylow descent.
No nauty, repository atlas, or external package.
"""
import itertools,json
from collections import Counter
N=10; PAIRS=list(itertools.combinations(range(N),2)); IDX={e:i for i,e in enumerate(PAIRS)}; ID=tuple(range(N))
def need(c,s):
    if not c: raise RuntimeError(s)
def comp(p,q): return tuple(p[q[i]] for i in range(N))
def inv(p): return tuple(p.index(i) for i in range(N))
def prep(p):
    ims=[]; rev=destrev=0
    for k,(i,j) in enumerate(PAIRS):
        a,b=p[i],p[j]; d=IDX[tuple(sorted((a,b)))]; ims.append(1<<d)
        if a>b: rev|=1<<k; destrev|=1<<d
    return p,ims,rev,destrev
def move(x,r,t=False):
    y=0
    while x:
        b=x&-x; y|=r[1][b.bit_length()-1]; x-=b
    return y^(r[3] if t else 0)
def sgn(x,r): return -1 if (x&r[2]).bit_count()%2 else 1
K=[tuple([(i+b)%7 for i in range(7)]+[7,8,9]) for b in range(7)]
GAMMA=[prep(tuple([(a*i)%7 for i in range(7)]+list(r))) for a in range(1,7) for r in itertools.permutations([7,8,9])]
GP={r[0]:r for r in GAMMA}
RC=(8,9,7); A=tuple([(2*i)%7 for i in range(7)]+[7,8,9]); B=tuple(list(range(7))+list(RC))
C3ROOTS=[(7,8,9),RC,(9,7,8)]
Q={tuple([(a*i)%7 for i in range(7)]+list(r)) for a in [1,2,4] for r in C3ROOTS}

def invariant_masks(t):
    seen=set(); oms=[]; base=0
    for i,j in PAIRS:
        if (i,j) in seen:continue
        directed={(p[i],p[j]) for p in K}; orb={tuple(sorted(e)) for e in directed}; seen|=orb
        need(not any((b,a) in directed for a,b in directed),'odd group orientation')
        oms.append(sum(1<<IDX[e] for e in orb))
        base|=sum(1<<IDX[tuple(sorted((a,b)))] for a,b in directed if a>b)
    vals=[]
    for sel in range(1<<len(oms)):
        x=base if t else 0
        for i,o in enumerate(oms):
            if sel>>i&1:x^=o
        vals.append(x)
    return vals,len(oms)
def canonical_group(P):
    return min(tuple(sorted(comp(comp(g[0],p),inv(g[0])) for p in P)) for g in GAMMA)
def records(t):
    vals,d=invariant_masks(t); seen=set(); out=[]
    for x in sorted(vals):
        if x in seen:continue
        copies={}; consistent=True
        for r in GAMMA:
            y=move(x,r,t); sign=1 if t else sgn(x,r)
            if y in copies and copies[y]!=sign:consistent=False
            copies[y]=sign
        seen.update(copies)
        if not consistent:continue
        stab=[r for r in GAMMA if move(x,r,t)==x]
        P={r[0] for r in stab}&Q
        need(len(P) in [1,3,9],'quotient Sylow size')
        a=len(stab)
        while a%3==0:a//=3
        need(len(stab)//a==len(P),'normal Q intersection is Sylow')
        out.append({'mask':x,'stab':len(stab),'copies':copies,'sylow':canonical_group(P)})
    return out,d

def determinant(A):
    A=[r[:] for r in A]; n=len(A); sign=1; prev=1
    for k in range(n-1):
        r=next((r for r in range(k,n) if A[r][k]),None)
        if r is None:return 0
        if r!=k:A[k],A[r]=A[r],A[k];sign=-sign
        pivot=A[k][k]
        for i in range(k+1,n):
            for j in range(k+1,n):
                num=A[i][j]*pivot-A[i][k]*A[k][j]
                need(num%prev==0,'Bareiss exact division');A[i][j]=num//prev
        for i in range(k+1,n):A[i][k]=0
        prev=pivot
    return sign*A[-1][-1] if n else 1

def match(M,p):
    owner={}
    def aug(i,seen):
        for j,v in enumerate(M[i]):
            if not v%p or j in seen:continue
            seen.add(j)
            if j not in owner or aug(owner[j],seen): owner[j]=i;return True
        return False
    need(all(aug(i,set()) for i in range(len(M))),'p-unit perfect matching')
    return sorted((i,j) for j,i in owner.items())

def predicate(x,t,gens):
    return any(all(move(move(x,r,t),GP[g],t)==move(x,r,t) for g in gens) for r in GAMMA)

def closure(gens):
    found={ID}; todo=[ID]
    while todo:
        x=todo.pop()
        for g in gens:
            y=comp(x,g)
            if y not in found:found.add(y);todo.append(y)
    return found

def main():
    # Independent finite witness that the diagonal F21 itself is nonpronormal.
    diag=closure([K[1],comp(A,B)])
    swap=tuple(list(range(7))+[8,7,9])
    other={comp(comp(swap,h),inv(swap)) for h in diag}
    generated=closure(list(diag|other))
    need(len(diag)==len(other)==21 and len(generated)==63,'nonpronormal witness group orders')
    need(not any({comp(comp(g,h),inv(g)) for h in diag}==other for g in generated),'nonpronormal: no conjugator in generated group')
    T,d=records(True);E,e=records(False); need(d==e==9,'fixed-cube dimension')
    need(len(T)==len(E),'Walsh square')
    M=[[sum(eps*(-1 if (r['mask']&y).bit_count()%2 else 1) for y,eps in c['copies'].items()) for c in E] for r in T]
    for i,r in enumerate(T):
        for j,s in enumerate(T):
            need(sum(M[i][k]*c['stab']*M[j][k] for k,c in enumerate(E))==(2**d*r['stab'] if i==j else 0),'weighted orthogonality')
    det=determinant(M); need(det and abs(det)&(abs(det)-1)==0,'dyadic determinant')
    matching=match(M,3)
    for i,j in matching:need(T[i]['sylow']==E[j]['sylow'],'whole quotient Sylow action preserved')
    tests=[('C7',[]),('F21 fixing three roots',[A]),('C7 x C3',[B]),('diagonal F21',[comp(A,B)]),('F21 x C3',[A,B])]
    counts=[]
    for name,gens in tests:
        tp=[predicate(r['mask'],True,gens) for r in T]; ep=[predicate(c['mask'],False,gens) for c in E]
        need(sum(tp)==sum(ep),'mixed-prime predicate counts')
        need(all(tp[i]==ep[j] for i,j in matching),'predicate respected by one matching')
        counts.append({'group':name,'tournament_types':sum(tp),'even_graph_types':sum(ep)})
    # Retain subgroup labels by the presence of the three quotient line classes.
    labels=[]
    for arr,t in [(T,True),(E,False)]:
        hist=Counter()
        for r in arr:
            bits=tuple(predicate(r['mask'],t,g) for _,g in tests[1:4]);hist[str(bits)]+=1
        labels.append(dict(sorted(hist.items())))
    need(labels[0]==labels[1],'complete joint one-quotient-prime profile')
    clean=lambda r:{k:v for k,v in r.items() if k not in ['copies','sylow']}
    print(json.dumps({'status':'PASS','n':10,'base':'C7 on 7 vertices with 3 fixed roots','normalizer_quotient':'C6 x S3','quotient_order':36,'diagonal_F21_nonpronormal_witness':{'subgroup_order':len(diag),'generated_group_order':len(generated),'conjugator_exists_in_generated_group':False,'external_conjugator':list(swap)},'fixed_cube_dimension':d,'fixed_objects_each':2**d,'orbit_types_each':len(T),'walsh_determinant':det,'walsh_determinant_abs_exponent':abs(det).bit_length()-1,'mixed_prime_counts':counts,'joint_line_predicate_histogram':labels[0],'tournament_representatives':list(map(clean,T)),'graph_representatives':list(map(clean,E)),'walsh_matrix':M,'matching_q3':matching},indent=2))
if __name__=='__main__':main()

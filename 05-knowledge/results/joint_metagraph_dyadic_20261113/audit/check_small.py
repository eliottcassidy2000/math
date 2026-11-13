#!/usr/bin/env python3
"""Independent labelled-orbit census and physical-path line enumeration, n<=6."""
from collections import Counter,defaultdict,deque
from fractions import Fraction as F
from itertools import combinations,permutations,product
from math import factorial
import json

gates=0
def need(ok,msg):
    global gates
    gates+=1
    if not ok:raise RuntimeError(msg)

def column_rank(columns):
    basis={}
    for col in columns:
        v=list(map(F,col))
        for p,b in sorted(basis.items()):
            q=v[p]
            if q:v=[x-q*y for x,y in zip(v,b)]
        p=next((i for i,x in enumerate(v) if x),None)
        if p is not None:
            q=v[p];basis[p]=[x/q for x in v]
    return len(basis)

def census(n):
    pairs=list(combinations(range(n),2));idx={e:k for k,e in enumerate(pairs)}
    top=(1<<len(pairs))-1
    def arc(t,a,b):
        return ((t>>idx[min(a,b),max(a,b)])&1)==(a<b)
    def transport(t,p):
        s=0
        for k,(a,b) in enumerate(pairs):
            if arc(t,p[a],p[b]):s|=1<<k
        return s
    generators=[]
    for j in range(n-1):
        p=list(range(n));p[j],p[j+1]=p[j+1],p[j];generators.append(p)
    owner={};orbits={}
    for t in range(top+1):
        if t in owner:continue
        queue=deque([t]);owner[t]=t;orbit=[t]
        while queue:
            u=queue.popleft()
            for p in generators:
                v=transport(u,p)
                if v not in owner:
                    owner[v]=t;queue.append(v);orbit.append(v)
                else:need(owner[v]==t,'orbit collision')
        orbits[t]=orbit
    reps=sorted(orbits);merged=sorted({min(t,owner[top^t]) for t in reps})
    dual={t:owner[top^t] for t in reps};aut={t:factorial(n)//len(orbits[t]) for t in reps}
    mid={t:min(t,dual[t]) for t in reps};midx={t:i for i,t in enumerate(merged)}
    A=[[F(0) for _ in merged] for _ in merged];B=[[F(0) for _ in merged] for _ in merged]
    hist=Counter();fixed=Counter();H={}
    allp=list(permutations(range(n)))
    for t in reps:
        paths=[p for p in allp if all(arc(t,p[j],p[j+1]) for j in range(n-1))]
        H[t]=len(paths);hist[len(paths)]+=1
        need(H[t]%2==1 and aut[t]%2==1,'odd count or stabilizer')
        if t==dual[t]:fixed[len(paths)]+=1
        i=midx[mid[t]]
        for p in paths:
            flip=sum(1<<idx[min(p[j],p[j+1]),max(p[j],p[j+1])] for j in range(n-1))
            u=owner[t^flip];j=midx[mid[u]];A[i][j]+=F(1,aut[t])
            reverse=[None]*n
            for a,b in zip(p,p[::-1]):reverse[a]=b
            if transport(t,reverse)==(top^t):B[i][j]+=F(1,aut[t])
    d=[sum(row) for row in A]
    for i,t in enumerate(merged):
        need(d[i]==(1 if dual[t]==t else 2)*F(H[t],aut[t]),'dart mass')
        for j in range(len(merged)):
            need(A[i][j]==A[j][i] and B[i][j]==B[j][i],'undirected matrix')
            need(A[i][j].denominator==B[i][j].denominator==1,'nonintegral incidence')
        need(A[i][i]%2==B[i][i]%2==0,'loop doubling')
    need(sum(d)==1<<((n-1)*(n-2)//2),'total tiling mass')
    profiles=[]
    for t in merged:
        s=tuple(sorted(sum(arc(t,a,b) for b in range(n) if a!=b) for a in range(n)))
        profiles.append(min(s,tuple(n-1-v for v in s[::-1])))
    pp=sorted(set(profiles));blocks=[[i for i,s in enumerate(profiles) if s==p] for p in pp]
    masses=[sum(d[i] for i in b) for b in blocks]
    P=[[A[i][j]/d[i] for j in range(len(d))] for i in range(len(d))]
    Q=[[sum(A[i][j] for i in a for j in b)/masses[ia] for b in blocks] for ia,a in enumerate(blocks)]
    rows=[[sum(P[i][j] for j in b) for b in blocks] for i in range(len(d))]
    lumpable=all(all(rows[i]==rows[b[0]] for i in b) for b in blocks)
    # Compare actual two-step block observations against the memoryless quotient.
    gram=[]
    for a,ba in enumerate(blocks):
        row=[]
        for b,bb in enumerate(blocks):
            true=sum(d[i]*P[i][j]*P[j][k] for i in ba for j in range(len(d)) for k in bb)/masses[a]
            shortcut=sum(Q[a][j]*Q[j][b] for j in range(len(blocks)))
            defect=true-shortcut
            # A weighted covariance gives exactly the same second-step defect.
            cov=F(0)
            for c,bc in enumerate(blocks):
                for j in bc:cov+=d[j]*(rows[j][a]-Q[c][a])*(rows[j][b]-Q[c][b])
            need(masses[a]*defect==cov,'Gram identity')
            row.append(defect)
        need(row[a]>=0,'negative return defect')
        gram.append(row)
    need(all(gram[a][a]==0 for a in range(len(blocks)))==lumpable,'lumpability iff')
    def act(v):return [sum(P[i][j]*v[j] for j in range(len(d))) for i in range(len(d))]
    def hidden(v):
        means=[sum(d[i]*v[i] for i in b)/masses[a] for a,b in enumerate(blocks)]
        return [v[i]-means[pp.index(profiles[i])] for i in range(len(d))]
    fresh=[[F(i in b) for i in range(len(d))] for b in blocks]
    columns=list(fresh);ranks=[column_rank(columns)]
    for _ in range(len(d)):
        fresh=[act(v) for v in fresh];columns.extend(fresh)
        ranks.append(column_rank(columns))
        if ranks[-1]==ranks[-2]:break
    need(ranks[-1]==len(d),'full small-order walk realization rank')
    hidden_char=None;memory=None;memory_hankel=None
    if n==5:
        coords=[];basis=[]
        for b in blocks:
            for i in b[:-1]:
                v=[F(0) for _ in d];v[i]=1;v[b[-1]]=-d[i]/d[b[-1]]
                coords.append(i);basis.append(v)
        image=[hidden(act(v)) for v in basis]
        C=[[v[i] for v in image] for i in coords]
        need(len(C)==3,'hidden dimension')
        tr=sum(C[i][i] for i in range(3))
        c2=sum(C[i][i]*C[j][j]-C[i][j]*C[j][i] for i,j in combinations(range(3),2))
        det=sum((1 if sum(p[i]>p[j] for i in range(3) for j in range(i+1,3))%2==0 else -1)*C[0][p[0]]*C[1][p[1]]*C[2][p[2]] for p in permutations(range(3)))
        hidden_char=[F(1),-tr,c2,-det]
        need(hidden_char==[F(1),F(-6590,186615),F(-16714,186615),F(2126,186615)],'independent cubic denominator')
    trans=owner[top];ti=midx[mid[trans]];ta=pp.index(profiles[ti])
    ret=sum(P[ti][j]*P[j][ti] for j in range(len(d)))
    coarse=sum(Q[ta][j]*Q[j][ta] for j in range(len(blocks)))
    need(len(blocks[ta])==1,'transitive score block')
    if n==5:
        v=hidden(act([F(i==ti) for i in range(len(d))]));memory=[]
        for _ in range(10):
            memory.append(act(v)[ti]);v=hidden(act(v))
        need(memory[0]==F(1,18),'memory constant')
        for j in range(7):
            need(186615*memory[j+3]==6590*memory[j+2]+16714*memory[j+1]-2126*memory[j],'memory recurrence')
        memory_hankel=sum((1 if sum(p[i]>p[j] for i in range(3) for j in range(i+1,3))%2==0 else -1)*memory[p[0]]*memory[1+p[1]]*memory[2+p[2]] for p in permutations(range(3)))
        need(memory_hankel==F(1,2059408867230),'minimal scalar memory rank')
    witness=None
    for b in blocks:
        for i,j in combinations(b,2):
            if rows[i]!=rows[j]:
                witness={'codes':[merged[i],merged[j]],'H':[H[merged[i]],H[merged[j]]],
                         'aut':[aut[merged[i]],aut[merged[j]]], 'degrees':[int(d[i]),int(d[j])],
                         'profile':profiles[i],'rows':[[str(x) for x in rows[k]] for k in [i,j]]}
                break
        if witness:break
    pairing_guard=None;vertex_automorphisms=None
    if n<=5:
        hgroups=[[t for t in reps if H[t]==h] for h in sorted(hist)]
        abstract=Counter();wired=Counter();colored=Counter()
        for images in product(*(list(permutations(g)) for g in hgroups)):
            J={t:u for g,im in zip(hgroups,images) for t,u in zip(g,im)}
            if any(J[J[t]]!=t or J[dual[t]]!=dual[J[t]] for t in reps):continue
            f=sum(J[t]==t for t in reps);abstract[f]+=1
            bar=[midx[mid[J[t]]] for t in merged]
            if all(A[i][j]==A[bar[i]][bar[j]] for i in range(len(d)) for j in range(len(d))):
                wired[f]+=1
                if all(B[i][j]==B[bar[i]][bar[j]] for i in range(len(d)) for j in range(len(d))):colored[f]+=1
        pairing_guard={'abstract_fixed_histogram':dict(sorted(abstract.items())),
                       'wired_fixed_histogram':dict(sorted(wired.items())),
                       'colored_fixed_histogram':dict(sorted(colored.items()))}
        degree_groups=[[i for i in range(len(d)) if d[i]==q] for q in sorted(set(d))]
        vertex_automorphisms=[]
        for images in product(*(list(permutations(g)) for g in degree_groups)):
            perm={i:j for g,im in zip(degree_groups,images) for i,j in zip(g,im)}
            if all(A[i][j]==A[perm[i]][perm[j]] for i in range(len(d)) for j in range(len(d))):
                vertex_automorphisms.append([merged[perm[i]] for i in range(len(d))])
    return {'n':n,'ordinary':len(reps),'merged':len(merged),'self_converse':sum(fixed.values()),
            'joint':dict(sorted(hist.items())),'fixed_joint':dict(sorted(fixed.items())),
            'score_profiles':pp,'lumpable':lumpable,'transitive_return2':str(ret),
            'score_return2':str(coarse),'first_row_witness':witness,
            'krylov_ranks':ranks,'hidden_characteristic':[str(v) for v in hidden_char] if hidden_char else None,
            'transitive_memory':[str(v) for v in memory] if memory else None,
            'memory_hankel_determinant':str(memory_hankel) if memory_hankel is not None else None,
            'pairing_guard':pairing_guard,'weighted_vertex_automorphisms':vertex_automorphisms,
            'node_metadata':[{'code':t,'H':H[t],'aut':aut[t],'degree':int(d[i]),'self_converse':dual[t]==t,
                              'neighbors':{str(merged[j]):int(A[i][j]) for j in range(len(d)) if A[i][j]}}
                             for i,t in enumerate(merged)] if n<=5 else None,
            'degree_histogram':dict(sorted(Counter(map(int,d)).items()))}

def main():
    rows=[census(n) for n in range(3,7)]
    need([r['ordinary'] for r in rows]==[2,4,12,56],'A000568 controls')
    need([r['self_converse'] for r in rows]==[2,2,8,12],'converse controls')
    need([r['lumpable'] for r in rows]==[True,True,False,False],'first lumpability boundary')
    need(rows[2]['transitive_return2']=='1/9' and rows[2]['score_return2']=='1/18','n5 return witness')
    need(rows[2]['krylov_ranks']==[7,10,10] and rows[3]['krylov_ranks']==[14,26,34,34],'minimal realization dimensions')
    need([len(r['weighted_vertex_automorphisms']) for r in rows[:3]]==[2,1,1],'small weighted vertex automorphism groups')
    need([min(r['pairing_guard']['abstract_fixed_histogram']) for r in rows[:3]]==[2,2,4],'abstract pairing minima')
    need([min(r['pairing_guard']['wired_fixed_histogram']) for r in rows[:3]]==[2,2,8],'wired pairing minima')
    need(all(r['pairing_guard']['wired_fixed_histogram']==r['pairing_guard']['colored_fixed_histogram'] for r in rows[:3]),'color guard at small orders')
    print(json.dumps({'status':'PASS','gates':gates,'universe':'all labelled tournaments n3..6, partitioned by adjacent-transposition orbits; every physical Hamilton path', 'rows':rows},indent=2,sort_keys=True))
if __name__=='__main__':main()

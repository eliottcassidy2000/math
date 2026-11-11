"""Exact two-parameter repeated-gadget family and five-dimensional interface proof."""
from pathlib import Path
from math import isqrt
from itertools import combinations,permutations
import json
from classify_modules import H,strong,need,inherited as a,two_birth as b,converse
HERE=Path(__file__).resolve().parent

def family(k,t):
    rows=a.ready_source((2,1),k)
    for _ in range(t):
        n=len(rows)
        rows=a.attach(a.attach(a.attach(rows,0),0),(1<<k)|(1<<n))
    return rows

def fixed_kernel(t):
    rows=a.induced(a.ready_source((2,1),1),[1,2,3,4])
    for _ in range(t):
        n=len(rows)
        rows=a.attach(a.attach(a.attach(rows,0),0),1|(1<<n))
    return rows

def mm(A,B):return [[sum(A[i][z]*B[z][j] for z in range(3)) for j in range(3)] for i in range(3)]
def mp(A,t):
    out=[[int(i==j) for j in range(3)] for i in range(3)]
    for _ in range(t):out=mm(out,A)
    return out

def marked_runs(rows,module,cut):
    """Brute placeholder enumeration, independent of subset path DP.

    A separator distinguishes two ordered bins. The first starts in cut;
    the second may be empty. Also count one-path words ending in cut.
    """
    others=[v for v in range(len(rows)) if v!=module];n=len(others)
    cover=[];endpoint=[]
    for r in range(1,n+3):
        cv=ep=0;length=n+r
        for positions in combinations(range(length),r):
            ps=set(positions)
            for perm in permutations(others):
                it=iter(perm);word=[module if i in ps else next(it) for i in range(length)]
                bad=[i for i in range(length-1) if not(rows[word[i]]>>word[i+1]&1)]
                if not bad and cut>>word[-1]&1:ep+=1
                if not(cut>>word[0]&1):continue
                for split in range(1,length+1):
                    if all(i==split-1 for i in bad):cv+=1
        cover.append(cv);endpoint.append(ep)
    return cover,endpoint

def formulas(k,t):
    D=1<<k
    if not t:return {'H':14*D+3,'A':14,'B':3,'q':5,'C':2,'residual_H':1}
    A=(19*12**t+(4*t+22)*3**t)//3
    B=(8*12**t+(2*t+17)*3**t)//9
    q=(24*12**t+(6*t+21)*3**t)//9
    C=(16*12**t+4*(t+1)*3**t)//9
    residual=(5*12**t+(2*t+5)*3**t)//9
    return {'H':A*D+B,'A':A,'B':B,'q':q,'C':C,'residual_H':residual}

def main():
    BA=a.ready_source((2,1),1);ids=[v for v in range(len(BA)) if v!=1]
    X=a.induced(BA,ids);cut=sum((BA[1]>>v&1)<<j for j,v in enumerate(ids))
    cover,endpoint=marked_runs(X,0,cut)
    need(cover==[32,15,0,0,0,0] and endpoint==[6,3,0,0,0,0],'all-k marked initial interface')
    need(b.module_runs(X)[1:]==[13,6,0,0,0],'initial X count interface')
    plus=[[3,6,2],[0,12,3],[0,0,3]]
    minus=[[3,-2,0],[0,0,-1],[0,0,3]]
    matrices=[]
    for t in range(1,21):
        P,N=mp(plus,t),mp(minus,t)
        need(P[1][2]==(12**t-3**t)//3 and P[0][2]==2*(12**t-3**t)//9,'positive cell powers')
        need(N[1][2]==-3**(t-1),'negative cell endpoint')
        need(N[0][2]==(0 if t==1 else 2*(t-1)*3**(t-2)),'negative joint channel')
        for k in range(1,11):
            D=1<<k;d=6*D+1;op=15*D+2;om=-3*D;h0=14*D+3
            value=h0*3**t+d*(P[0][2]+N[0][2])+op*P[1][2]+om*N[1][2]
            need(value==formulas(k,t)['H'],'matrix composition closed form')
        matrices.append({'t':t,'plus':P,'minus':N})
    kernels=[]
    for t in range(5):
        W=fixed_kernel(t);hp,en=a.dp(W);_,st=a.dp(converse(W));full=(1<<len(W))-1;J=(2,3)
        E=sum(sum(en[S][r] for r in J)*hp[full^S] for S in range(1,full+1))
        C=sum(en[S][i]*st[full^S][j] for S in range(1,full) for i in J for j in J if i!=j)
        residual=a.induced(W,[v for v in range(len(W)) if v not in J])
        f=formulas(1,t)
        need((hp[-1],E,C,H(residual))==(f['q'],f['A'],f['C'],f['residual_H']),'independent kernel DP')
        kernels.append({'t':t,'rows':W,'q':hp[-1],'A':E,'C':C,'residual_H':H(residual)})
    literal=[]
    for k in range(1,5):
        for t in range(3):
            rows=family(k,t);f=formulas(k,t);need(H(rows)==f['H'] and strong(rows),'literal family H and strongness')
            prof=a.root_profiles(rows)[k+4]
            need((prof['s'],prof['q'],prof['gap'],prof['delta'])==(f['q'],f['q'],f['C'],f['C']),'literal root profile')
            candidates=[]
            if t:
                for side,R in [('forward',rows),('converse',converse(rows))]:
                    profiles=a.root_profiles(R)
                    need(not any(p['gap'] in(1,2) for p in profiles),'no ready root after arbitrary re-marking/converse')
                    candidates.extend({'side':side,**p} for p in profiles if len(R)-1-R[p['root']].bit_count()<=2)
            # Inflate one already free prefix vertex: preparation and gadget commute.
            need(b.inflate(rows,0,2)==family(k+1,t),'commuting preparation and cell growth')
            literal.append({'k':k,'t':t,'H':f['H'],'profile':prof,'degree_eligible_roots':candidates})
    for k in range(1,5):
        hs=[formulas(k,t)['H'] for t in range(1,13)]
        need(all(hs[i+3]==18*hs[i+2]-81*hs[i+1]+108*hs[i] for i in range(len(hs)-3)),
             'order-three recurrence')
    for k in range(1,41):
        need((formulas(k,0)['H']-1)%4==2,'t0 nonsquare residue')
        need(formulas(k,1)['H']%17==0,'t1 fixed prime factor')
        z=formulas(k,2)['H']-1
        need(isqrt(z)**2!=z,'t2 explicit square obstruction control')
        for t in range(3,21):need(formulas(k,t)['H']%3==0,'higher-cell fixed factor')
    report={'status':'PASS','initial_X_rows':X,'cut':cut,'initial_cover_runs':cover,
            'initial_endpoint_runs':endpoint,'cell_plus':plus,'cell_minus':minus,
            'matrix_controls':matrices,'fixed_kernels':kernels,'literal_controls':literal,
            'coefficient_sequences':[{'t':t,**formulas(1,t)} for t in range(13)]}
    (HERE/'composition.json').write_text(json.dumps(report,sort_keys=True,indent=2)+'\n')
    print('PASS: repeated-gadget family, marked initial interface, independent kernels and commuting growth')
if __name__=='__main__':main()

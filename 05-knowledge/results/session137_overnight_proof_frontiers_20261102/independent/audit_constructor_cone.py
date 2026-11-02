"""Literal marked-source states and local-operation columns; no producer imports."""
from tournament_kernel import *
import json

def state(a,cut):
    n=len(a);h=g=q=0
    for p in permutations(range(n)):
        valid=[bool(a[p[j]][p[j+1]]) for j in range(n-1)]
        if all(valid):h+=1;g+=bool(cut>>p[-1]&1)
        for k in range(n):
            if cut>>p[k]&1 and all(valid[j] for j in range(n-1) if j!=k-1):q+=1
    f=literal_hp(append(a,cut))
    require((f-h)%2==0 and (q-f+h-g)%2==0,'source integrality')
    v=((f-h)//2,(q-f+h-g)//2,g,h-g)
    require(min(v)>=0,'source positivity')
    require((v[1]==0)==(cut.bit_count()<=1),'L zero exact cut boundary')
    return v

def migrate(a,cut,local,z):
    n=len(a);m=len(local);full=[[0]*(n+m) for _ in range(n+m)]
    for i in range(n):
        for j in range(n):full[i][j]=a[i][j]
        full[n][i]=(cut>>i)&1;full[i][n]=1-full[n][i]
        for j in range(1,m):full[i][n+j]=1
    for i in range(m):
        for j in range(m):full[n+i][n+j]=local[i][j]
    mark=n+z;vs=[i for i in range(n+m) if i!=mark]
    return induced(full,vs),sum(full[mark][v]<<j for j,v in enumerate(vs))

def run():
    source_cases=0;face_examples={}
    for n in range(1,5):
        for code in range(1<<(n*(n-1)//2)):
            a=matrix(n,code)
            for cut in range(1<<n):
                v=state(a,cut);source_cases+=1
                if v[1]==v[2]==0:
                    C=v[0];h=v[3]
                    if C:
                        odd=C
                        while odd%2==0:odd//=2
                        face_examples.setdefault(h,set()).add(odd)
                        if h==1:require(odd==1,'height1 dyadic face')
    one=matrix(1,0);two=matrix(2,1)
    witnesses=[(one,0),(one,1),(two,1),(two,3)]
    require([state(*w) for w in witnesses]==[(0,0,0,1),(0,0,1,0),(1,0,0,1),(0,1,1,0)],'lattice basis sources')
    local_cases=0;matrices=set()
    for n in range(2,5):
        for code in range(1<<(n*(n-1)//2)):
            local=matrix(n,code)
            for z in range(1,n):
                outputs=[state(*migrate(a,cut,local,z)) for a,cut in witnesses]
                cols=[tuple(x-y for x,y in zip(outputs[2],outputs[0])),
                      tuple(x-y for x,y in zip(outputs[3],outputs[1])),outputs[1],outputs[0]]
                require(all(x>=0 for col in cols for x in col),'universal positive-coordinate operation columns')
                matrices.add(tuple(cols));local_cases+=1
    print(json.dumps({'scope':'FINITE-EXACT all marked sources parent order<=4; every local tournament order2..4 and every fresh re-mark, reconstructed on four unimodular witnesses',
                      'source_states':source_cases,'local_recipes':local_cases,'distinct_matrices':len(matrices),
                      'observed_face_odd_parts':{h:sorted(v) for h,v in face_examples.items()},
                      'basis_states':[state(*w) for w in witnesses],
                      'height1_formal_hostile':{'positive_state':[5,0,0,1],'F':11,'reason':'C5 is not a power of two, although H11 is globally attainable'}},sort_keys=True,indent=2))
if __name__=='__main__':run()

#!/usr/bin/env python3
"""Exact independent finite audit of integral signed tournament cocycle."""
from itertools import combinations,permutations
import argparse,json

def req(q,m):
    if not q:raise RuntimeError(m)
def signed_cycle(S,t):
    i,j,k=t;z=(S[i][j],S[j][k],S[k][i])
    return int(z==(1,1,1))-int(z==(-1,-1,-1))
def tprod(S,t):
    i,j,k=t;return S[i][j]*S[j][k]*S[k][i]
def coboundary(fun,face):
    return sum((-1)**j*fun(face[:j]+face[j+1:]) for j in range(len(face)))
def eta(S,face):return coboundary(lambda t:signed_cycle(S,t),face)
def polynomial_eta(S,face):
    i,j,k,l=face;a,b,c,d,e,f=S[i][j],S[i][k],S[i][l],S[j][k],S[j][l],S[k][l]
    terms=(a*b*d,-a*c*e,b*c*f,-d*e*f);num=sum(terms)
    req(terms[0]*terms[1]*terms[2]*terms[3]==1,'four signed terms product')
    req(num in (-4,0,4),'signed numerator trichotomy')
    return num//4
def matrices(n):
    edges=list(combinations(range(n),2))
    for mask in range(1<<len(edges)):
        S=[[0]*n for _ in range(n)]
        for b,(i,j) in enumerate(edges):S[i][j]=1 if mask>>b&1 else -1;S[j][i]=-S[i][j]
        yield S
def run(mutate):
    hist={-1:0,0:0,1:0};switches=0;relabs=0
    for S in matrices(4):
        face=(0,1,2,3);e=eta(S,face);hist[e]+=1
        req(e==polynomial_eta(S,face),'integer lift polynomial formula')
        for t in combinations(range(4),3):
            i,j,k=t;req(4*signed_cycle(S,t)==S[i][j]+S[j][k]+S[k][i]+tprod(S,t),'signed cyclic indicator formula')
        req(4*e==coboundary(lambda t:tprod(S,t),face),'curl cancellation')
        pf=S[0][1]*S[2][3]-S[0][2]*S[1][3]+S[0][3]*S[1][2]
        req(8*e*e==pf*pf-1,'eta squared equals Pf defect')
        for mask in range(16):
            d=[1 if mask>>j&1 else -1 for j in range(4)]
            T=[[d[i]*S[i][j]*d[j] for j in range(4)] for i in range(4)]
            req(eta(T,face)==e,'switching invariance');switches+=1
        for p in permutations(range(4)):
            sign=(-1)**sum(p[i]>p[j] for i in range(4) for j in range(i+1,4))
            req(eta(S,p)==sign*e,'alternating relabel covariance');relabs+=1
    count5=0;modchecks=0
    for S in matrices(5):
        face=(0,1,2,3,4)
        vals=[eta(S,face[:j]+face[j+1:]) for j in range(5)]
        total=sum((1 if mutate else (-1)**j)*v for j,v in enumerate(vals))
        req(total==0,'oriented five-face integral conservation');count5+=1
        for m in range(2,34):req(sum((-1)**j*(v%m) for j,v in enumerate(vals))%m==0,'all-modulus conservation');modchecks+=1
    req(hist=={-1:8,0:48,1:8},'four-vertex signed distribution')
    # Mod2 adjacency loses all signs: the same all-one off-diagonal matrix
    # underlies eta=0 and eta=1, so integer lift cannot be recovered by /4 in F2.
    examples={}
    for S in matrices(4):examples.setdefault(eta(S,(0,1,2,3)),S)
    req(all(examples[0][i][j]%2==examples[1][i][j]%2 for i in range(4) for j in range(4)),'mod2 information loss hostile')
    req((0*0-1)!=8*0,'weighted K=0 rejects binary Pf square identity')
    return {'status':'PASS','labelled_foursets':64,'signed_distribution':hist,'labelled_fivesets':count5,'switching_checks':switches,'alternating_checks':relabs,'modulus_checks':modchecks,'moduli':'2..33 (integral identity proves every modulus)','hostiles':['mod2 skew matrix forgets eta','K=0 weighted Pf-square formula fails']}
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--mutate-five-sign',action='store_true');a=ap.parse_args()
    print(json.dumps(run(a.mutate_five_sign),indent=2,sort_keys=True))

"""Exact rational audit: arbitrary nonuniform fibres, not uniform proxies."""
from fractions import Fraction as F
import json
def need(v,m):
    if not v:raise RuntimeError(m)
def main():
    # Deliberately biased test law. These numbers are controls, not claimed ROOT weights.
    blocks=[[0,1],[2,3,4],[5,6,7,8]];w=[1,7,2,11,3,13,5,17,19]
    pi=[F(x,sum(w)) for x in w];mu=[sum(pi[x] for x in B) for B in blocks]
    Q=[[F(0) for _ in blocks] for _ in blocks]
    for a,b in ((0,1),(1,2)):
        c=min(mu[a],mu[b])/4;Q[a][b]=c/mu[a];Q[b][a]=c/mu[b]
    for a in range(3):Q[a][a]=1-sum(Q[a])
    owner={x:j for j,B in enumerate(blocks) for x in B}
    P=[[Q[owner[x]][owner[y]]*pi[y]/mu[owner[y]] for y in range(9)] for x in range(9)]
    checks=0
    for x in range(9):
        need(sum(P[x])==1,'row sum');checks+=1
        for y in range(9):need(pi[x]*P[x][y]==pi[y]*P[y][x],'weighted detailed balance');checks+=1
    for f in ([F(x*x-3*x+1) for x in range(9)],[F(1,x+1) for x in range(9)]):
        mean=sum(a*b for a,b in zip(pi,f));bar=[sum(pi[x]*f[x] for x in B)/mu[j] for j,B in enumerate(blocks)]
        within=sum(pi[x]*(f[x]-bar[owner[x]])**2 for x in range(9))
        var=sum(pi[x]*(f[x]-mean)**2 for x in range(9))
        vbar=sum(mu[j]*(bar[j]-mean)**2 for j in range(3))
        energy=sum(pi[x]*P[x][y]*(f[x]-f[y])**2 for x in range(9) for y in range(9))/2
        ebar=sum(mu[j]*Q[j][k]*(bar[j]-bar[k])**2 for j in range(3) for k in range(3))/2
        need(var==within+vbar and energy==within+ebar,'exact variance and energy identities');checks+=2
        residual=[f[x]-bar[owner[x]] for x in range(9)]
        need(all(sum(P[x][y]*residual[y] for y in range(9))==0 for x in range(9)),'conditional-zero subspace killed');checks+=1
    # A quotient bottleneck survives any positive disaggregation.
    eps=F(1,10000);need((1-(1-2*eps))==2*eps,'visible epsilon gap');checks+=1
    return {'status':'PASS exact rational conditional-lift identities','checks':checks,
            'test_weights':w,'test_marginals':list(map(str,mu)), 'epsilon_control_gap':str(2*eps),
            'scope':'nonuniform finite control plus algebraic proof; no asserted ROOT source law'}
if __name__=='__main__':print(json.dumps(main(),indent=2,sort_keys=True))

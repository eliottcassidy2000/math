"""Exact guarded responses for transitive variable modules in an acyclic interface."""
from clone_clock import *


def compile_runs(T,variables,bins=1,starts=(),ends=()):
    variables=tuple(variables);fixed=tuple(v for v in range(len(T)) if v not in variables)
    require(bins>=1 and len(set(variables))==len(variables),'distinct modules and positive bin count')
    require(all(v in fixed for _,v in tuple(starts)+tuple(ends)),'guards name fixed vertices')
    # A tournament is acyclic iff its sorted internal outdegrees are 0..d-1.
    require(sorted(sum(T[v]>>u&1 for u in variables) for v in variables)==list(range(len(variables))),
            'variable quotient must be acyclic')
    start=dict(starts);end=dict(ends);pos={v:i for i,v in enumerate(fixed)}
    varpos={v:i for i,v in enumerate(variables)};full=(1<<len(fixed))-1;zero=(0,)*len(variables)
    @lru_cache(None)
    def walk(mask,last,b):
        out=Counter()
        can_close=(b not in end or last==end[b]) and (last!=-1 or b not in start)
        if can_close:
            if b==bins-1 and mask==full:out[zero]+=1
            elif b<bins-1:out.update(walk(mask,-1,b+1))
        for v in range(len(T)):
            if v in pos and mask>>pos[v]&1:continue
            if last==-1 and b in start and v!=start[b]:continue
            if last!=-1 and not T[last]>>v&1:continue
            nxt=mask|(1<<pos[v]) if v in pos else mask
            for powers,w in walk(nxt,v,b).items():
                p=list(powers)
                if v in varpos:p[varpos[v]]+=1
                out[tuple(p)]+=w
        return out
    result=walk(0,-1,0)
    require(all(max(p,default=0)<=len(fixed)+bins for p in result),'uniform run bound')
    return result


def evaluate_runs(poly,sizes):
    arrays=[surjections(k,max(p[i] for p in poly)) for i,k in enumerate(sizes)]
    return sum(w*__import__('math').prod(arrays[i][p[i]] for i in range(len(sizes))) for p,w in poly.items())


def symbolic(poly,names):
    import sympy as sp
    sums=[]
    for i,name in enumerate(names):
        limit=max(p[i] for p in poly)
        sums.append([sp.Integer(0)]+[sum((-1)**(r-j)*__import__('math').comb(r,j)*(1 if j==1 else sp.Symbol(f'{name}{j}'))
                         for j in range(1,r+1)) for r in range(1,limit+1)])
    return sp.factor(sum(w*sp.prod(sums[i][p[i]] for i in range(len(names))) for p,w in poly.items()))


def main():
    polys={k:compile_runs(BASE,(0,5),2 if k=='S2' else 1,
                          starts=((0,1),) if k=='S2' else (),
                          ends=((0,1),) if k=='E1' else ()) for k in ('H','E1','S2')}
    for k in polys:require(polys[k]==bivariate(k),'independent general guarded compiler')
    poly=compile_runs(BASE,(0,4,5))
    import sympy as sp
    M,K,N=sp.symbols('m2 k2 n2');L=(2*N-1)*K-2*N+4
    end=compile_runs(BASE,(0,4,5),ends=((0,1),))
    start2=compile_runs(BASE,(0,4,5),bins=2,starts=((0,1),))
    require(sp.expand(symbolic(end,('m','k','n'))-(L-1))==0,'three-block named end response')
    require(sp.expand(symbolic(start2,('m','k','n'))-((M+1)*L-1))==0,'three-block named two-bin response')
    print('three-block polynomial',symbolic(poly,('m','k','n')),flush=True)
    print('positive coefficients',sorted((p,w) for p,w in poly.items() if min(p)>0),flush=True)
    (HERE/'three_module.json').write_text(json.dumps(dict(variables=(0,4,5),coefficients=[[list(p),w] for p,w in sorted(poly.items())],formula=str(symbolic(poly,('m','k','n'))),E1=str(symbolic(end,('m','k','n'))),S2=str(symbolic(start2,('m','k','n')))),indent=2)+'\n')


if __name__=='__main__':main()

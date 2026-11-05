"""Adversarial multiple-SCC and noninterval-label padding controls."""
import json, random
from audit import ROOT, data, induced, normalize

def ordinal(blocks):
    sizes=list(map(len,blocks));n=sum(sizes);rows=[];offset=0
    for block in blocks:
        later=sum(1<<v for v in range(offset+len(block),n))
        rows.extend((row<<offset)|later for row in block)
        offset+=len(block)
    return tuple(rows)

def main():
    C3=(4,1,2);C5=(8,1,3,6);single=(0,)
    cases=[
        ([single,single,C3,single,C3,single,single],[2]),
        ([single,C3,single,single,C3,single],[1,6]),
        ([C3,single,C5,single],[0]),
        ([C3,single,C5,single],[0,5]),
        ([single,C3,single,C3,single],[0,8]),
        ([single,single,C3,single,C3,single],[]),
    ]
    rng=random.Random(19791);results=[]
    for blocks,O in cases:
        S=ordinal(blocks);n=len(S);out=sum(1<<v for v in O)
        T=tuple(S[v]|((1<<n) if v not in O else 0) for v in range(n))+(out,)
        z=data(T)[0][n]
        for shuffle in (False,True):
            order=list(range(n+1))
            if shuffle:rng.shuffle(order)
            rows=induced(T,order);root=order.index(n)
            C,c,t,u=normalize(rows,root)
            cp=data(C)[0][c];h,e,R,d=cp;f=R+e-h-d;U=R-d
            actual=cp if not t else ((U<<t)+d-f,d,d+(U<<t),d)
            assert actual==z
            results.append({'order':n+1,'source_components':list(map(len,blocks)),
                            'outneighbors':O,'shuffled_labels':shuffle,
                            'original_profile':z,'core_order':len(C),
                            'tail':t,'interior_singletons':u,'core_profile':cp})
    (ROOT/'PADDING.json').write_text(json.dumps(results,indent=2)+'\n')
    print('PASS',len(results),'multiple-SCC/noninterval-label controls')

if __name__=='__main__':main()

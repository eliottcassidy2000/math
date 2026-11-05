"""Independent checks. Run from repo root with Python 3; no third-party dependencies."""
import itertools,json,random,hashlib
from transfer import automaton,advance,graph,held_karp,count
from explore import closure
from pathlib import Path

def valid_words(w,n):
    alphabets=[range(0,1<<w,1<<max(0,w-j)) for j in range(n)]
    return itertools.product(*alphabets)

def exhaustive(w,n):
    states,mat=automaton(w); final=[i for i,s in enumerate(states) if len(s)==1]
    checked=0
    def rec(v,word):
        nonlocal checked
        if len(word)==n:
            h=sum(v[i] for i in final)
            assert h==held_karp(graph(word,w)),(w,word,h)
            checked+=1;return
        j=len(word)
        for mask in range(0,1<<w,1<<max(0,w-j)):
            rec(advance(v,mat,mask),word+[mask])
    rec((1,)+(0,)*(len(states)-1),[])
    return checked

if __name__=='__main__':
    results={}
    for w,n in [(1,7),(2,7),(3,6)]:
        results[f'exhaustive_span{w}_order{n}']=exhaustive(w,n)
    rng=random.Random(197)
    for w in range(1,5):
        for _ in range(50):
            word=[rng.randrange(1<<w) for i in range(rng.randrange(1,11))]
            assert count(word,w)==held_karp(graph(word,w))
    results['random_held_karp_checks']=200
    h=[1,1,1]
    for n in range(3,41):h.append(sum(h[-3:]))
    for n in range(1,41):assert count([1]*n,1)==h[n]
    results['width1_tribonacci_verified_through_order']=40
    for w,M,expected in [(1,101,None),(2,47,[7,21,47]),(2,1001,[7,21,47]),(3,101,[7,21])]:
        r=closure(w,M)
        assert r['complete']
        gaps=[i for i in range(1,M+1,2) if i not in r['counts']]
        if expected is not None: assert gaps==expected
        assert all(i%2 for i in r['counts'])
        for h,word in r['counts'].items():assert count(list(map(int,word)),w)==h
        saved_path=Path(__file__).with_name(f'closure-{w}-{M}.json')
        if saved_path.exists():
            saved=json.loads(saved_path.read_text())
            assert saved['complete'] and saved['states']==r['states']
            assert {int(h):word for h,word in saved['counts'].items()}==r['counts']
            for h,word in saved['counts'].items():assert count(list(map(int,word)),w)==int(h)
        results[f'closure_span{w}_cap{M}']={'states':r['states'],'gaps':gaps}
    word=list(map(int,'0001321'))
    assert held_karp(graph(word,3))==47
    a=graph(word,3)
    reach=[{i}|{j for j in range(len(a)) if a[i][j]} for i in range(len(a))]
    for k in range(len(a)):
        for i in range(len(a)):
            if k in reach[i]:reach[i]|=reach[k]
    assert all(len(r)==len(a) for r in reach)
    results['47_witness']={'word':'0001321','order':7,'strong':True,'backward_edges_1based':[[i+1,j+1] for j in range(len(a)) for i in range(j) if a[j][i]]}
    print(json.dumps(results,indent=2))
    with open('session197/ribbon/VALIDATION.json','w') as f:json.dump(results,f,indent=2)

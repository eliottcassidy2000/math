from itertools import combinations,permutations
from pathlib import Path
import json,hashlib
from independent_tree_capacity import trees_to,canon,automorphisms

def check(x,msg):
    if not x:raise RuntimeError(msg)

def induced(n,edges,subset):
    at={v:i for i,v in enumerate(subset)}
    es=tuple((at[a],at[b]) for a,b in edges if a in at and b in at)
    return es if len(es)==len(subset)-1 else None

def has_asymmetric_subtree(n,edges,asym):
    for k in range(2,n+1):
        for subset in combinations(range(n),k):
            ee=induced(n,edges,subset)
            if ee is not None and canon(k,ee) in asym[k]:return True
    return False

def count(edges,mask,host):
    arcs=[(b,a) if mask>>i&1 else (a,b) for i,(a,b) in enumerate(edges)]
    return sum(all(host[p[a]]>>p[b]&1 for a,b in arcs) for p in permutations(range(len(host))))

def main():
    source=Path(__file__).with_name('minimality_witnesses.json')
    data=json.loads(source.read_text());trees=trees_to(8)
    asym={n:{c for c,e in ts.items() if len(automorphisms(n,e))==1} for n,ts in trees.items()}
    needs={n:{c for c,e in ts.items() if has_asymmetric_subtree(n,e,asym)} for n,ts in trees.items()}
    check(all(not needs[n] for n in range(1,7)),'unexpected small obstruction')
    result=[]
    for row in data['rows']:
        n=row['n'];have=set()
        for w in row['witnesses']:
            e=tuple(tuple(x) for x in w['edges']);m=w['orientation_mask'];h=w['host_out_bits']
            check(all(0<=x<(1<<n) for x in h),'host bits range')
            check(all(not (h[i]>>i&1) for i in range(n)),'host loops')
            check(all((h[i]>>j&1)+(h[j]>>i&1)==1 for i in range(n) for j in range(i+1,n)),'host tournament')
            trans=[sum(1<<j for j in range(i+1,n)) for i in range(n)]
            a=count(e,m,trans);b=count(e,m,h)
            check(a==w['transitive_count'] and b==w['host_count'],'count mismatch')
            check((a-b)%2,'no parity difference');have.add(canon(n,e))
            result.append({'n':n,'transitive_embeddings':a,'other_embeddings':b})
        check(have==needs[n],'incomplete nonhereditary coverage')
    p5=((0,1),(1,2),(2,3),(3,4));arcs=((1,0),(1,2),(2,3),(3,4))
    aut=automorphisms(5,p5)
    check(sum(set((p[a],p[b]) for a,b in arcs)==set(arcs) for p in aut)==1,'P5 oriented asymmetry')
    spectrum={}
    pairs=list(combinations(range(5),2))
    for mask in range(1<<10):
        host=[0]*5
        for i,(a,b) in enumerate(pairs):
            if mask>>i&1:host[a]|=1<<b
            else:host[b]|=1<<a
        c=count(p5,1,host);check(c%2==0,'P5 parity failed')
        spectrum[c]=spectrum.get(c,0)+1
    selfhost=[0]*5
    for a,b in arcs:selfhost[a]|=1<<b
    check(count(p5,1,selfhost)==1,'P5 arbitrary-host boundary')
    print(json.dumps({'status':'PASS','minimality_witnesses':result,'nonhereditary_trees_by_order':{n:len(v) for n,v in needs.items()},'P5_tournament_embedding_spectrum':spectrum,'P5_selfhost_embeddings':1,'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'scope':'Minimality certificate complete through order eight; P5 universal tournament-even conclusion is finite-exact at its fixed order five and does not assert evenness in arbitrary directed hosts.'},indent=2))
if __name__=='__main__':main()

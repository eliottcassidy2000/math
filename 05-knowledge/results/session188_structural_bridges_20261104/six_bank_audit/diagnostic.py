"""Independent exact degree classes and exhaustive original-module subset audit."""
import hashlib
import itertools
import json
from collections import Counter
from pathlib import Path

SOURCE=Path(__file__).resolve().parent.parent/'covered_laminar'
def need(p,msg):
    if not p:raise RuntimeError(msg)
def main():
    raw=(SOURCE/'mod4_six_vectors.bin').read_bytes()
    need(len(raw)==15859712 and max(raw)<4,'canonical byte bank')
    N=32768
    basis=[a+(c<<6) for c in range(64) for a in range(64) if a.bit_count()<=2 and c.bit_count()<=2]
    expected=json.loads((SOURCE/'degree_diagnostic.json').read_text())
    summaries=[];pairs=[]
    for d in range(1,5):
        indices=[i for i,m in enumerate(basis) if m.bit_count()<=d]
        # Sort full immutable vectors, rather than relying on hash signatures.
        ordered=sorted((bytes(raw[i*N+code] for i in indices),code) for code in range(N))
        groups=[sorted(code for _,code in g) for _,g in itertools.groupby(ordered,key=lambda item:item[0])]
        summaries.append(dict(maximum_degree=d,coordinates=len(indices),classes=len(groups),class_size_histogram={str(k):v for k,v in sorted(Counter(map(len,groups)).items())}))
        if d==2:
            pairs=sorted(g for g in groups if len(g)>1)
            need(all(len(g)==2 for g in pairs),'quadratic class cardinality')
    need(summaries==expected['degree_classes'],'every degree class statistic')
    need(pairs==sorted(expected['quadratic_collision_pairs']),'all exact quadratic pairs')
    edges=list(itertools.combinations(range(6),2));categories=Counter();records=[]
    for a,b in pairs:
        def arc(u,v):
            i,j=sorted((u,v));bit=edges.index((i,j));return bool(a>>bit&1)==(u<v)
        modules=[]
        for M in range(1,64):
            names=[i for i in range(6) if M>>i&1]
            if len(names)<2:continue
            if all(len({arc(v,u) for u in names})==1 for v in range(6) if not M>>v&1):modules.append(M)
        effects=[sum(1<<i for i,(u,v) in enumerate(edges) if M>>u&1 and M>>v&1) for M in modules]
        solutions=[]
        for bits in range(1<<len(modules)):
            effect=0
            for i,e in enumerate(effects):
                if bits>>i&1:effect^=e
            if effect==a^b:solutions.append([M for i,M in enumerate(modules) if bits>>i&1])
        need(len(solutions)<=1,'unique original-module representation')
        if not solutions:categories['outside_original_module_span']+=1;continue
        F=solutions[0]
        need(all(not(M&L) or M&L in (M,L) for M,L in itertools.combinations(F,2)),'laminar representation')
        maxima=[M for M in F if not any(M!=L and M&L==M for L in F)]
        for M in maxima:
            smaller=[L for L in F if L!=M and L&M==L]
            children=[L for L in smaller if not any(L!=K and L&K==L for K in smaller)]
            cover=0
            for L in children:cover|=L
            need(cover==M,'fully covered maximal set')
        categories['fully_covered_laminar']+=1;records.append([a,b,F])
    need(dict(categories)==expected['original_module_classification'],'module-pair categories')
    need(sorted(records)==sorted(expected['module_family_records']),'every exact family')
    a,b=pairs[0];difference=[[m&63,m>>6,(raw[i*N+a]-raw[i*N+b])%4] for i,m in enumerate(basis) if raw[i*N+a]!=raw[i*N+b]]
    need(difference==expected['first_pair']['complete_nonzero_mod4_difference'],'complete higher-degree contrast')
    need(Counter((x.bit_count()+y.bit_count()) for x,y,z in difference)=={3:48,4:48},'first pair degree census')
    return dict(status='PASS',method='Exact sorted vector tuples; exhaustive subsets of modules constructed from physical arc definition. No producer imports.',degree_classes=summaries,quadratic_pairs=len(pairs),module_categories=dict(categories),first_pair=[a,b],first_pair_degree_counts={3:48,4:48},bank_sha256=hashlib.sha256(raw).hexdigest(),source_receipt_sha256=hashlib.sha256((SOURCE/'degree_diagnostic.json').read_bytes()).hexdigest(),checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
if __name__=='__main__':print(json.dumps(main(),indent=2,sort_keys=True))

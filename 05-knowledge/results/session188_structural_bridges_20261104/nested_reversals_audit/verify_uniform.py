"""Independent support-switch and literal uniform-attachment audit.
Imports only this audit packet's independent tournament primitives.
"""
from itertools import combinations
from functools import lru_cache
from pathlib import Path
import hashlib,json,random
from verify_independent import rows,flip,extend,hp,op,is_module,substitute,need


def paths(T):
    out=[]
    def visit(P,used):
        out.append(tuple(P))
        for w in range(len(T)):
            if not used>>w&1 and T[P[-1]]>>w&1:visit(P+[w],used|1<<w)
    for u in range(len(T)):visit([u],1<<u)
    return out


def switch(pair,A):
    P,Q=pair
    p=next((i for i,x in enumerate(P) if x not in A),len(P))
    q=next((i for i,x in enumerate(Q) if x not in A),len(Q))
    return Q[:q]+P[p:],P[:p]+Q[q:]


def switch_checks():
    rng=random.Random(188813);cases=pairs=empty_p=reentries=0
    corpus=[]
    for n in range(3,6):
        codes=range(1<<(n*(n-1)//2)) if n<=4 else rng.sample(range(1024),32)
        corpus.extend(rows(n,c) for c in codes)
    for n in range(6,9):
        for m in (2,n-1):
            A=rows(m,rng.randrange(1<<(m*(m-1)//2)))
            q=n-m+1;Q=rows(q,rng.randrange(1<<(q*(q-1)//2)))
            corpus.append(substitute(Q,A))
    # Deliberate reentry control: quotient 0->1->2->0 and 0,1,2->3.
    corpus.append(substitute([10,12,9,0],[6,4,0]))
    for T in corpus:
        n=len(T);Ps=paths(T)
        for mask in range(1,1<<n):
            A={u for u in range(n) if mask>>u&1}
            if len(A)<2 or len(A)==n or not is_module(T,A):continue
            starts=[P for P in Ps if P[0] in A]
            for v in range(n):
                if v in A:continue
                S={(P,Q) for P in starts for Q in starts if Q[-1]==v and not set(P)&set(Q)}
                support={}
                for P,Q in S:
                    R,W=switch((P,Q),A)
                    need((R,W) in S,'switched pair invalid')
                    need((P,Q)!=(R,W),'fixed point')
                    need(switch((R,W),A)==(P,Q),'not involution')
                    key=frozenset(P+Q)
                    need(frozenset(R+W)==key,'support changed')
                    support[key]=support.get(key,0)^1
                    pairs+=1
                    if set(P)<=A:empty_p+=1
                    for Z in (P,Q):
                        exited=False
                        for u in Z:
                            if u not in A:exited=True
                            elif exited:reentries+=1;break
                need(not any(support.values()),'odd coefficient by full support')
                cases+=1
    need(reentries>0,'missing later-reentry control')
    return {'source_presentations':len(corpus),'module_anchor_cases':cases,'ordered_path_pairs':pairs,
            'pairs_with_empty_first_remainder':empty_p,'later_reentry_paths':reentries,
            'universe':'All labelled tournaments of orders 3 and 4; 32 seeded order-5 sources; six substituted sources of orders 6 through 8, plus an explicit reentry control. Every proper nontrivial module and exterior anchor.'}


def all_packings(Ms):
    yield ()
    for j,M in enumerate(Ms):
        for rest in all_packings([N for N in Ms[j+1:] if not M&N]):yield (M,)+rest


def battery_checks():
    rng=random.Random(188814);sources=[]
    for n in (3,4):
        for code in range(1<<(n*(n-1)//2)):
            T=rows(n,code);Ms=[M for M in range(1,1<<n) if 1<M.bit_count()<n and is_module(T,{u for u in range(n) if M>>u&1})]
            sources.extend((T,P) for P in all_packings(Ms))
    # Larger substituted inner modules, including compressed q=2 cases.
    for n in range(5,9):
        for m in (2,n-1):
            q=n-m+1;A=rows(m,rng.randrange(1<<(m*(m-1)//2)));Q=rows(q,rng.randrange(1<<(q*(q-1)//2)))
            sources.append((substitute(Q,A),((1<<m)-1,)))
    cases=q2=queries=0;examples=[]
    for T,pack in sources:
        n=len(T);N=(1<<n)-1;covered=0
        for M in pack:covered|=M
        free=N^covered
        if not free:continue
        blocks=list(pack)+[1<<v for v in range(n) if free>>v&1];q=len(blocks);O=op(T);U=O
        for M in pack:U=flip(U,{u for u in range(n) if M>>u&1})
        @lru_cache(None)
        def counts(a,c):return hp(extend(T,a,c)),hp(extend(O,a,c)),hp(extend(U,a,c))
        def delta(a,c):
            x,y,z=counts(a,c);need(x%2==y%2,'half contrast not integral');return ((x-y)//2)%2
        for v in range(n):
            if not free>>v&1:continue
            b=(hp(T,v)+hp(O,v))%2;d=0;As=[]
            for C in blocks:
                if C==1<<v:continue
                rep=(C&-C).bit_length()-1
                if T[v]>>rep&1:As.append(C|(1<<v));d+=1
                else:As.append(C)
            if d%2:As.append(1<<v)
            if (q-1+d)%2:As.append(0)
            need(len(As)==len(set(As)),'duplicate queries')
            need(len(As)==(q if q%2==0 else q-1+2*(d%2)),'battery length')
            need(sum(delta(a,1<<v) for a in As)%2==1-b,'compressed E1 battery')
            need(delta(1<<v,0)==b,'endpoint parity')
            data=[(a,1<<v,*counts(a,1<<v)) for a in As]+[(1<<v,0,*counts(1<<v,0))]
            need(all(y==z for a,c,x,y,z in data),'uniform restoration changed exact count')
            need(any((x-z)%4==2 for a,c,x,y,z in data),'no exact named separator')
            queries+=len(data);cases+=1;q2+=q==2
            if n>=5:examples.append({'rows':T,'blocks':blocks,'anchor':v,'compressed_order':q,'queries':[list(r) for r in data]})
    return {'batteries':cases,'q2_batteries':q2,'actual_query_checks':queries,'large_examples':examples,
            'universe':'Every pairwise disjoint proper-module packing in every labelled order-3/4 source, plus eight seeded substituted sources of orders 5 through 8. Counts are full actual Hamilton counts, including all polynomial degrees.'}


def main():
    root=Path(__file__).resolve().parent;source=root.parent/'returns_next'/'UNIFORM.md'
    print(json.dumps({'status':'PASS','switch':switch_checks(),'batteries':battery_checks(),
                     'proof_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
                     'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                     'helper_sha256':hashlib.sha256((root/'verify_independent.py').read_bytes()).hexdigest()},indent=2,sort_keys=True))
if __name__=='__main__':main()

"""Native repair to the older intrinsic pair-halving interface."""
from itertools import permutations,combinations
from support import need,save,junction
from ternary_family import build,first_unit
def pair_edges(T):
    return [(a,b) for a,b in combinations(range(len(T)),2) if not ((T[a]^T[b])&~((1<<a)|(1<<b)))]
def halve(T):
    edges=pair_edges(T);live=set(range(len(T)));matching=[]
    while live:
        selected=next(((a,b) for a in sorted(live) for b in live if a!=b and (min(a,b),max(a,b)) in edges
                       and sum((min(a,c),max(a,c)) in edges for c in live if c!=a)==1),None)
        need(selected is not None,'intrinsic path leaf supplies forced pair')
        a,b=selected;matching.append(tuple(sorted((a,b))));live.remove(a);live.remove(b)
    matching.sort();representatives=[p[0] for p in matching]
    return junction.sub(T,tuple(representatives)),matching
def main():
    out=[]
    for k in range(1,4):
        for flipped in (False,True):
            N=3**k;B=N//3;T=build(k,flipped);initial=T;n=len(T);trace=[]
            blocks=[[b,b+B,b+2*B,N+b] for b in range(B)]
            need(not pair_edges(T),'original cyclic-bottom family has no pair module')
            for b,block in enumerate(blocks):
                a,c=block[:2]
                if not T[a]>>c&1:a,c=c,a
                P=[]
                for j,V in enumerate(blocks):
                    good=[]
                    for word in permutations(V):
                        if not all(T[x]>>y&1 for x,y in zip(word,word[1:])):continue
                        if j==b and not any(word[i:i+2]==(a,c) for i in range(2)):continue
                        good.append(word)
                    need(bool(good),'local block Hamilton path and repair arc')
                    P.extend(min(good))
                P=tuple(P);i=P.index(a)
                need(P[i+1]==c and i<len(P)-2 and all(T[x]>>y&1 for x,y in zip(P,P[1:])),'full native guard from retained blocks')
                X=junction.ear(T,((1<<n)-1)^(1<<a));p1=P+(n,);p2=P[:i+1]+(n,)+P[i+1:]
                Y=junction.native(X,p1);need(all(Y[x]>>y&1 for x,y in zip(p2,p2[1:])),'second native guard')
                Z=junction.native(Y,p2);W=tuple(x&((1<<n)-1) for x in Z[:n])
                changed=[(x,y) for x,y in combinations(range(n),2) if ((T[x]>>y)^(W[x]>>y))&1]
                need(changed==[(min(a,c),max(a,c))],'one literal bottom arc reversal')
                trace.append({'helper_name':n+b,'arc':[a,c],'paths':[[x if x<n else n+b for x in Q] for Q in (p1,p2)]});T=W
            need(len(pair_edges(T))==3*B,'exact path-of-four pair modules')
            U,pairs1=halve(T);V,pairs2=halve(U)
            fibres=[tuple(sorted(x for j in pair for x in pairs1[j])) for pair in pairs2]
            need(set(fibres)=={tuple(sorted(block)) for block in blocks},'old names reconstruct every bottom fibre')
            base_names=[next(b for b,block in enumerate(blocks) if set(block)==set(fibre)) for fibre in fibres]
            need(all(((V[i]>>j)&1)==(first_unit(base_names[j]-base_names[i])==1) for i,j in combinations(range(B),2)),'actual odd-core orientation and original block addresses')
            need(len(U)==2*B and len(V)==B and not pair_edges(V),'two intrinsic halvings reach a pair-free odd core')
            need(all(row.bit_count()==(B-1)//2 for row in V),'regular odd core')
            out.append({'k':k,'flipped':flipped,'source_rows':initial,'final_rows':T,'native_moves':2*B,
                'helpers':B,'peak_order':n+1,'trace':trace,'first_pairs':pairs1,'first_quotient':U,'second_pairs':pairs2,
                'odd_core':V,'odd_core_original_fibres':fibres,'odd_core_block_addresses':base_names})
    save('halving_normalization.json',{'status':'FINITE_EXACT_NATIVE_CONTROLS_FOR_ALL_LEVEL_NORMALIZATION','programs':out,
        'scope':'order halving in this framed graph encoding, not division of Hamilton-path counts or a Collatz convergence proof'})
    print('PASS source-preserving native normalization and two intrinsic halvings through36 vertices',flush=True)
if __name__=='__main__':main()

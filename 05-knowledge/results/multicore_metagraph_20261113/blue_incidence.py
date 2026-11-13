"""Exact complement-line incidence, kept separate from arc-flip graphs."""
from rooted_growth import *

def opposite(T):return tuple(((1<<len(T))-1)^(1<<i)^row for i,row in enumerate(T))
def merged(T):return min(ordinary(T),ordinary(opposite(T)))
def blue(T):return T==opposite(relabel(T,tuple(reversed(range(len(T))))))

def sandwich(T,q):
    blocks=[tiling(q,(1<<((q-1)*(q-2)//2))-1) if q>1 else (0,)] if q else []
    blocks+=[T]
    if q:blocks+=[tuple(sum(1<<j for j in range(i+1,q)) for i in range(q))]
    n=sum(map(len,blocks));out=[];offset=0
    for block in blocks:
        size=len(block);later=((1<<n)-1)^((1<<(offset+size))-1)
        out.extend((row<<offset)|later for row in block);offset+=size
    return tuple(out)

def cyclic_sandwich(T,q):
    rows=list(sandwich(T,q));n=len(rows)
    for a in range(q):
        for c in range(q+len(T),n):
            rows[a]&=~(1<<c);rows[c]|=1<<a
    return tuple(rows)

def run(N=7):
    levels=[];graphs=[];first_partner_hostile=None
    for n in range(3,N+1):
        m=(n-1)*(n-2)//2;full=(1<<m)-1;nodes=Counter();blue_degree=Counter();edges=Counter()
        states=defaultdict(dict);partners=defaultdict(dict);node_counts={}
        for bits in range(1<<m):
            T=tiling(n,bits);node=merged(T);nodes[node]+=1
            isblue=reflection(n,bits)==bits;need(isblue==blue(T),'literal reversal anti-isomorphism')
            blue_degree[node]+=isblue
            partner=merged(tiling(n,bits^full));state=rooted(T,n-1)[0]
            states[state].setdefault(isblue,bits);partners[state].setdefault(partner,bits)
            if first_partner_hostile is None and len(partners[state])>1:
                first_partner_hostile=dict(n=n,rooted_rows=state,examples=[dict(partner=p,bits=b) for p,b in sorted(partners[state].items())])
            if bits<bits^full:edges[tuple(sorted((node,partner)))+(isblue,)]+=1
        for node,mass in nodes.items():
            T=cc.decode(n,node);aut=tc.automorphism_count(T);h=sum(endpoints(T));sc=ordinary(tuple(T))==ordinary(opposite(T))
            need(mass==h//aut*(1 if sc else 2),'merged orbit mass')
            need((blue_degree[node]%2==1) if sc else (blue_degree[node]==0),'blue activity exactly self-converse')
            node_counts[node]=dict(code=node,H=h,aut=aut,self_converse=sc,tilings=mass,blue_degree=blue_degree[node],black_degree=mass-blue_degree[node])
        for (a,b,isblue),count in edges.items():need(isblue or count%2==0,'mirror-paired black line multiplicity')
        degrees=Counter()
        for (a,b,isblue),count in edges.items():degrees[a]+=count;degrees[b]+=count
        need(degrees==nodes,'line half-edge mass including loops')
        level=dict(n=n,merged_nodes=len(nodes),tilings=sum(nodes.values()),blue_lines=sum(v for (a,b,c),v in edges.items() if c),black_lines=sum(v for (a,b,c),v in edges.items() if not c),rooted_color_splits=sum(len(s)>1 for s in states.values()),rooted_partner_splits=sum(len(s)>1 for s in partners.values()))
        levels.append(level);graphs.append(dict(n=n,nodes=list(node_counts.values()),edges=[dict(a=a,b=b,blue=c,multiplicity=v) for (a,b,c),v in sorted(edges.items())]))
        print(level,flush=True);rooted.cache_clear();ordinary.cache_clear()
    A=tiling(4,2);B=tiling(4,4)
    mapping=next(p for p in permutations(range(4)) if relabel(B,p)==A and p[-1]==3)
    padding=[]
    for q in range(9):
        U,V=sandwich(A,q),sandwich(B,q)
        p=tuple(range(q))+tuple(q+i for i in mapping)+tuple(range(q+4,2*q+4))
        need(relabel(V,p)==U,'boundary-fixed sandwich isomorphism')
        need(blue(U) and not blue(V),'all displayed sandwich colors')
        CU,CV=cyclic_sandwich(A,q),cyclic_sandwich(B,q)
        need(relabel(CV,p)==CU,'strong boundary-fixed sandwich isomorphism')
        need(blue(CU) and not blue(CV),'strong sandwich colors')
        from component_growth import components
        need(len(components(CU))==len(components(CV))==1,'strong sandwich controls')
        padding.append(dict(q=q,blue_rows=U,black_rows=V,strong_blue_rows=CU,strong_black_rows=CV,isomorphism=p))
    save('blue_incidence',dict(scope=f'All tilings and complement lines at orders3..{N}; no arc-flip edges; loops and multiplicities retained',levels=levels,graphs=graphs,first_rooted_partner_hostile=first_partner_hostile,boundary_window_controls=padding))
if __name__=='__main__':run(int(sys.argv[1]) if len(sys.argv)>1 else 7)

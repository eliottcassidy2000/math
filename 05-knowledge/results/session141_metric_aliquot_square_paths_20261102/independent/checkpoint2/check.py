from itertools import combinations
from collections import Counter
from math import isqrt, comb
import json
from fractions import Fraction as F
from itertools import permutations, product

def run(n):
    vertices=tuple(range(1,n+1))
    edges=tuple((i,j) for i in vertices for j in vertices if i<j and isqrt(i+j)**2==i+j)
    adj={i:tuple(j for j in vertices if j!=i and isqrt(i+j)**2==i+j) for i in vertices}
    paths=[]
    calls=0
    def visit(path,used):
        nonlocal calls
        calls+=1
        if len(path)==n: paths.append(tuple(path)); return
        for j in adj[path[-1]]:
            b=1<<(j-1)
            if not used&b: visit(path+[j],used|b)
    for v in vertices: visit([v],1<<(v-1))
    arcs=[frozenset(zip(p,p[1:])) for p in paths]
    compatible=lambda a,b: not any((v,u) in b for u,v in a)
    moments=Counter({0:2**len(edges)})
    profiles={}
    subsets=[]
    def walk(chosen,candidates,union):
        if chosen:
            k=len(chosen)
            profiles.setdefault(k,Counter())[len(union)]+=1
            moments[k]+=2**(len(edges)-len(union))
            subsets.append((chosen,len(union)))
        for pos,j in enumerate(candidates):
            walk(chosen+(j,),[k for k in candidates[pos+1:] if compatible(arcs[j],arcs[k])],union|arcs[j])
    walk((),list(range(len(paths))),frozenset())
    coeff=[sum((-1)**(j-k)*comb(j,k)*a for j,a in moments.items() if j>=k) for k in range(max(moments)+1)]
    out={'n':n,'guard_edges':len(edges),'directed_paths':len(paths),'dfs_calls':calls,'subset_union_profiles':profiles,'binomial_moments':moments,'guard_polynomial':coeff,'free_tournament_edges':comb(n,2)-len(edges),'paths':paths,'compatible_subsets':subsets}
    if sum(coeff)!=2**len(edges) or any(x<0 for x in coeff): raise RuntimeError('polynomial gate')
    return out


def require(c,label):
    if not c: raise RuntimeError(label)

def native_collatz():
    seen={5:0}; front={5}; counts=[1]
    for depth in range(1,27):
        nxt={2*m for m in front}|{(m-1)//3 for m in front if m%6==4}
        require(not(set(seen)&nxt),'root5 level disjointness')
        for n in nxt:seen[n]=depth
        front=nxt;counts.append(len(seen))
    for n,depth in seen.items():
        x=n
        for _ in range(depth):
            require(x!=5,'first hit clock premature')
            x=x//2 if x%2==0 else 3*x+1
        require(x==5,'first hit endpoint')
    require(counts[21]==704 and counts[26]==2263,'native counts')
    def bfs(kind,radius):
        origin=(0,0) if kind=='triangular' else (0,0,0)
        seen={origin};front={origin};sizes=[1]
        for r in range(1,radius+1):
            nxt=set()
            for v in front:
                if kind=='triangular':
                    a,b=v
                    nxt.update((a+da,b+db) for da,db in ((1,0),(-1,0),(0,1),(0,-1),(1,-1),(-1,1)))
                else:
                    a,b,t=v
                    nxt.update(((a,b,1),(a-1,b,1),(a,b-1,1)) if t==0 else ((a,b,0),(a+1,b,0),(a,b+1,0)))
            front=nxt-seen;seen|=front;sizes.append(len(seen))
            expected=1+3*r*(r+1)//(1 if kind=='triangular' else 2)
            require(len(seen)==expected,'lattice ball')
        return sizes
    honey=bfs('honeycomb',26);tri=bfs('triangular',26)
    require(next(i for i in range(27) if counts[i]>honey[i])==21,'honey first count violation')
    require(next(i for i in range(27) if counts[i]>tri[i])==26,'tri first count violation')
    return {'scope':'all root5 inverse states through26 and direct forward first-hit replay; lattice BFS through26','root5_counts':counts,'honeycomb_counts':honey,'triangular_counts':tri,'root21_excess':counts[21]-honey[21],'root26_excess':counts[26]-tri[26]}

def curvature():
    results=[]
    for q,W in product((F(1,3),F(2,5),F(49,100)),(F(1),F(1,2),F(1,10))):
        points={'':(W/2,F(1))}
        intervals={'':(F(0),W)}
        for k in range(6):
            for word in [w for w in points if len(w)==k]:
                a,b=intervals[word]; width=b-a
                for bit,lo,hi in [('0',a,a+q*width),('1',b-q*width,b)]:
                    intervals[word+bit]=(lo,hi); points[word+bit]=((lo+hi)/2,q**(k+1))
        g=W*(1-2*q); Lcosh=1+(1+W*W/4)*(1-q)**2/(2*q)
        mincosh=min((q+1/q)/2,1+(W*(1-q)/q)**2/2)
        for u,v in combinations(points,2):
            x,y=points[u];xx,yy=points[v]
            c=1+((x-xx)**2+(y-yy)**2)/(2*y*yy)
            l=0
            while l<min(len(u),len(v)) and u[l]==v[l]:l+=1
            D=len(u)+len(v)-2*l
            A=g*g*q**(-D);B=(W*W+2)*q**(-D)
            require(c<= (B+1/B)/2,'QI upper')
            if A>=1: require(c>=(A+1/A)/2,'QI lower')
            require(c>=mincosh,'uniform separation')
            if (len(v)==len(u)+1 and v[:-1]==u) or (len(u)==len(v)+1 and u[:-1]==v): require(c==Lcosh,'equal edge length')
        a=W*(1-q)/2
        require(a*a<1-q*q,'edge circle center strict monotone side')
        results.append({'q':str(q),'W':str(W),'vertices':len(points),'cosh_edge_length':str(Lcosh),'cosh_separation_lower_bound':str(mincosh)})
    # Exact dyadic extreme radius: root (1/2,1), outer x=2^(-k-1).
    extremes=[]
    for k in range(1,13):
        h=F(1,2**k);x=h/2
        c=1+((x-F(1,2))**2+(h-1)**2)/(2*h)
        expected=F(5,8)*2**k-F(1,4)+F(5,8)*h
        require(c==expected,'dyadic outer radius')
        extremes.append(str(c))
    return {'scope':'rational coordinates, every pair through depth6 in nine gapped models; analytic all-level proof audited separately','gapped_models':results,'dyadic_outer_cosh_radii_depth1_to12':extremes}


def common_block_tests():
    cases=0
    for n in range(1,6):
        paths=list(permutations(range(n)))
        for size in range(1,4 if n<=4 else 3):
            for bank in combinations(paths,size):
                arcsets=[set(zip(p,p[1:])) for p in bank]
                shared=set.intersection(*arcsets); union=set.union(*arcsets)
                succ=dict(shared); pred={v:u for u,v in shared}; blocks=[]
                for root in range(n):
                    if root in pred:continue
                    b=[root]
                    while b[-1] in succ:b.append(succ[b[-1]])
                    blocks.append(tuple(b))
                owner={v:i for i,b in enumerate(blocks) for v in b}
                joins=set()
                for p in bank:
                    order=[]
                    for v in p:
                        i=owner[v]
                        if not order or order[-1]!=i:order.append(i)
                    require(len(order)==len(blocks),'block revisited')
                    require(tuple(v for i in order for v in blocks[i])==p,'block direction or contiguity')
                    joins.update(zip(order,order[1:]))
                require(len(union)==n-len(blocks)+len(joins),'union normal form')
                compatible=not any((v,u) in union for u,v in union)
                singleton_ok=not any((j,i) in joins and len(blocks[i])==len(blocks[j])==1 for i,j in joins)
                require(compatible==singleton_ok,'singleton conflict iff')
                cases+=1
    def hpaths(vertices,edges):return [p for p in permutations(vertices) if all(e in edges for e in zip(p,p[1:]))]
    blowups=0
    for mask in range(16):
        options=[(0,1),(0,2),(1,2),(2,1)];edges={e for j,e in enumerate(options) if mask>>j&1}
        blocks=[(0,),(1,2),(3,4,5)]
        expanded=set().union(*(set(zip(b,b[1:])) for b in blocks))|{(blocks[i][-1],blocks[j][0]) for i,j in edges}
        expected={tuple(v for i in p for v in blocks[i]) for p in hpaths(range(3),edges)}
        actual=set(hpaths(range(6),expanded));require(actual==expected,'source-forced blowup')
        blowups+=1
    require(len(hpaths(range(3),{(0,1),(1,2),(2,0)}))==3,'cycle hostile base')
    require(len(hpaths(range(4),{(0,1),(1,2),(2,3),(3,0)}))==4,'cycle hostile blowup')
    return {'common_block_banks':cases,'scope':'all one/two/three-path banks n<=4 and all one/two-path banks n5','source_zero_digraphs_on_three_vertices_checked':blowups,'source_free_cycle_hostile_counts':[3,4]}

def main():
    square={str(n):run(n) for n in range(1,31)}
    expected={25:[4294962432,4608,256],26:[17179857952,10240,928,64],27:[68719413784,54976,7360,488,120,8],28:[274877723488,158096,21840,2928,528,64],29:[1099511486464,129024,10240,2048],30:[4398046212096,274432,20480,4096]}
    for n,e in expected.items():require(square[str(n)]['guard_polynomial']==e,'square coefficients')
    print(json.dumps({'status':'PASS','assert_statements_used':False,'curvature':curvature(),'native_collatz':native_collatz(),'common_blocks':common_block_tests(),'square_paths':square},separators=(',',':'),sort_keys=True))
if __name__=='__main__':main()

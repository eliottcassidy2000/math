"""All minimum positive-integer connectors for the user's square-sum pieces."""
from math import isqrt
from itertools import permutations,product,combinations
from functools import lru_cache
from collections import defaultdict
from support import need,save

BASE=(527,2,7,44,5,4,13,3)
USER=(527,2,7,317,44,5,4,12,13,3)
def square(n):return isqrt(n)**2==n

def bridges(a,b,excluded):
    a,b=max(a,b),min(a,b);d=a-b;out=[]
    need(d>0,'distinct old endpoints')
    for u in range(1,isqrt(d)+1):
        if d%u:continue
        v=d//u
        if (v-u)%2:continue
        c=((v-u)//2)**2-b
        if c>0 and c not in excluded:
            need(square(a+c) and square(b+c),'divisor-pair bridge');out.append(c)
    return tuple(sorted(set(out)))

def components(V):
    unseen=set(V);out=[]
    while unseen:
        start=min(unseen);todo=[start];unseen.remove(start);C=[]
        while todo:
            v=todo.pop();C.append(v)
            for w in sorted(unseen.copy()):
                if square(v+w):unseen.remove(w);todo.append(w)
        out.append(tuple(sorted(C)))
    return tuple(out)

def path_count(V):
    V=tuple(V);n=len(V);adj=tuple(sum(square(a+b)<<j for j,b in enumerate(V) if i!=j) for i,a in enumerate(V));full=(1<<n)-1
    @lru_cache(None)
    def f(S,v):
        if S==full:return 1
        return sum(f(S|1<<j,j) for j in range(n) if adj[v]>>j&1 and not S>>j&1)
    return sum(f(1<<i,i) for i in range(n))

def classify(V):
    V=tuple(V);old=set(V);C=components(V)
    blocks=[tuple(p for p in permutations(c) if all(square(a+b) for a,b in zip(p,p[1:]))) for c in C]
    need(all(blocks),'traceable components')
    catalog={(a,b):bridges(a,b,old) for a,b in combinations(sorted({p[0] for B in blocks for p in B}),2)}
    paths=set()
    for order in permutations(range(len(C))):
        for chosen in product(*(blocks[i] for i in order)):
            menus=[catalog[tuple(sorted((A[-1],B[0])))] for A,B in zip(chosen,chosen[1:])]
            for connectors in product(*menus):
                if len(set(connectors))!=len(connectors):continue
                p=chosen[0]
                for c,B in zip(connectors,chosen[1:]):p+=(c,)+B
                need(len(set(p))==len(p) and all(square(a+b) for a,b in zip(p,p[1:])),'actual square-sum path')
                paths.add(min(p,p[::-1]))
    grouped=defaultdict(list)
    for p in sorted(paths):grouped[tuple(sorted(set(p)-old))].append(p)
    sets=sorted(grouped,key=lambda c:(max(c),sum(c),c))
    controls=[]
    user_connectors=tuple(sorted(set(USER+(1,8) if 8 in old else USER)-old))
    for c in dict.fromkeys(sets[:3]+[user_connectors]+sets[-2:]):
        need(c in grouped,'positive connector set in complete classification')
        count=path_count(V+c);need(count==2*len(grouped[c]),'independent undirected subset DP')
        controls.append({'connectors':c,'oriented_path_count':count})
    hostile=None
    if len(V)==8:
        candidates=sorted({x for B in catalog.values() for x in B})
        for c in sorted(combinations(candidates,2),key=lambda c:(max(c),sum(c),c)):
            if c not in grouped and len(components(V+c))==1:
                need(path_count(V+c)==0,'connected is insufficient for a minimum splice')
                hostile={'connectors':c,'connected':True,'Hamilton_paths':0};break
    return {'old_vertices':V,'components':C,'minimum_new_vertices':len(C)-1,'positive_integer_scope':'all, with no numerical cutoff',
      'catalog':[{'endpoints':p,'connectors':c} for p,c in sorted(catalog.items()) if c],
      'connector_sets':sets,'unoriented_path_count':len(paths),'paths_by_connector_set':[{'connectors':c,'paths':grouped[c]} for c in sets],
      'independent_controls':controls,'connectivity_hostile':hostile}

def main():
    A=classify(BASE);B=classify(BASE+(8,))
    need(len(A['components'])==3 and len(B['components'])==4,'three pieces then isolated8')
    need(min(USER,USER[::-1]) in {tuple(p) for d in A['paths_by_connector_set'] for p in d['paths']},'user supplied bridge')
    extended=USER+(1,8);need(min(extended,extended[::-1]) in {tuple(p) for d in B['paths_by_connector_set'] for p in d['paths']},'user1--8 extension')
    alternate=classify((527,2,7,317,5,4,13,3))
    save('square_splices.json',{'status':'FINITE_EXACT_COMPLETE_MINIMUM_CONNECTOR_CLASSIFICATION','base':A,'with8':B,'alternate_317_anchor':alternate,
      'user_path':USER,'user_squares':[isqrt(a+b) for a,b in zip(USER,USER[1:])],'extended_path':extended,
      'response_depths':{'base':2,'with8':3},'anchor_correction':'Incoming S173 recovers original44 anchor: connectors are12,317. Treating317 as an anchor instead gives a different valid thirteen-pair problem, retained separately.',
      'scope':'Hamilton paths through all old and added distinct positive integers; paths counted modulo reversal'})
    print('PASS square splices: connector sets / unoriented paths',[(len(d['connector_sets']),d['unoriented_path_count']) for d in (A,B)],'best',[d['connector_sets'][0] for d in (A,B)],'hostile',A['connectivity_hostile'],flush=True)
if __name__=='__main__':main()

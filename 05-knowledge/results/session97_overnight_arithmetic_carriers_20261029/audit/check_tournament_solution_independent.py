"""Independent synchronous-iteration and Tarjan-SCC tournament solution audit."""
from itertools import combinations,permutations
import json

def need(b,m):
    if not b:raise ArithmeticError(m)
def source_scc(V,E):
    adjacency={v:[] for v in V}
    for x,y in E:adjacency[x].append(y)
    indices={};low={};stack=[];active=set();parts=[]
    def visit(v):
        indices[v]=low[v]=len(indices);stack.append(v);active.add(v)
        for w in adjacency[v]:
            if w not in indices:visit(w);low[v]=min(low[v],low[w])
            elif w in active:low[v]=min(low[v],indices[w])
        if low[v]==indices[v]:
            group=set()
            while True:
                w=stack.pop();active.remove(w);group.add(w)
                if w==v:break
            parts.append(group)
    for v in V:
        if v not in indices:visit(v)
    return frozenset(x for S in parts if not any(y in S and x not in S for x,y in E) for x in S)

def audit(n,E):
    sets=sorted((frozenset(v for v in range(n) if mask>>v&1) for mask in range(1<<n)),key=len)
    depth={};fixed={}
    for W in sets:
        children={y:frozenset(x for x in W if (x,y) in E) for y in W}
        depth[W]=0 if len(W)<=1 else 1+max([0]+[depth[C] for C in children.values()])
        fixed[W]=source_scc(W,{(x,y) for y,C in children.items() for x in fixed[C]})
    states=[]
    for kind in [0,1]:
        current={}
        for W in sets:
            scores={x:sum((x,y) in E for y in W) for x in W}
            current[W]=W if kind==0 else frozenset(x for x in W if scores[x]==max(scores.values()))
        for k in range(n+1):
            for W in sets:
                if k>=depth[W]:need(current[W]==fixed[W],'all later iterations stabilized')
            current={W:source_scc(W,{(x,y) for y in W for x in current[frozenset(z for z in W if (z,y) in E)]}) for W in sets}
    need(depth[sets[-1]]<=n-1,'depth bound')

def main():
    total=0
    for n in range(1,6):
        pairs=list(combinations(range(n),2))
        for bits in range(1<<len(pairs)):
            E={(a,b) if bits>>i&1 else (b,a) for i,(a,b) in enumerate(pairs)}
            audit(n,E);total+=1
    E={(0,1),(2,0),(0,3),(4,0),(1,2),(1,3),(1,4),(2,3),(2,4),(3,4)}
    V=range(5);A=[2,0,1,0,3];B=[4,0,1,2,3]
    need(source_scc(V,E)==frozenset(V),'strong host')
    need(sum({(p[x],p[y]) for x,y in E}==E for p in permutations(V))==1,'rigid host')
    wins=[source_scc(V,{(v[y],y) for y in V}) for v in [A,B]]
    need(wins==[frozenset([0,1,2]),frozenset(V)],'different winner sets')
    for p in permutations(V):
        for v,W in zip([A,B],wins):
            mapped={(p[v[y]],p[y]) for y in V}
            need(source_scc(V,mapped)==frozenset(p[x] for x in W),'equivariant transport')
    print(json.dumps({'labelled_tournaments':total,'seeds':['all','maximum outdegree'],
        'host_rigid_and_strong':True,'winner_sets':[sorted(w) for w in wins],
        'host_relabellings_checked':120},sort_keys=True))
if __name__=='__main__':main()

"""Separate physical isomorphisms from exact continuation-equivalent cores."""
from context import *
from collections import defaultdict,Counter
from itertools import combinations
import networkx as nx

def graph(r,marked):
    G=nx.DiGraph();G.add_nodes_from(range(len(r['rows'])))
    for i,row in enumerate(r['rows']):
        G.add_edges_from((i,j) for j in range(len(r['rows'])) if row>>j&1)
        G.nodes[i]['color']=(r['roots'].index(i)+1 if i in r['roots'] else 3 if i==r['clock'] else 0) if marked else 0
    return G

def run():
    data=load('endpoint_profiles.json.gz',HERE)['records'];groups=defaultdict(list);stats=Counter();pairs=[]
    for r in data:groups[tuple(map(tuple,r['observer']))].append(r)
    for group in groups.values():
        for a,b in combinations(group,2):
            gm=nx.algorithms.isomorphism.DiGraphMatcher(graph(a,True),graph(b,True),node_match=lambda x,y:x['color']==y['color'])
            marked=gm.is_isomorphic()
            unmarked=marked or nx.is_isomorphic(graph(a,False),graph(b,False))
            stats['pairs']+=1;stats['marked_isomorphic' if marked else 'unmarked_only' if unmarked else 'different_unmarked_cores']+=1
            pairs.append({'left':{k:v for k,v in a.items() if k!='observer'},'right':{k:v for k,v in b.items() if k!='observer'},'named_isomorphic':marked,'unmarked_isomorphic':unmarked,'common_observer':a['observer']})
    widths=Counter(len(r['observer'][0]) for r in data)
    out={'scope':'All78 retained source/free-clock presentations; named roots ordered, clock distinguished. Equality of complete17-query exponential profiles proves equality of every future12-command H response at every positive clock length. Original seed marking and arbitrary operations are not quotiented.','classes':len(groups),'widths':dict(widths),'statistics':dict(stats),'pairs':pairs}
    save('observer_rewrites.json.gz',out);print({k:v for k,v in out.items() if k!='pairs'})
if __name__=='__main__':run()

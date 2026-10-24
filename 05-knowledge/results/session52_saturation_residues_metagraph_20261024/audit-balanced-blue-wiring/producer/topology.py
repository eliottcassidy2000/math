#!/usr/bin/env python3
"""Construct verifiable spanning attachment to the inherited connected core."""
from pathlib import Path
from collections import deque,Counter
import json
HERE=Path(__file__).resolve().parent
def need(ok,msg):
    if not ok:raise RuntimeError(msg)
def main():
    d=json.loads((HERE/'normal.json').read_text());vertices={r['atlas_index'] for r in d['nodes']};core={r['atlas_index'] for r in d['nodes'] if r['refinement_error']>0}
    adj={v:set() for v in vertices};witness={}
    for w,a,p,c,kind,b,q in d['ports']:
        if kind=='internal' and a!=b:
            adj[a].add(b);witness.setdefault(tuple(sorted((a,b))),min(w,c))
    def bfs(roots,allowed):
        parent={v:None for v in sorted(roots)};dist={v:0 for v in roots};todo=deque(sorted(roots))
        while todo:
            v=todo.popleft()
            for w in sorted(adj[v]&allowed):
                if w not in parent:parent[w]=v;dist[w]=dist[v]+1;todo.append(w)
        return parent,dist
    cp,cd=bfs({min(core)},core);need(set(cp)==core,'positive core connected')
    parent,dist=bfs(core,vertices);need(set(parent)==vertices,'every node attaches to core')
    tree=[(v,p) for v,p in cp.items() if p is not None]+[(v,p) for v,p in parent.items() if p is not None]
    need(len(tree)==len(vertices)-1,'spanning tree edge count')
    # Explicit edge witnesses are actual C-paired words in normal.json.
    tree=[[a,b,witness[tuple(sorted((a,b)))]] for a,b in tree]
    paths=[]
    for v in sorted(vertices-core):
        path=[v]
        while parent[path[-1]] is not None:path.append(parent[path[-1]])
        need(path[-1] in core and len(path)-1==dist[v],'shortest core attachment')
        paths.append({'node':v,'distance':dist[v],'path_to_core':path})
    core_edges=sum(m for a,b,m in d['adjacency_multiplicities'] if a in core and b in core)
    need(core_edges==90,'inherited induced core edges')
    crossing=sum(m for a,b,m in d['adjacency_multiplicities'] if (a in core)!=(b in core))
    zero=sum(m for a,b,m in d['adjacency_multiplicities'] if a not in core and b not in core)
    out={'status':'PASS','positive_core_nodes':sorted(core),'positive_core_component_sizes':[len(core)],'positive_core_edges':core_edges,'zero_error_nodes':len(vertices-core),'core_to_zero_edges':crossing,'zero_to_zero_edges':zero,'attachment_distance_histogram':sorted(Counter(dist[v] for v in vertices-core).items()),'spanning_tree_schema':['source_atlas_index','target_atlas_index','smaller_C_pair_word'],'spanning_tree':tree,'shortest_attachments':paths,'excluded_nonbalanced_old_indices':d['known_positive_control']['excluded_nonbalanced_old_indices'],'interpretation':'The inherited twenty-node balanced positive-error core was already connected. The four old isolated positive-error nodes are nonbalanced and excluded, not joined by this graph.'}
    need(90+crossing+zero==11452,'edge partition')
    (HERE/'topology.json').write_text(json.dumps(out,sort_keys=True,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k not in ('spanning_tree','shortest_attachments','positive_core_nodes')},sort_keys=True))
if __name__=='__main__':main()

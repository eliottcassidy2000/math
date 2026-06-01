"""
Exact planarity test via Demoucron-Malgrange-Pertuiset (DMP) face-tracing algorithm.
Works on any simple graph. Returns True/False. Self-contained, no deps.
We test each biconnected component? For tiny graphs (<=7 vertices) we just run DMP on the
whole graph after checking it's connected (planarity is per-component; a graph is planar iff
every connected component is planar; DMP requires connected, ideally biconnected, but we
handle general connected graphs by adding components independently).
"""
from itertools import combinations

def _planar_connected(verts, edge_set):
    # edge_set: set of frozenset pairs among verts (a connected graph)
    verts = list(verts)
    E = set(edge_set)
    n = len(verts)
    m = len(E)
    if n <= 4:
        return True
    if m > 3*n - 6:
        return False
    adj = {v:set() for v in verts}
    for e in E:
        a,b = tuple(e)
        adj[a].add(b); adj[b].add(a)

    # Find an initial cycle. DFS to get a cycle.
    def find_cycle():
        parent={}; visited=set()
        for s in verts:
            if s in visited: continue
            stack=[(s,None)]
            while stack:
                u,p=stack.pop()
                if u in visited:
                    continue
                visited.add(u); parent[u]=p
                for w in adj[u]:
                    if w==p: continue
                    if w in visited:
                        # found back edge u-w : reconstruct cycle
                        cyc=[u]; x=u
                        # walk parents of u until reach w
                        path_u=[]; cur=u
                        anc=set()
                        cur2=u
                        chain={}
                        # build ancestor chain of u
                        au=[]; c=u
                        while c is not None:
                            au.append(c); c=parent[c]
                        if w in au:
                            idx=au.index(w)
                            cyc=au[:idx+1]
                            return cyc
                    else:
                        stack.append((w,u))
        return None

    cyc = find_cycle()
    if cyc is None:
        return True  # forest
    cycle_edges = set()
    L=len(cyc)
    for i in range(L):
        cycle_edges.add(frozenset((cyc[i],cyc[(i+1)%L])))

    # Embedding represented as set of faces; each face is a set of vertices on its boundary,
    # plus we need boundary as cyclic for path insertion. We track faces as frozenset of vertices
    # AND keep the actual edge boundary to test if a path's endpoints lie on the face.
    # Represent face boundary as a list (cyclic). Start: two faces = the cycle (inside & outside).
    placed_edges = set(cycle_edges)
    faces = [list(cyc), list(cyc)]  # inner and outer, same boundary vertex cycle

    def boundary_vertices(face):
        return set(face)

    # Bridges: connected pieces of G - placed_edges that attach to placed vertices,
    # plus chords (single edges between placed vertices not yet placed).
    def compute_bridges():
        placed_v = set()
        for e in placed_edges:
            placed_v |= set(e)
        bridges=[]
        # chords: edges between two placed vertices, not yet placed
        for e in E - placed_edges:
            a,b=tuple(e)
            if a in placed_v and b in placed_v:
                bridges.append({'edges':{e},'verts':{a,b},'attach':{a,b}})
        # components of unplaced vertices
        unplaced = set(verts)-placed_v
        seen=set()
        for s in unplaced:
            if s in seen: continue
            comp=set(); stack=[s]; battach=set(); bedges=set()
            while stack:
                u=stack.pop()
                if u in comp: continue
                comp.add(u); seen.add(u)
                for w in adj[u]:
                    e=frozenset((u,w))
                    if e in placed_edges: continue
                    bedges.add(e)
                    if w in placed_v:
                        battach.add(w)
                    elif w not in comp:
                        stack.append(w)
            bridges.append({'edges':bedges,'verts':comp|battach,'attach':battach})
        return bridges

    while True:
        if placed_edges == E:
            return True
        bridges = compute_bridges()
        if not bridges:
            return True
        # for each bridge, find faces that contain all its attachments
        def admissible(bridge):
            att=bridge['attach']
            res=[]
            for fi,f in enumerate(faces):
                fv=set(f)
                if att <= fv:
                    res.append(fi)
            return res
        # pick bridge with fewest admissible faces
        chosen=None; chosen_faces=None
        for b in bridges:
            af=admissible(b)
            if len(af)==0:
                return False
            if chosen is None or len(af)<len(chosen_faces):
                chosen=b; chosen_faces=af
        # find a path through the chosen bridge between two attachment points (subset of attach)
        b=chosen
        att=list(b['attach'])
        # build subgraph of bridge edges
        badj={}
        for e in b['edges']:
            a,c=tuple(e)
            badj.setdefault(a,set()).add(c); badj.setdefault(c,set()).add(a)
        # find a path between two attachment vertices using only bridge edges, passing through interior
        path=None
        if len(att)==1:
            # bridge attaches at single vertex -> it's a "leaf bridge"; place all its edges,
            # it doesn't split a face. Just add edges, keep face.
            placed_edges |= b['edges']
            continue
        # BFS path between att[0] and some other attachment
        import collections as _c
        start=att[0]
        targets=set(att[1:])
        prev={start:None}; q=_c.deque([start])
        end=None
        while q:
            u=q.popleft()
            if u in targets:
                end=u; break
            for w in badj.get(u,()):
                if w not in prev:
                    prev[w]=u; q.append(w)
        if end is None:
            # attachments not connected within bridge (shouldn't happen for a true bridge comp);
            # just place edges
            placed_edges |= b['edges']
            continue
        # reconstruct path
        path=[]; cur=end
        while cur is not None:
            path.append(cur); cur=prev[cur]
        path.reverse()
        u0,u1=path[0],path[-1]
        # place this path into one admissible face that contains BOTH endpoints
        target_face=None
        for fi in chosen_faces:
            fv=set(faces[fi])
            if u0 in fv and u1 in fv:
                target_face=fi; break
        if target_face is None:
            return False
        # split face: the face boundary is a cyclic list; path endpoints split it into two arcs;
        # new faces = arc1 + path , arc2 + reversed path
        fb=faces[target_face]
        # positions of u0,u1 on boundary
        try:
            i0=fb.index(u0); i1=fb.index(u1)
        except ValueError:
            return False
        if i0>i1: i0,i1=i1,i0; pth=path[::-1]
        else: pth=path
        arc1=fb[i0:i1+1]
        arc2=fb[i1:]+fb[:i0+1]
        interior=pth[1:-1]
        face_a = arc1 + interior[::-1]
        face_b = arc2 + interior
        faces.pop(target_face)
        faces.append(face_a)
        faces.append(face_b)
        # place path edges and vertices
        for i in range(len(path)-1):
            placed_edges.add(frozenset((path[i],path[i+1])))

def is_planar(N, edges):
    """N vertices labelled 0..N-1, edges = iterable of frozenset/tuple pairs."""
    E=set(frozenset(e) for e in edges)
    adj={v:set() for v in range(N)}
    for e in E:
        a,b=tuple(e); adj[a].add(b); adj[b].add(a)
    # components
    seen=set()
    for s in range(N):
        if s in seen: continue
        comp=set(); stack=[s]
        while stack:
            u=stack.pop()
            if u in comp: continue
            comp.add(u); seen.add(u)
            for w in adj[u]: stack.append(w)
        ce=set(e for e in E if set(e)<=comp)
        if not _planar_connected(comp, ce):
            return False
    return True

if __name__=="__main__":
    # tests
    def K(n):
        return n, [frozenset((i,j)) for i,j in combinations(range(n),2)]
    for n in range(1,8):
        N,e=K(n)
        print("K%d planar?"%n, is_planar(N,e), "(expect", n<=4,")")
    # K3,3
    e=[frozenset((i,3+j)) for i in range(3) for j in range(3)]
    print("K3,3 planar?", is_planar(6,e), "(expect False)")
    # a cycle
    print("C6 planar?", is_planar(6,[frozenset((i,(i+1)%6)) for i in range(6)]),"(expect True)")
    # K5 minus an edge (planar)
    N,e=K(5); e=[x for x in e if x!=frozenset((0,1))]
    print("K5-e planar?", is_planar(5,e),"(expect True)")
    # Petersen graph (nonplanar)
    pet_edges=[(0,1),(1,2),(2,3),(3,4),(4,0),(0,5),(1,6),(2,7),(3,8),(4,9),(5,7),(7,9),(9,6),(6,8),(8,5)]
    print("Petersen planar?", is_planar(10,[frozenset(e) for e in pet_edges]),"(expect False)")

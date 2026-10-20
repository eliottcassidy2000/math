#!/usr/bin/env python3
"""Cyclic-five snark poles: exact graphs, labelled tensors, rooted certificates."""
from collections import Counter
from itertools import combinations, chain, product, permutations
from pathlib import Path
import json
import random
import runpy
import networkx as nx

ROOT = Path(__file__).resolve().parents[1]
PRE = runpy.run_path(str(ROOT / "04-computation/blue_quotient_20261020b_snarks.py"))
P = PRE["PETERSEN"]
edge = PRE["edge"]
adjacency = PRE["adjacency"]
col = PRE["pole_colourings"]
check = PRE["check"]
suppress = PRE["suppress_interior"]


def deletion_pole(g, removed):
    removed = set(removed)
    cut = [e for e in sorted(g) if len(set(e) & removed) == 1]
    live = [next(v for v in e if v not in removed) for e in cut]
    pole = {e for e in g if not set(e) & removed}
    pole.update(edge(-i-1, v) for i, v in enumerate(live))
    return pole, cut


def negator(g=P, path=(0, 1, 2)):
    u, w, v = path
    adj = adjacency(g)
    check(w in adj[u] and v in adj[w], "negator path")
    ports = sorted(adj[u]-{w})+sorted(adj[v]-{w})+sorted(adj[w]-{u, v})
    check(len(ports) == 5, "negator boundary")
    result = {e for e in g if not set(e) & set(path)}
    result.update(edge(-i-1, a) for i, a in enumerate(ports))
    return result


def loupekine(twist=0, factors=None):
    """Three Petersen negators in a ring, residuals at one cubic vertex."""
    if factors is None:
        factors = [negator() for _ in range(3)]
    g = set()
    connectors = []
    for i, n in enumerate(factors):
        rename = {v: len(connectors)*100+v+1 for v in adjacency(n) if v >= 0}
        g.update(edge(rename[a], rename[b]) for a, b in n if min(a,b) >= 0)
        ends = [rename[next(iter(adjacency(n)[-j-1]))] for j in range(5)]
        connectors.append(ends)
    for i in range(3):
        a, b = connectors[i][2:4]
        c, d = connectors[(i+1)%3][:2]
        if i == 2 and twist:
            c, d = d, c
        g.update((edge(a,c), edge(b,d), edge(0,connectors[i][4])))
    return g, connectors


def double_star():
    # Sage official DoubleStarSnark adjacency constructor, read 2026-10-20.
    d = {0:[1,14,15],1:[0,2,11],2:[1,3,7],3:[2,4,18],4:[3,5,14],
         5:[10,4,6],6:[5,21,7],7:[8,2,6],8:[9,13,7],9:[24,8,10],
         10:[9,11,5],11:[1,10,12],12:[11,27,13],13:[8,12,14],
         14:[0,4,13],15:[0,16,29],16:[15,20,23],17:[25,18,28],
         18:[3,17,19],19:[18,26,23],20:[16,28,21],21:[20,6,22],
         22:[26,21,29],23:[16,24,19],24:[25,9,23],25:[24,17,29],
         26:[27,19,22],27:[12,26,28],28:[17,27,20],29:[25,22,15]}
    return {edge(a,b) for a, ns in d.items() for b in ns}


def cyclic_cut_at_most(g, k=4):
    """Enumerate every deletion set of size 0..k using bitsets."""
    verts = sorted(adjacency(g))
    ix = {v:i for i,v in enumerate(verts)}
    edges = [(ix[a],ix[b]) for a,b in sorted(g)]
    n = len(verts)
    full = (1<<n)-1
    adj = [0]*n
    for a,b in edges:
        adj[a] |= 1<<b
        adj[b] |= 1<<a
    checked = 0
    by_size = [0]*(k+1)
    for removed in chain.from_iterable(combinations(range(len(edges)),size)
                                       for size in range(k+1)):
        by_size[len(removed)] += 1
        if len(removed)==k:
            checked += 1
        work = adj[:]
        for j in removed:
            a,b = edges[j]
            work[a] ^= 1<<b
            work[b] ^= 1<<a
        rest = full
        cyclic = []
        while rest:
            found = rest & -rest
            frontier = found
            while frontier:
                reach = 0
                while frontier:
                    bit = frontier & -frontier
                    frontier ^= bit
                    reach |= work[bit.bit_length()-1]
                frontier = reach & ~found
                found |= frontier
            rest &= ~found
            degree_sum = 0
            bits = found
            while bits:
                bit = bits & -bits
                bits ^= bit
                degree_sum += (work[bit.bit_length()-1] & found).bit_count()
            if degree_sum//2 >= found.bit_count():
                cyclic.append(found)
        if len(cyclic) >= 2:
            return False, checked, [sorted(g)[j] for j in removed],by_size
    return True, checked, None,by_size


def five_cut_witness(g,side):
    side = set(side)
    other = set(adjacency(g))-side
    cut = sorted(e for e in g if len(set(e)&side)==1)
    check(len(cut)==5,"five-cut size")
    for part in (side,other):
        h = nx.Graph([e for e in g if set(e)<=part])
        check(nx.is_connected(h) and h.number_of_edges()>=h.number_of_nodes(),
              "five-cut cyclic sides")
    return {"side":sorted(side),"edges":cut}


def canonical_word(w):
    order = {}
    return tuple(order.setdefault(c,len(order)) for c in w)


def signature(counts):
    orbit = {w:n for w,n in counts.items() if canonical_word(w) == w}
    check(sum(orbit.values())*6 == sum(counts.values()), "free colour action")
    defects = [0]*len(next(iter(counts)))
    for w,n in orbit.items():
        sizes = Counter(w)
        check(sorted(sizes.values()) in ([1,1,5],[1,3,3]), "seven-port parity")
        sign = 1 if sorted(sizes.values()) == [1,1,5] else -1
        for i,c in enumerate(w):
            if sizes[c] == 1:
                defects[i] += sign*n
    return {"support":len(counts), "mass":sum(counts.values()),
            "weights":sorted(set(counts.values())), "defects":defects}


def nn_pole(twist=0):
    g, connectors = loupekine(twist)
    removed = set(range(201, 211))
    result, cut = deletion_pole(g, removed)
    check(len(cut) == 5, "NN five-cut")
    return result, cut


def matching_deletions(edges, size):
    edges = tuple(edges)
    def visit(start, chosen, used):
        if len(chosen) == size:
            yield tuple(chosen)
            return
        if len(edges)-start < size-len(chosen):
            return
        for j in range(start,len(edges)):
            a,b = edges[j]
            if a not in used and b not in used:
                yield from visit(j+1,chosen+[edges[j]],used|{a,b})
    yield from visit(0,[],set())


def nn_minor():
    pole,cut = nn_pole()
    counts,_ = col(pole,5)
    tested = 0
    for deleted in matching_deletions([e for e in sorted(pole) if min(e)>=0],6):
        tested += 1
        retained = pole-set(deleted)
        if not nx.is_connected(nx.Graph(list(retained))):
            continue
        reduced = suppress(retained)
        smaller,_ = col(reduced,5)
        if set(smaller) == set(counts):
            check(all(counts[w] == 3*smaller[w] for w in counts), "NN scale")
            return {"deleted":deleted,"retained":sorted(retained),
                    "reduced":sorted(reduced),"tested":tested}
    raise RuntimeError("NN minor search failed")


def nn_cyclic_gate():
    pole,_ = nn_pole()
    adj = adjacency(pole)
    ports = [next(iter(adj[-i-1])) for i in range(5)]
    check(len(set(ports))==5,"distinct NN port owners")
    check(not any(edge(a,b) in pole for a,b in combinations(ports,2)),
          "independent NN port owners")
    vertices = sorted(v for v in adj if v>=0)
    ix = {v:i for i,v in enumerate(vertices)}
    masks = [sum(1<<ix[v] for v in e if v>=0) for e in pole]
    minima = [None]*16
    checked = 0
    for mask in range(1<<15):
        size = mask.bit_count()
        if size<3:
            continue
        checked += 1
        boundary = sum((mask&e).bit_count()==1 for e in masks)
        check(boundary>=5,"NN internal cyclic gate")
        minima[size] = boundary if minima[size] is None else min(minima[size],boundary)
    check(checked==32647,"all NN internal subsets size at least3")
    return {"port_owners":ports,"subsets_checked":checked,"minimum_by_size":minima}


def triple_orbits(g):
    graph = nx.Graph(list(g))
    automorphisms = list(nx.algorithms.isomorphism.GraphMatcher(graph,graph).
                         isomorphisms_iter())
    pending = {vs for vs in combinations(sorted(adjacency(g)),3)
               if sum(set(e)<=set(vs) for e in g) == 1}
    result = []
    while pending:
        triple = min(pending)
        orbit = {tuple(sorted(a[v] for v in triple)) for a in automorphisms}
        check(orbit <= pending, "triple orbit overlap")
        pending -= orbit
        result.append((triple,len(orbit)))
    return result


def triple_profiles(g):
    rows = []
    for triple,multiplicity in triple_orbits(g):
        pole,cut = deletion_pole(g,triple)
        counts,_ = col(pole,7)
        rows.append({"triple":triple,"orbit":multiplicity,"cut":cut,
                     **signature(counts)})
    return rows


def petersen_support_candidate_audit():
    words = sorted({canonical_word(w) for w in product(range(3),repeat=7)
                    if sorted(Counter(w).values()) in ([1,1,5],[1,3,3])})
    index = {w:i for i,w in enumerate(words)}
    check(len(index)==91,"seven-port orbit universe")
    base,_ = deletion_pole(P,(0,1,3))
    counts,_ = col(base,7)
    support = {canonical_word(w) for w in counts}
    masks = {sum(1<<index[canonical_word(tuple(w[i] for i in perm))] for w in support)
             for perm in permutations(range(7))}
    check(len(masks)==630,"all distinct port-permuted Petersen supports")
    rows = []
    for triple,multiplicity in triple_orbits(double_star()):
        pole,_ = deletion_pole(double_star(),triple)
        counts,_ = col(pole,7)
        mask = sum(1<<index[w] for w in counts if canonical_word(w)==w)
        eligible = sum(candidate&mask==candidate for candidate in masks)
        rows.append({"triple":triple,"multiplicity":multiplicity,"eligible":eligible,
                     "signature":signature(counts)})
    check(sum(r["multiplicity"] for r in rows if r["eligible"]==0)==100,
          "strict domination hostile triple census")
    check(next(r["eligible"] for r in rows if r["triple"]==(0,1,3))==2,
          "successful interface has two candidate supports")
    return {"candidate_supports":len(masks),"rows":rows}


def random_rooted_petersen(g, triple, attempts=10000, seed=1729):
    pole,cut = deletion_pole(g,triple)
    original,_ = col(pole,7)
    candidate = sorted(e for e in pole if min(e)>=0)
    required = (len(adjacency(g))-10)//2
    rng = random.Random(seed)
    stats = Counter()
    for attempt in range(attempts):
        rng.shuffle(candidate)
        chosen = []
        used = set()
        for e in candidate:
            if not used & set(e):
                chosen.append(e)
                used.update(e)
                if len(chosen) == required:
                    break
        if len(chosen) != required:
            continue
        stats["matching"] += 1
        retained = pole-set(chosen)
        if not nx.is_connected(nx.Graph(list(retained))):
            continue
        stats["connected"] += 1
        try:
            smaller = suppress(retained)
            closed = suppress(g-set(chosen))
        except RuntimeError:
            continue
        stats["simple"] += 1
        if not nx.is_isomorphic(nx.Graph(list(closed)),nx.Graph(list(P))):
            continue
        stats["Petersen"] += 1
        counts,_ = col(smaller,7)
        if set(counts) <= set(original):
            return {"triple":triple,"attempt":attempt+1,"seed":seed,
                    "deleted":sorted(chosen),"pole":sorted(pole),
                    "retained":sorted(retained),"reduced":sorted(smaller),
                    "closed":sorted(closed),"cut":cut,"stats":dict(stats),
                    "old_signature":signature(original),
                    "new_signature":signature(counts)}
    return {"triple":triple,"attempts":attempts,"stats":dict(stats)}


def trace_suppression_paths(retained):
    adj = adjacency(retained)
    branch = {v for v in adj if v < 0 or len(adj[v]) != 2}
    used = set()
    paths = []
    for a in sorted(branch):
        for b in sorted(adj[a]):
            if edge(a,b) in used:
                continue
            path = [a,b]
            used.add(edge(a,b))
            while path[-1] not in branch:
                prev,now = path[-2:]
                nxt = next(iter(adj[now]-{prev}))
                used.add(edge(now,nxt))
                path.append(nxt)
            paths.append(path)
    check(used == set(retained), "untraced cycle or edge")
    check(len({edge(p[0],p[-1]) for p in paths}) == len(paths), "parallel paths")
    check({edge(p[0],p[-1]) for p in paths} == suppress(retained), "path reduction")
    return paths


def nn_seven_lift(twist=0):
    g,_ = loupekine(twist)
    deleted = set(nn_minor()["deleted"])
    closed = suppress(g-deleted)
    check(nx.is_isomorphic(nx.Graph(list(closed)),nx.Graph(list(P))), "NN Petersen")
    # A genuine nonempty seven-port interface inside the third negator.
    triples = [t for t in combinations(sorted(v for v in adjacency(g) if v>=201),3)
               if sum(set(e)<=set(t) for e in g)==1]
    rows = []
    for triple in triples:
        original,cut = deletion_pole(g,triple)
        reduced = suppress(original-deleted)
        old,_ = col(original,7)
        new,_ = col(reduced,7)
        check(set(old)==set(new), "NN seven-port support")
        check(all(old[w]==3*new[w] for w in old), "NN seven-port scale")
        rows.append({"triple":triple,"cut":cut,"old":signature(old),
                     "new":signature(new)})
    check(bool(rows), "no seven-port lift")
    return {"twist":twist,"closed":sorted(closed),"rows":rows,
            "paths":trace_suppression_paths(nn_pole()[0]-deleted)}


def flower_nn_family(k):
    flower = PRE["OLD"]["flower"](k)
    g,_ = loupekine(0,[negator(),negator(),negator(flower,(1,0,2))])
    deleted = set(nn_minor()["deleted"])
    reduced = suppress(g-deleted)
    rename = {10:2,108:1,110:0}
    rename.update({v:v-201 for v in adjacency(reduced) if v>=201})
    check(set(rename)==set(adjacency(reduced)), "family restoration vertices")
    restored = {edge(rename[a],rename[b]) for a,b in reduced}
    check(restored == flower, "family exact flower restoration")
    return g


def cubic_multigraph_count(edges):
    """Independent perfect-matching/two-factor count; keeps edge multiplicity."""
    if any(a==b for a,b in edges):
        return 0
    adj = {}
    for j,(a,b) in enumerate(edges):
        adj.setdefault(a,[]).append((j,b))
        adj.setdefault(b,[]).append((j,a))
    check(all(len(ns)==3 for ns in adj.values()), "cubic closure degree")
    total = 0
    def visit(unmatched,chosen):
        nonlocal total
        if not unmatched:
            seen = set()
            cycles = 0
            for v in adj:
                if v in seen:
                    continue
                component = {v}
                stack = [v]
                while stack:
                    x = stack.pop()
                    for e,y in adj[x]:
                        if e not in chosen and y not in component:
                            component.add(y)
                            stack.append(y)
                if len(component)%2:
                    return
                cycles += 1
                seen.update(component)
            total += 2**cycles
            return
        v = min(unmatched,key=lambda x:sum(y in unmatched for _,y in adj[x]))
        for e,w in adj[v]:
            if w in unmatched:
                visit(unmatched-{v,w},chosen|{e})
    visit(set(adj),set())
    return total


def closure_bank():
    kernels = []
    for triple in combinations(range(7),3):
        rest = sorted(set(range(7))-set(triple))
        a = rest[0]
        for b in rest[1:]:
            c,d = sorted(set(rest)-{a,b})
            pole = {edge(0,-i-1) for i in triple}
            pole |= {edge(-a-1,-b-1),edge(-c-1,-d-1)}
            kernels.append(("J",pole))
    for i in range(7):
        owner = [0,0,1,2,3,4,4]
        owner[i],owner[3] = owner[3],owner[i]
        pole = {(0,1),(1,2),(2,3),(3,4)}
        pole.update(edge(-j-1,v) for j,v in enumerate(owner))
        kernels.append(("K"+str(i),pole))
    check(len(kernels)==112, "112 closure kernels")
    return kernels


def glue_poles(first,second):
    left = adjacency(first)
    right = adjacency(second)
    offset = max(v for v in left if v>=0)+1
    rename = {v:offset+i for i,v in enumerate(sorted(v for v in right if v>=0))}
    result = [e for e in first if min(e)>=0]
    for a,b in second:
        x = rename[a] if a>=0 else next(iter(left[a]))
        y = rename[b] if b>=0 else next(iter(left[b]))
        result.append(edge(x,y))
    return result


def independent_closure_audit(old,new):
    first,_ = col(old,7)
    second,_ = col(new,7)
    rows = []
    for label,kernel in closure_bank():
        values,_ = col(kernel,7)
        a = sum(n*values[w] for w,n in first.items())
        b = sum(n*values[w] for w,n in second.items())
        check(a==cubic_multigraph_count(glue_poles(old,kernel)), "old closure count")
        check(b==cubic_multigraph_count(glue_poles(new,kernel)), "new closure count")
        rows.append([label,a,b])
    return rows


def physical_kempe_audit(pole):
    counts_n = [0]*7
    def visit(word,colours):
        for i,colour in enumerate(word):
            if word.count(colour)!=1:
                continue
            h = nx.Graph()
            h.add_edges_from(e for e,c in colours.items() if c!=colour)
            ends = [sorted(v for v in component if v<0)
                    for component in nx.connected_components(h)
                    if any(v<0 for v in component)]
            check(len(ends)==3 and all(len(e)==2 for e in ends), "three Kempe paths")
            if all(word[-a-1]!=word[-b-1] for a,b in ends):
                counts_n[i] += 1
    counts,_ = col(pole,7,visit)
    defects = signature(counts)["defects"]
    check(counts_n==[-6*d for d in defects], "physical Kempe identity")
    return counts_n


DS_DELETED = {(4,5),(7,8),(9,24),(12,13),(15,29),(17,18),(19,23),
              (20,28),(21,22),(26,27)}


def verify_double_star_certificate():
    g = double_star()
    triple = (0,1,3)
    pole,cut = deletion_pole(g,triple)
    check(DS_DELETED <= pole, "certificate deletes existing internal edges")
    check(len(set(v for e in DS_DELETED for v in e))==20, "certificate matching")
    retained = pole-DS_DELETED
    reduced = suppress(retained)
    closed = suppress(g-DS_DELETED)
    iso = nx.algorithms.isomorphism.GraphMatcher(nx.Graph(list(closed)),nx.Graph(list(P)))
    check(iso.is_isomorphic(), "double-star target Petersen")
    original,_ = col(pole,7)
    smaller,_ = col(reduced,7)
    check(set(smaller) < set(original), "strict support domination")
    check(len(original)==192 and len(smaller)==120, "double-star support sizes")
    check(set(original.values())=={4,6} and set(smaller.values())=={1},
          "double-star multiplicities")
    remainder = Counter({w:original[w]-4*smaller[w] for w in original
                         if original[w]-4*smaller[w]})
    check(all(n>0 for n in remainder.values()), "nonnegative tensor remainder")
    check(signature(remainder)["defects"]==[-16,8,0,-8,0,0,0],
          "nonrealizable remainder hostile")
    return {"original_graph":sorted(g),"triple":triple,"cut":cut,
            "exterior_edges":[(0,1)],"deleted":sorted(DS_DELETED),
            "pole":sorted(pole),"retained":sorted(retained),"reduced":sorted(reduced),
            "paths":trace_suppression_paths(retained),"closed":sorted(closed),
            "petersen_vertex_map":iso.mapping,
            "old_signature":signature(original),"new_signature":signature(smaller),
            "old_physical_N":physical_kempe_audit(pole),
            "new_physical_N":physical_kempe_audit(reduced),
            "remainder4_signature":signature(remainder),
            "independent_closures":independent_closure_audit(pole,reduced)}


def main():
    destination = ROOT / "05-knowledge/results/black_descent_20261020c"
    destination.mkdir(exist_ok=True)
    print("NetworkX",nx.__version__,"exact labelled supports; independent cubic closures")
    print("Definitions: labelled ports -1..-7; proper edge colours 0,1,2; orbit factor6")
    check(cubic_multigraph_count([(0,1)]*3)==6,"parallel-edge positive control")
    check(cubic_multigraph_count([(0,0),(0,1),(1,1)])==0,"looped dumbbell hostile")
    a,b = loupekine(0)[0],loupekine(1)[0]
    check(not nx.is_isomorphic(nx.Graph(list(a)),nx.Graph(list(b))), "two distinct Loupekine")
    graph_audit = []
    for name,g in [("Loupekine0",loupekine(0)[0]),
                   ("Loupekine1",loupekine(1)[0]),("DoubleStar",double_star())]:
        cut_audit = cyclic_cut_at_most(g)
        check(cut_audit[0], "cyclic-five input")
        colourings,matchings = PRE["OLD"]["matching_colour_count"](g)
        check(colourings==0, "input snark")
        check(all(len(ns)==3 for ns in adjacency(g).values()),"input cubic")
        witness = five_cut_witness(g,range(15) if name=="DoubleStar" else range(4,11))
        graph_audit.append({"name":name,"edges":sorted(g),
                            "girth":PRE["girth"](nx.Graph(list(g))),
                            "C":colourings,"PM":matchings,
                            "cut_counts_sizes0to4":cut_audit[3],
                            "cyclic_five_cut":witness})
        print(name,"n",len(adjacency(g)),"girth",PRE["girth"](nx.Graph(list(g))),
              "C",colourings,"PM",matchings,"all_deletions_sizes0to4",cut_audit[3],flush=True)
    (destination/"snarks_graph_audit.json").write_text(json.dumps(graph_audit,indent=2)+"\n")
    dyad_counts,_ = col(negator(),5)
    check(len(dyad_counts)==36 and set(dyad_counts.values())=={1}, "Petersen dyad")
    print("Petersen negator: 36 boundary words, unique extensions")
    minor = nn_minor()
    (destination/"snarks_nn_minor.json").write_text(json.dumps(minor,indent=2)+"\n")
    print("NN rooted minor: f_NN=3 f_P3; deleted",minor["deleted"],flush=True)
    gate = nn_cyclic_gate()
    (destination/"snarks_nn_cyclic_gate.json").write_text(json.dumps(gate,indent=2)+"\n")
    print("NN universal cyclic-five gate",gate,flush=True)
    hosts = []
    for k in (5,7):
        host = PRE["OLD"]["flower"](k)
        cut_audit = cyclic_cut_at_most(host)
        check(cut_audit[0], "small flower host cyclic-five")
        hosts.append({"k":k,"cut_counts_sizes0to4":cut_audit[3]})
        print("Flower host J"+str(k),"all_deletions_sizes0to4",cut_audit[3],flush=True)
    (destination/"snarks_flower_hosts.json").write_text(json.dumps(hosts,indent=2)+"\n")
    for twist in range(2):
        lift = nn_seven_lift(twist)
        (destination/("snarks_loupekine"+str(twist)+"_lift.json")).write_text(
            json.dumps(lift,indent=2)+"\n")
        print("NN seven lift",twist,"cases",len(lift["rows"]),
              "first",lift["rows"][0],flush=True)
    for name,g in [("Loupekine0",loupekine(0)[0]),
                   ("Loupekine1",loupekine(1)[0]),("DoubleStar",double_star())]:
        rows = triple_profiles(g)
        labelled = sum(r["orbit"] for r in rows)
        check(labelled==len(g)*(len(adjacency(g))-6), "complete one-edge triple universe")
        print("SEVEN",name,"orbits",len(rows),"triples",labelled,
              "support_range",(min(r["support"] for r in rows),max(r["support"] for r in rows)),
              "support120_triples",sum(r["orbit"] for r in rows if r["support"]==120),flush=True)
        (destination / ("snarks_"+name+"_triples.json")).write_text(
            json.dumps(rows,indent=2)+"\n")
    # Every such Petersen pole is in one automorphism orbit and has support120.
    p_rows = triple_profiles(P)
    check(sum(r["orbit"] for r in p_rows)==60 and len(p_rows)==1 and
          p_rows[0]["support"]==120, "Petersen seven-pole benchmark")
    certificate = verify_double_star_certificate()
    print("DS certified minor: old",certificate["old_signature"],
          "new",certificate["new_signature"],flush=True)
    print("DS physical N:",certificate["old_physical_N"],"->",
          certificate["new_physical_N"],flush=True)
    print("DS independent closure counts: 112/112; extra seven",
          certificate["independent_closures"][-7:],flush=True)
    print("DS f_old-4f_new nonrealizable:",certificate["remainder4_signature"],flush=True)
    (destination/"snarks_double_star_minor.json").write_text(
        json.dumps(certificate,indent=2)+"\n")
    candidates = petersen_support_candidate_audit()
    (destination/"snarks_double_star_support_candidates.json").write_text(
        json.dumps(candidates,indent=2)+"\n")
    print("Petersen support candidates",candidates["candidate_supports"],
          "double-star triples with no dominated candidate",
          [(r["triple"],r["multiplicity"]) for r in candidates["rows"] if r["eligible"]==0],
          flush=True)
    family = []
    for k in (5,7,9,11):
        g = flower_nn_family(k)
        cut_audit = cyclic_cut_at_most(g)
        check(cut_audit[0], "NN-flower cyclic-five finite sample")
        colouring,matching = PRE["OLD"]["matching_colour_count"](g)
        check(colouring==0, "NN-flower noncolourable")
        row = {"k":k,"n":len(adjacency(g)),"edges":sorted(g),
               "C":colouring,"PM":matching,"all_4deletions":cut_audit[1],
               "cut_counts_sizes0to4":cut_audit[3],
               "cyclic_five_cut":five_cut_witness(g,{v for v in adjacency(g) if v<201})}
        family.append(row)
        print("NN-flower",k,"n",row["n"],"PM",matching,
              "all_deletions_sizes0to4",cut_audit[3],"exact minor restores J"+str(k),flush=True)
    (destination/"snarks_nn_flower_family.json").write_text(json.dumps(family,indent=2)+"\n")
    print("PASS: graph constructors, cyclic cuts, ordered supports, rooted paths, "
          "physical Kempe defects, independent 112 closure counts")


if __name__ == "__main__":
    main()

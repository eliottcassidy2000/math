#!/usr/bin/env python3
"""Exact observer, tree-gauge, deficient-slice, and path-incidence controls.
Standard library only. Every check remains live under python -O.
"""
from itertools import combinations, permutations, product
from collections import Counter, defaultdict
from functools import lru_cache
import hashlib, json
from pathlib import Path

CHECKS = 0

def require(condition, message):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(message)

class Space:
    def __init__(self, n):
        self.n = n
        self.edges = tuple(combinations(range(n), 2))
        self.idx = {e: i for i, e in enumerate(self.edges)}
        self.N = len(self.edges)
        self.full = (1 << self.N) - 1
        self.perms = tuple(permutations(range(n)))
        self.cuts = tuple(self.cut(s) for s in range(1 << (n-1)))
        self.class_of = tuple(min(x ^ c for c in self.cuts) for x in range(1 << self.N))
        self.classes = defaultdict(list)
        for x, k in enumerate(self.class_of):
            self.classes[k].append(x)
    def bit(self, u, v):
        return 1 << self.idx[tuple(sorted((u,v)))]
    def cut(self, vertex_mask):
        return sum(1 << i for i,(u,v) in enumerate(self.edges) if ((vertex_mask>>u) ^ (vertex_mask>>v)) & 1)
    def arc(self, x, u, v):
        return bool((u < v) ^ bool(x & self.bit(u,v)))
    def components(self, h):
        roots = list(range(self.n))
        def root(u):
            while roots[u] != u:
                u = roots[u]
            return u
        for i,(u,v) in enumerate(self.edges):
            if h >> i & 1:
                roots[root(u)] = root(v)
        return len({root(u) for u in range(self.n)})
    def order_mask(self, order):
        pos = {v:i for i,v in enumerate(order)}
        return sum(1 << i for i,(u,v) in enumerate(self.edges) if pos[u] > pos[v])
    def relabel_graph(self, x, perm):
        return sum(self.bit(perm[u],perm[v]) for i,(u,v) in enumerate(self.edges) if x >> i & 1)
    def relabel_tournament(self, x, perm):
        return self.relabel_graph(x, perm) ^ self.order_mask(perm)
    def path_mask(self, p):
        return sum(self.bit(u,v) for u,v in zip(p,p[1:]))
    def paths(self, x):
        return tuple(p for p in self.perms if all(self.arc(x,u,v) for u,v in zip(p,p[1:])))
    def normalize(self, x, tree):
        adjacency = [[] for _ in range(self.n)]
        for i,(u,v) in enumerate(self.edges):
            if tree >> i & 1:
                adjacency[u].append((v,i))
                adjacency[v].append((u,i))
        potentials = {0:0}
        stack = [0]
        while stack:
            u = stack.pop()
            for v,i in adjacency[u]:
                if v not in potentials:
                    potentials[v] = potentials[u] ^ ((x>>i)&1)
                    stack.append(v)
        if len(potentials) != self.n:
            raise ValueError('Normalization requested on disconnected graph')
        vertex_mask = sum(z << v for v,z in potentials.items())
        return x ^ self.cut(vertex_mask)
    def shape(self, h):
        deg = sorted(sum(bool(h & self.bit(u,v)) for v in range(self.n) if v != u) for u in range(self.n))
        return {(1,1,2,2):'P4',(1,1,1,3):'star',(0,2,2,2):'triangle+isolate'}.get(tuple(deg),'other')

summary = {'status':'PASS', 'scope':'all finite claims in note; finite controls as stated', 'orders':[]}
for n in range(2,6):
    s = Space(n)
    qdim = s.N-n+1
    require(len(s.classes) == 1 << qdim, 'wrong number of switching classes')
    require(all(len(v)==1<<(n-1) for v in s.classes.values()), 'wrong orbit size')
    trees = [h for h in range(1<<s.N) if h.bit_count()==n-1 and s.components(h)==1]
    path = s.path_mask(tuple(range(n)))
    star = sum(s.bit(0,v) for v in range(1,n))
    constraints_checked = 0
    for h in range(1<<s.N):
        m,c = h.bit_count(),s.components(h)
        beta = m-n+c
        # Every pattern for n<=4, three deterministic patterns at n=5.
        if n<=4:
            indices = [i for i in range(s.N) if h>>i&1]
            patterns = [sum(((a>>j)&1)<<i for j,i in enumerate(indices)) for a in range(1<<m)]
        else:
            patterns = sorted({0,h,h & 0x155})
        for a in patterns:
            observed = Counter(s.class_of[x] for x in range(1<<s.N) if x & h == a)
            require(all(v==1<<(c-1) for v in observed.values()), f'intersection multiplicity n={n},h={h},a={a}')
            require(len(observed)==1<<(s.N-m-c+1), 'image class count')
            require(len(observed)*(1<<beta)==len(s.classes), 'cycle codimension count')
            constraints_checked += 1
    normalized = {}
    for tree in trees:
        for representative in s.classes:
            z = s.normalize(representative,tree)
            require(z & tree == 0, 'tree not normalized')
            require(s.class_of[z] == representative, 'normalization changed switching class')
            require(all(s.normalize(x,tree)==z for x in s.classes[representative]), 'orbit representative nonunique')
            normalized[tree,representative] = z
    pairs = product(trees,trees) if n<=4 else ((tree,target) for tree in trees for target in {path,star})
    transition_checks = 0
    for tree,target in pairs:
        for representative in s.classes:
            x = normalized[tree,representative]
            y = s.normalize(x,target)
            require(y==normalized[target,representative], 'transition class mismatch')
            require(s.normalize(y,tree)==x, 'transition not invertible')
            # Translation of free-host complement through change of tree.
            shifted = s.normalize(x ^ (s.full ^ tree),target)
            delta = s.normalize(s.full ^ tree,target)
            require(shifted == y ^ delta, 'complement conjugation not linear')
            transition_checks += 1
    all_paths = {}
    ordered_incidence = set()
    for x in range(1<<s.N):
        paths = s.paths(x)
        all_paths[x]=paths
        for p in paths:
            y=x ^ (s.full ^ s.path_mask(p))
            require((x,y) not in ordered_incidence, 'path witness not unique')
            ordered_incidence.add((x,y))
            require(p in s.paths(y), 'path lost under offpath complement')
        # Exact cyclic-triangle polynomial.
        polynomial=0
        direct=0
        for i,j,k in combinations(range(n),3):
            a=bool(x&s.bit(i,j)); b=bool(x&s.bit(j,k)); c=bool(x&s.bit(i,k))
            polynomial += c+a*b-a*c-b*c
            direct += ((s.arc(x,i,j) and s.arc(x,j,k) and s.arc(x,k,i)) or
                       (s.arc(x,j,i) and s.arc(x,k,j) and s.arc(x,i,k)))
        require(polynomial==direct, 'cyclic-triangle formula failed')
        for perm in s.perms:
            y=s.relabel_tournament(x,perm)
            require(all(s.arc(y,perm[u],perm[v])==s.arc(x,u,v) for u,v in s.edges), 'affine relabelling formula')
            require((x.bit_count()^y.bit_count())&1 == s.order_mask(perm).bit_count()&1, 'orientation parity cocycle')
    require(all((y,x) in ordered_incidence for x,y in ordered_incidence), 'path incidence asymmetric')
    # Independently compare the literal incidence quotient with the generator law.
    actual_rows = defaultdict(Counter)
    for x,y in ordered_incidence:
        actual_rows[s.class_of[x]][s.class_of[y]] += 1
    generators = Counter(s.class_of[s.full ^ s.path_mask(p)] for p in s.perms)
    for representative in s.classes:
        expected = Counter()
        for generator,multiplicity in generators.items():
            expected[s.class_of[representative ^ generator]] += multiplicity
        require(actual_rows[representative] == expected, 'weighted Cayley incidence quotient')
    unseen=set(s.classes)
    quotient_components=[]
    while unseen:
        root=unseen.pop()
        reached={root}
        stack=[root]
        while stack:
            x=stack.pop()
            for generator in generators:
                y=s.class_of[x^generator]
                if y not in reached:
                    reached.add(y)
                    unseen.discard(y)
                    stack.append(y)
        quotient_components.append(len(reached))
    expected_components=2 if n>=3 and n%4==1 else 1
    require(len(quotient_components)==expected_components, 'quotient connectivity boundary')
    if n>=3:
        require(0 not in generators, 'path-complement generator collapsed to a cut')
    require(len(ordered_incidence)==len(s.perms)*(1<<qdim), 'global path incidence count')
    # Spanning-tree switching sum, applied to the whole directed path family.
    require(all(sum(len(all_paths[x]) for x in members)==len(s.perms) for members in s.classes.values()), 'switching-class H sum')
    if n<=4:
        for p in s.perms:
            bp=s.order_mask(p)
            for q in s.perms:
                bq=s.order_mask(q)
                inv=bp^bq
                for x in range(1<<s.N):
                    require((x^bp)^inv == x^bq, 'observer cocycle failed')
    row={'n':n,'labelled_tournaments':1<<s.N,'labelled_switching_classes':len(s.classes),'trees':len(trees),
         'constraint_assignments':constraints_checked,'tree_transition_states':transition_checks,
         'directed_path_complement_incidences':len(ordered_incidence),'H_histogram':dict(sorted(Counter(map(len,all_paths.values())).items())),
         'cut_quotient_components':sorted(quotient_components),'cut_quotient_generator_multiplicities':dict(sorted(Counter(generators.values()).items()))}
    if n==3:
        switched=s.cut(1<<1)
        require(len(all_paths[0])==1 and len(all_paths[switched])==3, 'H-changing switching hostile')
        row['H_switching_hostile']={'source_mask':0,'switch_vertex':1,'target_mask':switched,'H':[1,3]}
    if n==4:
        require(Counter(s.shape(h) for h in range(1<<s.N) if h.bit_count()==3)=={'P4':12,'star':4,'triangle+isolate':4}, 'three shape census')
        hosts={kind:next(h for h in range(1<<s.N) if h.bit_count()==3 and s.shape(h)==kind) for kind in ['P4','star','triangle+isolate']}
        row['host_slice_coverage']={kind:dict(sorted(Counter(Counter(s.class_of[x] for x in range(1<<s.N) if x & (s.full^host)==0).values()).items())) for kind,host in hosts.items()}
        row['path_to_star_table']=[]
        for a,b,c in product(range(2),repeat=3):
            x=a*s.bit(0,2)+b*s.bit(0,3)+c*s.bit(1,3)
            y=s.normalize(x,star)
            bits=tuple(int(bool(y&s.bit(*edge))) for edge in [(1,2),(1,3),(2,3)])
            require(bits==(a,b^c,a^b), 'literal path-to-star map')
            shifted=s.normalize(x^(s.full^path),star)
            require(shifted==y^s.bit(1,2), 'path complement not one triangle-edge flip')
            row['path_to_star_table'].append({'path_bits':[a,b,c],'star_bits':bits})
        require(all(s.shape(x^y)=='P4' for x,y in ordered_incidence), 'non-P4 path complement support')
        require(len(generators)==6 and set(generators.values())=={4}, 'n4 quotient not six generators of weight four')
        require(set(generators)==set(s.classes)-{0,s.class_of[s.full]}, 'n4 omitted matching generator wrong')
        shape_classes={shape:Counter(s.class_of[h] for h in range(1<<s.N) if h.bit_count()==3 and s.shape(h)==shape) for shape in ['P4','star','triangle+isolate']}
        require(set(shape_classes['star'])=={0} and set(shape_classes['triangle+isolate'])=={s.class_of[s.full]}, 'star/triangle cut classes wrong')
        require(set(shape_classes['P4'])==set(generators) and set(shape_classes['P4'].values())=={2}, 'path cut-class multiplicity wrong')
        row['three_edge_switching_classes']={shape:dict(sorted(count.items())) for shape,count in shape_classes.items()}
        row['path_incidence_difference_shapes']=dict(Counter(s.shape(x^y) for x,y in ordered_incidence))
        row['three_edge_graph_tournament_H']=dict((shape,dict(sorted(Counter(len(all_paths[h]) for h in range(1<<s.N) if h.bit_count()==3 and s.shape(h)==shape).items()))) for shape in ['P4','star','triangle+isolate'])
    summary['orders'].append(row)
# Larger-order rank controls enumerate path generators only, not all tournaments.
summary['generator_rank_controls']=[]
for n in range(3,9):
    edges=tuple(combinations(range(n),2))
    idx={e:i for i,e in enumerate(edges)}
    full=(1<<len(edges))-1
    basis={}
    difference_basis={}
    first_generator=None
    weights=Counter()
    def bit(u,v): return 1<<idx[tuple(sorted((u,v)))]
    for p in permutations(range(n)):
        if p[0]>p[-1]: continue  # one representative per undirected path
        support=full ^ sum(bit(u,v) for u,v in zip(p,p[1:]))
        # A direct fundamental-triangle encoder, independent of normalize().
        q=sum(1<<k for k,(i,j) in enumerate(combinations(range(1,n),2))
              if (bool(support&bit(i,j)) ^ bool(support&bit(0,i)) ^ bool(support&bit(0,j))))
        require(q!=0, 'zero generator in larger-order check')
        weights[q]+=2
        if first_generator is None: first_generator=q
        z=q^first_generator
        while z:
            pivot=z.bit_length()-1
            if pivot in difference_basis: z ^= difference_basis[pivot]
            else:
                difference_basis[pivot]=z
                break
        z=q
        while z:
            pivot=z.bit_length()-1
            if pivot in basis: z ^= basis[pivot]
            else:
                basis[pivot]=z
                break
    dim=(n-1)*(n-2)//2
    expected_rank=dim-(n%4==1)
    require(len(basis)==expected_rank, 'all-order generator-span rank check')
    expected_difference_rank=dim-(n%2==1)
    require(len(difference_basis)==expected_difference_rank, 'all-order difference-span rank check')
    require((len(basis)==len(difference_basis)+1)==(n%4==3), 'bipartiteness rank condition')
    summary['generator_rank_controls'].append({'n':n,'quotient_dimension':dim,'rank':len(basis),'difference_rank':len(difference_basis),
         'distinct_generators':len(weights),'generator_multiplicity_histogram':dict(sorted(Counter(weights.values()).items()))})
summary['checks']=CHECKS
summary['source_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
print(json.dumps(summary,indent=2,sort_keys=True))

"""Exact native controls for KERNELS.md and MOVES.md; standard library only."""
from collections import Counter
from functools import lru_cache
from itertools import combinations, product
from pathlib import Path
import json
import random
import sys

HERE = Path(__file__).resolve().parent
CHECKS = 0


def need(value, message):
    global CHECKS
    CHECKS += 1
    if not value:
        raise ValueError(message)


def words(d):
    return tuple(product(range(3), repeat=d))


def arc(x, y):
    for a, b in zip(x, y):
        if a != b:
            return (b-a) % 3 == 1
    return False


def tree_perm(d, turns):
    bank = words(d)
    index = {x: i for i, x in enumerate(bank)}
    return tuple(index[tuple((x[k]+turns[x[:k]]) % 3 for k in range(d))]
                 for x in bank)


def preserves(A, p):
    return all(A[i][j] == A[p[i]][p[j]]
               for i in range(len(A)) for j in range(len(A)))


def components(m, edges):
    parent = list(range(m))
    def root(i):
        while parent[i] != i:
            i = parent[i]
        return i
    for i, j in edges:
        parent[root(i)] = root(j)
    return len({root(i) for i in range(m)})


def hierarchy_controls():
    records = []
    rng = random.Random(3127)
    for d in range(1, 5):
        bank = words(d)
        A = tuple(tuple(arc(x, y) for y in bank) for x in bank)
        marks = tuple(i for i, x in enumerate(bank) if x[-1] == 0)
        profiles = [tuple((A[i][j], A[j][i]) for j in marks)
                    for i in range(len(bank)) if i not in marks]
        need(len(profiles) == len(set(profiles)), 'internal marks distinguish all unmarked vertices')
        hist = Counter()
        for i, j in combinations(range(len(bank)), 2):
            r = next(k for k in range(d) if bank[i][k] != bank[j][k])
            common = sum(A[i][k] and A[j][k] for k in range(len(bank)))
            need(common == (3**d-2*3**(d-r-1)-1)//2, 'common-out hierarchy')
            hist[common] += 1
        prefixes = tuple(q for r in range(d) for q in words(r))
        for _ in range(30):
            p = tree_perm(d, {q: rng.randrange(3) for q in prefixes})
            need(len(set(p)) == 3**d and preserves(A, p), 'tree rotation is actual tournament symmetry')
        records.append(dict(depth=d, vertices=len(bank), marks=len(marks),
                            common_out_histogram=sorted(hist.items())))
    return records


def nine_networks():
    d = 2
    bank = words(d)
    prefixes = tuple(q for r in range(d) for q in words(r))
    G = tuple(tree_perm(d, dict(zip(prefixes, t))) for t in product(range(3), repeat=4))
    need(len(set(G)) == 81, 'all 81 source-prefix normal forms distinct')
    identity = G.index(tuple(range(9)))
    # Actual full block test, independent of the prefix-parameter formula.
    relations = []
    for level in range(3):
        B = tuple(tuple(x[:level] == y[:level] for y in bank) for x in bank)
        H = {(i, j) for i, g in enumerate(G) for j, h in enumerate(G)
             if all(B[u][v] == B[g[u]][h[v]] for u in range(9) for v in range(9))}
        relations.append(H)
        need(len(H) == (6561, 2187, 81)[level], 'actual block stabilizer order')
    maximum = {(i, j): max(k for k in range(3) if (i, j) in relations[k])
               for i in range(81) for j in range(81)}
    buckets = Counter((maximum[i,j], maximum[i,k], maximum[j,k])
                      for i, j, k in product(range(81), repeat=3))
    edges = ((0,1),(0,2),(1,2))
    census = []
    for labels in product(range(3), repeat=3):
        count = sum(v for key, v in buckets.items() if all(x >= y for x,y in zip(key, labels)))
        cs = [components(3, [e for e,l in zip(edges,labels) if l>=k]) for k in (1,2)]
        need(count == 3**(cs[0]+3*cs[1]), 'nested network exact count')
        census.append(dict(labels=list(labels), components=cs, frames=count, owner_base=3*cs[1]))
    # Mixed path has a nonfunctional first edge and a functional second edge.
    mixed = [(i,j,k) for i,j in relations[1] for k in range(81) if (j,k) in relations[2]]
    kernels = []
    for cell in range(3):
        kernels.append({g[cell] for g in mixed if all(g[j] == identity for j in range(3) if j != cell)})
    need(list(map(len,kernels)) == [27,1,1], 'intersect incident one-sided kernels')
    # Right cosets from literal permutation multiplication.
    indices = {g:i for i,g in enumerate(G)}
    qmaps=[]
    for K in kernels:
        mapping={}
        for gi,g in enumerate(G):
            if gi not in mapping:
                coset={indices[tuple(g[G[ki][u]] for u in range(9))] for ki in K}
                label=min(coset)
                mapping.update({x:label for x in coset})
        qmaps.append(mapping)
    fibres=Counter(tuple(qmaps[j][g[j]] for j in range(3)) for g in mixed)
    need(len(fibres)==81 and set(fibres.values())=={27}, 'every accepted quotient tuple has 27 lifts')
    marked=((0,0),(0,3),(0,6),(1,0),(1,3),(1,6))
    residual=[g for g in mixed if all(G[g[c]][u]==u for c,u in marked)]
    need(len(residual)==1, 'six actual owners kill mixed path group')
    # Lower-bound witnesses: six disjoint supports of independent bottom rotations.
    witnesses=[]
    for cells in ((0,), (1,2)):
        for a in range(3):
            turns={q:0 for q in prefixes}; turns[(a,)]=1
            gi=indices[tree_perm(2,turns)]
            g=tuple(gi if c in cells else identity for c in range(3))
            need(g in mixed, 'independent lower-bound symmetry')
            support={(c,u) for c in range(3) for u in range(9) if G[g[c]][u]!=u}
            need(all(not support.intersection(s) for s in witnesses), 'disjoint support lower bounds')
            witnesses.append(support)
    return dict(networks=census, mixed=dict(frames=len(mixed), kernel_sizes=list(map(len,kernels)),
                                           quotient_frames=len(fibres), owner_base=6))


def matched_chain(m):
    def edge(x,y):
        i,a=x; j,b=y
        if i==j: return (b-a)%3==1
        if i>j: return not edge(y,x)
        if j==i+1: return a==b
        return True
    paths=[]
    for t in range(3):
        path=[(i,(t-i+k)%3) for i in range(m) for k in range(3)]
        need(len(set(path))==3*m and all(edge(x,y) for x,y in zip(path,path[1:])), 'actual Hamiltonian path')
        paths.append(path)
    for a,b in product(range(3),repeat=2):
        if a==b: continue
        move=tuple((int(k==b)-int(k==a)) for i in range(m) for k in range(3))
        need(sum(v!=0 for v in move)==2*m and sum(map(abs,move))==2*m, 'wide primitive move')
    need(len({p[0] for p in paths})==3, 'indexed owner distinguishes paths')
    if m<=8:
        feasible=[t for t in product(range(3),repeat=m) if all(a==b for a,b in zip(t,t[1:]))]
        need(feasible==[tuple([a]*m) for a in range(3)], 'one-hot assignment exhaustive')
    # Conformal primitive bank in a root-coordinate box: all zero-sum triples.
    bank=[q for q in product(range(-2,3),repeat=3) if sum(q)==0 and any(q)]
    primitives=[]
    for q in bank:
        below=lambda h: all(hv*qv>=0 and abs(hv)<=abs(qv) for hv,qv in zip(h,q))
        if not any(h!=q and below(h) for h in bank): primitives.append(q)
    need(len(primitives)==6 and all(sum(map(abs,q))==2 for q in primitives), 'all primitive root transfers')
    return dict(cells=m, native_paths=3, atom_support=2*m, root_states=3)


RULE={'C':'CA','A':'B','B':'CA'}
LETTER={'C':(1,(1,2,0),(0,1,2)), 'A':(1,(1,0,2),(2,0,-1)), 'B':(1,(0,1,2),(1,1,-2))}
EMPTY=(0,(0,1,2),(0,0,0))


def concatenate(U,V):
    n,p,w=U; m,q,v=V
    return n+m, tuple(q[p[t]] for t in range(3)), tuple(w[t]+v[p[t]] for t in range(3))


@lru_cache(None)
def summary(letter,depth):
    if depth==0: return LETTER[letter]
    out=EMPTY
    for c in RULE[letter]: out=concatenate(out,summary(c,depth-1))
    return out


def prefix(letter,depth,k):
    if not 0<=k<=summary(letter,depth)[0]: raise ValueError('invalid address')
    if k==0: return EMPTY
    if k==summary(letter,depth)[0]: return summary(letter,depth)
    if depth==0: return LETTER[letter]
    out=EMPTY
    for c in RULE[letter]:
        size=summary(c,depth-1)[0]
        take=min(k,size)
        out=concatenate(out,prefix(c,depth-1,take))
        k-=take
        if not k: break
    return out


def address(letter,depth,k):
    if depth==0: return letter
    for c in RULE[letter]:
        size=summary(c,depth-1)[0]
        if k<size: return address(c,depth-1,k)
        k-=size
    raise ValueError('invalid interior slot')


def fibonacci_controls():
    word='C'
    for d in range(13):
        value=EMPTY
        need(prefix('C',d,0)==EMPTY,'empty prefix')
        for k,c in enumerate(word,1):
            value=concatenate(value,LETTER[c])
            need(prefix('C',d,k)==value,'literal versus compressed owned prefix')
        need(summary('C',d)==value,'whole summary')
        word=''.join(RULE[c] for c in word)
    depth=100; n=summary('C',depth)[0]
    for k in (0,1,2,n//3,n//2,n-2,n-1):
        need(prefix('C',depth,k+1)==concatenate(prefix('C',depth,k),LETTER[address('C',depth,k)]),
             'large exact address continuation')
    zero=(0,(1,0,2),(0,0,0))
    need(zero[0]==EMPTY[0] and zero[1]!=EMPTY[1], 'zero metric length does not erase phase event')
    return dict(literal_depth=12, compressed_depth=depth, compressed_slots=n,
                summary=summary('C',depth), zero_length_before_after=[EMPTY[1],zero[1]])


def main():
    out=dict(status='FINITE-EXACT scoped controls', hierarchy=hierarchy_controls(),
             nine=nine_networks(), chains=[matched_chain(m) for m in list(range(1,9))+[16,31,64]],
             fibonacci=fibonacci_controls())
    out['checks_without_receipt']=CHECKS
    text=json.dumps(out,indent=2,sort_keys=True)+'\n'
    target=HERE/'kernels_moves.json'
    if '--write' in sys.argv: target.write_text(text,encoding='utf8')
    else: need(json.loads(target.read_text(encoding='utf8'))==json.loads(text),'frozen receipt')
    print(json.dumps(dict(status='PASS',checks=CHECKS,nine_networks=len(out['nine']['networks']),
                         mixed=out['nine']['mixed'],compressed_slots=out['fibonacci']['compressed_slots'])))


if __name__=='__main__': main()

"""Independent literal permutation controls for the scoped decoder."""
import itertools as it
import json
import random
from decoder import decode, stabilizer

checks=0
def require(test, message):
    global checks
    checks+=1
    if not test: raise RuntimeError(message)

def preserves(rows,g):
    return all(((rows[i]>>j)&1)==((rows[g[i]]>>g[j])&1)
               for i in range(len(rows)) for j in range(i))

def construct(sizes,crossings):
    cells=[]; start=0
    for p in sizes:
        cells.append(list(range(start,start+p)));start+=p
    rows=[0]*start
    for C in cells:
        p=len(C)
        for x in range(p):
            for d in range(1,(p-1)//2+1):rows[C[x]]|=1<<C[(x+d)%p]
    for i,X in enumerate(cells):
        for j in range(i+1,len(cells)):
            Y=cells[j];B=crossings[i,j]
            for x,a in enumerate(X):
                for y,b in enumerate(Y):
                    if B[x][y]: rows[a]|=1<<b
                    else: rows[b]|=1<<a
    return rows,cells

def literal_rotations(rows,cells):
    answer=[]
    for shifts in it.product(*(range(len(C)) for C in cells)):
        g=list(range(len(rows)))
        for C,k in zip(cells,shifts):
            for i,v in enumerate(C):g[v]=C[(i+k)%len(C)]
        if preserves(rows,g): answer.append(g)
    return answer

def verify(rows,cells):
    d=decode(rows)
    require({frozenset(c) for c in d['cells']}=={frozenset(c) for c in cells},'stable partition')
    autos=literal_rotations(rows,cells)
    require(d['group_order']==len(autos),'group count')
    require(all(preserves(rows,g) for g in d['generators']),'generator arcs')
    generated={tuple(range(len(rows)))}
    for g,(p,_) in zip(d['generators'],d['free_components']):
        expanded=set();power=list(range(len(rows)))
        for _ in range(p):
            expanded.update(tuple(power[v] for v in h) for h in generated)
            power=[g[v] for v in power]
        generated=expanded
    require(generated==set(map(tuple,autos)),'full generated permutation group')
    orbits=[{g[v] for g in autos} for v in range(len(rows))]
    require(d['orbit_sizes']==list(map(len,orbits)),'literal vertex orbits')
    distinct={frozenset(o) for o in orbits}
    e=sum(next(len(C) for C in cells if next(iter(o)) in C)-len(o) for o in distinct)
    require(e==d['error'],'error')
    return d

def matrix(p,perm):return [[int(perm[x]==y) for y in range(p)] for x in range(p)]

def main():
    for p in (3,5,7):
        rows,cells=construct([p],{})
        auto=[g for g in it.permutations(range(p)) if preserves(rows,g)]
        require(len(auto)==p,'full internal permutation action')
        verify(rows,cells)
    five=[]
    for degree in (1,2):
        choices=list(it.combinations(range(5),degree))
        for rows in it.product(choices,repeat=5):
            if any(sum(j in r for r in rows)!=degree for j in range(5)):continue
            B=[[int(j in r) for j in range(5)] for r in rows]
            five.append(B)
    require(len(five)==2160,'complete labelled five crossing count')
    hist={}
    for B in five:
        rows,cells=construct([5,5],{(0,1):B});d=verify(rows,cells)
        hist[d['error']]=hist.get(d['error'],0)+1
        converse=[((1<<len(rows))-1)^r^(1<<i) for i,r in enumerate(rows)]
        require(decode(converse)['error']==d['error'],'converse')
    require(hist=={0:60,40:2100},'five crossing error census')
    perms=list(it.permutations(range(3)))
    triangles={}
    for pp in it.product(perms,repeat=3):
        rows,cells=construct([3]*3,{e:matrix(3,p) for e,p in zip(((0,1),(0,2),(1,2)),pp)})
        d=verify(rows,cells);triangles[d['error']]=triangles.get(d['error'],0)+1
    require(triangles=={0:108,18:108},'triangle census')
    rng=random.Random(510051)
    degree2=[B for B in five if sum(B[0])==2]
    for _ in range(50):
        rows,cells=construct([5]*3,{e:rng.choice(degree2) for e in ((0,1),(0,2),(1,2))})
        verify(rows,cells)
    rows,cells=construct([3,5],{(0,1):[[1]*5 for _ in range(3)]})
    require(verify(rows,cells)['group_order']==15,'mixed prime product')
    # A bad regular-seven cell is prime and regular, but NOT consecutive cyclic.
    bad=[74,84,97,22,37,11,56]
    try:decode(bad)
    except ValueError:pass
    else:raise RuntimeError('accepted bad regular seven outside domain')
    require(True,'domain rejection')
    # A single rigid crossing pins both cells despite acyclic interaction.
    B=matrix(5,[0,1,2,4,3]);rows,cells=construct([5,5],{(0,1):B})
    d=verify(rows,cells)
    require(d['error']==40 and any(e[2]=='rigid' for e in d['edges']),'tree rigidity')
    for _ in range(30):
        permutation=list(range(len(rows)));rng.shuffle(permutation)
        relabelled=[0]*len(rows)
        for i in range(len(rows)):
            for j in range(len(rows)):
                if (rows[i]>>j)&1:relabelled[permutation[i]]|=1<<permutation[j]
        dd=decode(relabelled)
        require(dd['error']==d['error'] and dd['group_order']==d['group_order'],'relabel invariance')
    print(json.dumps({'status':'PASS','checks':checks,'five_labelled_error_histogram':hist,
                      'triangle_labelled_error_histogram':triangles},sort_keys=True))

if __name__=='__main__':main()

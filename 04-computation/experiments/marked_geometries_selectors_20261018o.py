#!/usr/bin/env python3
"""Marked bridges at 5, 27, 21, 9; boundary-selector obstruction and list gate.

Finite geometries are not declared planar unit-distance realizations. Every
incidence bijection is output and replayed, with original labels retained.
"""
import argparse
from collections import Counter
from itertools import combinations, permutations, product
import json
from pathlib import Path
import networkx as nx

CHECKS=0


def check(ok, why):
    global CHECKS
    CHECKS+=1
    if not ok:raise RuntimeError(why)


def span(rows):
    s={0}
    for r in rows:s|={x^r for x in s}
    return s


def q(x):
    b=[(x>>i)&1 for i in range(6)]
    return (b[0]*b[1]+b[2]*b[3]+b[4]+b[4]*b[5]+b[5])%2


def polar(x,y):return q(x^y)^q(x)^q(y)


def graph(points, lines):
    g=nx.Graph()
    index={p:i for i,p in enumerate(points)}
    for i in range(len(points)):g.add_node(i,part=0)
    for i,line in enumerate(lines):
        g.add_node(len(points)+i,part=1)
        for p in line:g.add_edge(index[p],len(points)+i)
    return g


def incidence_map(P,L,P2,L2):
    gm=nx.algorithms.isomorphism.GraphMatcher(graph(P,L),graph(P2,L2),
        node_match=lambda a,b:a['part']==b['part'])
    check(gm.is_isomorphic(),'incidence isomorphism exists')
    m=gm.mapping; f={p:P2[m[i]] for i,p in enumerate(P)}
    check({frozenset(f[p] for p in line) for line in L}==set(map(frozenset,L2)),
          'replay all labelled incidence blocks')
    return f


def fm(a,b):
    v=0
    while b:
        if b&1:v^=a
        a<<=1
        if a&4:a^=7
        b>>=1
    return v


def normalize(v):
    a=next(x for x in v if x)
    ai=next(b for b in range(1,4) if fm(a,b)==1)
    return tuple(fm(ai,x) for x in v)


def dot(a,b):
    r=0
    for x,y in zip(a,b):r^=fm(x,y)
    return r


def mobius(a,b,c,d):
    def val(x):
        if x==11:return a*pow(c,-1,11)%11 if c else 11
        y=(c*x+d)%11
        return (a*x+b)*pow(y,-1,11)%11 if y else 11
    return tuple(val(x) for x in range(12))


def witt10():
    gens=[mobius(1,1,0,1),mobius(3,0,0,1),mobius(0,10,1,0)]
    g=list(range(12))
    for a,b in ((2,10),(3,4),(5,9),(6,7)):g[a],g[b]=b,a
    gens.append(g)
    def orbit(base):
        seen={frozenset(base)}; pending=list(seen)
        while pending:
            s=pending.pop()
            for g in gens:
                t=frozenset(g[x] for x in s)
                if t not in seen:seen.add(t);pending.append(t)
        return seen
    for base in combinations(range(12),6):
        h=orbit(base)
        if len(h)==132:break
    check(len(h)==132,'inherited Witt orbit')
    blocks=sorted({tuple(sorted(s-{10,11})) for s in h if {10,11}<=s})
    check(len(blocks)==30,'W10 size')
    for a in combinations(range(10),3):
        check(sum(set(a)<=set(b) for b in blocks)==1,'W10 triple unique')
    return blocks


def nae(a,b,c):return not(a==b==c)


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True)
    args=ap.parse_args()
    P=[x for x in range(1,64) if q(x)==0]
    L=sorted({tuple(sorted((x,y,x^y))) for x,y in combinations(P,2) if polar(x,y)==0})
    check((len(P),len(L))==(27,45),'elliptic quadric sizes')
    neigh={p:{s for s in P if s!=p and polar(p,s)==0} for p in P}
    local={}
    for p in P:
        pairs=sorted({tuple(sorted((s,s^p))) for s in neigh[p]})
        check(len(neigh[p])==10 and len(pairs)==5,'marked neighbourhood is five pairs')
        check(all((polar(a,b)==0)==(a^b==p) for a,b in combinations(neigh[p],2)),
              'induced neighbour graph exactly matching')
        local[p]=pairs
    for p in P:
        for line in L:
            if p not in line:check(sum(x in neigh[p] for x in line)==1,'GQ axiom')
    # A nonsingular polar hyperplane retains 15 original points and 15 lines.
    sections={}
    for h in range(1,64):
        if q(h)!=1:continue
        S={p for p in P if polar(h,p)==0}; T=[l for l in L if set(l)<=S]
        check((len(S),len(T))==(15,15),'doily section')
        sections[h]=sorted(S)
    check(len({tuple(v) for v in sections.values()})==36,'36 distinct sections')
    h=min(sections); S=sections[h]; T=[l for l in L if set(l)<=set(S)]
    duads=list(combinations(range(6),2))
    def matchings(xs):
        if not xs:yield ();return
        for b in xs[1:]:
            rest=tuple(x for x in xs[1:] if x!=b)
            for m in matchings(rest):yield ((xs[0],b),)+m
    ML=list(matchings(tuple(range(6))))
    doily_map=incidence_map(S,T,duads,ML)
    # Independent classical 27-class model in the Picard lattice of the
    # plane blown up at six points: H^2=1, E_i^2=-1, other products zero.
    classes=[];names=[]
    for i in range(6):
        classes.append((0,)+tuple(int(j==i) for j in range(6)));names.append(f'E{i}')
    for i,j in combinations(range(6),2):
        classes.append((1,)+tuple(-int(k in (i,j)) for k in range(6)));names.append(f'L{i}{j}')
    for i in range(6):
        classes.append((2,)+tuple(-int(j!=i) for j in range(6)));names.append(f'Q{i}')
    ip=lambda a,b:a[0]*b[0]-sum(x*y for x,y in zip(a[1:],b[1:]))
    meets=lambda i,j:ip(classes[i],classes[j])==1
    triples=[t for t in combinations(range(27),3) if all(meets(i,j) for i,j in combinations(t,2))]
    cubic_map=incidence_map(P,L,list(range(27)),triples)
    check(len(triples)==45,'Picard model tritangent triples')
    for i,j in combinations(range(27),2):
        dif=tuple(a-b for a,b in zip(classes[i],classes[j]))
        check((ip(dif,dif)==-2)==(not meets(i,j)), 'root difference means skew lines, not meeting')
    # Nine points of the Hermitian curve inside the 21-point plane over F4.
    PG=sorted({normalize(v) for v in product(range(4),repeat=3) if any(v)})
    lines=[tuple(i for i,p in enumerate(PG) if dot(a,p)==0) for a in PG]
    U=[i for i,p in enumerate(PG) if sum(x!=0 for x in p)%2==0]
    intersections=[tuple(i for i in l if i in U) for l in lines]
    UL=sorted(l for l in intersections if len(l)==3)
    check((len(PG),len(U),len(UL))==(21,9,12),'Hermitian unital sizes')
    check(Counter(map(len,intersections))=={1:9,3:12},'tangent/secant counts')
    check(all(len(l)==5 for l in lines),'projective line size five')
    A=list(product(range(3),repeat=2))
    AL=sorted({tuple(sorted(((x+t*a)%3,(y+t*b)%3) for t in range(3)))
               for x,y in A for a,b in A if (a,b)!=(0,0)})
    u_to_affine=incidence_map(U,UL,A,AL)
    W=witt10(); WL=[tuple(x for x in b if x!=9) for b in W if 9 in b]
    w_to_affine=incidence_map(list(range(9)),WL,A,AL)
    # Mark a point in the derived affine plane. Its eight remaining points
    # come in four opposite pairs, exactly the matching derived from W10.
    origin=next(p for p,a in w_to_affine.items() if a==(0,0))
    affine_to_w={a:p for p,a in w_to_affine.items()}
    reps=[(1,0),(0,1),(1,1),(1,2)]
    ordered=[(sign*x%3,sign*y%3) for sign in (1,2) for x,y in reps]
    opposite_pairs={frozenset((affine_to_w[r],affine_to_w[(2*r[0]%3,2*r[1]%3)])) for r in reps}
    derived_pairs={frozenset(set(b)-{9,origin}) for b in W if {9,origin}<=set(b)}
    check(opposite_pairs==derived_pairs,'same marked W10 opposite pairs')
    pair_masks=[(1<<i)|(1<<(i+4)) for i in range(4)]
    even_pairs=span([pair_masks[0]^p for p in pair_masks[1:]])
    extensions=[even_pairs|{w^t for w in even_pairs} for t in (15,30)]
    check(all(len(c)==16 and Counter(x.bit_count() for x in c)=={0:1,4:14,8:1} for c in extensions),
          'oriented affine residue gives exactly two H8 extensions')
    GL=[]; atlas_orbit=set(); atlas_stabilizer=0; code_stabilizer=0
    free=span([255,240]); original_charts={frozenset(r^v for v in free) for r in (0,170,204,102)}
    # Code E+<15> equals the standard H8 in the retained order.
    check(extensions[0]==span([170,204,240,255]),'same literal H8 port order')
    for a,b,c,d in product(range(3),repeat=4):
        det=(a*d-b*c)%3
        if det==0:continue
        perm=[ordered.index(((a*x+b*y)%3,(c*x+d*y)%3)) for x,y in ordered]
        relabel=lambda w:sum(1<<perm[i] for i in range(8) if w>>i&1)
        image=set(map(relabel,extensions[0])); parity=sum(perm[i]>=4 for i in range(4))%2
        check(image==extensions[0 if det==1 else 1],'determinant controls oriented H8 extension')
        check(parity==(det==2),'orientation parity equals quadratic determinant character')
        if det==1:
            code_stabilizer+=1
            new=frozenset(frozenset(map(relabel,ch)) for ch in original_charts)
            atlas_orbit.add(new);atlas_stabilizer+=new==original_charts
        GL.append(dict(matrix=[a,b,c,d],determinant=det,port_permutation=perm))
    check(len(GL)==48 and code_stabilizer==24,'GL and SL local actions')
    ordered_atlases=sorted(atlas_orbit,key=lambda a:sorted(sorted(c) for c in a))
    atlas_actions=set()
    for g in GL:
        if g['determinant']!=1:continue
        perm=g['port_permutation']
        relabel=lambda w:sum(1<<perm[i] for i in range(8) if w>>i&1)
        action=tuple(ordered_atlases.index(frozenset(frozenset(map(relabel,c)) for c in atlas))
                     for atlas in ordered_atlases)
        check(sum(action[i]>action[j] for i,j in combinations(range(4),2))%2==0,
              'four-atlas action is even')
        atlas_actions.add(action)
    check(len(ordered_atlases)==4 and len(atlas_actions)==12,'tetrahedral A4 action on chart atlases')
    # The finite geometry does not carry planar unit lengths. Authenticate
    # a cheap hostile: any two unital points share a secant => full graph K9.
    check(all(sum(a in l and b in l for l in UL)==1 for a,b in combinations(U,2)),
          'unital collinearity graph K9, cannot be planar unit-distance')
    # Four Hamming charts as common discrete derivatives on two squares.
    C=span([170,204,240,255,768]); charts={}
    for sel in product(range(2),repeat=2):
        words=[]
        for starts in product(range(2),repeat=3):
            p,q0,k=starts; alpha,beta=sel
            bits=[p,p^alpha,p^beta,p^alpha^beta,
                  q0,q0^alpha,q0^beta,q0^alpha^beta,k,k]
            words.append(sum(b<<i for i,b in enumerate(bits)))
        charts[str(sel)]=sorted(words)
    check(set().union(*(set(x) for x in charts.values()))==C,'four-selector relation exact')
    check(sum(map(len,charts.values()))==32,'charts disjoint')
    hostile=(0,170,204); majority=(hostile[0]&hostile[1])|(hostile[0]&hostile[2])|(hostile[1]&hostile[2])
    check(all(x in C for x in hostile) and majority==136 and majority not in C,
          'majority hostile forbids arbitrary two-list selector')
    # Standalone triangle plus three pendant inputs, with internal lists
    # {0,1,2}, exactly implements NAE. Without lists, 6-colour support is full.
    naetable=[]
    for ports in product(range(2),repeat=3):
        accepted=sum(all(a!=b for a,b in zip(inner,ports)) for inner in permutations(range(3)))
        six=sum(all(a!=b for a,b in zip(inner,ports)) for inner in permutations(range(6),3))
        check(bool(accepted)==nae(*ports),'triangle implements listed NAE')
        check(six>0,'unlisted six-colour triangle does not implement NAE')
        naetable.append(dict(ports=ports,three_list_extensions=accepted,six_extensions=six))
    # NAE4 split provides a composition-ready OR3 gadget, and four OR3
    # clauses implement a three-bit even-parity gate. Check hidden choices.
    for bits in product(range(2),repeat=4):
        possible=any(nae(bits[0],bits[1],w) and nae(1-w,bits[2],bits[3]) for w in range(2))
        check(possible==(len(set(bits))>1),'NAE4 composition with complemented port')
    for bits in product(range(2),repeat=3):
        possible=True
        for forbidden in product(range(2),repeat=3):
            if sum(forbidden)%2!=1:continue
            literals=[a^b for a,b in zip(bits,forbidden)]
            possible &= any(nae(literals[0],literals[1],w) and nae(1-w,literals[2],0) for w in range(2))
        check(possible==(sum(bits)%2==0),'composed even-parity list gate')
    # Exact planar drawing: coordinates (a+b sqrt(3), c+d sqrt(3))/2.
    # Interior equilateral triangle 0,(1,0),(1/2,sqrt3/2); outer pendant
    # ports (-1,0),(2,0),(1/2,1+sqrt3/2).
    pts=[(0,0,0,0),(2,0,0,0),(1,0,0,1),(-2,0,0,0),(4,0,0,0),(1,0,2,1)]
    def dist2(v,w):
        a,b,c,d=[x-y for x,y in zip(v,w)]
        return (a*a+3*b*b+c*c+3*d*d,2*a*b+2*c*d)
    E=[(i,j) for i,j in combinations(range(6),2) if dist2(pts[i],pts[j])==(4,0)]
    check(E==[(0,1),(0,2),(0,3),(1,2),(1,4),(2,5)],'all unit distances authenticated')
    out=dict(status='FINITE-EXACT; accompanying proofs scoped in note',checks=CHECKS,
        quadric=dict(points=P,lines=L,local_matchings=local,doily_section_h=h,
                     section_to_duads=[[p,doily_map[p]] for p in S],section_count=len(sections),
                     to_cubic_line=[[p,names[cubic_map[p]]] for p in P],picard_classes=dict(zip(names,classes))),
        unital=dict(PG_points=PG,PG_lines=lines,points=U,secants=UL,
                    to_affine=[[p,u_to_affine[p]] for p in U]),
        witt=dict(blocks=W,derived_at=9,to_affine=[[p,w_to_affine[p]] for p in range(9)],
                  second_mark=origin,ordered_nonzero_vectors=ordered,
                  hamming_extensions=[sorted(c) for c in extensions],GL2_3=GL,
                  SL2_chart_atlas_orbit=len(atlas_orbit),SL2_chart_atlas_stabilizer=atlas_stabilizer,
                  four_atlases=[[sorted(c) for c in sorted(a,key=lambda s:sorted(s))] for a in ordered_atlases],
                  atlas_actions=sorted(atlas_actions)),
        selectors=dict(charts=charts,majority_inputs=hostile,majority_output=majority,
                       nae_table=naetable,geometric_coordinates=pts,unit_edges=E,
                       unpaid='three-colour list forcing and compatible gluing at fixed physical ports'))
    args.out.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(dict(checks=CHECKS,quadric=[len(P),len(L)],unital=[len(PG),len(U),len(UL)],
                         selectors=4,unit_edges=len(E)),indent=2))


if __name__=='__main__':main()

"""Exact labelled excursion forests, affine sorting, and carry controls.

No floating comparison decides a result. Run with Python 3.10+.
"""
from collections import Counter
from fractions import Fraction as Q
from itertools import product
import json


def require(ok, label):
    if not ok:
        raise RuntimeError(label)


def shortcut(n):
    return (3*n+1)//2 if n & 1 else n//2


def compose(later, earlier):
    a, b = later
    c, d = earlier
    return a*c, a*d+b


def block(word):
    A, C, den = 1, 0, 1
    for bit in word:
        A, C, den = (3*A, 3*C+den, 2*den) if bit else (A, C, 2*den)
    return A, C, den


def affine(word):
    a, c, d = block(word)
    return Q(a,d), Q(c,d)


def residue(word):
    a,c,d=block(word)
    return (-c*pow(a,-1,d)) % d


def curvature(pair):
    a,b=pair
    require(b>0,'curvature needs an odd letter')
    return (a-1)/b


def carry_stream(n):
    # r_s=n mod 2^s; carries in the digit multiplication 3n+1.
    out=[]
    for s in range(1,n.bit_length()+3):
        den=1<<s
        out.append((3*(n % den)+1)//den)
    while out and out[-1]==0:
        out.pop()
    return out


def tiles_for_step(n):
    p=n&1
    c=p
    tiles=[]
    digits=[]
    for s in range(n.bit_length()+3):
        b=(n>>s)&1
        value=(1+2*p)*b+c
        q,c_next=value%2,value//2
        tile=(p,b,c,q,c_next)
        require(c in (0,1,2) and c_next in (0,1,2),'tile carries')
        if p==0:
            require(c==c_next==0,'even row has no carries')
        tiles.append(tile)
        digits.append(q)
        c=c_next
    require(digits[0]==0 and c==0,'left parity and finite right boundary')
    result=sum(bit<<s for s,bit in enumerate(digits[1:]))
    require(result==shortcut(n),'eight-tile local computation')
    return tiles


def digit_poly(n,z):
    return sum(((n>>s)&1)*z**s for s in range(n.bit_length()))


def forest_polynomial(states,z):
    # Independent block identity retaining positions of every carry.
    B=Q(0)
    A=Q(0)
    power3=1
    for t,n in enumerate(states[:-1]):
        if n&1:
            C_over_z=sum(c*z**s for s,c in enumerate(carry_stream(n)))
            B=3*B+z**t
            A=3*A+z**t*C_over_z
            power3*=3
    require(z**(len(states)-1)*digit_poly(states[-1],z)==
        power3*digit_poly(states[0],z)+B-(2-z)*A,'weighted forest polynomial')
    return B,A


def path(n, cap):
    states=[n]
    for _ in range(cap):
        if n==1:
            break
        n=shortcut(n)
        states.append(n)
    return states


def forest(states):
    """Nodes are matched up/down intervals; children are ordered in time.

    Flat edges and unmatched boundary edges remain explicit in `atoms`.
    The root forest includes open nodes at truncated endpoints.
    """
    nodes=[]
    stack=[]
    boundary=[]
    for i,(x,y) in enumerate(zip(states,states[1:])):
        step=y.bit_length()-x.bit_length()
        require(step in (-1,0,1),'height step')
        if step==1:
            node={'start':i,'end':None,'parent':stack[-1] if stack else None,'children':[]}
            nodes.append(node)
            node_id=len(nodes)-1
            if stack:
                nodes[stack[-1]]['children'].append(node_id)
            stack.append(node_id)
        elif step==-1:
            if stack:
                nodes[stack.pop()]['end']=i+1
            else:
                boundary.append(i)
    bits=[n&1 for n in states[:-1]]
    for node in reversed(nodes):
        if node['end'] is None:
            continue
        i,j=node['start'],node['end']
        word=bits[i:j]
        A,C,D=block(word)
        x,y=states[i],states[j]
        require(D*y==A*x+C,'unwrapped node identity')
        require(x.bit_length()==y.bit_length(),'matched node same band')
        require(all(z.bit_length()>x.bit_length() for z in states[i+1:j]),'primitive height excursion')
        require(C>0,'matched node has upward odd step')
        k=Q(A-D,C)
        require(Q(y-x)==Q(C,D)*(1+k*x),'curvature displacement')
        require((y<x)==(k<Q(-1,x)),'source-dependent descent threshold')
        # Independently rebuild with immediate subtrees as macro atoms.
        macro=(Q(1),Q(0))
        cursor=i
        atoms=[]
        for child_id in node['children']:
            child=nodes[child_id]
            require(child['end'] is not None,'a closed node has no open child')
            for t in range(cursor,child['start']):
                macro=compose(affine([bits[t]]),macro)
                atoms.append(['edge',t])
            child_map=tuple(Q(v) for v in child['map'])
            macro=compose(child_map,macro)
            atoms.append(['node',child_id])
            cursor=child['end']
        for t in range(cursor,j):
            macro=compose(affine([bits[t]]),macro)
            atoms.append(['edge',t])
        require(macro==(Q(A,D),Q(C,D)),'tree composition retains flats')
        node.update(source=x,target=y,word=''.join(map(str,word)),length=j-i,
                    odd=sum(word),carry=C,denominator=D,map=[str(z) for z in macro],
                    curvature=str(k),atoms=atoms,
                    direction='down' if y<x else 'up' if y>x else 'equal')
    return nodes,boundary,stack


def word_census():
    rows=[]
    maps=[]
    for length in range(1,11):
        residues=set()
        for word in product((0,1),repeat=length):
            n=residue(word)
            residues.add(n)
            x=n+(1<<length)
            y=x
            for bit in word:
                require((y&1)==bit,'source cylinder realizes word')
                y=shortcut(y)
            a,c,d=block(word)
            require(d*y==a*x+c,'word identity independent replay')
            if length<=5 and sum(word):
                maps.append((word,affine(word)))
        require(len(residues)==1<<length,'source residues biject parity words')
        rows.append([length,len(residues)])
    count=ties=0
    for w,f in maps:
        for v,g in maps:
            A,B=f
            C,D=g
            fg,gf=compose(f,g),compose(g,f)
            kf,kg=curvature(f),curvature(g)
            require(fg[0]==gf[0],'same slope under swap')
            require(fg[1]-gf[1]==B*D*(kf-kg),'intrinsic comparator')
            require(curvature(fg)==(B*kf+A*D*kg)/(B+A*D),'curvature weighted mean')
            require(min(kf,kg)<=curvature(fg)<=max(kf,kg),'mean bounds')
            if w+v!=v+w:
                require(residue(w+v)!=residue(v+w),'nontrivial swap changes source')
            if kf==kg:
                ties+=1
            count+=1
    return rows,count,ties


def main():
    rows,comparisons,ties=word_census()
    seeds=list(range(2,502))+[703,871,6171,77031,837799,63728127]
    seeds += [(1<<r)-1 for r in (16,32,64,128)]
    counts=Counter()
    example=None
    forest27=None
    tile_types=set()
    for n in seeds:
        states=path(n,4000)
        nodes,boundary,open_nodes=forest(states)
        counts['seeds']+=1
        counts['steps']+=len(states)-1
        counts['open_nodes']+=len(open_nodes)
        counts['unmatched_down']+=len(boundary)
        for node in nodes:
            if node['end'] is None:
                continue
            counts['closed_nodes']+=1
            counts['return_'+node['direction']]+=1
            x,y=node['source'],node['target']
            if y>x and example is None:
                example={k:node[k] for k in ('source','target','word','length','odd','carry','denominator','curvature')}
        for x in states:
            tile_types.update(tiles_for_step(x))
            cs=carry_stream(x)
            K=sum(cs)
            require(K==1+3*x.bit_count()-(3*x+1).bit_count(),'binary carry identity')
            require(all(c in (0,1,2) for c in cs),'ternary carry alphabet')
            if x&1:
                require(K>=3,'odd positive carry minimum')
        if n==27:
            forest27={k:sum(node.get('direction')==k for node in nodes) for k in ('up','down','equal')}
            # One full tree sample retains exact integer and chronology.
            forest27['first_node']={k:v for k,v in nodes[0].items() if k!='atoms'}
    # Prefixes keep the live boundary; no fictitious completed return.
    for horizon in (1,2,5,10,20,40):
        forest(path(27,horizon))
    polynomial_controls=0
    for source in (1,3,27,51,703,6171,(1<<32)-1):
        for horizon in (1,3,8,20,40):
            states=path(source,horizon)
            for z in (Q(1),Q(3,2),Q(2)):
                forest_polynomial(states,z)
                polynomial_controls+=1
    w=(1,0)
    v=(1,1,0)
    require(residue(w+v)==25 and residue(v+w)==11,'minimal formal block swap control')
    require(len(tile_types)==8,'exact number of tile types')
    for L in range(9,70):
        n=15*(1<<L)+51
        y=45*(1<<(L-1))+77
        require(shortcut(n)==y,'padded polynomial hostile')
        # Both sides have degree <= L+4. Distinct exact integer evaluations
        # exceeding that degree also certify the polynomial identity.
        for z in range(L+6):
            expected=z**(L-1)*(z-1)*(z**4-1)+z*(z-1)*(z**4-z*z+1)
            require(digit_poly(y,z)-digit_poly(n,z)==expected,'weighted depth hostile identity')
            if z>1:
                require(expected>0,'weighted depth grows')
    print(json.dumps({'word_residue_census':rows,'affine_pair_checks':comparisons,
        'comparator_ties':ties,'forest_counts':dict(counts),'growing_return':example,
        'forest27':forest27,'swap':{'words':['10110','11010'],'sources_mod32':[25,11]},
        'local_tiles_p_b_c_q_cnext':sorted(tile_types),
        'digit_polynomial_hostile':'61 padded pairs, exact polynomial identity; positive for every real z>1',
        'weighted_forest_polynomial_controls':polynomial_controls,
        'scope':'exact finite forests and algebra; no universal descent or stopping bound'},indent=2))


if __name__=='__main__':
    main()

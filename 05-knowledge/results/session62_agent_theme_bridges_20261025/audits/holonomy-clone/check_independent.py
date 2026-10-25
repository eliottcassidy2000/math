#!/usr/bin/env python3
"""Independent exact audit; no imports from the candidate implementations."""
from collections import Counter
from itertools import combinations, permutations, product
from math import prod, lcm
import json

CHECKS = 0

def need(condition, message):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(message)

BITS = tuple(product(range(2), repeat=3))
PERMS = tuple(permutations(range(3)))
EVEN = tuple(x for x in BITS if sum(x) % 2 == 0)

def transport(x, p):
    # Pullback convention: internally independent of candidate's pushforward.
    return tuple(x[p[i]] for i in range(3))

def xor(x, y):
    return tuple(a ^ b for a, b in zip(x, y))

def holonomy_audit():
    census = Counter()
    gain_cases = 0
    for beta in range(4):
        for linear in product(PERMS, repeat=beta):
            fixed_even = [x for x in EVEN if all(transport(x, p) == x for p in linear)]
            f = {1:0, 2:1, 4:2}[len(fixed_even)]
            image_fibers = Counter(tuple(xor(x, transport(x, p)) for p in linear) for x in BITS)
            lost = 2 * beta - 2 + f
            need(lost >= 0, 'realizable holonomy has nonnegative obstruction')
            need(len(image_fibers) == 2**(2-f), 'full compatible gain count')
            need(all(m == 2**(1+f) for m in image_fibers.values()), 'all affine solution fibers')
            parity_compatible = 0
            for gains in product(BITS, repeat=beta):
                gain_cases += 1
                scalar_ok = all(sum(w) % 2 == 0 for w in gains)
                if scalar_ok:
                    parity_compatible += 1
                roots = [x for x in BITS if all(xor(x, transport(x,p)) == w for p,w in zip(linear,gains))]
                need(len(roots) == image_fibers.get(gains, 0), 'literal root equations versus fiber census')
                if roots:
                    need(scalar_ok, 'endpoint section projects to completion section')
                    need(Counter(sum(x)%2 for x in roots) == {0:2**f,1:2**f}, 'each completion has exactly 2^f lifts')
            need(parity_compatible == len(image_fibers) * 2**lost, 'exact forgotten obstruction index')
            census[(beta,f,lost)] += 1
    rot = (1,2,0)
    need([x for x in BITS if transport(x,rot)==x] == [(0,0,0),(1,1,1)], 'nontrivial holonomy with two sections')
    need(any(transport(x,rot)!=x for x in BITS), 'same transport not trivialized')
    return {'gain_cases':gain_cases, 'linear_cases':sum(census.values()),
            'census':[{'beta':b,'fixed_even_dimension':f,'lost_dimension':d,'cases':v}
                      for (b,f,d),v in sorted(census.items())]}

def fano_audit():
    forced = tuple(frozenset((2*i,2*i+1,6)) for i in range(3))
    transversal = [frozenset(2*i+x[i] for i in range(3)) for x in BITS]
    valid = []
    all_pairs = set(combinations(range(7),2))
    for indices in combinations(range(8),4):
        lines = forced + tuple(transversal[i] for i in indices)
        multiplicity = Counter(tuple(sorted(p)) for line in lines for p in combinations(line,2))
        if set(multiplicity) == all_pairs and set(multiplicity.values()) == {1}:
            valid.append(frozenset(lines))
    need(len(valid)==2, 'exactly two Fano completions from literal pair incidence')
    halves = [frozenset(forced + tuple(transversal[i] for i,x in enumerate(BITS) if sum(x)%2==p)) for p in range(2)]
    need(set(valid)==set(halves), 'the two completions are cube parity halves')
    for p in PERMS:
        for w in BITS:
            # Slot map implements the same pullback coordinate action.
            pointmap = {6:6}
            for i in range(3):
                for b in range(2):
                    pointmap[2*p[i]+b] = 2*i+(b^w[i])
            for eta in range(2):
                mapped = frozenset(frozenset(pointmap[x] for x in line) for line in halves[eta])
                need(mapped == halves[eta ^ (sum(w)%2)], 'actual seven-point Fano incidence action')
    return {'candidate_four_line_sets':70,'fano_completions':len(valid),'affine_incidence_maps':48}

def cycles(p):
    unused = set(range(len(p)))
    ans = []
    while unused:
        x = next(iter(unused)); cycle = []
        while x in unused:
            unused.remove(x); cycle.append(x); x = p[x]
        ans.append(tuple(cycle))
    return ans

def mask_action(mask,p):
    out = 0
    while mask:
        bit = mask & -mask
        out |= 1 << p[bit.bit_length()-1]
        mask -= bit
    return out

def value_minus_one(p):
    return sum(c * (-1)**i for i,c in enumerate(p))

def divide_plus_one(p):
    q = [0]*(len(p)-1)
    q[-1] = p[-1]
    for i in range(len(p)-2,0,-1):
        q[i-1] = p[i]-q[i]
    need(p[0] == q[0], 'exact division by t+1')
    return q

def multiplicity(p):
    d=0
    while len(p)>1 and value_minus_one(p)==0:
        p=divide_plus_one(p); d+=1
    return d,value_minus_one(p),p

def group_orbit_audit(group, odd=True):
    n=len(group[0]); unseen=set(range(1<<n)); orbit_of={}; reps=[]
    while unseen:
        m=min(unseen); orb={mask_action(m,p) for p in group}
        need(orb <= unseen, 'subset orbit partition')
        j=len(reps); reps.append(m)
        for x in orb: orbit_of[x]=j
        unseen -= orb
    P=[0]*(n+1)
    for m in reps:P[m.bit_count()]+=1
    d,lead,Q=multiplicity(P)
    types=[cycles(p) for p in group]
    mincycles=min(map(len,types))
    point_orbits=[]; unseen=set(range(n))
    while unseen:
        i=min(unseen); orb={p[i] for p in group}
        point_orbits.append(orb); unseen-=orb
    r=len(point_orbits)
    for F in range(1<<n):
        invariant = all(mask_action(F,p)==F for p in group)
        singleton_witness = all(orbit_of[F ^ mask_action(F,p)]==orbit_of[0] for p in group)
        need(invariant == singleton_witness, 'fixed-XOR descent necessity at the orbit of F')
    if odd:
        need(all(len(c)%2 for cs in types for c in cs), 'all cycles odd')
        need(d==mincycles, 'subset root order equals minimum cycles')
        need(lead*len(group)==sum(prod(map(len,cs)) for cs in types if len(cs)==d), 'positive leading constant')
        need(d>=r and (d-r)%2==0, 'common invariant lower bound and parity')
        flags=[sum(1<<x for x in o) for o in point_orbits]
        toggles=[sum(flags[i] for i in range(r) if m>>i&1) for m in range(1<<r)]
        cube_unseen=set(range(len(reps))); cube_leads=[]
        while cube_unseen:
            j=min(cube_unseen); orb={orbit_of[reps[j]^f] for f in toggles}
            need(len(orb)==1<<r, 'free common-invariant toggle cube on subset orbits')
            need(orb<=cube_unseen, 'toggle cubes partition subset orbits')
            cube_unseen-=orb
            cp=[0]*(n+1)
            for q in orb:cp[reps[q].bit_count()]+=1
            cd,cl,_=multiplicity(cp)
            need(cd==r, 'individual Boolean cube vanishing order exactly r')
            cube_leads.append(cl)
        need((sum(cube_leads)==0)==(d>r), 'extra cancellation is between common-invariant cubes')
        need(len(reps)%(2**d)==0, 'integral extra divisibility')
    return {'points':n,'group_order':len(group),'point_orbits':r,'root_order':d,
            'minimum_cycles':mincycles,'leading_coefficient':lead,'polynomial':P,'residual_polynomial':Q}

def cyclic_parts(parts):
    shifts=[]
    for a in range(lcm(*parts)):
        p=[];offset=0
        for m in parts:
            p.extend(offset+(j+a)%m for j in range(m));offset+=m
        shifts.append(tuple(p))
    return shifts

def projections(forms):
    return [tuple(3*j+(z+f[0]*a+f[1]*b)%3 for j,f in enumerate(forms) for z in range(3))
            for a,b in product(range(3),repeat=2)]

def regular(moduli):
    points=list(product(*(range(q) for q in moduli))); ix={v:i for i,v in enumerate(points)}
    return [tuple(ix[tuple((x[j]+a[j])%moduli[j] for j in range(len(moduli)))] for x in points) for a in points]

def group_audit():
    cases={
      'empty':[()],
      'trivial3':[tuple(range(3))],
      'C3_regular':regular((3,)),
      'C9_regular':regular((9,)),
      'C3xC3_regular':regular((3,3)),
      'C3xC3_four_projections':projections(((1,0),(0,1),(1,1),(1,2))),
      'Frobenius21_on7':[tuple((pow(2,j,7)*x+b)%7 for x in range(7)) for j in range(3) for b in range(7)],
      'C15_three_plus_five':cyclic_parts((3,5)),
      'C3_diagonal_two_orbits':cyclic_parts((3,3)),
      'C3_independent_two_orbits':projections(((1,0),(0,1))),
    }
    out={name:group_orbit_audit(g) for name,g in cases.items()}
    out['C2_even_hostile']=group_orbit_audit([(0,1),(1,0)],odd=False)
    need(out['C2_even_hostile']['root_order']==0 and out['C2_even_hostile']['minimum_cycles']==1,'even group requires qualification')
    need(out['C3xC3_regular']['root_order']==3 and out['C3xC3_regular']['point_orbits']==1,'strict regular gap example')
    need(out['C3xC3_four_projections']['root_order']==6 and out['C3xC3_four_projections']['point_orbits']==4,'synchronization hostile')
    return out

def build_base(k,forms):
    points=list(product(range(3),repeat=k)); labels=[('x',a) for a in points]+[('y',i,c) for i in range(len(forms)) for c in range(3)]
    def beats(u,v):
        if u[0]=='y' and v[0]=='x':return not beats(v,u)
        if u[0]=='x' and v[0]=='y':return sum(a*b for a,b in zip(u[1],forms[v[1]]))%3==v[2]
        if u[0]=='y':return u[1]<v[1] if u[1]!=v[1] else (v[2]-u[2])%3==1
        difference=next(b-a for a,b in zip(u[1],v[1]) if a!=b)
        return difference%3==1
    matrix=[[False if u==v else beats(u,v) for v in labels] for u in labels]
    return labels,matrix

def homogeneous(A):
    return [(i,j) for i in range(len(A)) for j in range(i+1,len(A)) if all(A[i][h]==A[j][h] for h in range(len(A)) if h not in (i,j))]

def tournament_audit():
    out={}
    for k in range(1,5):
        forms=[tuple(int(i==j) for j in range(k)) for i in range(k)]
        labels,A=build_base(k,forms)
        pairs=homogeneous(A)
        if k==1:
            need(len(pairs)==3,'k1 boundary pair witness');out['k1_pairs']=pairs;continue
        need(not pairs,'literal base pair-free')
        scores=[sum(row) for row in A]
        need(all(scores[i]==(3**k-1)//2+k for i,u in enumerate(labels) if u[0]=='x'),'X score formula')
        need(all(scores[i]==1+2*3**(k-1)+3*(k-1-u[1]) for i,u in enumerate(labels) if u[0]=='y'),'observer score formula')
        need(len(set(scores))==k+1,'intrinsic score classes separate all observer parts and X')
        # Exhaust all observer rotations. Incidence forces at most these maps;
        # every map is checked on every arc, without a graph-search algorithm.
        index={u:i for i,u in enumerate(labels)}; candidates=[]
        for shift in product(range(3),repeat=k):
            p=[index[('x',tuple((a+b)%3 for a,b in zip(u[1],shift)))] if u[0]=='x'
               else index[('y',u[1],(u[2]+shift[u[1]])%3)] for u in labels]
            need(all(A[i][j]==A[p[i]][p[j]] for i in range(len(A)) for j in range(len(A))),'forced observer map preserves every arc')
            candidates.append(p)
        expanded=[(i,a) for i,u in enumerate(labels) for a in range(2 if u[0]=='x' else 1)]
        B=[[False if z==w else (a<b if i==j else A[i][j]) for w,(j,b) in enumerate(expanded)] for z,(i,a) in enumerate(expanded)]
        actual=homogeneous(B)
        expected=[(i,i+1) for i,(j,a) in enumerate(expanded) if labels[j][0]=='x' and a==0]
        need(actual==expected,'literal clone pairs are precisely intended blocks')
        out['k'+str(k)]={'base_vertices':len(labels),'inflated_vertices':len(B),'clone_pairs':len(actual),'automorphisms_proved_by_observers':len(candidates),'gap_minimum_cycles':3**(k-1)}
    forms=((1,0),(0,1),(1,1),(1,2)); labels,A=build_base(2,forms)
    need(not homogeneous(A),'four-form base pair-free')
    need(set(map(sum,A))=={8,16,13,10,7},'four-form intrinsic degree classes')
    ix={u:i for i,u in enumerate(labels)}; admitted=[]
    for shifts in product(range(3),repeat=4):
        a,b=shifts[:2]
        if shifts[2]!=(a+b)%3 or shifts[3]!=(a+2*b)%3:continue
        p=[ix[('x',((u[1][0]+a)%3,(u[1][1]+b)%3))] if u[0]=='x' else ix[('y',u[1],(u[2]+shifts[u[1]])%3)] for u in labels]
        need(all(A[i][j]==A[p[i]][p[j]] for i in range(len(A)) for j in range(len(A))),'four-form forced automorphisms')
        admitted.append(p)
    need(len(admitted)==9,'four-form full observer-compatible group')
    expanded=[(i,a) for i,u in enumerate(labels) for a in range(2 if u[0]=='y' else 1)]
    B=[[False if z==w else (a<b if i==j else A[i][j]) for w,(j,b) in enumerate(expanded)] for z,(i,a) in enumerate(expanded)]
    need(len(homogeneous(B))==12,'T33 literal clone pairs')
    out['T33']={'vertices':len(B),'clone_pairs':12,'automorphisms_proved_by_observers':9}
    return out

def common_action_audit():
    group=regular((3,3)); points=list(product(range(3),repeat=2)); index={x:i for i,x in enumerate(points)}
    generators=[group[index[(1,0)]],group[index[(0,1)]]]
    image={}; even_image=set()
    for mask in range(1<<9):
        value=tuple(mask ^ mask_action(mask,p) for p in generators)
        image.setdefault(value,[]).append(mask)
        if mask.bit_count()%2==0:even_image.add(value)
    need(len(image)==256 and set(map(len,image.values()))=={2},'nine-gap full coboundary image by enumeration')
    need(len(even_image)==256,'nine-gap even coboundary injectivity by enumeration')
    hv=18-(len(image).bit_length()-1)
    he=16-(len(even_image).bit_length()-1)
    need((hv,he)==(10,8),'same action has full H1 ten and forgotten H1 eight')
    return {'gaps':9,'vertex_assignments':512,'common_fixed_dimension':1,
            'minimum_individual_fixed_dimension':min(len(cycles(p)) for p in group),
            'full_coboundary_image':len(image),'even_coboundary_image':len(even_image),'full_H1':hv,'forgotten_H1':he}

def exterior_audit():
    # Independent exact determinant and rank over Q; standard-library only.
    from fractions import Fraction
    def determinant_poly(M):
        n=len(M); ans=[0]*(n+1)
        for p in permutations(range(n)):
            sign=(-1)**sum(p[i]>p[j] for i in range(n) for j in range(i+1,n))
            poly=[sign]
            for i,j in enumerate(p):
                q=[0]*(len(poly)+1)
                for d,c in enumerate(poly):q[d]+=(i==j)*c;q[d+1]+=M[i][j]*c
                poly=q
            ans=[a+b for a,b in zip(ans,poly)]
        return ans
    def fixed_dimension(M):
        n=len(M);A=[[Fraction(M[i][j]-(i==j)) for j in range(n)] for i in range(n)];rank=0
        for j in range(n):
            pivot=next((i for i in range(rank,n) if A[i][j]),None)
            if pivot is None:continue
            A[rank],A[pivot]=A[pivot],A[rank];c=A[rank][j];A[rank]=[x/c for x in A[rank]]
            for i in range(n):
                if i!=rank:
                    c=A[i][j];A[i]=[a-c*b for a,b in zip(A[i],A[rank])]
            rank+=1
        return n-rank
    def mm(A,B):return [[sum(a*b for a,b in zip(row,col)) for col in zip(*B)] for row in A]
    I=[[1,0],[0,1]];R=[[0,-1],[1,-1]];S=[[0,1],[1,0]]
    C3=[I,R,mm(R,R)];D3=C3+[mm(S,r) for r in C3]
    groups={'C2_sign':[[[1]],[[-1]]], 'C3_real_rotation':C3,'D3_real_standard':D3,
            'C3_rotation_plus_trivial':[[row+[0] for row in M]+[[0,0,1]] for M in C3]}
    out={}
    for name,G in groups.items():
        polys=[determinant_poly(g) for g in G]
        P=[sum(v)//len(G) for v in zip(*polys)]
        need(all(sum(v)%len(G)==0 for v in zip(*polys)),'exterior averaged coefficients integral')
        d,lead,_=multiplicity(P);minimum=min(fixed_dimension(g) for g in G)
        need(d==minimum,'real exterior valuation equals minimum fixed dimension')
        need(lead>0,'real complementary determinant positivity')
        out[name]={'polynomial':P,'root_order':d,'minimum_fixed_dimension':minimum}
    # C3 acting as scalar primitive omega on C^3: wedge degree j has
    # character omega^j, so only j=0,3 are invariant.
    d,_,_=multiplicity([1,0,0,1])
    need(d==1,'complex scalar root order')
    out['complex_C3_scalar_hostile']={'polynomial':[1,0,0,1],'root_order':d,'minimum_fixed_dimension':0}
    return out

if __name__=='__main__':
    result={'holonomy':holonomy_audit(),'fano':fano_audit(),'odd_actions':group_audit(),
            'literal_tournaments':tournament_audit(),'exterior':exterior_audit(),'common_action_bridge':common_action_audit()}
    result.update(status='PASS',checks=CHECKS)
    print(json.dumps(result,sort_keys=True,indent=2))

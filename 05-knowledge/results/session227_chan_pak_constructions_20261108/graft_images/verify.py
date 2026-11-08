#!/usr/bin/env python3
"""Portable exact controls for state hulls and graph-level graft images.
Only the standard library. Every guard is active under python -O.
Run: python3 -I -B verify.py (and with -O).
"""
from collections import Counter
from fractions import Fraction as F
import importlib.util
import itertools
import json
from math import gcd
from pathlib import Path

spec = importlib.util.spec_from_file_location('inherited', Path(__file__).with_name('inherited_verify.py'))
v = importlib.util.module_from_spec(spec)
spec.loader.exec_module(v)

def check(ok, message):
    if not ok:
        raise RuntimeError(message)

def facets(x):
    D,E,U,d=x
    return [E,d-E,D-d,U-D-d+E]

def lattice(x):
    D,E,U,d=x
    return all(F(z).denominator==1 for z in x) and D%2==d%2==1 and (U-E)%2==0

def c3_block(m):
    check(isinstance(m,int) and m>=1,'family domain')
    u, w=m,m+1
    rows=[sum(1<<j for j in range(i+1,m))|(1<<u) for i in range(m)]
    return rows+[1<<w,(1<<m)-1]

def transitive(n):
    return [sum(1<<j for j in range(i+1,n)) for i in range(n)]

def induced(rows, vertices):
    return [sum(1<<j for j,y in enumerate(vertices) if rows[x]>>y&1) for x in vertices]

def relabel(rows, perm):
    """Old vertex i receives new name perm[i]."""
    n=len(rows);check(sorted(perm)==list(range(n)),'permutation domain')
    out=[0]*n
    for i in range(n):
        out[perm[i]]=sum(1<<perm[j] for j in range(n) if rows[i]>>j&1)
    return out

def rooted_iso(a, ar, b, br):
    if len(a)!=len(b):
        return False
    other_a=[i for i in range(len(a)) if i!=ar]
    other_b=[i for i in range(len(b)) if i!=br]
    for tail in itertools.permutations(other_b):
        perm=[0]*len(a);perm[ar]=br
        for x,y in zip(other_a,tail):perm[x]=y
        if relabel(a,perm)==b:return True
    return False

def modules_at_root(rows, root, k):
    n=len(rows);out=[]
    if not 1<=k<n:return out
    for tail in itertools.combinations([i for i in range(n) if i!=root],k-1):
        vertices=[root]+list(tail);mask=sum(1<<i for i in vertices)
        ext=[i for i in range(n) if not mask>>i&1]
        if any((rows[i]&mask) not in (0,mask) for i in ext):continue
        losses=[i for i in ext if not rows[i]&mask]
        out.append((vertices,losses))
    return out

def image_witnesses(rows, root, b, br):
    out=[]
    for vertices,losses in modules_at_root(rows,root,len(b)):
        if len(losses)==1 and rooted_iso(induced(rows,vertices),0,b,br):
            out.append((vertices,losses[0]))
    return out

def inverse_b4(y):
    A,B,C,F0=map(F,y)
    return [C+5*B-(A+13*F0)/2,
            C+4*B-(5*A+11*F0)/4,
            F0-B,(A+3*F0-8*B)/8]

def labelled_images(source_order):
    """Build every source and every relabeling independently of module tests."""
    out=set()
    for a in v.tournaments(source_order):
        rows,root=v.graft(a,0,v.B4)
        old_other=[i for i in range(len(rows)) if i!=root]
        for tail in itertools.permutations(range(1,len(rows))):
            p=[0]*len(rows);p[root]=0
            for i,j in zip(old_other,tail):p[i]=j
            out.add(tuple(relabel(rows,p)))
    return out

def main():
    # Independent literal permutation profiles, no strongness or graph filters.
    counts=Counter();finite_max={}
    for n in range(1,6):
        c=0
        for a in v.tournaments(n):
            for root in range(n):
                x=v.profile(a,root)
                check(x==v.permutation_profile(a,root),'subset/permutation profile disagreement')
                check(all(z>=0 for z in facets(x)),'universal cone')
                check(lattice(x),'universal affine lattice')
                if x[0]==1:
                    check(x[2]==1<<a[root].bit_count(),'transitive H1 profile')
                counts['rooted_profiles']+=1;c+=1
        finite_max[str(n)]=c
    rays=[[1,1,1,1],[1,0,2,1],[1,0,1,0],[0,0,1,0]]
    check(abs(v.determinant(rays))==1,'unimodular cone coordinate matrix')
    for x in [[-3,2,5,1],[9,2,21,5],[1,0,6,1]]:
        coeff=facets(x)
        check([sum(coeff[i]*rays[i][j] for i in range(4)) for j in range(4)]==x,'cone decomposition identity')
    basis_graphs=[([0],0),([2,0],0),([6,4,0],0),([2,4,1],0),(c3_block(2),0)]
    states=[v.profile(a,r) for a,r in basis_graphs]
    expected=[[1,1,1,1],[1,0,2,1],[1,0,4,1],[3,1,3,1],[5,1,9,3]]
    check(states==expected,'affine-lattice physical witnesses')
    differences=[[s[j]-states[0][j] for j in range(4)] for s in states[1:]]
    check(abs(v.determinant(differences))==8,'affine-lattice index')
    z=[1,0,6,1]
    check(all(x>=0 for x in facets(z)) and lattice(z),'nonphysical hostile passes both hulls')
    check(z[2]&(z[2]-1)!=0,'nonphysical hostile is not a power of two')

    family=[]
    for m in range(1,11):
        b=c3_block(m);H=(1<<m)+1;Q=2*3**m+(1<<(m+2))-2
        check(v.strong(b),'cyclic family strong')
        for root in range(m+2):
            d=v.stats(b,root)
            check(d['H']==H and d['Q']==Q,'cyclic family H,Q')
            if m<=4:
                check(d==v.permutation_stats(b,root),'cyclic family permutation statistics')
                counts['family_permutation_statistics']+=1
            Delta=d['q']*d['Q']-d['H']*d['b']
            if root>=m:
                check(d['s']==d['e']==d['q']==1 and d['b']==1<<(m+1),'singleton root statistics')
                check(Delta==2*(3**m+2**m-1-4**m),'singleton Delta formula')
                if root==m:
                    check(v.profile(b,root)==[H,1,H,1],'cone third-ray sequence')
                if m>=2:check(Delta<0,'negative Delta')
            elif m>=2:
                check(d['q']==(1<<(m-1))+1,'block deletion H')
                check(d['b']==2*3**(m-1)+(1<<(m+1))-2,'block deletion Q')
                check(Delta==(2**m+4)*3**(m-1)+3*2**m,'positive Delta formula')
                check(d['s']==1<<(m-1-root) and d['e']==1<<root,'block root endpoints')
                check(Delta>0,'positive Delta')
            else:check(Delta==0,'triangle singular boundary')
            K=v.operator(b,root)
            check(v.determinant(K)==d['s']*d['e']*Delta,'family graft determinant')
            if m>=3 and m%2:
                check(all(c%3==0 for c in K[0]),'odd-block full-rank scalar factor three')
            counts['family_rooted_stats']+=1
        family.append(dict(m=m,order=m+2,H=H,Q=Q,singleton_delta=2*(3**m+2**m-1-4**m),
                           block_delta=((2**m+4)*3**(m-1)+3*2**m if m>=2 else 0)))
    for n in range(2,13):
        check(v.profile(transitive(n),0)==[1,0,1<<(n-1),1],'cone fourth-ray sequence')
        counts['transitive_limit_controls']+=1

    check(v.operator(v.B4,0)==v.K and v.determinant(v.K)==-8,'inherited B4 matrix')
    for y in itertools.product(range(-4,5),repeat=4):
        x=inverse_b4(y)
        integral=all(t.denominator==1 for t in x)
        check(integral==((y[0]+3*y[3])%8==0),'B4 exact integer image lattice')
        check(v.mv(v.K,x)==list(y),'B4 inverse right identity')
        counts['integer_image_grid']+=1
    image_differences=[v.mv(v.K,row) for row in differences]
    check(abs(v.determinant(image_differences))==64,'B4 physical affine-lattice image index')
    nonphysical_output=v.mv(v.K,z)
    check(nonphysical_output==[145,31,159,37] and inverse_b4(nonphysical_output)==z,'transported hostile')

    bank=[a for n in range(1,5) for a in v.tournaments(n)]
    affine_generators=[[2,0,0,0],[0,-1,1,0],[0,0,2,0],[0,0,0,2]]
    scalar_moduli=Counter()
    for b in bank:
        d=v.stats(b,0);K=v.operator(b,0)
        increments=[v.mv(K,z)[0] for z in affine_generators]
        H,Q=d['H'],d['Q']
        check(increments==[2*H,Q,2*(Q-H),2*(2*H-Q)],'scalar lattice increments')
        g=gcd(*increments)
        check(g==gcd(2*H,Q),'exact scalar affine lattice modulus')
        check(v.mv(K,states[0])[0]==H,'scalar affine lattice offset')
        scalar_moduli[str(g)]+=1
        counts['scalar_affine_lattice_gadgets']+=1
    for a in bank:
        x=v.profile(a,0)
        for b in bank:
            g,r=v.graft(a,0,b)
            check(v.profile(g,r)==v.mv(v.operator(b,0),x),'literal small graft closure')
            ws=image_witnesses(g,r,b,0)
            check(any(set(M)==set(range(len(a),len(g))) and old==0 for M,old in ws),'literal graph image witness')
            counts['literal_graft_pairs']+=1
        g,r=v.graft(a,0,v.B4)
        check(inverse_b4(v.profile(g,r))==x,'physical B4 inverse')
        # Full graph for a prescribed two-invertible-graft chain.
        h,rr=v.graft(g,r,v.B4)
        y=v.profile(h,rr)
        check(inverse_b4(inverse_b4(y))==x,'two-step inverse physical chain')
        check(v.mv(v.mm(v.K,v.K),x)==y,'two-step image map')
        counts['literal_two_step_chains']+=1

    image_census={}
    for n in [5,6]:
        expected_images=labelled_images(n-4)
        found=0;total=0
        for rows in v.tournaments(n):
            actual=bool(image_witnesses(rows,0,v.B4,0))
            check(actual==(tuple(rows) in expected_images),'exact graph recognition census')
            found+=actual;total+=1
        check(found==len(expected_images),'constructed graph image count')
        image_census[str(n)]={'all_labelled_root0_tournaments':total,'rooted_B4_images':found}

    good=[18,32,3,7,14,29];bad=[18,40,3,5,14,29]
    check(v.strong(good) and v.strong(bad),'collision strongness')
    for rows in (good,bad):
        check(v.profile(rows,0)==[37,7,43,9],'physical graph-image collision')
        check(v.permutation_profile(rows,0)==[37,7,43,9],'collision independent permutation counts')
    check(modules_at_root(bad,0,4)==[],'nongraft has no root-containing four-module')
    check(image_witnesses(good,0,v.B4,0)==[([0,2,3,4],1)],'good graph exact image witness')
    check(induced(good,[0,2,3,4])==v.B4,'good rooted module exact order')
    check(inverse_b4([37,7,43,9])==[1,0,2,1],'collision physical inverse')
    check(v.profile([2,0],0)==[1,0,2,1],'inverse physical edge source')

    print(json.dumps(dict(status='PASS',scope='Finite exact controls; all-order claims require PROOF.md',
        controls=dict(sorted(counts.items())),rooted_profile_universe=finite_max,
        filtering='None in the stated tournament universes; no strongness or image prefilter',
        family=family,cone_rays=rays,affine_lattice_generators=states,
        affine_lattice_index=8,B4_integer_image_index=8,B4_affine_image_index=64,
        small_gadget_scalar_affine_moduli=dict(sorted(scalar_moduli.items())),
        graph_image_census=image_census,
        cone_lattice_nonphysical_source=z,transported_nonphysical_source=nonphysical_output,
        physical_state_collision={'good':good,'bad':bad,'root':0,'state':[37,7,43,9],
             'common_inverse':[1,0,2,1],'good_witness':[[0,2,3,4],1],
             'bad_root_containing_four_modules':[]},
        no_all_odd_coverage_claim=True),indent=2,sort_keys=True))

if __name__=='__main__':
    main()

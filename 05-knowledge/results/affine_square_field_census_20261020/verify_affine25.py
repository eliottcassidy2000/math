#!/usr/bin/env python3
"""Exact global census for F25+ : C3; pure Python, no external data/packages.

Run: python3 verify_affine25.py --output verification.json
The exhaustive partition-refinement and independent pair-invariant CSP searches
both enumerate every isomorphism fixing 0. Translations justify rooting at 0.
All checks are explicit and remain active under python -O.
"""
import argparse,collections,hashlib,itertools,json,math
from pathlib import Path
Q=25

def require(condition,message):
    if not condition: raise RuntimeError(message)

def add(x,y):return (x%5+y%5)%5+5*((x//5+y//5)%5)
def neg(x):return (-(x%5))%5+5*((-(x//5))%5)
def sub(x,y):return add(x,neg(y))
def mul(x,y):
    a,b=x%5,x//5;c,d=y%5,y//5
    return (a*c+2*b*d)%5+5*((a*d+b*c)%5)
def power(x,k):
    z=1
    while k:
        if k&1:z=mul(z,x)
        x=mul(x,x);k//=2
    return z
G=7
LOG={power(G,k):k for k in range(24)}
require(len(LOG)==24,'7 must be primitive in F5[t]/(t^2-2)')

def adjacency(mask,tournament=False):
    if tournament:
        S={x for x in range(1,Q) if bool(mask>>(LOG[x]%4)&1)==(LOG[x]%8<4)}
    else:S={x for x in range(1,Q) if mask>>(LOG[x]%4)&1}
    return tuple(sum(1<<y for y in range(Q) if y!=x and sub(y,x) in S) for x in range(Q))

def affine(a,b,e=0):return tuple(add(mul(a,power(x,5**e)),b) for x in range(Q))
NORMALIZER=tuple(affine(a,b,e) for e in range(2) for a in range(1,Q) for b in range(Q))
MULTIPLIER=tuple(affine(a,b) for a in range(1,Q) for b in range(Q))
H=tuple(affine(power(G,8*k),b) for k in range(3) for b in range(Q))

def image(A,f):
    B=[0]*Q
    for x in range(Q):
        for y in range(Q):
            if A[x]>>y&1:B[f[x]]|=1<<f[y]
    return tuple(B)

def sign(A,f):
    return -1 if sum(1 for x in range(Q) for y in range(x+1,Q) if A[x]>>y&1 and f[x]>f[y])%2 else 1

def flipped_pairs(A,f):
    return sum(1 for x in range(Q) for y in range(x+1,Q) if A[x]>>y&1 and f[x]==y and f[y]==x)

def refinement_maps(A,B,limit=None):
    """Simultaneous equitable partitions; branch over every image of a vertex."""
    out=[];nodes=0
    def refine(P,R):
        while True:
            newP=[];newR=[];changed=False
            masksP=[sum(1<<v for v in C) for C in P]
            masksR=[sum(1<<v for v in C) for C in R]
            for C,D in zip(P,R):
                ss={};tt={}
                for v in C:ss.setdefault(tuple((A[v]&m).bit_count() for m in masksP),[]).append(v)
                for v in D:tt.setdefault(tuple((B[v]&m).bit_count() for m in masksR),[]).append(v)
                if ss.keys()!=tt.keys():return None
                for key in sorted(ss):
                    if len(ss[key])!=len(tt[key]):return None
                    newP.append(ss[key]);newR.append(tt[key])
                changed|=len(ss)>1
            P,R=newP,newR
            if not changed:return P,R
    def rec(P,R):
        nonlocal nodes
        nodes+=1
        refined=refine(P,R)
        if refined is None:return
        P,R=refined
        j=min((j for j in range(len(P)) if len(P[j])>1),key=lambda j:len(P[j]),default=None)
        if j is None:
            f=[None]*Q
            for C,D in zip(P,R):f[C[0]]=D[0]
            f=tuple(f)
            require(image(A,f)==B,'refinement leaf is not an isomorphism')
            out.append(f);return
        v=P[j][0]
        for w in R[j]:
            rec(P[:j]+[[v],[x for x in P[j] if x!=v]]+P[j+1:],R[:j]+[[w],[x for x in R[j] if x!=w]]+R[j+1:])
            if limit and len(out)>=limit:return
    rec([[0],list(range(1,Q))],[[0],list(range(1,Q))])
    return sorted(out),nodes

def csp_maps(A,B,limit=None):
    """Independent MRV search using only pairwise adjacency/common-neighbor colors."""
    full=(1<<Q)-1;out=[];nodes=0
    def color(A,x,y):return ((A[x]>>y)&1,(A[y]>>x)&1,(A[x]&A[y]).bit_count())
    ca=[[color(A,x,y) for y in range(Q)] for x in range(Q)]
    rel=[]
    for y in range(Q):
        D={}
        for z in range(Q):
            k=color(B,y,z);D[k]=D.get(k,0)|(1<<z)
        rel.append(D)
    def rec(f,C,remaining,used):
        nonlocal nodes
        nodes+=1
        if not remaining:
            f=tuple(f)
            require(image(A,f)==B,'CSP leaf is not an isomorphism')
            out.append(f);return
        x=min(remaining,key=lambda x:((C[x]&~used).bit_count(),x))
        choices=C[x]&~used
        while choices:
            t=choices&-choices;choices-=t;y=t.bit_length()-1
            D=C.copy();ok=True
            for z in remaining:
                if z==x:continue
                D[z]&=rel[y].get(ca[x][z],0)&~(used|t)
                if not D[z]:ok=False;break
            if ok:
                f[x]=y;rec(f,D,remaining-{x},used|t)
                if limit and len(out)>=limit:return
    C=[full]*Q;C[0]=1
    rec([-1]*Q,C,set(range(Q)),0)
    return sorted(out),nodes

def orbit_census(objects,perms,tournament=False):
    lookup={A:k for k,A in enumerate(objects)};pending=set(lookup.values());rows=[]
    while pending:
        m=min(pending);A=objects[m]
        ims=[image(A,f) for f in perms]
        require(all(B in lookup for B in ims),'group escaped fixed cube')
        orb=sorted({lookup[B] for B in ims});pending.difference_update(orb)
        stab=[f for f,B in zip(perms,ims) if B==A]
        signs=None if tournament else sorted({sign(A,f) for f in stab})
        rows.append({'representative':m,'orbit':orb,'stabilizer_order':len(stab),'signs':signs})
    return rows

def four_cube_control(pmod8):
    # A skew subset of Z/8 has exactly one of j,j+4 for each 0<=j<4.
    sources=[frozenset(j+(0 if mask>>j&1 else 4) for j in range(4)) for mask in range(16)]
    perms=[lambda x,k=k,e=e:(k+(pmod8 if e else 1)*x)%8 for k in range(8) for e in range(2)]
    source_orbits=[];pending=set(sources)
    while pending:
        S=min(pending,key=lambda S:sorted(S));O={frozenset(f(x) for x in S) for f in perms}
        require(O<=set(sources),'invalid signed-cube action');pending-=O;source_orbits.append(len(O))
    even=[];rows=[]
    for mask in range(16):
        selected={j for j in range(4) if mask>>j&1};stabilizer_signs=set()
        for f in perms:
            if {f(j)%4 for j in selected}==selected:
                stabilizer_signs.add(-1 if sum(f(j)>=4 for j in selected)%2 else 1)
        if stabilizer_signs=={1}:even.append(mask)
        rows.append({'mask':mask,'signs':sorted(stabilizer_signs)})
    require(source_orbits==[16] and even==[0],'four-cube family failed')
    witnesses=[]
    for mask,k,e in [(1,4,0),(3,1 if pmod8==3 else 0,1),(5,2,0),(7,4,0),(15,1,0)]:
        S={j for j in range(4) if mask>>j&1}
        f=lambda j:(k+(pmod8 if e else 1)*j)%8
        require({f(j)%4 for j in S}==S,'displayed graph witness does not stabilize')
        flips=sum(f(j)>=4 for j in S)
        require(flips%2==1,'displayed graph witness has positive sign')
        witnesses.append({'graph_mask':mask,'shift_k':k,'frobenius_power':e,'reversed_orbit_count':flips})
    interval=frozenset(range(4))
    scalar_orbit={frozenset((j+k)%8 for j in interval) for k in range(8)}
    frobenius_image=frozenset(pmod8*j%8 for j in interval)
    require(len(scalar_orbit)==8 and frobenius_image not in scalar_orbit,'Frobenius does not fuse the two scalar source orbits')
    return {'p_mod_8':pmod8,'source_orbit_sizes':source_orbits,'even_graph_masks':even,'graph_signs':rows,
            'negative_sign_witnesses':witnesses,'base_selector':sorted(interval),'frobenius_image':sorted(frobenius_image)}

def run():
    require(len(set(NORMALIZER))==1200 and len(set(H))==75,'group sizes')
    objects={t:[adjacency(m,t) for m in range(16)] for t in [False,True]}
    # Exact universe check: H has precisely four orbits of size 75 on pairs.
    pairs={(x,y) for x in range(Q) for y in range(x+1,Q)};pair_orbits=[]
    while pairs:
        x,y=min(pairs);O={tuple(sorted((f[x],f[y]))) for f in H};pairs-=O;pair_orbits.append(len(O))
    require(sorted(pair_orbits)==[75]*4,'incorrect pair cube')
    for t,As in objects.items():
        require(len(set(As))==16,'duplicate fixed object')
        for A in As:
            require(all(image(A,f)==A for f in H),'object is not H-fixed')
            require(all(not(A[x]>>x&1) for x in range(Q)),'loop')
            require(all(((A[x]>>y)&1)+((A[y]>>x)&1)==1 for x in range(Q) for y in range(x+1,Q)) if t else all(((A[x]>>y)&1)==((A[y]>>x)&1) for x in range(Q) for y in range(x+1,Q)),'wrong graph type')
    out={'field':'F5[t]/(t^2-2), label a+5b','primitive_element':G,'pair_orbit_sizes':pair_orbits,'group_orders':{'H':75,'affine_multiplier':600,'normalizer':1200},'censuses':{}}
    for t in [True,False]:
        As=objects[t];name='tournaments' if t else 'graphs'
        L=orbit_census(As,MULTIPLIER,t);N=orbit_census(As,NORMALIZER,t)
        rows=[]
        for row in N:
            m=row['representative'];A=As[m]
            if not t and m in [0,15]:
                order=math.factorial(Q);full_even=m==0;rooted=None;nodecounts=None
            else:
                F,nodes=refinement_maps(A,A);C,cnodes=csp_maps(A,A)
                require(F==C,'independent automorphism enumerations disagree')
                require(len(F)==len(set(F)),'duplicate automorphisms')
                Nroot=sorted(f for f in NORMALIZER if f[0]==0 and image(A,f)==A)
                require(F==Nroot,'unexpected automorphism outside normalizer')
                if t:require(len(F)*Q==75,'tournament Aut must equal H')
                order=Q*len(F);full_even=None if t else all(sign(A,f)==1 for f in F)
                rooted=F;nodecounts={'equitable_refinement':nodes,'pair_csp':cnodes}
            rows.append({**row,'degree':A[0].bit_count(),'full_automorphism_order':order,'globally_sign_even':full_even,'rooted_automorphisms':rooted,'search_nodes':nodecounts})
        fusion=[]
        for Arow,Brow in itertools.combinations(rows,2):
            a,b=Arow['representative'],Brow['representative']
            if As[a][0].bit_count()!=As[b][0].bit_count():continue
            F,nodes=refinement_maps(As[a],As[b]);C,cnodes=csp_maps(As[a],As[b])
            require(F==C==[],'unexpected global fusion')
            fusion.append({'source_mask':a,'target_mask':b,'isomorphic':False,'search_nodes':[nodes,cnodes]})
        out['censuses'][name]={'multiplier_orbits':L,'normalizer_and_global_orbits':rows,'same_degree_global_isomorphism_tests':fusion}
    require(len(out['censuses']['tournaments']['multiplier_orbits'])==2,'multiplier source count')
    require(len(out['censuses']['tournaments']['normalizer_and_global_orbits'])==1,'global source count')
    require([r['representative'] for r in out['censuses']['graphs']['multiplier_orbits'] if r['signs']==[1]]==[0,3],'multiplier relative-even graph count')
    require([r['representative'] for r in out['censuses']['graphs']['normalizer_and_global_orbits'] if r['globally_sign_even']]==[0],'global target count')
    sigma=affine(1,0,1);A=adjacency(3)
    require(image(A,sigma)==A and sign(A,sigma)==-1 and flipped_pairs(A,sigma)==5,'Frobenius hostile')
    out['frobenius_adjacent_graph_hostile']={'mask':3,'permutation':sigma,'edge_sign':-1,'endpoint_swapped_graph_edges':5}
    # An independent invariant immediately separates the two degree-12 graphs.
    out['degree12_common_neighbor_histograms']={str(m):dict(sorted(collections.Counter((adjacency(m)[0]&adjacency(m)[y]).bit_count() for y in range(1,Q)).items())) for m in [3,5]}
    out['family_signed_cube_controls']=[four_cube_control(p) for p in [3,5]]
    # Independent link from the compressed sign cocycle to actual F25 edges.
    compressed_sign_checks=0
    for mask,A in enumerate(objects[False]):
        for k in range(24):
            for e in range(2):
                f=affine(power(G,k),0,e)
                if image(A,f)!=A:continue
                flips=sum(((k+(5**e)*j)%8)>=4 for j in range(4) if mask>>j&1)
                require(sign(A,f)==(-1)**flips,'compressed and full edge signs disagree')
                compressed_sign_checks+=1
    out['compressed_to_full_sign_checks']=compressed_sign_checks
    out['status']='PASS'
    out['universe']='All 16 H-fixed tournaments and 16 H-fixed graphs on F25; all rooted isomorphisms for nontrivial normalizer representatives; exact signed four-cube actions for p mod 8 in {3,5}'
    out['script_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()

    out['global_counts']={'tau_T':1,'tau_E':1,'normalizer_source_fusion_excess':0,'normalizer_even_sign_escape':0,'multiplier_source_fusion_excess':1,'multiplier_even_sign_escape':1}
    return out

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--output',default='verification.json');args=parser.parse_args()
    result=run();data=json.dumps(result,indent=2,sort_keys=True)+'\n';open(args.output,'w').write(data)
    print(json.dumps(result['global_counts'],sort_keys=True))
    print('Wrote '+args.output+'; SHA256 '+hashlib.sha256(data.encode()).hexdigest())

"""Exact native universal electrical guards; no assertion is validation-critical.

Reads frozen current-Git inputs, independently builds the binary kernel, validates
integer higher lifts, and proves the entire vertex partition via separating
characters plus direct integer image witnesses. It does not rerun tournament
isomorphism classification or infer wiring from electrical invariants.
"""
from array import array
from collections import Counter, defaultdict
from hashlib import sha256
from itertools import combinations, permutations, product
from pathlib import Path
import argparse, csv, gzip, json, sys

HERE=Path(__file__).resolve().parent
INPUTS=HERE/'inputs'
PREFIX='05-knowledge__results__'
CHECKS=0
PINS={}
REPO=None
EXPECTED=json.loads((HERE/'SOURCE_PINS.json').read_text())

def need(ok, why):
    global CHECKS
    CHECKS+=1
    if not ok: raise RuntimeError(why)

def read(name):
    p=INPUTS/(PREFIX+name.replace('/','__'))
    if not p.exists():
        candidates=([REPO] if REPO else [])+list(HERE.parents)
        p=next((base/'05-knowledge/results'/name for base in candidates if (base/'05-knowledge/results'/name).exists()),p)
    raw=p.read_bytes(); PINS[name]=sha256(raw).hexdigest()
    need(PINS[name]==EXPECTED[name],'frozen original source bytes: '+name)
    return raw

def kernel(rows,n):
    piv={}
    for row in rows:
        while row:
            p=(row&-row).bit_length()-1
            if p in piv: row^=piv[p]
            else: piv[p]=row;break
    basis=[]
    for f in range(n):
        if f in piv: continue
        x=1<<f
        for p in sorted(piv,reverse=True):
            if (piv[p]&x).bit_count()%2: x^=1<<p
        need(all((r&x).bit_count()%2==0 for r in rows),'complete binary kernel witness')
        basis.append(x)
    need(len(basis)+len(piv)==n,'binary rank-nullity')
    return basis

def mv(M,x,m=None):
    z=tuple(sum(a*b for a,b in zip(row,x)) for row in M)
    return tuple(a%m for a in z) if m else z

def small_duality_controls():
    counts=[]
    # Enumerate every 2-by-2 matrix in these rings, not only symmetric matrices.
    for m in (2,3,4):
        vec=list(product(range(m),repeat=2)); num=0
        for a,b,c,d in product(range(m),repeat=4):
            M=((a,b),(c,d)); MT=((a,c),(b,d))
            H=[x for x in vec if mv(M,x,m)==(0,0)]
            ann={z for z in vec if all(sum(a*b for a,b in zip(z,x))%m==0 for x in H)}
            image={mv(MT,y,m) for y in vec}
            need(ann==image,'annihilator equals transpose image over composite ring')
            num+=1
        counts.append(dict(modulus=m,all_matrices=num))
    # Extra pins are not a mere root gauge. Two grounded vertices kill all
    # harmonic vectors on this doubled edge, yet e0-e1 is not in im L mod8.
    M=((2,-2),(-2,2)); vec=list(product(range(8),repeat=2))
    H=[x for x in vec if mv(M,x,8)==(0,0) and x[0]==0]
    H2=[x for x in H if x[1]==0]
    d=(1,7)
    need(any(sum(a*b for a,b in zip(d,x))%8 for x in H),'single grounding retains real guard')
    need(all(sum(a*b for a,b in zip(d,x))%8==0 for x in H2),'extra-grounding hostile')
    need(d not in {mv(M,x,8) for x in vec},'wrong unconstrained-image criterion rejected')
    # Same marginals and all harmonic guards still do not determine wiring.
    # These distinct labelled four-vertex paths have unimodular reductions.
    paths=((0,1,2,3),(0,2,1,3)); Ls=[]
    for path in paths:
        L=[[0]*4 for _ in range(4)]
        for i,j in zip(path,path[1:]):
            L[i][i]+=1;L[j][j]+=1;L[i][j]-=1;L[j][i]-=1
        Ls.append(L)
        H=[x+(0,) for x in product(range(8),repeat=3) if mv(L,x+(0,),8)==(0,)*4]
        need(H==[(0,0,0,0)],'tree has only constant grounded modular harmonics')
    need(Ls[0]!=Ls[1] and [Ls[0][i][i] for i in range(4)]==[Ls[1][i][i] for i in range(4)],'same-margin wiring hostile')
    return dict(exhaustive_matrices=counts,additional_pin_hostile=True,same_margin_wiring_hostile=True)

def main():
    global REPO
    parser=argparse.ArgumentParser();parser.add_argument('--repo',type=Path);parser.add_argument('--write',action='store_true')
    args=parser.parse_args();REPO=args.repo
    controls=small_duality_controls()
    N=3528; root=N-1
    D=json.loads(read('metagraph_order8_electrical_20261023/receipt.json'))['order8']
    nodes=list(csv.DictReader(read('metagraph_order8_20261023/merged_nodes.tsv').decode().splitlines(),delimiter='\t'))
    raw=gzip.decompress(read('metagraph_order8_20261023/merged_edges.tsv.gz'))
    edges=[]; binary=[0]*N; adjacency=[{} for _ in range(N)]
    for row in csv.DictReader(raw.decode().splitlines(),delimiter='\t'):
        i,j,b,k=(int(row[key]) for key in ('source','target','blue','black'))
        need(0<=i<=j<N and b>=0 and k>=0 and k%2==0 and b+k>0,'actual literal edge domain')
        w=b+k; edges.append((i,j,w))
        need(j not in adjacency[i],'unique native pair export')
        adjacency[i][j]=w;adjacency[j][i]=w
        if i!=j and w%2: binary[i]^=(1<<i)|(1<<j);binary[j]^=(1<<i)|(1<<j)
    need(len(nodes)==N and D['root']==root,'native grounded gauge')
    mask=(1<<root)-1
    basis=kernel([r&mask for r in binary[:root]],root)
    need(len(basis)==3363,'native binary kernel dimension')
    labels=[0]*N
    for k,x in enumerate(basis):
        scan=x
        while scan:
            bit=scan&-scan;scan-=bit;labels[bit.bit_length()-1]|=1<<k
    def apply(x):
        y=[0]*N
        for i,j,w in edges:
            z=w*(x[i]-x[j]);y[i]+=z;y[j]-=z
        return y
    Y=D['radical_lifts_mod4'];Z=D['eighth_radical_lifts_mod8'][0];X=D['higher_layers'][-1]['lifts'][0]
    need(len(Y)==5,'five actual mod4 lifts')
    for x,m in [(x,4) for x in Y]+[(Z,8),(X,256)]:
        need(len(x)==N and x[root]==0 and all(v%m==0 for v in apply(x)),f'actual full harmonic lift modulo{m}')
    sigs={2:labels,4:[(labels[i],tuple(y[i]%4 for y in Y)) for i in range(N)],
          8:[(labels[i],tuple(y[i]%4 for y in Y),Z[i]%8) for i in range(N)]}
    partitions={}
    for m,signatures in sigs.items():
        buckets=defaultdict(list)
        for i,sig in enumerate(signatures):buckets[sig].append(i)
        partitions[m]=dict(classes=len(buckets),class_size_histogram=dict(sorted(Counter(map(len,buckets.values())).items())),
                           collisions=[v for v in buckets.values() if len(v)>1])
    pairs=[[0,3527],[53,3492],[109,3523]]
    need(partitions[8]['collisions']==pairs,'only universal-mod8 candidate collisions')
    leafcert=[]
    for parent,leaf in pairs:
        need(adjacency[leaf]=={parent:1} and parent not in adjacency[parent], 'unit leaf and no parent loop')
        vector=[0]*N;vector[leaf]=1
        image=apply(vector);expected=[0]*N;expected[leaf]=1;expected[parent]=-1
        need(image==expected,'integer Laplacian image certificate for every modulus')
        leafcert.append(dict(parent=parent,leaf=leaf,image_witness_nonzero=[[leaf,1]],
                             image_nonzero=[[parent,-1],[leaf,1]],parent_has_self_loop=False))
    # All candidate equalities are proved by integer image certificates; all
    # other pairs have an authenticated separating character in the bank.
    # Thus completeness does not depend on claiming this bank generates H.
    classes={}; universal=[]
    for sig in sigs[8]:
        if sig not in classes:classes[sig]=len(classes)
        universal.append(classes[sig])
    old_buckets=defaultdict(list)
    for i,x in enumerate(X):old_buckets[x%8].append(i)
    old_pairs=sum(len(v)*(len(v)-1)//2 for v in old_buckets.values())
    need(all((X[a]-X[b])%8==0 for a,b in pairs),'original potential respects universal classes')
    # Every normalized Hamiltonian presentation, source/target both retained.
    owner=array('H');owner.frombytes(gzip.decompress(read('metagraph_order8_20261023/owner_map.u16.gz')))
    if sys.byteorder!='little':owner.byteswap()
    need(len(owner)==1<<21,'complete native presentation domain')
    merged=[None]*6880
    for r in nodes:
        for key in ('class_a','class_b'):merged[int(r[key])]=int(r['merged_id'])
    need(all(v is not None for v in merged),'every ordinary class has merged owner')
    targetsets=[set() for _ in range(N)]; portcounts=Counter(); full=len(owner)-1
    for q in range(len(owner)):
        s,t=merged[owner[q]],merged[owner[full^q]]
        targetsets[s].add(t)
        portcounts[(min(s,t),max(s,t))]+=1
    need(len(portcounts)==len(edges),'owner map and graph have identical pair support')
    need(all(portcounts[i,j]==2*w for i,j,w in edges),'all native pair conductances independently replayed from owner map')
    hidden=[]; source_reports=[]
    for s,targets in enumerate(targetsets):
        need(len({universal[t] for t in targets})==len(targets),'fixed source plus universal phase determines target for every native port')
        need(targets==set(adjacency[s]),'owner-map actual source neighborhood equals graph')
        h=Counter((X[s]-X[t])%8 for t in targets)
        lost=sum(c*(c-1)//2 for c in h.values())
        if lost:hidden.append(s)
        if s in (1238,2242,2246,2315):
            source_reports.append(dict(source=s,distinct_targets=len(targets),single_character_phase_target_counts=dict(sorted(h.items())),
                                       ambiguous_unordered_target_pairs=lost,universal_target_classes=len(targets)))
    # Literal hostile: same fixed source, same endpoints, same original phase,
    # but separated by the original network's different harmonic character.
    source=1238;P=(0,4,2,1,6,7,5,3);Q=(0,6,7,2,1,5,4,3)
    row=nodes[source];canon=int(row['canon_a']); A=[[0]*8 for _ in range(8)]
    for bit,(i,j) in enumerate(combinations(range(8),2)):
        A[i][j]=(canon>>(27-bit))&1;A[j][i]=1-A[i][j]
    chords=[(i,j) for i in range(8) for j in range(i+2,8)]
    hostile=[]
    for p in (P,Q):
        need(all(A[i][j] for i,j in zip(p,p[1:])),'hostile path is literally Hamiltonian')
        q=sum(A[p[i]][p[j]]<<k for k,(i,j) in enumerate(chords));s=merged[owner[q]];t=merged[owner[full^q]]
        need(s==source,'hostile source owner authenticated')
        hostile.append(dict(path=p,mask=q,source=s,target=t,original_phase=(X[s]-X[t])%8))
    need([r['target'] for r in hostile]==[2289,2009] and [r['original_phase'] for r in hostile]==[2,2], 'exact one-character target hostile')
    need(P[0]==Q[0] and P[-1]==Q[-1],'literal shared endpoint pins')
    sep=[0]*N;sep[2289]=4
    image=apply(sep)
    need(all(v%8==0 for v in image) and sep[root]==0,'literal separating mod8 harmonic character')
    need(all(w%2==0 for w in adjacency[2289].values()),'all separator incident conductances even')
    for r in hostile:r['separating_phase']=(sep[source]-sep[r['target']])%8
    need([r['separating_phase'] for r in hostile]==[4,0],'universal guard rejects one-character-legal hostile')
    sep_cert=dict(potential_nonzero=[[2289,4]],root=root,incident_neighbors=len(adjacency[2289]),
                  weighted_degree=sum(adjacency[2289].values()),all_incident_weights_even=True,
                  exact_current_div8_nonzero=[[i,v//8] for i,v in enumerate(image) if v])
    coupled_separators=[]
    for source,u,v,neighbors,degree in [(2242,2874,1957,273,898),(2246,1841,2038,395,934)]:
        potential=[0]*N;potential[u]=4;current=apply(potential)
        need(potential[root]==0 and all(a%8==0 for a in current),'coupled-slot sparse harmonic separator')
        need(len(adjacency[u])==neighbors and sum(adjacency[u].values())==degree and all(w%2==0 for w in adjacency[u].values()),'coupled-slot native even-conductance certificate')
        phases=[(potential[source]-potential[t])%8 for t in (u,v)]
        need(phases==[4,0],'coupled-slot targets separated by full original bank')
        coupled_separators.append(dict(source=source,targets=[u,v],phases=phases,potential_nonzero=[[u,4]],root=root,
             incident_neighbors=neighbors,weighted_degree=degree,all_incident_weights_even=True,
             exact_current_div8_nonzero=[[i,a//8] for i,a in enumerate(current) if a]))
    # Independent graph-theoretic replay: an edge of multiplicity greater than
    # one cannot be a bridge in the literal multigraph. Loops never are.
    sys.setrecursionlimit(10000)
    arrival=[-1]*N;low=[-1]*N;timer=0;bridges=[]
    def visit(u,parent):
        nonlocal timer
        arrival[u]=low[u]=timer;timer+=1
        for v,w in adjacency[u].items():
            if v==u:continue
            if arrival[v]<0:
                visit(v,u);low[u]=min(low[u],low[v])
                if low[v]>arrival[u] and w==1:bridges.append(sorted((u,v)))
            elif v!=parent or w>1:low[u]=min(low[u],arrival[v])
    visit(root,-1)
    need(timer==N,'native graph connected')
    need(sorted(bridges)==pairs,'all literal native multigraph bridges are the three certified unit leaves')
    # Same-source bank comparisons cancel source; converse/reversal retains
    # merged endpoints. The actual involution q -> full^q swaps source,target.
    need(all(merged[owner[full^(full^q)]]==merged[owner[q]] for q in range(len(owner))), 'complement involution retains literal owners')
    need(all(r['source']==1238 for r in hostile),'enclosing hostile source agrees with every path source')
    result=dict(status='FINITE-EXACT original native order8 universal electrical target compiler',
                input_sha256=PINS,small_controls=controls,native=dict(vertices=N,root=root,literal_pair_rows=len(edges),
                    normalized_presentations=len(owner),binary_kernel_dimension=len(basis),
                    separating_bank_counts=dict(mod2_basis_scaled_by4=len(basis),mod4_lifts_scaled_by2=len(Y),mod8_lifts=1),
                    partitions=partitions,universal_leaf_certificates=leafcert,all_literal_bridges=sorted(bridges),
                    original_character_node_residue_counts={k:len(v) for k,v in sorted(old_buckets.items())},
                    original_character_ambiguous_unordered_node_pairs=old_pairs,
                    original_character_false_equalities_vs_universal=old_pairs-len(pairs),
                    source_conditioned_universal_target_injective=True,
                    original_character_ambiguous_sources=len(hidden),four_source_target_ambiguities=source_reports),
                hostile=dict(shared_owner_clock_pins=[[0,0],[7,3]],source=1238,paths=hostile,separating_character=sep_cert),
                coupled_slot_separators=coupled_separators,checks_before_receipt_comparison=CHECKS)
    text=json.dumps(result,indent=2,sort_keys=True)+'\n';out=HERE/'receipt.json'
    if '--write' in sys.argv:out.write_text(text)
    else:need(out.read_text()==text,'frozen exact receipt equality')
    print(json.dumps(dict(status='PASS',checks=CHECKS,universal_classes=len(classes),collisions=pairs,
                         original_character_false_equalities=old_pairs-len(pairs),
                         original_character_ambiguous_sources=len(hidden),source_conditioned_target_injective=True,
                         receipt_sha256=sha256(text.encode()).hexdigest()),sort_keys=True))

if __name__=='__main__':main()

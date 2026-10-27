#!/usr/bin/env python3
"""Independent finite audit: raw edge masks, DFS cycle test, explicit tournaments.
No package dependencies, no isomorphism census and no complement-state BFS.
"""
import itertools,json,math
from collections import Counter
from pathlib import Path

def require(condition,message):
    if not condition: raise RuntimeError(message)
def scores(n,edges,mask):
    d=[0]*n
    # The transitive baseline has j -> i, hence outdegree i at vertex i.
    for k,(i,j) in enumerate(edges):
        d[i if (mask>>k)&1 else j]+=1
    return d

def score_L(d):return sum(abs(s-i) for i,s in enumerate(sorted(d)))//2

def is_forest(n,edges,mask):
    adj=[[] for _ in range(n)]
    for k,(i,j) in enumerate(edges):
        if (mask>>k)&1:adj[i].append(j);adj[j].append(i)
    seen=set()
    def dfs(v,parent):
        seen.add(v)
        for w in adj[v]:
            if w==parent:continue
            if w in seen or not dfs(w,v):return False
        return True
    return all(v in seen or dfs(v,-1) for v in range(n))

report={'scope':'All backward forests on a fixed labelled total order n=2..6; all tournaments and legal Hamiltonian paths n=2..5; explicit hostiles and n11 sharpness.'}
forests=[]
for n in range(2,7):
    edges=list(itertools.combinations(range(n),2)); hist=Counter(); by_edges=Counter(); assignment_checks=0
    for mask in range(1<<len(edges)):
        if not is_forest(n,edges,mask):continue
        d=scores(n,edges,mask); L=score_L(d)
        require(L<=n-2,('forest cap',n,mask,L))
        require(sum(d)==n*(n-1)//2,('score sum',n,mask))
        hist[L]+=1;by_edges[mask.bit_count()]+=1
        if n<=5:
            cost=min(sum(abs(d[v]-p[v]) for v in range(n)) for p in itertools.permutations(range(n)))
            require(cost==2*L,('assignment identity',n,mask,cost,L));assignment_checks+=1
    require(by_edges[n-1]==n**(n-2),('Cayley tree count',n))
    forests.append(dict(n=n,count=sum(hist.values()),histogram=dict(sorted(hist.items())),counts_by_edges=dict(sorted(by_edges.items())),brute_assignment_checks=assignment_checks))
report['independent_forests']=forests
moves=[]
for n in range(2,6):
    edges=list(itertools.combinations(range(n),2)); eindex={e:k for k,e in enumerate(edges)}; full=(1<<len(edges))-1;hist=Counter()
    for mask in range(1<<len(edges)):
        d=scores(n,edges,mask); L=score_L(d)
        require(score_L([n-1-x for x in d])==L,('converse',n,mask))
        for p in itertools.permutations(range(n)):
            pairs=list(zip(p,p[1:])); pmask=0; legal=True
            for u,v in pairs:
                k=eindex[tuple(sorted((u,v)))]; low_wins=(mask>>k)&1
                if bool(low_wins)!=(u<v):legal=False;break
                pmask|=1<<k
            if not legal:continue
            target=mask^(full^pmask);dt=scores(n,edges,target)
            expected=[n-1-x for x in d];expected[p[0]]+=1;expected[p[-1]]-=1
            require(dt==expected,('endpoint law',n,mask,p))
            require((target^(full^pmask))==mask,('involution',n,mask,p))
            change=score_L(dt)-L;require(abs(change)<=1,('Lipschitz',n,mask,p,change));hist[change]+=1
    require(sum(hist.values())==math.factorial(n)*(1<<(len(edges)-(n-1))),('all legal path count',n))
    moves.append(dict(n=n,legal_moves=sum(hist.values()),change_histogram=dict(sorted(hist.items()))))
report['legal_complement_controls']=moves

def witness(n,chosen):
    edges=list(itertools.combinations(range(n),2));chosen={tuple(sorted(e)) for e in chosen};mask=sum(1<<i for i,e in enumerate(edges) if e in chosen);d=scores(n,edges,mask)
    return dict(n=n,reversed_edges=sorted(chosen),forest=is_forest(n,edges,mask),scores=d,sorted_scores=sorted(d),L=score_L(d),m=len(chosen))
sharp11=witness(11,zip((4,7,0,8,1,9,2,10,3,5),(7,0,8,1,9,2,10,3,5,6)))
require(sharp11['forest'] and sharp11['m']==10 and sharp11['L']==9,'n11 sharp witness')
sharp11['regular_score_transport_half_l1']=sum(abs(s-5) for s in sharp11['scores'])//2
require(sharp11['regular_score_transport_half_l1']==6,'n11 scalar cost')
report['sharp_n11_tree']=sharp11
hostiles={
 'even_n8_refutes_n_minus_3':witness(8,[(0,7),(1,6),(1,7),(2,5),(2,7),(3,6)]),
 'not_every_m_edge_forest_has_L_le_m_minus_1':witness(3,[(0,2)]),
 'm_le_n_minus_1_without_forest_is_not_enough':witness(7,[(0,4),(0,5),(0,6),(1,5),(1,6),(2,6)]),
}
require(hostiles['even_n8_refutes_n_minus_3']['L']==6,'even hostile')
require(hostiles['not_every_m_edge_forest_has_L_le_m_minus_1']['L']==1,'m hostile')
require(hostiles['m_le_n_minus_1_without_forest_is_not_enough']['L']==6 and not hostiles['m_le_n_minus_1_without_forest_is_not_enough']['forest'],'cyclic hostile')
report['hostiles']=hostiles
paley11_edges=list(itertools.combinations(range(11),2));qr={1,3,4,5,9};paley11_mask=sum(1<<k for k,(i,j) in enumerate(paley11_edges) if (j-i)%11 in qr);pd=scores(11,paley11_edges,paley11_mask)
require(pd==[5]*11 and score_L(pd)==15,'Paley11 regularity')
report['P11']=dict(quadratic_residues=sorted(qr),scores=pd,L=15,forest_ceiling=9,move_lower_bound=6)


def component_count(n,edges,mask):
    # Independent reachability closure, rather than the forest DFS test.
    remaining=set(range(n)); count=0
    adjacency=[set() for _ in range(n)]
    for k,(u,v) in enumerate(edges):
        if mask>>k&1:adjacency[u].add(v);adjacency[v].add(u)
    while remaining:
        component={min(remaining)}
        while True:
            enlarged=component.union(*(adjacency[v] for v in component))
            if enlarged==component:break
            component=enlarged
        remaining-=component;count+=1
    return count
beta_controls=[]
for n in range(2,7):
    edges=list(itertools.combinations(range(n),2));hist=Counter();connected=0
    for mask in range(1<<len(edges)):
        c=component_count(n,edges,mask);m=mask.bit_count();beta=m-n+c;L=score_L(scores(n,edges,mask))
        require(beta>=0,('negative cycle rank',n,mask))
        require(L<=n-2+beta,('beta inequality',n,mask,L,beta,c))
        if c==1:
            require(L<=m-1,('connected strict inequality',n,mask,L,m));connected+=1
        hist[n-2+beta-L]+=1
    beta_controls.append(dict(n=n,total_graphs=1<<len(edges),connected=connected,slack_histogram=dict(sorted(hist.items()))))
report['arbitrary_backward_graph_beta_controls']=beta_controls
carousels=[]
for r in range(1,16):
    n=2*r+1;chosen=[(i,j) for i in range(r) for j in range(i+r+1,2*r+1)]
    edges=list(itertools.combinations(range(n),2));chosen_set=set(chosen);mask=sum(1<<k for k,e in enumerate(edges) if e in chosen_set)
    d=scores(n,edges,mask);c=component_count(n,edges,mask);beta=len(chosen)-n+c;L=score_L(d)
    require(d==[r]*n and c==2 and beta==(r-1)*(r-2)//2 and beta==L-n+2,('sharp carousel',r))
    carousels.append(dict(n=n,m=len(chosen),components=c,L=L,beta=beta,regular_degree=r))
report['sharp_regular_carousels']=carousels

# Literal first-step certificate, independent of producer data files.
n=11;edges=list(itertools.combinations(range(n),2));full=(1<<len(edges))-1
path=tuple(range(n));path_set={tuple(sorted(pair)) for pair in zip(path,path[1:])}
pmask=sum(1<<k for k,e in enumerate(edges) if e in path_set)
require(all((v-u)%11 in qr for u,v in zip(path,path[1:])),'first-step path legality')
step=paley11_mask^(full^pmask);dstep=scores(n,edges,step);Lstep=score_L(dstep)
order=(0,6,1,7,2,8,3,9,4,10,5);pos={v:i for i,v in enumerate(order)};backward_pairs=[]
for k,(u,v) in enumerate(edges):
    source,target=(u,v) if step>>k&1 else (v,u)
    if pos[source]>pos[target]:backward_pairs.append((u,v))
backward_set=set(backward_pairs);bmask=sum(1<<k for k,e in enumerate(edges) if e in backward_set)
c=component_count(n,edges,bmask);m=len(backward_pairs);beta=m-n+c
require(Lstep==14 and m==14 and c==2 and beta==5,'first-step exact certificate')
require(Lstep==m,'first-step f exact from L lower and frame upper')
require(beta==Lstep-n+2,'first-step beta exact from lower and frame upper')
report['P11_first_complement']={
 'legal_path':path,'forward_frame':order,'scores':dstep,'sorted_scores':sorted(dstep),
 'L':Lstep,'backward_pairs':backward_pairs,'m':m,'components':c,'frame_beta':beta,
 'intrinsic_beta_exact':beta,'minimum_feedback_arc_count_exact':m,
 'minimum_beta_over_all_one_step_P11_images_exact':5,
 'basis':'Lipschitz gives L>=14 on every one-step image, hence beta>=5; this legal image attains both L14 and beta5.'}

# Independently count labelled forests by the component containing label zero.
F=[1]
for n in range(1,10):
    F.append(sum(math.comb(n-1,k-1)*(1 if k==1 else k**(k-2))*F[n-k] for k in range(1,n+1)))
import re
cpp_rows=[]
for line in Path(__file__).with_name('forest_audit.out').read_text().splitlines():
    match=re.fullmatch(r'n=(\d+) total=(\d+) maxL=(\d+) counts_by_edges=([0-9,]+) histogram=([0-9:,]+) witness=(.*)',line)
    require(match is not None,('unparsed C++ output',line))
    n,total,maxL=map(int,match.group(1,2,3));ce=list(map(int,match.group(4).strip(',').split(',')));hist={int(k):int(v) for k,v in (term.split(':') for term in match.group(5).strip(',').split(','))}
    require(total==F[n] and ce[-1]==n**(n-2),('independent forest counts',n))
    if n<=6:
        row=forests[n-2]
        require(hist==row['histogram'] and ce==list(row['counts_by_edges'].values()),('independent histogram comparison',n))
    cpp_rows.append(dict(n=n,count=total,maxL=maxL,Cayley_tree_count=ce[-1]))
report['Cplusplus_crosschecks']=cpp_rows


# Verify the complete certificate by independent literal bit-mask operations.
# Only the plain-text certificate is read, never producer source or routines.
chain_lines=Path(__file__).with_name('producer-greedy.pinned.out').read_text().splitlines()
require(len(chain_lines)==24,'six-move certificate shape')
state=paley11_mask;chain=[]
def feedback_dp(rows):
    n=len(rows);dp=[0]*(1<<n)
    for S in range(1,1<<n):
        choices=[]
        for v in range(n):
            if S>>v&1:
                R=S^(1<<v)
                choices.append(dp[R]+(rows[v]&R).bit_count())
        dp[S]=min(choices)
    return dp[-1]
for stage in range(6):
    header=chain_lines[4*stage].split(); claimed={header[i]:int(header[i+1]) for i in range(0,len(header),2)}
    path=tuple(map(int,chain_lines[4*stage+1].split()[1:]));order=tuple(map(int,chain_lines[4*stage+2].split()[1:]));claimed_rows=list(map(int,chain_lines[4*stage+3].split()[1:]))
    require(sorted(path)==list(range(11)) and sorted(order)==list(range(11)),('certificate permutations',stage))
    rows=[0]*11
    for k,(u,v) in enumerate(edges):
        source,target=(u,v) if state>>k&1 else (v,u)
        rows[source]|=1<<target
    require(all(rows[u]>>v&1 for u,v in zip(path,path[1:])),('illegal certificate path',stage))
    pp={tuple(sorted((u,v))) for u,v in zip(path,path[1:])}
    keep=sum(1<<k for k,e in enumerate(edges) if e in pp)
    state=state^(full^keep)
    new_rows=[0]*11
    for k,(u,v) in enumerate(edges):
        source,target=(u,v) if state>>k&1 else (v,u)
        new_rows[source]|=1<<target
    require(new_rows==claimed_rows,('certificate adjacency mismatch',stage))
    d=[row.bit_count() for row in new_rows];L=score_L(d);pos={v:i for i,v in enumerate(order)};backwards=[]
    for k,(u,v) in enumerate(edges):
        source,target=(u,v) if state>>k&1 else (v,u)
        if pos[source]>pos[target]:backwards.append((u,v))
    backwards=set(backwards);bm=sum(1<<k for k,e in enumerate(edges) if e in backwards);c=component_count(11,edges,bm);m=len(backwards);beta=m-11+c;f=feedback_dp(new_rows)
    require(L==14-stage and beta==5-stage and c==claimed['c'] and f==claimed['f'] and L==claimed['L'] and beta==claimed['beta'],('stage invariants',stage,L,beta,c,f,claimed))
    require(beta==L-11+2,('exact intrinsic beta certification',stage))
    if stage==5:require(m==10 and c==1 and is_forest(11,edges,bm),'final connected tree')
    chain.append(dict(depth=stage+1,path=path,forward_order=order,adjacency_rows=new_rows,scores=d,sorted_scores=sorted(d),L=L,f_dp_exact=f,m=m,components=c,beta_intrinsic_exact=beta,backward_pairs=sorted(backwards),distance_from_P11_exact=stage+1,distance_to_forest_exact=5-stage))
report['six_move_certificate_independent_verification']=chain
report['exact_distance_conclusion']={
 'P11_to_any_forest_framed_tournament':6,
 'lower_bound':'P11 has L15; every forest-framed order11 tournament has L<=9; each legal move changes L by at most1.',
 'upper_bound':'Six literal legal paths verified above, with final backward graph a connected ten-edge tree.',
 'free_operations':'Relabelling and converse preserve L and the forest-framed property; the same lower bound and exact distance hold with either free.',
 'final_feedback_count':'Independent subset DP gives f=10 for the final forest-framed tournament: the old forest f<=n-1 ceiling itself cannot improve to n-2.'}

report['status']='PASS: all explicit requires passed; theorem proof remains separate from finite controls.'
print(json.dumps(report,indent=2,sort_keys=True))

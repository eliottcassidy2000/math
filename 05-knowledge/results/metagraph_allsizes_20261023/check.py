"""A realized all-depth obstruction to vertex-incidence refinement.

Two explicit 13-vertex circulants; independent Hamilton DP directions;
all vertex tags, pair tags and root stabilizers; exact extension controls.
This is not a census of all thirteen-vertex tournaments.
"""
from collections import Counter
from functools import lru_cache
from hashlib import sha256
from itertools import combinations, permutations
from pathlib import Path
import argparse
import json

HERE = Path(__file__).resolve().parent
CHECKS = 0


def need(ok, why):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise ValueError(why)


def enc(obj):
    return (json.dumps(obj,sort_keys=True,separators=(',',':'))+'\n').encode()


def circulant(n, differences):
    need(type(n) is int and n >= 3 and n % 2 == 1, 'odd literal order')
    need(type(differences) is tuple and all(type(x) is int and 0 < x < n for x in differences), 'literal difference set')
    need(len(set(differences)) == len(differences) == (n-1)//2, 'one choice per opposite pair')
    need(set(differences).isdisjoint({n-x for x in differences}), 'no opposite pair')
    A = tuple(tuple(int((v-u)%n in differences) for v in range(n)) for u in range(n))
    audit(A)
    return A


def audit(A):
    n = len(A)
    need(type(A) is tuple and n > 0 and all(type(r) is tuple and len(r) == n for r in A), 'literal square tuple')
    need(all(type(x) is int and x in (0,1) for r in A for x in r), 'literal binary adjacency')
    need(all(A[i][i] == 0 for i in range(n)), 'zero diagonal')
    need(all(A[i][j]+A[j][i] == 1 for i,j in combinations(range(n),2)), 'one arc per pair')


def tags(A):
    n=len(A); cyc=[]; v3=[0]*n; v4=[0]*n; ce=Counter()
    for i,j,k in combinations(range(n),3):
        if A[i][j] == A[j][k] == A[k][i]:
            cyc.append((i,j,k))
            for v in (i,j,k):v3[v]+=1
            for edge in combinations((i,j,k),2):ce[edge]+=1
    strong4=[]
    for vertices in combinations(range(n),4):
        if sorted(sum(A[u][v] for v in vertices) for u in vertices) == [1,1,2,2]:
            strong4.append(vertices)
            for v in vertices:v4[v]+=1
    need(sum(c*(c-1)//2 for c in ce.values()) == len(strong4), 'independent edge-overlap / strong-four identity')
    return tuple((sum(A[v]),v3[v],v4[v]) for v in range(n)), ce, len(cyc), len(strong4)


def edge_profiles(A):
    n=len(A); full=Counter(); local=[]
    for u in range(n):
        row=Counter()
        for v in range(n):
            if not A[u][v]:continue
            p=(sum(A[u][w]*A[v][w] for w in range(n)),
               sum(A[w][u]*A[w][v] for w in range(n)),
               sum(A[v][w]*A[w][u] for w in range(n)))
            row[p]+=1; full[p]+=1
        local.append(tuple(sorted(row.items())))
    return tuple(sorted(full.items())),tuple(local)


def h_forward(A):
    n=len(A); rows=[sum(A[u][v]<<v for v in range(n)) for u in range(n)]
    dp=[[0]*n for _ in range(1<<n)]
    for u in range(n):dp[1<<u][u]=1
    for mask in range(1,1<<n):
        for u,value in enumerate(dp[mask]):
            if not value:continue
            left=rows[u]&~mask
            while left:
                bit=left&-left;left-=bit
                dp[mask|bit][bit.bit_length()-1]+=value
    return sum(dp[-1]),tuple(dp[-1])


def h_backward_start(A, start):
    """Independent remaining-set recursion, for one literal initial vertex."""
    n=len(A)
    @lru_cache(None)
    def suffix(u, left):
        if not left:return 1
        return sum(suffix(v,left^(1<<v)) for v in range(n) if left>>v&1 and A[u][v])
    return suffix(start,((1<<n)-1)^(1<<start)),suffix.cache_info().currsize


def stabilizer(A):
    """Exhaust all adjacency-preserving permutations fixing vertex zero."""
    n=len(A); image={0:0}; used={0}; solutions=[]; visited=0
    def search():
        nonlocal visited
        visited+=1
        if len(image)==n:
            solutions.append(tuple(image[u] for u in range(n)))
            return
        possibilities=[]
        for u in range(n):
            if u in image:continue
            candidates=[v for v in range(n) if v not in used and all(A[u][x] == A[v][y] for x,y in image.items())]
            if not candidates:return
            possibilities.append((len(candidates),u,candidates))
        _,u,candidates=min(possibilities)
        for v in candidates:
            image[u]=v;used.add(v);search();used.remove(v);del image[u]
    search()
    for p in solutions:
        need(all(A[i][j]==A[p[i]][p[j]] for i in range(n) for j in range(n)), 'literal stabilizer permutation')
    return solutions,visited


def extend_sources(A,k):
    need(type(k) is int and k >= 0, 'literal nonnegative prefix length')
    n=len(A)
    B=tuple(tuple(int(i<j) if i<k else (0 if j<k else A[i-k][j-k]) for j in range(n+k)) for i in range(n+k))
    audit(B)
    return B


def stable_paired_partition(A,B):
    """Actual refinement to stability on the disjoint carrier, no fixed-round cap.

    Inclusion of the old color makes each round a refinement. Equal color
    sets/multiplicities and identical stable quotient prove all later rounds.
    Incoming as well as outgoing colors are retained in this stronger control.
    """
    n=len(A);need(len(B)==n,'equal size'); colors=[list(tags(T)[0]) for T in (A,B)]
    def compress(records):
        keys=sorted(set(records[0]+records[1])); lookup={x:i for i,x in enumerate(keys)}
        return [[lookup[x] for x in row] for row in records]
    colors=compress(colors); rounds=0
    while True:
        need(Counter(colors[0])==Counter(colors[1]),'refinement signatures still equal')
        records=[]
        for T,row in zip((A,B),colors):
            records.append([(row[u],tuple(sorted(row[v] for v in range(n) if T[u][v])),
                            tuple(sorted(row[v] for v in range(n) if T[v][u]))) for u in range(n)])
        new=compress(records);rounds+=1
        oldclasses=len(set(colors[0]+colors[1])); newclasses=len(set(new[0]+new[1]))
        need(newclasses>=oldclasses,'old color retained')
        colors=new
        if newclasses==oldclasses:break
        need(rounds<=2*n,'finite partition refinement bound')
    need(Counter(colors[0])==Counter(colors[1]),'identical stable quotient multiplicities')
    return {'rounds_to_stability':rounds,'classes_per_tournament':len(set(colors[0])),
            'class_sizes':sorted(Counter(colors[0]).values())}


def controls():
    # Exhaustive literal permutation controls, independent of either DP.
    for n in range(1,7):
        T=tuple(tuple(int(i<j) for j in range(n)) for i in range(n));audit(T)
        direct=sum(all(T[p[i]][p[i+1]] for i in range(n-1)) for p in permutations(range(n)))
        need(h_forward(T)[0]==sum(h_backward_start(T,u)[0] for u in range(n))==direct==1,'transitive positive control')
    T=circulant(5,(1,2))
    direct=sum(all(T[p[i]][p[i+1]] for i in range(4)) for p in permutations(range(5)))
    need(h_forward(T)[0]==5*h_backward_start(T,0)[0]==direct,'nontrivial brute-force control')
    # A degree mismatch really is detected; the refinement is not vacuous.
    other=tuple(tuple(int(i<j) for j in range(5)) for i in range(5))
    try:stable_paired_partition(T,other)
    except ValueError:need(True,'hostile unequal tags rejected')
    else:need(False,'hostile accepted')


def run():
    controls();specs=((2,7,8,9,10,12),(1,2,4,7,8,10)); cases=[]; mats=[]
    for S in specs:
        A=circulant(13,S);mats.append(A)
        need(all(A[i][j]==A[(i+1)%13][(j+1)%13] for i in range(13) for j in range(13)), 'translation transitivity')
        need(all(A[(2*i)%13][(2*(i+1))%13] for i in range(13)), 'literal Hamilton cycle, hence strong')
        need(all(A[i][j]==A[(-j)%13][(-i)%13] for i in range(13) for j in range(13)), 'negation is converse isomorphism')
        vertex,ce,c3,s4=tags(A)
        need(set(vertex)=={(6,21,112)},'all actual vertex tags')
        need(c3==91 and s4==364,'global incidence counts')
        H,endcounts=h_forward(A);one,states=h_backward_start(A,0)
        need(H==13*one and len(set(endcounts))==1,'independent forward/all-start and backward/fixed-start H')
        autos,visits=stabilizer(A)
        need(autos==[tuple(range(13))],'complete rooted stabilizer is trivial')
        aut=13*len(autos)
        full,local=edge_profiles(A)
        need(len(set(local))==1,'vertexwise pair profile is uniform')
        moments=[sum(c for c in ce.values()),sum(c*(c-1)//2 for c in ce.values()),sum(c*(c-1)*(c-2)//6 for c in ce.values())]
        need(sum(count*p[2] for p,count in full)==moments[0],'oriented and unoriented edge incidence agree')
        cases.append({'n':13,'differences':S,'outgoing_rows':[sum(A[u][v]<<v for v in range(13)) for u in range(13)],
                      'uniform_vertex_tag':vertex[0],'cyclic_triangles':c3,'strong_four':s4,'H':H,'fixed_start_H':one,
                      'backward_states':states,'Aut':aut,'rooted_aut_search_nodes':visits,'ordinary_tilings':H//aut,
                      'self_converse':True,'edge_profile':full,'outedge_profile_per_vertex':local[0],
                      'edge_binomial_moments_orders_1_2_3':moments})
    need(cases[0]['H']==3683797 and cases[1]['H']==3704857,'reproduced discovered H values')
    need(cases[0]['H']!=cases[1]['H'],'realized nonisomorphism witness')
    need(cases[0]['edge_binomial_moments_orders_1_2_3']==[273,364,234] and cases[1]['edge_binomial_moments_orders_1_2_3']==[273,364,221], 'next overlap moment repairs pair')
    need(cases[0]['edge_profile']!=cases[1]['edge_profile'],'pair incidence sidecar separates')
    extensions=[]
    for k in range(9):
        A,B=[extend_sources(T,k) for T in mats]
        stable=stable_paired_partition(A,B)
        need(stable['class_sizes']==([1]*k+[13]),'source singletons and one unresolved core cell')
        for T in (A,B):
            v,ce,c3,s4=tags(T)
            need(v[k:]==((6,21,112),)*13,'core tags unchanged by source prefix')
            need(all(v[i]==(13+k-1-i,0,0) for i in range(k)),'prefix tags exact')
        if k<=2:
            need(h_forward(A)[0]==cases[0]['H'] and h_forward(B)[0]==cases[1]['H'],'finite direct extension-H control')
        extensions.append({'prefix_length':k,'n':13+k,**stable})
    return {'status':'PROVED all-depth incidence-refinement obstruction at every n>=13 via source extension; FINITE-EXACT base pair and bounded controls',
            'cases':cases,'extension_controls':extensions,'checks':CHECKS}


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true');args=parser.parse_args()
    result=run();raw=enc(result);target=HERE/'receipt.json'
    if args.write:target.write_bytes(raw)
    else:need(target.read_bytes().replace(b'\r\n',b'\n')==raw,'frozen exact receipt')
    print(json.dumps({'checks':CHECKS,'sha256_lf':sha256(raw).hexdigest(),'cases':[[r['H'],r['Aut'],r['ordinary_tilings'],r['edge_binomial_moments_orders_1_2_3']] for r in result['cases']]}))

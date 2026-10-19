#!/usr/bin/env python3
"""Standard-library fixed controls for five bounded-degree port-word seams.
No poset census; all-length conclusions require the accompanying proofs.
"""
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations, product
import json

ALPHA = [(), (6,), (8,), (6,8), (8,6)]
BETA = [(), (0,), (3,), (0,3), (3,0)]
GRAPHS = [
    ((6,0),(8,3)),
    ((6,3),(8,0)),
    ((8,0),(8,3)),
    ((6,0),(8,0),(8,3)),
    ((6,3),(8,0),(8,3)),
]
GRAPH_NAMES = ['matching_identity','matching_swap','star_at_8','L_at_0','L_at_3']
CORES = {
    'B9': {'pred':(0,0,1,0,8,24,3,5,189),'labels':tuple(range(9))},
    'B8': {'pred':(0,1,0,4,12,1,3,95),'labels':(0,2,3,4,5,6,7,8)},
}


def matmul(A,B):
    assert len(A[0]) == len(B)
    return [[sum(A[i][k]*B[k][j] for k in range(len(B)))
             for j in range(len(B[0]))] for i in range(len(A))]


def scalar(left,M,right):
    return matmul(matmul([left],M),[[x] for x in right])[0][0]


def shuffle_words(a,b):
    for positions in combinations(range(len(a)+len(b)),len(a)):
        selected=set(positions); ia=ib=0; w=[]
        for k in range(len(a)+len(b)):
            if k in selected:
                w.append(a[ia]); ia+=1
            else:
                w.append(b[ib]); ib+=1
        yield tuple(w)


def shuffle_kernel(graph):
    allowed=set(graph)
    return [[sum(all((u,v) in allowed or w.index(u)<w.index(v)
                     for u in a for v in b) for w in shuffle_words(a,b))
             for b in BETA] for a in ALPHA]


def explicit_extensions(p):
    def rec(s,w):
        if len(w)==len(p):
            yield w
        else:
            for x,ds in enumerate(p):
                if not(s>>x&1) and ds&s==ds:
                    yield from rec(s|(1<<x),w+(x,))
    return rec(0,())


def ecount(p,pair=None):
    p=list(p)
    if pair is not None:
        p[pair[1]] |= 1<<pair[0]
    full=(1<<len(p))-1
    @lru_cache(None)
    def rec(s):
        if s==full:return 1
        return sum(rec(s|(1<<x)) for x,ds in enumerate(p)
                   if not(s>>x&1) and ds&s==ds)
    return rec(0)


def skeleton_extensions(n):
    """Independent explicit 20-skeleton insertion construction, in original labels."""
    for six in shuffle_words((0,2,7),(3,4,5)):
        skeleton=six+(8,)
        if n==8:
            for gap in range(skeleton.index(0)+1,8):
                yield skeleton[:gap]+(6,)+skeleton[gap:]
        else:
            for q1,q6 in combinations(range(9),2):
                w=[];j=0
                for k in range(9):
                    if k==q1:w.append(1)
                    elif k==q6:w.append(6)
                    else:w.append(skeleton[j]);j+=1
                if w.index(0)<w.index(6):yield tuple(w)


def profiles(w):
    beta=();alpha=()
    for x in w:
        if x not in (0,3):break
        beta+=(x,)
    for x in reversed(w):
        if x not in (6,8):break
        alpha=(x,)+alpha
    assert alpha
    return BETA.index(beta), ALPHA.index(alpha)-1


def profile_tables(core):
    p,labels=core['pred'],core['labels']
    words={tuple(labels[x] for x in w) for w in explicit_extensions(p)}
    independent=list(skeleton_extensions(len(p)))
    assert len(independent)==len(set(independent)) and set(independent)==words
    N=[[0]*4 for _ in range(5)];C=[[0]*4 for _ in range(5)]
    for w in words:
        i,j=profiles(w)
        N[i][j]+=1;C[i][j]+=w.index(5)<w.index(7)
    pair=(labels.index(5),labels.index(7))
    assert sum(map(sum,N))==ecount(p)
    assert sum(map(sum,C))==ecount(p,pair)
    return N,C


def validate(p):
    for j,ds in enumerate(p):
        assert ds < (1<<j)
        for i in range(j):
            if ds>>i&1:assert p[i]&~ds==0
    return p


def assemble(names,graphs):
    p=[]; offsets=[]
    for name in names:
        offset=len(p);offsets.append(offset)
        p.extend((ds<<offset)|((1<<offset)-1) for ds in CORES[name]['pred'])
    for k,g in enumerate(graphs):
        A,B=CORES[names[k]],CORES[names[k+1]]
        for u,v in GRAPHS[g]:
            x=offsets[k]+A['labels'].index(u)
            y=offsets[k+1]+B['labels'].index(v)
            p[y]&=~(1<<x)
    return validate(tuple(p)),offsets


def graph_data(p):
    adj=[[] for _ in p];edges=[]
    for i,j in combinations(range(len(p)),2):
        if not(p[j]>>i&1):
            edges.append((i,j));adj[i].append(j);adj[j].append(i)
    timer=0;tin=[-1]*len(p);low=[0]*len(p);bridges=[]
    def dfs(v,parent):
        nonlocal timer
        tin[v]=low[v]=timer;timer+=1
        for w in adj[v]:
            if w==parent:continue
            if tin[w]>=0:low[v]=min(low[v],tin[w])
            else:
                dfs(w,v);low[v]=min(low[v],low[w])
                if low[w]>tin[v]:bridges.append(tuple(sorted((v,w))))
    dfs(0,-1)
    assert all(t>=0 for t in tin)
    return {'max_degree':max(map(len,adj)), 'bridges':sorted(bridges),
            'incomparability_edges':len(edges), 'max_reference_span':max(j-i for i,j in edges)}


def product_data(names,graphs,kernels):
    """Exact all-marked products via independent prefix and suffix propagation."""
    m=len(names);prefix=[[1]*5];after=[]
    for k,name in enumerate(names):
        row=matmul([prefix[k]],CORES[name]['N'])[0];after.append(row)
        if k+1<m:prefix.append(matmul([row],kernels[graphs[k]][1:])[0])
    suffix=[None]*m;suffix[-1]=[1]*4
    for k in range(m-2,-1,-1):
        right=matmul(CORES[names[k+1]]['N'],[[x] for x in suffix[k+1]])
        suffix[k]=[row[0] for row in matmul(kernels[graphs[k]][1:],right)]
    total=sum(after[-1]);marked=[]
    for k,name in enumerate(names):
        assert scalar(prefix[k],CORES[name]['N'],suffix[k])==total
        count=scalar(prefix[k],CORES[name]['C'],suffix[k]);marked.append(count)
        assert 24*total<=55*count<=31*total
    return total,marked


def run():
    K=[shuffle_kernel(g) for g in GRAPHS]
    left=sorted({tuple(row) for kernel in K for row in kernel})
    right=sorted({tuple(row[j] for row in kernel) for kernel in K for j in range(5)})
    assert len(left)==8 and len(right)==6
    tables={}
    for name,core in CORES.items():
        N,C=profile_tables(core);core['N']=N;core['C']=C
        corners=[]
        for i,l in enumerate(left):
            for j,r in enumerate(right):
                den=scalar(l,N,r[1:]);num=scalar(l,C,r[1:])
                lower=55*num-24*den;upper=31*den-55*num
                assert den>0 and lower>=0 and upper>=0
                corners.append({'left_generator':i,'right_generator':j,
                                'denominator':den,'numerator':num,'probability':str(F(num,den)),
                                'lower_deficit_55C_minus_24N':lower,
                                'upper_deficit_31N_minus_55C':upper})
        probabilities=[F(x['numerator'],x['denominator']) for x in corners]
        expected=(F(1375,3072),F(265,494)) if name=='B9' else (F(24,55),F(146,275))
        assert (min(probabilities),max(probabilities))==expected
        tables[name]={'prefix_by_nonempty_suffix_counts':N,'marked_counts':C,
                      'probability_range':list(map(str,expected)),
                      'minimum_lower_deficit':min(x['lower_deficit_55C_minus_24N'] for x in corners),
                      'minimum_upper_deficit':min(x['upper_deficit_31N_minus_55C'] for x in corners),
                      'all_48_generator_checks':corners}
    examples=[(('B9','B9'),(g,)) for g in range(5)]
    examples += [(('B8','B9'),(2,)),(('B8','B9','B8','B9'),(4,2,3)),
                 (('B9','B8','B9'),(3,4))]
    fixed=[]
    for names,graphs in examples:
        p,offsets=assemble(names,graphs);total,marked=product_data(names,graphs,K)
        assert total==ecount(p)
        for k,name in enumerate(names):
            labels=CORES[name]['labels']
            pair=tuple(offsets[k]+labels.index(x) for x in (5,7))
            assert marked[k]==ecount(p,pair)
        graph=graph_data(p);assert graph['max_degree']<=7
        guaranteed=names[-1]=='B9' and all(name!='B8' or graphs[k]!=2 for k,name in enumerate(names[:-1]))
        if guaranteed:assert graph['bridges']==[]
        if names==('B8','B9') and graphs==(2,):assert graph['bridges']==[(5,7)]
        fixed.append({'word':list(names),'seam_indices':list(graphs),'vertices':len(p),
                      'extensions':total,'marked_probabilities':[str(F(z,total)) for z in marked],
                      'graph':graph})
    # 48 blocks from the Thue-Morse sequence, ending with a B9 sentinel.
    # Seam choices are deterministic and mixed; B8's outgoing star is avoided.
    names=['B8' if i.bit_count()%2 else 'B9' for i in range(48)];names[-1]='B9'
    graphs=[]
    for i,name in enumerate(names[:-1]):
        g=(i.bit_count()+2*(i+1).bit_count())%5
        if name=='B8' and g==2:g=3
        graphs.append(g)
    total,marked=product_data(names,graphs,K);p,_=assemble(names,graphs);graph=graph_data(p)
    assert 200<=len(p)<=500 and graph['max_degree']<=7 and not graph['bridges']
    probabilities=[F(z,total) for z in marked]
    return {'status':'all fixed exact controls passed; not a census or an all-length computational proof',
            'alpha_profiles':ALPHA,'beta_profiles':BETA,'graph_names':GRAPH_NAMES,
            'deleted_comparison_graphs':GRAPHS,'full_5_by_5_shuffle_kernels':K,
            'left_row_generators':left,'right_column_generators':right,'cores':tables,
            'fixed_poset_controls':fixed,
            'long_product_control':{'blocks':len(names),'vertices':len(p),'word':names,
                'seam_indices':graphs,'extensions':total,
                'all_marked_probabilities':[str(z) for z in probabilities],
                'minimum_probability':str(min(probabilities)),'maximum_probability':str(max(probabilities)),
                'graph':graph,'full_poset_ideal_DP_performed':False}}


if __name__=='__main__':
    print(json.dumps(run(),indent=2))

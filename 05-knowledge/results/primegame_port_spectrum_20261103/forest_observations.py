"""Exact support bound and finite ranks for source-oblivious H observations."""
from itertools import combinations,product
from math import comb,factorial,gcd
from pathlib import Path
from fractions import Fraction
import argparse,json
from port_scout import add,rows
from growth_identity import sink
from verify import count

HERE=Path(__file__).resolve().parent

def linear_forest(n,edges):
    parent=list(range(n));degree=[0]*n
    def root(x):
        while parent[x]!=x:x=parent[x]
        return x
    for a,b in edges:
        degree[a]+=1;degree[b]+=1
        if degree[a]>2 or degree[b]>2:return False
        x,y=root(a),root(b)
        if x==y:return False
        parent[x]=y
    return True

def component_recurrence(limit):
    seq=[1]
    for n in range(1,limit+1):
        seq.append(sum(comb(n-1,k-1)*(1 if k==1 else factorial(k)//2)*seq[n-k] for k in range(1,n+1)))
    return seq

def endpoint_recurrence(limit):
    seq=[1]
    for n in range(1,limit+1):
        seq.append(sum(comb(n-1,k-1)*(1 if k==1 else comb(k,2))*seq[n-k] for k in range(1,n+1)))
    return seq

def endpoint_groups(n,edges):
    groups={}
    for mask in range(1<<len(edges)):
        chosen=[e for i,e in enumerate(edges) if mask>>i&1]
        if not linear_forest(n,chosen):continue
        neighbors=[[] for _ in range(n)]
        for a,b in chosen:neighbors[a].append(b);neighbors[b].append(a)
        used=set();blocks=[];sign=1
        for v in range(n):
            if v in used or len(neighbors[v])!=1:continue
            path=[v];previous=-1;current=v
            while True:
                onward=[w for w in neighbors[current] if w!=previous]
                if not onward:break
                previous,current=current,onward[0];path.append(current)
            used.update(path)
            if path[0]>path[-1]:path.reverse()
            for a,b in zip(path,path[1:]):sign*=1 if a<b else -1
            blocks.append((sum(1<<v for v in path),path[0],path[-1]))
        key=tuple(sorted(blocks))
        groups.setdefault(key,[]).append((mask,sign))
    return groups

def insert_mod(pivots,vector,p):
    a=[x%p for x in vector]
    for j in range(len(a)):
        if not a[j]:continue
        if j in pivots:
            b=pivots[j];q=a[j]
            a=[(x-q*y)%p for x,y in zip(a,b)]
        else:
            inv=pow(a[j],-1,p);pivots[j]=[(x*inv)%p for x in a]
            break

def rank_mod(rows_in,p):
    pivots={}
    for vector in rows_in:insert_mod(pivots,vector,p)
    return len(pivots)

def rank_integer(rows_in):
    pivots={}
    for vector in rows_in:
        a=list(vector)
        for j in range(len(a)):
            if not a[j]:continue
            if j not in pivots:
                pivots[j]=a
                break
            b=pivots[j];q=a[j];t=b[j]
            a=[t*x-q*y for x,y in zip(a,b)]
            common=0
            for x in a:common=gcd(common,x)
            if common>1:a=[x//common for x in a]
    return len(pivots)

def attach_cut(a,mask):
    n=len(a)
    return [row if mask>>v&1 else row|(1<<n) for v,row in enumerate(a)]+[mask]

def walsh(values):
    a=list(values);step=1
    while step<len(a):
        for start in range(0,len(a),2*step):
            for j in range(start,start+step):
                u,v=a[j],a[j+step];a[j]=u+v;a[j+step]=u-v
        step*=2
    return a

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',default='forest_observations.json');parser.add_argument('--five',action='store_true');args=parser.parse_args()
    seq=component_recurrence(12);blocks=endpoint_recurrence(12);censuses=[];rank_rows=[];closure=[];checks=0;path_law=None
    def need(ok,message):
        nonlocal checks
        checks+=1
        if not ok:raise RuntimeError(message)
    for n in range(1,7):
        es=list(combinations(range(n),2))
        forests=[mask for mask in range(1<<len(es)) if linear_forest(n,[e for i,e in enumerate(es) if mask>>i&1])]
        need(len(forests)==seq[n],'labelled linear-forest recurrence')
        groups=endpoint_groups(n,es)
        need(len(groups)==blocks[n],'endpoint-block recurrence')
        need(sum(map(len,groups.values()))==seq[n],'endpoint groups partition every forest')
        censuses.append({'n':n,'all_edge_sets':1<<len(es),'linear_forests':len(forests),'endpoint_blocks':len(groups)})
    for n in (3,4):
        old_edges=list(combinations(range(n),2));new_edges=list(combinations(range(n+1),2))
        old_groups=endpoint_groups(n,old_edges);new_groups=endpoint_groups(n+1,new_edges)
        valid={mask for group in old_groups.values() for mask,_ in group}
        forbidden=set(range(1<<len(old_edges)))-valid
        maps=[new_edges.index(e) for e in old_edges];tested=0
        for cut in range(1<<n):
            constant=sum(1<<new_edges.index((v,n)) for v in range(n) if not (cut>>v&1))
            extended=[constant|sum(((mask>>i)&1)<<j for i,j in enumerate(maps)) for mask in range(1<<len(old_edges))]
            for group in new_groups.values():
                # Use z=2a-1 to clear all powers of two in the basis values.
                degree=group[0][0].bit_count()
                values=[sum(sign*(-1)**(degree+(a&forest).bit_count()) for forest,sign in group) for a in extended]
                spectrum=walsh(values)
                ok=all(spectrum[s]==0 for s in forbidden)
                ok=ok and all(len({sign*spectrum[mask] for mask,sign in g})==1 for g in old_groups.values())
                need(ok,'every basis polynomial pulls back into the old endpoint space')
                tested+=1
        closure.append({'old_n':n,'new_basis_size':len(new_groups),'cuts':1<<n,'basis_pullbacks_verified':tested})
    for n,depth in [(3,3),(4,4)]+([(5,3)] if args.five else []):
        es=list(combinations(range(n),2));m=len(es)
        groups=endpoint_groups(n,es)
        forbidden=[mask for mask in range(1<<m) if not linear_forest(n,[e for i,e in enumerate(es) if mask>>i&1])]
        observations=[];by_depth=[]
        for k in range(depth+1):
            for word in product(range(n+1),repeat=k):
                values=[]
                for mask in range(1<<m):
                    a=rows(n,mask)
                    for op in word:a=sink(a) if op==n else add(a,op)
                    values.append(count(tuple(a)))
                spectrum=walsh(values)
                need(all(spectrum[s]==0 for s in forbidden),'no nonlinear-forest Fourier support')
                need(all(len({sign*spectrum[mask] for mask,sign in group})==1 for group in groups.values()),'internal path order is summed with its exact sign')
                if n==4 and path_law is None:
                    delta=[(-1)**((a&41).bit_count())+(-1)**((a&26).bit_count()) for a in range(64)]
                    need(sum(d*v for d,v in zip(delta,values))==0,'same-endpoint path-law annihilator')
                    correction=sum(d*v*v for d,v in zip(delta,values))
                    if correction:
                        moments=[[str(sum(Fraction(4+sg*d,256)*v**r for d,v in zip(delta,values))) for r in (1,2)] for sg in (1,-1)]
                        path_law={'old_path_masks':[41,26],'word':word,'moments':moments,'density':'mu_plus/minus(a)=(1 +/- (chi_41(a)+chi_26(a))/4)/64'}
                # Once the support/sign identities are checked, one Fourier
                # coefficient per endpoint block is a rank-preserving encoding.
                observations.append(values if n<5 else [sign*spectrum[mask] for mask,sign in (g[0] for g in groups.values())])
            ranks=[rank_mod(observations,p) for p in (101,103)]
            exact_rank=rank_integer(observations)
            need(all(r<=exact_rank<=blocks[n] for r in ranks),'exact rational rank and endpoint-block bounds')
            by_depth.append({'depth':k,'contexts':len(observations),'rank_Q':exact_rank,'ranks_mod_101_103':ranks})
        if n==5:
            pivots={}
            for vector in observations:insert_mod(pivots,vector,101)
            sampled=[]
            for index in range(6**4):
                if len(pivots)==blocks[n]:break
                encoded=(137*index)%(6**4);word=[]
                for _ in range(4):word.append(encoded%6);encoded//=6
                values=[]
                for mask in range(1<<m):
                    a=rows(n,mask)
                    for op in word:a=sink(a) if op==n else add(a,op)
                    values.append(count(tuple(a)))
                spectrum=walsh(values)
                need(all(spectrum[s]==0 for s in forbidden),'depth-four forest support')
                need(all(len({sign*spectrum[mask] for mask,sign in group})==1 for group in groups.values()),'depth-four endpoint blocks')
                vector=[sign*spectrum[mask] for mask,sign in (g[0] for g in groups.values())]
                observations.append(vector);insert_mod(pivots,vector,101);sampled.append(word)
            exact=rank_integer(observations);modular=[len(pivots),rank_mod(observations,103)]
            need(all(r<=exact<=blocks[n] for r in modular),'depth-four rational and modular rank')
            rank_rows.append({'n':n,'linear_forest_upper_bound':seq[n],'endpoint_block_upper_bound':blocks[n],'layers':by_depth,'rank_coordinates':'one verified Fourier coefficient per endpoint block','depth4_probe':{'selection':'base-six words from 137*i mod 6^4; stop at upper rank or exhaustion','tested_words':sampled,'rank_Q':exact,'ranks_mod_101_103':modular}})
            continue
        all_cuts=[]
        for first,second in product(range(1<<n),range(1<<(n+1))):
            vector=[]
            for mask in range(1<<m):
                a=attach_cut(attach_cut(rows(n,mask),first),second)
                vector.append(count(tuple(a)))
            spectrum=walsh(vector)
            need(all(spectrum[s]==0 for s in forbidden),'arbitrary two-birth forest support')
            need(all(len({sign*spectrum[mask] for mask,sign in group})==1 for group in groups.values()),'arbitrary births respect endpoint blocks')
            all_cuts.append(vector)
        cut_ranks=[rank_mod(all_cuts,p) for p in (101,103)]
        exact=rank_integer(all_cuts)
        need(all(r<=exact<=blocks[n] for r in cut_ranks),'two-birth rank check')
        rank_rows.append({'n':n,'linear_forest_upper_bound':seq[n],'endpoint_block_upper_bound':blocks[n],'layers':by_depth,'arbitrary_two_birth_contexts':len(all_cuts),'two_birth_rank_Q':exact,'two_birth_ranks_mod_101_103':cut_ranks})
    result={'checks':checks,'sequence_L0_to_L12':seq,'sequence_B0_to_B12':blocks,'complete_edge_set_censuses':censuses,'closed_basis_pullbacks':closure,'actual_growth_ranks':rank_rows,'hidden_path_law':path_law,'rank_scope':'Fraction-free integer elimination gives exact rational ranks; two prime-field computations are independent lower-bound controls. Attaining the all-context endpoint-block bound proves the all-context rank at that order.'}
    (HERE/args.output).write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(result))

if __name__=='__main__':main()

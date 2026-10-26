#!/usr/bin/env python3
"""Independent signed-component full closure certificate at orders eight to ten.

The producer's n10 diagnostic reused the earlier independent Gaussian.
This second implementation uses no Gaussian methods for n10 equations.
"""
from core import old_words,endpair_subcells,path,need,literal_source,literal_ledger
from itertools import combinations,product
from collections import Counter
import json

class Signs:
    def __init__(self,n):
        self.n=n
        self.pairs=list(combinations(range(n),2))
        self.index={p:i+1 for i,p in enumerate(self.pairs)}
        self.parent=list(range(len(self.pairs)+1));self.sign=[0]*len(self.parent)
        self.size=[1]*len(self.parent);self.ok=True
    def copy(self):
        t=Signs(self.n);t.parent=self.parent[:];t.sign=self.sign[:];t.size=self.size[:];t.ok=self.ok
        return t
    def root(self,v):
        if self.parent[v]!=v:
            p,d=self.root(self.parent[v]);self.sign[v]^=d;self.parent[v]=p
        return self.parent[v],self.sign[v]
    def equal(self,a,b,d):
        a,x=self.root(a);b,y=self.root(b);d^=x^y
        if a==b:
            if d:self.ok=False
            return
        if self.size[a]<self.size[b]:a,b=b,a
        self.parent[b]=a;self.sign[b]=d;self.size[a]+=self.size[b]
    def literal(self,a,b,flips=frozenset()):
        p=(min(a,b),max(a,b))
        return self.index[p],int(a>b)^int(p in flips)
    def force(self,a,b,value):
        v,d=self.literal(a,b);self.equal(v,0,value^d)
    def anti(self,q,flips=frozenset()):
        r=dict(zip(q,reversed(q)))
        for a,b in self.pairs:
            v,d=self.literal(a,b,flips);w,e=self.literal(r[b],r[a],flips)
            self.equal(v,w,d^e)
    def forced(self,a,b,flips):
        v,d=self.literal(a,b,flips);r,e=self.root(v);z,f=self.root(0)
        return r==z and (d^e^f)==1
    def dimension(self):
        need(self.ok,'dimension of inconsistent signed family')
        return len({self.root(v)[0] for v in self.parent})-1
    def singleton(self):
        need(self.dimension()==0,'singleton decoder only')
        z,c=self.root(0);rows=[0]*self.n
        for (a,b),v in self.index.items():
            r,d=self.root(v);need(r==z,'all variables forced')
            u,w=(a,b) if d^c else (b,a);rows[u]|=1<<w
        return tuple(rows)

def source(n,I,opened,chord):
    t=Signs(n);t.anti(list(range(n)))
    for v in range(n-1):t.force(v,v+1,1)
    t.force(n-1,0,1)
    for v in range(1,n-1):t.force(0,v,v not in I)
    t.force(1,n-2,opened);t.force(2,n-2,chord)
    need(t.ok,'nonempty raw source')
    return t

def models(t):
    need(t.ok,'enumeration requires consistent signed family')
    z,c=t.root(0)
    free=sorted({t.root(v)[0] for v in t.parent}-{z})
    for choice in range(1<<len(free)):
        values={z:c}|{r:(choice>>j)&1 for j,r in enumerate(free)}
        rows=[0]*t.n
        for (a,b),v in t.index.items():
            r,d=t.root(v);u,w=(a,b) if values[r]^d else (b,a)
            rows[u]|=1<<w
        yield tuple(rows)

def clone_partition(rows):
    """Exhaust every vertex subset; do not infer maximal clones from weights.

    A clone block is a transitive module. Return its actual internal order.
    Explicit maximality, disjointness, and full arc reconstruction are checked.
    """
    from core import arc
    n=len(rows);modules=[];orders={}
    for mask in range(1,1<<n):
        block=[v for v in range(n) if mask>>v&1]
        outside=[v for v in range(n) if not(mask>>v&1)]
        if any(any(arc(rows,u,v)!=arc(rows,block[0],v) for u in block[1:]) for v in outside):continue
        order=sorted(block,key=lambda u:-sum(arc(rows,u,v) for v in block if v!=u))
        if any(not arc(rows,u,v) for i,u in enumerate(order) for v in order[i+1:]):continue
        modules.append(mask);orders[mask]=order
    maximal=[mask for mask in modules if not any(mask!=other and mask&other==mask for other in modules)]
    need(sum(mask.bit_count() for mask in maximal)==n and sum(maximal)==(1<<n)-1,
         'maximal transitive modules partition every literal vertex')
    need(all(not(a&b) for a,b in combinations(maximal,2)),'maximal clone blocks disjoint')
    blocks=sorted((orders[mask] for mask in maximal),key=min)
    owner={v:(i,j) for i,block in enumerate(blocks) for j,v in enumerate(block)}
    for u,v in combinations(range(n),2):
        i,a=owner[u];j,b=owner[v]
        reconstructed=int(a<b) if i==j else arc(rows,blocks[i][0],blocks[j][0])
        need(arc(rows,u,v)==reconstructed,'exact arcs reconstruct from ordered blocks and literal quotient')
    return blocks

def main():
    from pathlib import Path
    from hashlib import sha256
    from itertools import permutations
    from core import arc,move,mismatch,twins
    data_bytes=Path(__file__).with_name('found.json').read_bytes()
    data=json.loads(data_bytes);table={tuple(item['source_rows']):item for item in data}
    need(len(data)==len(table)==30,'thirty distinct lower-order proposed sources')
    totals={};failures=[];replay=[];exhausted=set();partitions=[]
    for n in (8,9,10):
        counts=Counter();dimensions=Counter()
        for ports in combinations(range(2,n-1),3):
            I=set(ports)
            for opened,chord in product((False,True),repeat=2):
                s=source(n,I,opened,chord);dimensions[s.dimension()]+=1;counts['all_source_cells']+=1
                if not opened and {n-3,n-2}<=I:
                    h=next(iter(I-{n-3,n-2}));leaves=endpair_subcells(s,h,chord);counts['endpair_cells']+=1
                    kind='endpair'
                else:
                    label,qs=old_words(n,I,opened,chord);leaves=[(label,[],qs)];counts['base_cells']+=1
                    kind='base'
                refined=[]
                for label,guards,qs in leaves:
                    t=s.copy()
                    for a,b,value in guards:t.force(a,b,value)
                    need(t.ok,'nonempty refined source leaf')
                    refined.append({'label':label,'guards':guards,'dimension':t.dimension(),'words':qs})
                    counts['leaves']+=1;flips=set();blue=[];illegal=[]
                    for step,q in enumerate(qs,1):
                        need(sorted(q)==list(range(n)),'every old word is a complete permutation')
                        for a,b in zip(q,q[1:]):
                            if not t.forced(a,b,flips):
                                need(t.forced(b,a,flips),'failed join is forced false, not merely undetermined')
                                illegal.append({'step':step,'ordered_join':[a,b],'forced_arc_value':0})
                        trial=t.copy();trial.anti(q,flips)
                        if trial.ok:blue.append((step,trial))
                        flips^=set(t.pairs)^path(q)
                    need(flips==set(t.pairs)^path(range(n)),'every old formal ledger has exact symbolic target')
                    if illegal:
                        # Legality failure is not a blue-singleton condition.
                        # Cover the ENTIRE exact source leaf, without conditioning.
                        needed=set(models(t));reason='whole illegal source leaf'
                        need(n==8 and label=='h4_backward' and t.dimension()==3 and len(needed)==8,
                             'exact whole-leaf illegal boundary')
                        counts['illegal_leaves']+=1
                    else:
                        for step,trial in blue:
                            need(trial.dimension()==0,'all legal old-blue loci are singleton')
                        needed=set().union(*(set(models(trial)) for _,trial in blue))
                        reason='complete old-blue locus union'
                    if not needed:
                        counts['uniform_success_leaves']+=1
                        continue
                    counts['failed_leaves']+=1
                    failures.append({'n':n,'ports':ports,'middle_open':opened,'x_to_c':chord,
                        'kind':kind,'old_label':label,'guards':guards,'source_dimension':t.dimension(),
                        'illegal_joins':illegal,'blue_loci':[{'step':step,'dimension':q.dimension()} for step,q in blue],
                        'repair_universe_reason':reason,'exact_sources_requiring_replacement':len(needed)})
                    for rows in sorted(needed):
                        need(rows not in exhausted,'disjoint literal replacement source universe')
                        need(rows in table,('missing exact repair certificate',n,ports,rows))
                        item=table[rows]
                        need(item['n']==n and tuple(item['ports'])==ports
                             and item['middle_open']==opened and item['x_to_c']==chord,
                             'certificate exact source labels')
                        literal_source(rows,I,opened,chord)
                        need(all(arc(rows,a,b)==value for a,b,value in guards),'literal complementary guards')
                        if not illegal:
                            old_flags,_,_=literal_ledger(rows,qs,False)
                            need(any(old_flags),'literal old-blue hostile is real')
                        else:
                            formal=rows;old_flags=[];legal=[]
                            for q in qs:
                                legal.append(all(arc(formal,a,b) for a,b in zip(q,q[1:])))
                                old_flags.append(mismatch(formal,q) is None)
                                formal=move(formal,q)
                            need(legal==[False,False,True],'actual forced illegal boundary steps')
                        flags,trace,target=literal_ledger(rows,item['words'],True)
                        need(flags==[False]*3,'all three proposed steps are actually black')
                        for step in trace:
                            step['ordered_maximal_transitive_clone_blocks']=clone_partition(step['source_rows'])
                        replay.append({'n':n,'ports':ports,'source_rows':rows,'repair_reason':reason,
                            'old_words':qs,'old_formal_blue_flags':old_flags,'new_steps':trace,
                            'literal_target_rows':target,'target_twin_pairs':twins(target),
                            'target_ordered_maximal_transitive_clone_blocks':clone_partition(target)})
                        exhausted.add(rows);counts['literal_sources_repaired']+=1
                need(sum(1<<leaf['dimension'] for leaf in refined)==1<<s.dimension(),'full partition cardinality')
                for first,second in combinations(refined,2):
                    t=s.copy()
                    for a,b,value in first['guards']+second['guards']:t.force(a,b,value)
                    need(not t.ok,'pairwise disjoint leaf guards')
                if kind=='endpair':
                    partitions.append({'n':n,'ports':ports,'x_to_c':chord,'source_dimension':s.dimension(),'leaves':refined})
        expected={8:(40,34,6,47,45,2,16,4),9:(80,72,8,89,84,5,5,7),10:(140,130,10,151,142,9,9,10)}[n]
        actual=tuple(counts[k] for k in ('all_source_cells','base_cells','endpair_cells','leaves','uniform_success_leaves','failed_leaves','literal_sources_repaired'))
        need(actual==expected[:-1],('complete lower-order source partition',n,actual))
        need(dimensions=={expected[-1]:expected[0]},('unfiltered raw source dimension',n))
        totals[n]={'counts':dict(counts),'raw_dimensions':dict(dimensions),'literal_marked_sources':expected[0]*(1<<expected[-1])}
    need(exhausted==set(table),'supplied thirty sources equal exhaustive required replacement universe')
    hostile=(6,4,8,19,103,79,15);literal_source(hostile,{3,4,5},False,False)
    target=move(hostile,list(range(7)))
    HP=[q for q in permutations(range(7)) if all(arc(target,a,b) for a,b in zip(q,q[1:]))]
    need(len(HP)==9 and all(mismatch(target,q) is None for q in HP),'sharp order-seven pure-blue obstruction')
    print(json.dumps({'status':'PASS','method':'independent signed-component raw source and anti equations, complete illegal leaves and singleton-blue unions, literal arc replay',
       'input_sha256':sha256(data_bytes).hexdigest(),'orders':totals,'failed_old_leaves':failures,
       'endpair_partitions':partitions,'literal_repairs':replay,'total_repaired_sources':len(exhausted),
       'order7_sharpness_hostile':{'source':hostile,'target':target,'permutations_tested':5040,'directed_paths':9,'black_paths':0},
       'scope':'FINITE-EXACT complete closed marked endpoint-four closure at n8,n9,n10. Combined with sealed n>=11 theorem gives n>=8; order7 failure makes8 the optimal integer all-order threshold for this unfiltered local replacement statement. No whole black-core connectivity conclusion.'},indent=2,sort_keys=True))

if __name__=='__main__':main()


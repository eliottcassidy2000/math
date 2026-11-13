#!/usr/bin/env python3
"""Colored ancestor coverage and exact boundary-carry port selection."""
from dataclasses import dataclass, field
from collections import defaultdict
from itertools import product
from pathlib import Path
import json

def require(ok,message):
    if not ok: raise RuntimeError(message)

@dataclass
class Tree:
    children: tuple=()
    index: int=-1
    color: int=-1
    mask: int=0
    cls: tuple=()
    active: int=0
    weight: int=0

def make_tree(colors):
    n=len(colors)
    if n==1:
        i,c=colors[0]
        return Tree(index=i,color=c,mask=1<<i,cls=('L',c))
    require(n%3==0,'not a ternary level')
    m=n//3
    children=tuple(make_tree(colors[k*m:(k+1)*m]) for k in range(3))
    classes=tuple(x.cls for x in children)
    canonical=min(classes[k:]+classes[:k] for k in range(3))
    a=int(classes[0]==classes[1]==classes[2])
    return Tree(children=children,mask=sum(x.mask for x in children),
                cls=('N',canonical),active=a,weight=a)

def nodes(t):
    if t.children:
        yield t
        for c in t.children: yield from nodes(c)

def leaves(mask):
    while mask:
        bit=mask&-mask; mask-=bit; yield bit.bit_length()-1

def cover(t,selected):
    return sum(v.weight for v in nodes(t) if v.mask&selected)

def greedy_layers(t,eligible):
    reachable={0}; result=[reachable]
    for k in range(eligible.bit_count()):
        nxt=set()
        for S in reachable:
            choices=[(cover(t,S|1<<x)-cover(t,S),x) for x in leaves(eligible&~S)]
            best=max(g for g,_ in choices)
            nxt.update(S|1<<x for g,x in choices if g==best)
        result.append(nxt); reachable=nxt
    return result

@dataclass
class Table:
    kind: str
    mask: int
    counts: dict=field(default_factory=dict)
    children: tuple=()
    branches: dict=field(default_factory=dict)
    weight: int=0
    index: int=-1
    role: int=-1
    roles: int=0

    def state_of(self,S):
        S &= self.mask
        if self.kind=='leaf':
            k=int(bool(S)); q=[0]*self.roles
            if k and self.role>=0: q[self.role]=1
            return (k,k,k,0,tuple(q))
        if self.kind=='charge':
            k,l,r,c,q=self.children[0].state_of(S)
            return (k,l,r,c+self.weight*bool(k),q)
        a,b=(x.state_of(S) for x in self.children)
        return (a[0]+b[0],a[1],b[2],a[3]+b[3],tuple(x+y for x,y in zip(a[4],b[4])))

    def unrank(self,key,rank):
        require(0<=rank<self.counts[key],'rank out of range')
        if self.kind=='leaf': return self.mask if key[0] else 0
        if self.kind=='charge':
            k,l,r,c,q=key
            return self.children[0].unrank((k,l,r,c-self.weight*bool(k),q),rank)
        L,R=self.children
        for a,b in self.branches[key]:
            block=L.counts[a]*R.counts[b]
            if rank<block:
                i,j=divmod(rank,R.counts[b])
                return L.unrank(a,i)|R.unrank(b,j)
            rank-=block
        raise RuntimeError('unrank fell through')

    def rank(self,S):
        key=self.state_of(S); require(key in self.counts,'infeasible selected set')
        if self.kind=='leaf': return 0
        if self.kind=='charge': return self.children[0].rank(S)
        L,R=self.children; a,b=L.state_of(S),R.state_of(S); out=0
        for x,y in self.branches[key]:
            if (a,b)==(x,y): return out+L.rank(S)*R.counts[b]+R.rank(S)
            out+=L.counts[x]*R.counts[y]
        raise RuntimeError('rank fell through')

def merge(L,R,K,caps,hardcore):
    out=Table('merge',L.mask|R.mask,children=(L,R))
    out.counts=defaultdict(int); out.branches=defaultdict(list)
    for a in sorted(L.counts):
        for b in sorted(R.counts):
            if a[0]+b[0]>K or (hardcore and a[2] and b[1]): continue
            q=tuple(x+y for x,y in zip(a[4],b[4]))
            if any(x>y for x,y in zip(q,caps)): continue
            key=(a[0]+b[0],a[1],b[2],a[3]+b[3],q)
            out.counts[key]+=L.counts[a]*R.counts[b]
            out.branches[key].append((a,b))
    return out

def build_table(t,eligible,K,roles=None,caps=(),hardcore=False):
    roles={} if roles is None else roles
    if not t.children:
        q=(0,)*len(caps); counts={(0,0,0,0,q):1}
        role=roles.get(t.index,-1)
        if eligible&t.mask and K:
            z=list(q)
            if role>=0: z[role]=1
            if all(x<=y for x,y in zip(z,caps)): counts[(1,1,1,0,tuple(z))]=1
        return Table('leaf',t.mask,counts,index=t.index,role=role,roles=len(caps))
    tables=[build_table(c,eligible,K,roles,caps,hardcore) for c in t.children]
    merged=tables[0]
    for c in tables[1:]: merged=merge(merged,c,K,caps,hardcore)
    out=Table('charge',t.mask,children=(merged,),weight=t.weight)
    for (k,l,r,c,q),n in merged.counts.items(): out.counts[(k,l,r,c+t.weight*bool(k),q)]=n
    return out

def query(table,k,covered,quota=None):
    return [(key,n) for key,n in sorted(table.counts.items())
            if key[0]==k and key[3]==covered and (quota is None or key[4]==quota)]

def select(table,k,covered,rank=0,quota=None):
    for key,n in query(table,k,covered,quota):
        if rank<n: return table.unrank(key,rank)
        rank-=n
    raise ValueError('requested target/rank has no physical marking')

def valid(S,eligible,hardcore):
    return S&~eligible==0 and (not hardcore or not(S&(S<<1)))

def check_table(t,eligible,K,roles=None,caps=(),hardcore=False):
    roles={} if roles is None else roles
    table=build_table(t,eligible,K,roles,caps,hardcore)
    direct=defaultdict(set)
    for S in range(t.mask+1):
        if S.bit_count()>K or not valid(S,eligible,hardcore): continue
        key=table.state_of(S)
        if any(x>y for x,y in zip(key[4],caps)): continue
        direct[key].add(S)
    require(set(table.counts)==set(direct),'DP state universe mismatch')
    for key,sets in direct.items():
        require(table.counts[key]==len(sets),'DP count mismatch')
        decoded={table.unrank(key,i) for i in range(len(sets))}
        require(decoded==sets,'DP lift universe mismatch')
        for i in range(len(sets)): require(table.rank(table.unrank(key,i))==i,'rank inverse')
    return table

def group_permutations(depth):
    if depth==0: return [(0,)]
    child=group_permutations(depth-1); n=3**(depth-1); out=[]
    for shift in range(3):
        for a,b,c in product(child,repeat=3):
            maps=(a,b,c)
            out.append(tuple(((i//n+shift)%3)*n+maps[i//n][i%n] for i in range(3*n)))
    return out

def max_tree_weight(t,S):
    vertices=list(leaves(S)); parent={v:v for v in vertices}
    def root(v):
        while parent[v]!=v:
            parent[v]=parent[parent[v]]; v=parent[v]
        return v
    edges=[]
    for i,x in enumerate(vertices):
        for y in vertices[i+1:]:
            w=sum(v.weight for v in nodes(t) if v.mask>>x&1 and v.mask>>y&1)
            edges.append((w,x,y))
    total=0
    for w,x,y in sorted(edges,reverse=True):
        a,b=root(x),root(y)
        if a!=b: parent[a]=b; total+=w
    return total

def literal_graph(t):
    if not t.children:
        return ((2,80,3,7,173,143,189,15) if t.color==0 else
                (130,156,153,209,225,207,135,0))
    blocks=[literal_graph(c) for c in t.children]; offsets=[0,len(blocks[0]),len(blocks[0])+len(blocks[1])]
    rows=[]
    for i,R in enumerate(blocks):
        other=(i+1)%3; out=((1<<len(blocks[other]))-1)<<offsets[other]
        rows.extend((v<<offsets[i])|out for v in R)
    return tuple(rows)

def root_child(rows,r):
    n=len(rows); out=list(rows); out[r]|=1<<n; out.append(((1<<n)-1)^(1<<r))
    return tuple(out)

def is_auto(rows,p):
    return all(((rows[i]>>j)&1)==((rows[p[i]]>>p[j])&1)
               for i in range(len(rows)) for j in range(len(rows)))

def strong(rows):
    full=(1<<len(rows))-1
    for start in range(len(rows)):
        seen=0; todo=1<<start
        while todo:
            b=todo&-todo; todo-=b; seen|=b
            todo|=rows[b.bit_length()-1]&~seen
        if seen!=full:return False
    return True

def main():
    greedy_states=0
    for colors in product((0,1),repeat=9):
        t=make_tree(list(enumerate(colors))); optimum=[-1]*10
        for S in range(512): optimum[S.bit_count()]=max(optimum[S.bit_count()],cover(t,S))
        for k,layer in enumerate(greedy_layers(t,511)):
            for S in layer: require(cover(t,S)==optimum[k],'colored greedy failure')
            greedy_states+=len(layer)
    weighted_cases=0; mst_controls=0
    for weights in product(range(3),repeat=4):
        t=make_tree(list(enumerate((0,)*9)))
        for v,w in zip(nodes(t),weights):v.weight=w
        eligible=0b101010101; optimum=[-1]*6
        for S in range(512):
            if S&~eligible:continue
            optimum[S.bit_count()]=max(optimum[S.bit_count()],cover(t,S))
            if S:
                require(cover(t,S)==sum(cover(t,1<<x) for x in leaves(S))-max_tree_weight(t,S),
                        'weighted maximum spanning tree identity')
                mst_controls+=1
        for k,layer in enumerate(greedy_layers(t,eligible)):
            require(all(cover(t,S)==optimum[k] for S in layer),'weighted eligible greedy')
        weighted_cases+=1
    fixture=(0,0,0,1,1,1,0,0,1)
    t=make_tree(list(enumerate(fixture))); E=sum(v.active for v in nodes(t))
    group=[p for p in group_permutations(2) if all(fixture[i]==fixture[p[i]] for i in range(9))]
    for S in range(512):
        actual=sum(all(p[i]==i for i in leaves(S)) for p in group)
        require(actual==3**(E-cover(t,S)),'point stabilizer formula')
    tables=0
    for eligible in (511,0b111101111,0b100010101,0):
        for hardcore in (False,True):
            check_table(t,eligible,9,hardcore=hardcore); tables+=1
    # Literal roles from row9 of the authenticated three-color substitution.
    word='C'
    while len(word)<9: word=''.join({'C':'CA','A':'B','B':'CA'}[x] for x in word)
    row=word[:9][::-1]; role={i:'CAB'.index(c) for i,c in enumerate(row)}
    table=check_table(t,511,3,role,(1,1,1),True); tables+=1
    answers=query(table,3,E,(1,1,1)); total=sum(n for _,n in answers)
    require(total>0,'colored row target unexpectedly empty')
    chosen=select(table,3,E,0,(1,1,1))
    require(cover(t,1<<0|1<<6)==1 and cover(t,1<<0|1<<3)==2,'raw-distance hostile')
    tiny=make_tree(list(enumerate((0,0,0))))
    gap=check_table(tiny,0b101,2,hardcore=True)
    require(select(gap,2,1)==0b101,'forced-zero gap was reconnected')
    # Exact-target gap with two eligible leaves in a depth-three colored hierarchy.
    rigid=(0,0,1,0,1,1,0,0,1)
    big=make_tree(list(enumerate((0,)*9+rigid+rigid)))
    eligible=(1<<0)|(1<<9); bt=build_table(big,eligible,1)
    covered={key[3] for key in bt.counts if key[0]==1}
    require(covered=={0,2},'constrained exponent gap changed')
    physical=[]; old=literal_graph(t); rotations=group_permutations(2)
    for S in (1, (1<<0)|(1<<6), (1<<0)|(1<<3), chosen):
        rows=old
        for x in leaves(S): rows=root_child(rows,8*x)
        require(all((rows[i]&((1<<len(old))-1))==old[i] for i in range(len(old))),
                'physical birth changed an old arc')
        require(strong(rows),'physical extension lost strongness')
        aut=0
        for p in rotations:
            lifted=tuple(8*p[i//8]+i%8 for i in range(len(old)))+tuple(range(len(old),len(rows)))
            aut+=is_auto(rows,lifted)
        expected=3**(E-cover(t,S))
        require(aut==expected,'literal birth stabilizer control')
        physical.append({'selected':list(leaves(S)),'old_order':len(old),'new_order':len(rows),
                         'residual_group_order':expected,'old_arcs_preserved':True})
    receipt={'status':'FINITE-EXACT controls; the all-tree proofs are separate',
             'complete_colored_depth2_greedy':{'colorings':512,'all_tie_prefix_states':greedy_states},
             'weighted_eligibility_controls':{'weight_vectors':weighted_cases,'maximum_tree_checks':mst_controls},
             'independent_group_control':{'ambient_group':81,'color_stabilizer':len(group),'marking_masks':512},
             'DP_tables_against_full_mask_enumeration':tables,
             'raw_distance_hostile':{'colors':fixture,'first':0,'bad_second':6,'good_second':3,
                                     'remaining_group_orders':[3,1]},
             'three_color_row':{'row':row,'quota':[1,1,1],'marks':3,'target_group_order':1,
                                'valid_markings':total,'rank0':list(leaves(chosen))},
             'forced_zero_gap':{'eligible':[0,2],'selected':[0,2]},
             'exact_exponent_gap':{'initial_exponent':4,'eligible':[0,9],
                                    'budget':1,'attainable_residual_exponents':[2,4]},
             'literal_birth_controls':physical,
             'scope':'Selections control point stabilizers; root-child births provide the stated physical lift. No fixed Hamilton target claim.',
             'reproduce':'python3 -B 05-knowledge/results/zeckendorf_target_selection_20261113/charts/selector.py'}
    Path(__file__).with_name('receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt,indent=2))

if __name__=='__main__': main()

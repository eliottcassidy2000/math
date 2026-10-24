#!/usr/bin/env python3
"""Exact independently recomputed controls. Standard library only, no asserts."""
from collections import Counter, defaultdict
from functools import lru_cache
from itertools import combinations, permutations, product
from pathlib import Path
import argparse, hashlib, json

HERE=Path(__file__).resolve().parent
CHECKS=0

def need(ok,why):
    global CHECKS
    CHECKS+=1
    if not ok: raise RuntimeError(why)

def powerset(xs):
    xs=tuple(xs)
    return [frozenset(c) for k in range(len(xs)+1) for c in combinations(xs,k)]

def allocated(families):
    return [t for t in product(*families) if sum(map(len,t))==len(frozenset().union(*t))]

def occupancy(families):
    """Nonnegative commutative square-zero product, coefficients count tuples."""
    p={frozenset():1}
    for family in families:
        q=defaultdict(int)
        for used,n in p.items():
            for edge in family:
                if not used&edge:q[used|edge]+=n
        p=dict(q)
    return p

def hall(lists):
    return all(len(set().union(*(lists[i] for i in J)))>=len(J) for J in powerset(range(len(lists))))

@lru_cache(None)
def packing(family):
    return max((len(s) for s in powerset(family) if sum(map(len,s))==len(frozenset().union(*s))),default=0)

def greedy_certificate(families,universe,kind):
    K=0
    for f in families:
        if not f:return False
        if kind=='exact':
            good=all(any(not e&u for e in f) for u in powerset(universe) if len(u)<=K)
        elif kind=='pack':good=packing(f)>K
        elif kind=='degrees':
            dd=sorted((sum(u in e for e in f) for u in universe),reverse=True)
            good=len(f)>sum(dd[:K])
        else:raise ValueError(kind)
        if not good:return False
        K+=max(map(len,f))
    return True

def all_greedy_prefixes_extend(families):
    prefixes=[frozenset()]
    for f in families:
        nxt=[]
        for used in prefixes:
            extensions=[used|e for e in f if not used&e]
            if not extensions:return False
            nxt+=extensions
        prefixes=nxt
    return True

def hall_checks():
    cases=0;positive=0
    for m in range(1,4):
        for n in range(1,5):
            lists=powerset(range(n))
            for menus in product(lists,repeat=m):
                fs=[[frozenset([u]) for u in menu] for menu in menus]
                literal=allocated(fs);alg=occupancy(fs)
                need(bool(literal)==hall(menus),'Hall iff actual singleton allocation')
                need(sum(alg.values())==len(literal),'singleton occupancy multiplicities')
                cases+=1;positive+=bool(literal)
    # Distinct records over one resource are not three resources.
    full_lifts=[('a',k) for k in range(3)]
    need(len(allocated([[frozenset([a]) for a,k in full_lifts]]*2))==0,'kernel lifts cannot supply second capacity unit')
    need(hall([set(full_lifts),set(full_lifts)]),'hostile raw-lift Hall falsely passes')
    need(not hall([{'a'},{'a'}]),'correct normalized Hall fails')
    p=occupancy([[frozenset('a'),frozenset('b')]]*2)
    need(p=={frozenset('ab'):2},'positive coefficient 2 survives rather than exterior cancellation')
    return dict(bipartite_catalogue_cases=cases,feasible_cases=positive,range='demands 1..3, resources 1..4, all neighbor-list tuples',same_two_lists_coefficient=2)

def hypergraph_checks():
    U=frozenset(range(4));edges=tuple(frozenset(e) for e in combinations(U,2));allfamilies=[tuple(sorted(f,key=lambda x:tuple(sorted(x)))) for f in powerset(edges)]
    cases=0;feasible=0;union_pass_false=0;certs=Counter()
    for fs in product(allfamilies,repeat=2):
        a=allocated(fs);p=occupancy(fs)
        need(len(a)==sum(p.values()),'all two-colour rank-two occupancy products')
        if a:feasible+=1
        union_pass=all(packing(tuple(set().union(*(set(fs[i]) for i in J))))>=len(J) for J in powerset(range(2)))
        if union_pass and not a:union_pass_false+=1
        need(not a or union_pass,'uncoloured matching necessary only')
        for kind in ('exact','pack','degrees'):
            cert=greedy_certificate(fs,U,kind)
            if cert:
                certs[kind]+=1;need(all_greedy_prefixes_extend(fs),'depletion certificate extends every greedy prefix')
        cases+=1
    small=[f for f in allfamilies if len(f)<=2]
    triples=0
    for fs in product(small,repeat=3):
        p=occupancy(fs);a=allocated(fs)
        need(sum(p.values())==len(a),'three-colour rank-two allocation exactness')
        need(not a,'three 2-resource demands cannot fit in four resources')
        triples+=1
    A=(frozenset((0,1)),frozenset((2,3)));B=(frozenset((0,2)),frozenset((1,3)))
    need(packing(A)==packing(B)==2,'sharp hostile has two internally disjoint options each')
    need(not allocated([A,B]),'four-edge rainbow hostile')
    need(all(len(set().union(*(set().union(*([A,B][i])) for i in J)))>=2*len(J) for J in powerset(range(2))),'weighted resource neighborhood Hall passes hostile')
    B_one=B+(frozenset((0,4)),);B_pack=B+(frozenset((4,5)),)
    need(len(allocated([A,B_one]))==1,'one fresh native key repairs one branch')
    need(len(allocated([A,B_pack]))==2 and packing(B_pack)==3,'three-pack repair extends either first choice')
    need(greedy_certificate([A,B_pack],range(6),'pack'),'strict pack threshold positive')
    need(greedy_certificate([A,B_pack],range(6),'degrees'),'resource-degree criterion positive')
    independent=tuple(frozenset((2*i,2*i+1)) for i in range(5))
    for kind in ('exact','pack','degrees'):
        need(greedy_certificate([independent]*3,range(10),kind),'three-demand depletion positive '+kind)
        need(all_greedy_prefixes_extend([independent]*3),'all three-demand positive greedy prefixes extend')
    need(not greedy_certificate([A,B_one],range(5),'exact'),'all-prefix criterion is not necessary for some allocation')
    # Empty resource sets really consume no capacity in the abstract compiler.
    need(len(allocated([tuple([frozenset()])]*4))==1,'zero-resource witnesses abstractly valid')
    need(greedy_certificate([tuple([frozenset()])]*4,range(0),'exact'),'empty-edge exact depletion criterion')
    return dict(two_demand_catalogue_cases=cases,two_demand_feasible=feasible,uncoloured_matching_false_positives=union_pass_false,depletion_certificates=dict(certs),three_demand_catalogue_cases=triples,range='all subsets of six 2-edges on four resources; triples restricted to at most two options per demand',hostile=[[[0,1],[2,3]],[[0,2],[1,3]]])


def charged_product(menus,modulus):
    """Each tuple is (tag, full resource occurrence tuple, charge, grade)."""
    dp={(frozenset(),0,0,()):1}
    for menu in menus:
        nxt=defaultdict(int)
        for (used,q,g,tags),mult in dp.items():
            for tag,resources,charge,grade in menu:
                rr=frozenset(resources)
                if len(rr)!=len(resources) or used&rr:continue
                nxt[used|rr,(q+charge)%modulus,g+grade,tags+(tag,)]+=mult
        dp=dict(nxt)
    return dp

def charge_relation_checks():
    menu=[(f'{u}:{k}',(u,), (u+k)%3,k) for u in range(3) for k in range(3)]
    menu.append(('internal-repeat',(0,0),0,0))
    rows=charged_product([menu,menu],3)
    brute=Counter()
    for a,b in product(menu,repeat=2):
        if len(set(a[1]+b[1]))!=len(a[1])+len(b[1]):continue
        brute[frozenset(a[1]+b[1]),(a[2]+b[2])%3,a[3]+b[3],(a[0],b[0])]+=1
    need(rows==dict(brute),'charge grade full witness tags exact')
    for q,g in product(range(3),range(5)):
        need(sum(v for (u,c,h,t),v in rows.items() if (c,h)==(q,g))==sum(1 for a,b in product(menu,repeat=2) if len(set(a[1]+b[1]))==len(a[1])+len(b[1]) and ((a[2]+b[2])%3,a[3]+b[3])==(q,g)),'coefficient extraction exact fixed charge and grade')
    # Private resources: only a higher-arity factor empties the catalogue.
    priv=[[(k,(i,),0,0) for k in range(3)] for i in range(3)]
    full=charged_product(priv,3)
    need(sum(full.values())==27,'all private-resource phase tuples preserved before relation')
    survive=[z for z in full if z[3][0]==z[3][1] and z[3][1]==z[3][2] and z[3][2]==(z[3][0]+1)%3]
    need(not survive,'frustrated relation cannot be erased by occupancy projection')
    # A single representative per resource and charge can discard the only compatible lifts.
    keep=[z for z in charged_product(priv[:2],3) if z[3]==(1,1)]
    need(len(keep)==1,'retained full lift pair survives declared relation')
    need((0,0)!=(1,1),'arbitrary phase-zero representative loses compatible pair')
    return dict(charged_full_assignments=sum(rows.values()),charge_modulus=3,grades=list(range(5)),private_phase_tuples=27,frustrated_survivors=len(survive))

def repair_assignment_checks():
    # Abstract full-record binary banks. Native charged supply is not asserted.
    q=5;labels=(0,1,2)
    grade={r:{sum(labels[j] for j in I)%q for I in combinations(range(3),r)} for r in range(4)}
    h1=set(range(q))-grade[1];h2=set(range(q))-grade[2]
    forbidden={(a-b)%q for a in h2 for b in h1}
    allowed=[]
    for lam in range(q):
        support=grade[2]|{(lam+a)%q for a in grade[1]}
        need((support==set(range(q)))==(lam not in forbidden),'exact forbidden differences define supplied repair graph')
        if lam not in forbidden:allowed.append(lam)
    need(allowed==[3,4],'small two-candidate repair graph')
    full=0
    for assignment in permutations(allowed):
        for debts in product(range(q),repeat=2):
            occupied=[]
            for i,(lam,debt) in enumerate(zip(assignment,debts)):
                opts=[]
                for V in combinations(range(4),2):
                    if sum((labels+(lam,))[j] for j in V)%q!=debt:continue
                    resources=tuple(('base',i,j,int(j in V)) if j<3 else ('repair',lam,int(j in V)) for j in range(4))
                    opts.append(resources)
                need(bool(opts),'matched actual binary repair supplies each fixed-grade debt')
                occupied+=list(opts[0])
            need(len(occupied)==len(set(occupied))==8,'private full reservations make simultaneous repaired sections disjoint')
            full+=1
    candidates=((3,2,1),(4,4,0),(3,1,1))  # (difference label, actual U charge, frozen outside bit)
    translated=0
    for assignment in permutations(range(3),2):
        outside=next(p for p in range(3) if p not in assignment)
        lam0,u0,frozenbit=candidates[outside]
        outside_charge=(u0+lam0*frozenbit)%q
        for debts in product(range(q),repeat=2):
            occupied=[('repair',outside,frozenbit)];actual_charge=outside_charge;actual_grade=frozenbit
            for i,(candidate,debt) in enumerate(zip(assignment,debts)):
                lam,baseline,_=candidates[candidate]
                chosen=next(V for V in combinations(range(4),2) if (baseline+sum((labels+(lam,))[j] for j in V))%q==debt)
                occupied += [('base',i,j,int(j in chosen)) if j<3 else ('repair',candidate,int(j in chosen)) for j in range(4)]
                actual_charge+=baseline+sum((labels+(lam,))[j] for j in chosen);actual_grade+=len(chosen)
            need(actual_charge%q==(sum(debts)+outside_charge)%q,'candidate-specific U baselines and frozen outside charge counted once')
            need(actual_grade==4+frozenbit,'unassigned full slot retains exact population contribution')
            need(len(occupied)==len(set(occupied))==9,'two-option reservations are nine occupied slot units, not eighteen')
            translated+=1
    return dict(scope='abstract binary full-record banks, not imported native numerical charges',charge_modulus=5,receiver_labels=list(labels),requested_grade=2,forbidden_labels=sorted(forbidden),successful_actual_candidate_labels=allowed,matched_assignment_debt_vectors_checked=full,candidate_specific_baseline_and_unused_slot_cases=translated)


P=((3,0,6,7,1,5,4,2),(3,1,6,7,0,5,4,2));G=(4,0,2,3,1,5,6,7)

def seed():
    a=[[0]*8 for _ in range(8)]
    for t,(i,j) in enumerate(combinations(range(8),2)):
        a[i][j]=(29558785>>(27-t))&1;a[j][i]=1-a[i][j]
    return tuple(map(tuple,a))

def source(r):
    a=seed();q=8*r
    return tuple(tuple(0 if i==j else a[i%8][j%8] if i<q and j<q and i//8==j//8 else int(i<j) if i<q and j<q else int(i<q) if (i<q)!=(j<q) else int((j-i)%3==1) for j in range(q+3)) for i in range(q+3))

def path(r,b,k):return tuple(8*j+u for j in range(r) for u in P[(b>>j)&1])+tuple(8*r+(k+j)%3 for j in range(3))
def transpose(a):return tuple(zip(*a))
def partner(a,p):
    where={u:i for i,u in enumerate(p)}
    return tuple(tuple(0 if u==v else a[u][v] if abs(where[u]-where[v])==1 else a[v][u] for v in range(len(a))) for u in range(len(a)))
def word(a,p):return tuple(a[p[i]][p[j]] for i in range(len(p)) for j in range(i+2,len(p)))
def normalized_closure(a,p):
    w=word(a,p);v=word(transpose(a),p[::-1])
    return frozenset((w,v,tuple(1-x for x in w),tuple(1-x for x in v)))
def isom(a,b,p):return all(a[u][v]==b[p[u]][p[v]] for u in range(len(a)) for v in range(len(a)))
def strong(a):
    for root in range(len(a)):
        seen={root};todo=[root]
        while todo:
            u=todo.pop()
            for v in range(len(a)):
                if a[u][v] and v not in seen:seen.add(v);todo.append(v)
        if len(seen)!=len(a):return False
    return True

def native_checks():
    base=seed();deg=list(map(sum,base));autom=[]
    for p in permutations(range(8)):
        if all(deg[i]==deg[p[i]] for i in range(8)) and isom(base,base,p):autom.append(p)
    need(autom==[tuple(range(8))],'seed rigidity checked on all permutations')
    hp=[p for p in permutations(range(8)) if all(base[u][v] for u,v in zip(p,p[1:]))]
    need(len(hp)==409,'complete native seed Hamiltonian census')
    pp=[p for p in hp if all(p[t]==P[0][t] for t in (0,2,3,5,6,7))]
    need(set(pp)==set(P),'complete six-common-pin native seed menu')
    need(strong(base),'seed genuinely strongly connected')
    rows=[]
    for r in range(1,5):
        a=source(r);n=len(a);records={(b,k):path(r,b,k) for b in range(1<<r) for k in range(3)}
        targets={z:partner(a,p) for z,p in records.items()};closures={z:normalized_closure(a,p) for z,p in records.items()}
        for b in range(1<<r):
            need(len({closures[b,k] for k in range(3)})==1,'all three kernel lifts consume same complete resource closure')
            need(len(closures[b,0])==4,'one full normalized K/converse unit has four words')
        need(len(set().union(*(closures[b,0] for b in range(1<<r))))==4*(1<<r),'distinct native resources have disjoint entire closures')
        for (b,k),p in records.items():
            need(sorted(p)==list(range(n)) and all(a[u][v] for u,v in zip(p,p[1:])),'literal native Hamiltonian path')
            need(all(p[8*j+t]==8*j+P[0][t] for j in range(r) for t in (0,2,3,5,6,7)),'all actual common indexed pins retained')
            need((p[8*r]-8*r)%3==k,'canonical phase offset retained')
            need(partner(transpose(a),p[::-1])==transpose(targets[b,k]),'literal converse square')
            need(strong(targets[b,k]),'actual native partner is strong')
            for l in range(3):
                perm=tuple(u if u<8*r else 8*r+(u-8*r+l-k)%3 for u in range(n))
                need(isom(targets[b,k],targets[b,l],perm),'kernel rotation preserves actual target via named transport')
            for j in range(r):
                g=G if not (b>>j)&1 else tuple(G.index(u) for u in range(8))
                perm=tuple(8*j+g[u%8] if u//8==j else u for u in range(n))
                need(isom(targets[b,k],targets[b^(1<<j),k],perm),'endpoint-marked prefix gluing authenticates constant target class')
            for (c,l),t in records.items():
                need(sum(u!=v for u,v in zip(p,t))==2*(b^c).bit_count()+3*(k!=l),'exact all-pair native owner-clock distance')
                need(sum(u!=v for u,v in zip(p[::-1],t[::-1]))==sum(u!=v for u,v in zip(p,t)),'converse retains exact cost')
                need((p[8*r]==t[8*r])==(k==l),'shared suffix owner relation iff phases agree')
                need((p[::-1][2]==t[::-1][2])==(k==l),'reflected converse clock retains shared equality')
        rows.append(dict(prefix_blocks=r,source_order=n,full_records=len(records),normalized_resources=1<<r,closed_words=4*(1<<r)))
    return dict(seed_hamiltonian_paths=len(hp),seed_automorphisms=len(autom),rows=rows)


def literal_menu_states(r,menus):
    a=source(r);states=[]
    for keytuples in product(*menus):
        keys=tuple(u for e in keytuples for u in e)
        for phases in product(range(3),repeat=len(keys)):
            records=tuple(path(r,b,k) for b,k in zip(keys,phases))
            if len({p[8*r] for p in records})!=1:continue
            reservoirs=[normalized_closure(a,p) for p in records]
            if sum(map(len,reservoirs))!=len(set().union(*reservoirs)):continue
            states.append((keys,phases))
    return states

def components(states,dist,budget):
    unseen=set(states);cc=[]
    while unseen:
        root=unseen.pop();seen={root};todo=[root]
        while todo:
            x=todo.pop()
            nxt=[y for y in unseen if dist(x,y)<=budget]
            unseen.difference_update(nxt);seen.update(nxt);todo+=nxt
        cc.append(seen)
    return sorted(map(len,cc))

def native_allocation_controls():
    A=((0,1),(2,3));B=((0,2),(1,3))
    no=literal_menu_states(2,[A,B]);need(not no,'native complete full-record capacity hostile')
    no3=literal_menu_states(3,[A,B]);need(not no3,'same-carrier S3 hostile before releasing unused keys')
    repair_one=literal_menu_states(3,[A,B+((0,4),)])
    repair_pack=literal_menu_states(3,[A,B+((4,5),)])
    need(len(repair_one)==3,'one-new-key native repair has three full lifts')
    need(len(repair_pack)==6,'three-pack native repair has six full assemblies')
    r=3;a=source(r);F=(0,1,4,5);Gk=(2,3,4,5)
    states=tuple(product(range(2),range(3)))
    records={(m,k):tuple(path(r,b,k) for b in (F if m==0 else Gk))+(path(r,6,m),) for m,k in states}
    for (m,k),rr in records.items():
        need(len({p[8*r] for p in rr[:4]})==1,'old actual shared-owner equality')
        rrsets=[normalized_closure(a,p) for p in rr]
        need(sum(map(len,rrsets))==len(set().union(*rrsets))==20,'five disjoint complete normalized closures at actual arrival')
    def dist(x,y):return sum(sum(u!=v for u,v in zip(p,q)) for p,q in zip(records[x],records[y]))
    def support(x,y):return sum(p!=q for p,q in zip(records[x],records[y]))
    for x,y in product(states,repeat=2):
        expected=0 if x==y else 12 if x[0]==y[0] else 7 if x[1]==y[1] else 19
        need(dist(x,y)==expected,'six-state exact literal endpoint clock cost')
        es=0 if x==y else 4 if x[0]==y[0] else 3 if x[1]==y[1] else 5
        need(support(x,y)==es,'six-state exact slot support')
    cc7=components(states,dist,7);cc12=components(states,dist,12)
    need(cc7==[2,2,2] and cc12==[6],'paid boundary plus internal component threshold')
    # Paired endpoint occupancy uses separate alphabets, not one union test.
    paired=Counter();wrong_single=0
    for x,y in product(states,repeat=2):
        if x[0]!=0 or y[0]!=1:continue
        resource_before=[normalized_closure(a,p) for p in records[x]]
        resource_after=[normalized_closure(a,p) for p in records[y]]
        if len(set().union(*(resource_before+resource_after)))==sum(map(len,resource_before+resource_after)):wrong_single+=1
        # Two arrivals have already been separately authenticated.
        paired[support(x,y),dist(x,y)]+=1
    need(dict(paired)=={(3,7):3,(5,19):6},'two-alphabet boundary coefficients retain actual minimizers')
    need(wrong_single==0,'one occupancy alphabet falsely kills every legal mode crossing')
    # A fixed allocation phase change never releases or reacquires any resource.
    for m in range(2):
        r0=[normalized_closure(a,p) for p in records[m,0]]
        r1=[normalized_closure(a,p) for p in records[m,1]]
        need(r0==r1,'fixed-key phase change continuously holds each actual resource')
    return dict(hostile_source_order=19,hostile_local_full_witnesses_per_demand=6,hostile_global_allocations=0,one_new_key_repair_full_records=len(repair_one),three_pack_repair_full_records=len(repair_pack),boundary_source_order=27,boundary_full_states=6,boundary_minimum_slots=3,boundary_minimum_clocks=7,budget7_components=cc7,budget12_components=cc12,paired_boundary_coefficients=[{'slots':s,'clocks':c,'coefficient':v} for (s,c),v in sorted(paired.items())],one_alphabet_wrong_survivors=wrong_single)


def source_pins():
    pins=json.loads((HERE/'SOURCE_PINS.json').read_text())
    for row in pins['files']:
        b=(HERE/row['local']).read_bytes()
        need(len(b)==row['bytes'] and hashlib.sha256(b).hexdigest()==row['sha256'],'byte-preserved historical source '+row['local'])
    return dict(commit=pins['commit'],files=len(pins['files']))

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output');args=ap.parse_args()
    result=dict(status='PASS',source_pins=source_pins(),hall=hall_checks(),hypergraphs=hypergraph_checks(),charge_relations=charge_relation_checks(),uniform_repair_assignment=repair_assignment_checks(),native_authentication=native_checks(),native_allocation=native_allocation_controls())
    result['checks']=CHECKS
    text=json.dumps(result,indent=2,sort_keys=True)+'\n'
    if args.output:Path(args.output).write_text(text)
    else:print(text,end='')

if __name__=='__main__':main()

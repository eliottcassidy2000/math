#!/usr/bin/env python3
"""Portable exact labelled-source and raw-tree ledger controls, Python >=3.9.
No third-party libraries, floating point, assert, network, or repository imports.
"""
from collections import Counter, defaultdict
from fractions import Fraction as F
from itertools import permutations, product
from pathlib import Path
import hashlib
import json

HERE = Path(__file__).resolve().parent
PIN = '96e35642110104cc547516a63f37eebcf930f707'
VHASH = 'df88f03b576b13c98c0c5378fc1cf61685eeb6817073681266d2a781dc329b26'

def check(ok, label):
    if not ok:
        raise RuntimeError(label)

def add(a,b): return (a[0]+b[0],a[1]+b[1])
def scale(c,a): return (c*a[0],c*a[1])
def total(xs):
    ans=(0,0)
    for x in xs: ans=add(ans,x)
    return ans

def Q(x): return x[0]*x[0]-x[0]*x[1]+x[1]*x[1]
def beta(x,y): return 2*x[0]*y[0]-x[0]*y[1]-x[1]*y[0]+2*x[1]*y[1]
def energy(xs): return sum(Q(x) for x in xs)
def debt(es): return sum(max(x-y,0) for x,y in zip(es,es[1:]))
def dist(z,a,b): return max(min(a,b)-z,z-max(a,b),0)

def matchings(labels):
    labels=tuple(labels)
    if not labels:
        yield ()
        return
    a=labels[0]
    for k in range(1,len(labels)):
        b=labels[k]
        rest=labels[1:k]+labels[k+1:]
        for tail in matchings(rest): yield ((a,b),)+tail

def matching_sums(xs,matching): return [add(xs[i],xs[j]) for i,j in matching]
def matching_energies(xs): return [energy(matching_sums(xs,m)) for m in matchings(range(len(xs)))]

def greedy_tree(xs):
    xs=list(xs); initial=energy(xs); P=Q(total(xs)); n=len(xs)
    check(n>0 and n&(n-1)==0,'power-two leaves')
    es=[initial]; witnesses=[]
    blocks=[(j,) for j in range(n)]
    while len(xs)>1:
        m=len(xs); E=energy(xs)
        choices=list(matchings(range(m)))
        chosen=min(choices,key=lambda z: energy(matching_sums(xs,z)))
        ys=matching_sums(xs,chosen); new=energy(ys)
        check(F(new)<=F(P+(m-2)*E,m-1),'minimum no greater than matching mean')
        check(F(max(new-E,0))<=F(P,m),'per-coarsening selected debt')
        blocks=[blocks[i]+blocks[j] for i,j in chosen]
        witnesses.append([list(b) for b in blocks]);xs=ys;es.append(new)
    es.reverse()
    D=debt(es)
    check(F(D)<=F((n-1)*P,n),'compatible greedy tree bound')
    return {'energies':es,'debt':D,'tree_blocks_bottom_up':witnesses}

def algebra_controls():
    palette=[(-1,0),(0,0),(1,0),(0,1),(2,-1)]
    tested=0
    for xs in product(palette,repeat=4):
        P=Q(total(xs)); C=energy(xs); ms=matching_energies(xs)
        check(sum(ms)==P+2*C,'three perfect-matchings identity')
        check(all(F(P,2)<=x<=2*C for x in ms),'sharp middle energy range')
        ds=[debt([P,m,C]) for m in ms]
        check(all(d==max(P-C,0)+dist(m,P,C) for m,d in zip(ms,ds)),'insertion excess identity')
        best_bound=min(F(P,2),max(F(C)-F(P,2),0))
        check(min(ds)-max(P-C,0)<=best_bound,'sharp deterministic local excess bound')
        worst=(dist(2*C,P,C)+2*dist(F(P,2),P,C))/3
        check(F(sum(ds),3)-max(P-C,0)<=worst,'sharp expected local excess bound')
        tested+=1
    high_arity=[]
    for n in (2,4,6,8):
        xs=[((i*i+2*i)%7-3,(3*i+1)%5-2) for i in range(n)]
        P=Q(total(xs));C=energy(xs);ms=matching_energies(xs)
        mean=F(sum(ms),len(ms))
        check(mean==F(P+(n-2)*C,n-1),'higher arity matching mean')
        high_arity.append({'n':n,'perfect_matchings':len(ms),'P':P,'C':C,'mean':str(mean)})
    hostile=[(2,-1),(2,-1),(0,1),(0,1)]
    P=Q(total(hostile));C=energy(hostile);ms=matching_energies(hostile)
    check((P,C,ms)==(16,16,[32,8,8]),'sharp Eisenstein hostile')
    # Labelled unit-source realization; no source deleted, no weight changed.
    sources=[('a1',0,(1,0)),('a2',0,(1,0)),('a3',0,(0,-1)),
             ('b1',1,(1,0)),('b2',1,(1,0)),('b3',1,(0,-1)),
             ('c1',2,(0,1)),('d1',3,(0,1))]
    check([total([v for _,b,v in sources if b==j]) for j in range(4)]==hostile,'unit-source realization')
    check(sum(Q(v) for _,_,v in sources)==8,'unit-source full norm retained')
    # Scalar integer alternative, independent of Eisenstein geometry.
    scalar=[(5,0),(3,0),(-1,0),(-1,0)]
    check((Q(total(scalar)),energy(scalar),matching_energies(scalar))==(36,36,[68,20,20]),'integer scalar hostile')
    # Disjoint terminals still use ONE global positive part per level.
    aggregate=[]
    for copies in (1,2,3,4,5,6):
        choices=sorted({sum(v) for v in product(ms,repeat=copies)})
        p=c=16*copies
        aggregate.append({'copies':copies,'minimum_excess':min(dist(m,p,c) for m in choices)})
        check((0 in [dist(m,p,c) for m in choices])==(copies%3==0),'whole-level cancellation across terminals')
    # Pairwise positives, per-parent positives, and global positives differ.
    per_parent=[4,-4]
    check(max(sum(per_parent),0)==0 and sum(max(v,0) for v in per_parent)==4,'positive-part placement hostile')
    return {'four_tuple_universe':len(palette)**4,'four_tuples_checked':tested,
            'higher_arity':high_arity,'sharp_hostile':{'children':hostile,'P':P,'C':C,'middle_energies':ms,
            'local_debts':[debt([P,m,C]) for m in ms],'mean_middle':str(F(sum(ms),3)),
            'mean_local_debt':str(F(sum(debt([P,m,C]) for m in ms),3)),
            'labelled_unit_sources':sources,'full_source_norm':8},'terminal_aggregation':aggregate}

def tree_controls():
    xs=[(1,0)]*4+[(-1,0)]*4;n=len(xs);P=Q(total(xs));C=energy(xs)
    totals=[0]*4;sum_D=0;count=0
    for order in permutations(range(n)):
        es=[]
        for b in (8,4,2,1):
            es.append(sum(Q(total(xs[i] for i in order[j:j+b])) for j in range(0,n,b)))
        for k,E in enumerate(es): totals[k]+=E
        sum_D+=debt(es);count+=1
    means=[F(x,count) for x in totals]
    for b,mean in zip((8,4,2,1),means):
        check(mean==F((b-1)*P+(n-b)*C,n-1),'compatible uniform tree level mean')
    check(debt(means)==0 and sum_D>0,'Jensen obstruction in compatible ensemble')
    count_greedy=0
    for signs in product((-1,1),repeat=8):
        ys=[(s,0) for s in signs];greedy_tree(ys);count_greedy+=1
    sharp=[]
    for n in (1,2,4,8):
        r=greedy_tree([(1,0)]*n)
        check(r['debt']==n*n-n,'greedy bound sharp on coherent sources')
        sharp.append({'n':n,'debt':r['debt'],'P':n*n})
    chosen=greedy_tree([(2,-1),(2,-1),(0,1),(0,1)])
    hostile_family=[]
    for h in range(1,11):
        n=2**h
        es=[0]+[n*n//(2**k) for k in range(1,h+1)]
        D=debt(es)
        check(D==n*n//2-n,'fixed low-bit debt family')
        hostile_family.append({'n':n,'source_norm':n,'terminal_energy':0,'native_shaped_debt':D})
    # Independent marginal partitions can fail refinement compatibility.
    coarse=[set(range(4)),set(range(4,8))]
    fine=[{i,i+4} for i in range(4)]
    check(not all(any(b<=a for a in coarse) for b in fine),'nonnested independent partitions')
    return {'permutations_checked':count,'compatible_tree_mean_energies':list(map(str,means)),
            'debt_of_mean_energies':str(debt(means)),'mean_actual_debt':str(F(sum_D,count)),
            'greedy_sign_words_checked':count_greedy,'greedy_coherent_sharpness':sharp,
            'greedy_sharp_hostile_tree':chosen,'unbounded_fixed_tree_family':hostile_family}

LOWER=((1,2),(-2,-1),(1,-1))
UPPER=((0,1),(-1,-1),(1,0))
def native_charge(s,N):
    # Literal complete B-chain, preserving zero original charges.
    A=(0,0);rows=[0,0,0];x=s
    while x<4*N:
        A=add(A,(LOWER if x<2*N else UPPER)[x%3])
        w=(0,-3,3) if x<2*N else (-1,-1,2)
        rows=[rows[r]+w[(x+r)%3] for r in range(3)]
        x=4*x+1
    a,b=A
    check(rows==[2*a-b,-a-b,-a+2*b],'independent native real-row reconstruction')
    check(Q(A) in (0,1,3),'native palette')
    return A

def native_edge(s):
    a=2 if s%8==1 else 1;t=a;y=(3*s+1)//2**a
    while y%8==5:y=(y-1)//4;t+=2
    check(2**t*(3*y+1)==3*(3*s+1)+2**a,'native exact affine clock')
    return y,t

def native_records(N,V):
    records=[]
    for s in range(1,4*N,2):
        if s%8==5:continue
        x=s;t=d=0;p=0;seen=set()
        while x not in V:
            check(x not in seen and d<10000,'finite route unresolved, not dropped')
            seen.add(x);p=x;x,dt=native_edge(x);t+=dt;d+=1
        A=native_charge(s,N);v=scale(-1 if t%2 else 1,A)
        records.append({'s':s,'u':x,'p':p,'t':t,'d':d,'A':list(A),'v':list(v),'m':(s-1)//2})
    return records

def native_ledger(records,all_records,N,cutoff):
    S=sum(Q(r['A']) for r in all_records);Sarr=sum(Q(r['A']) for r in records)
    terminals=sorted({r['u'] for r in records});qs=[];q=1
    while q<2*N:qs.append(q);q*=2
    qs.append(q);bins=[];es=[];zero_counts=[]
    for q in qs:
        values=defaultdict(lambda:(0,0))
        for r in records:
            key=(r['u'],r['m']%q);values[key]=add(values[key],r['v'])
        bins.append(values);es.append(sum(Q(v) for v in values.values()))
        # Unoccupied bins are implicit zeros, not absent from the partition.
        zero_counts.append({'q':q,'conceptual_bins':q*len(terminals),'occupied_by_labels':len(values),
                            'unoccupied_bins':q*len(terminals)-len(values)})
    deltas=[a-b for a,b in zip(es,es[1:])]
    pair_positive=0;pair_signed=0;pairs=0;by_level=defaultdict(int)
    by_terminal=defaultdict(list)
    for r in records:by_terminal[r['u']].append(r)
    # Direct pair enumeration is independent of bin-sum polarization.
    for rows in by_terminal.values():
        for j,r in enumerate(rows):
            for z in rows[:j]:
                gap=abs(r['m']-z['m']);check(gap>0,'distinct source label/odd address')
                lca=gap&-gap;b=beta(r['v'],z['v'])
                by_level[lca]+=b;pair_positive+=max(b,0);pair_signed+=b;pairs+=1
    check([by_level[q] for q in qs[:-1]]==deltas,'each labelled pair allocated once at actual first differing bit')
    check(es[-1]==Sarr,'separated source norm exact')
    check(es[0]-Sarr==sum(deltas)==pair_signed,'native signed telescope')
    D=sum(max(d,0) for d in deltas)
    check(es[0]<=S+D,'pending sources retained in denominator')
    per_parent=0
    for j,q in enumerate(qs[:-1]):
        fine=bins[j+1];total_beta=0
        # Parents with no labelled descendants have two zero children; skip
        # their zero arithmetic only, while retaining conceptual bin counts.
        parents={(u,r%q) for u,r in fine}
        for u,r in parents:
            b=beta(fine.get((u,r),(0,0)),fine.get((u,r+q),(0,0)))
            total_beta+=b;per_parent+=max(b,0)
        check(total_beta==deltas[j],'whole-level direct bin polarization')
    check(D<=per_parent<=pair_positive,'positive parts after progressively finer sums')
    # Re-pair native four-grandchild blocks independently, retaining all labels.
    native_pairing_controls=[]
    for j in range(len(qs)-2):
        q=qs[j];fine=bins[j+2];middle_sums=[0,0,0]
        patterns=[((0,1),(2,3)),((0,2),(1,3)),((0,3),(1,2))]
        parents={(u,r%q) for u,r in fine}
        for u,r in parents:
            four=[fine.get((u,r+k*q),(0,0)) for k in range(4)]
            for z,pat in enumerate(patterns):middle_sums[z]+=energy(matching_sums(four,pat))
        check(middle_sums[1]==es[j+1],'actual low-bit pairing identified')
        check(sum(middle_sums)==es[j]+2*es[j+2],'native-labelled mean middle identity')
        native_pairing_controls.append({'q':q,'P':es[j],'C':es[j+2],
                                        'global_same_pattern_middle_energies':middle_sums,
                                        'mean_middle':str(F(sum(middle_sums),3)),
                                        'two_level_debts':[debt([es[j],m,es[j+2]]) for m in middle_sums]})
    return {'N':N,'depth_cutoff':cutoff,'sources_total':len(all_records),'arrived_sources':len(records),
            'zero_charge_arrivals':sum(Q(r['A'])==0 for r in records),'S_full':S,'S_arrived':Sarr,
            'terminal_energy':es[0],'moduli':qs,'energies':es,'level_deltas':deltas,'D_plus':D,
            'sum_positive_parent_bin_interactions':per_parent,'sum_positive_labelled_pair_interactions':pair_positive,
            'same_terminal_labelled_pairs':pairs,'empty_bin_ledger':zero_counts,
            'four_child_controls':native_pairing_controls}

def native_controls():
    raw=(HERE/'frontier.txt').read_bytes()
    check(hashlib.sha256(raw).hexdigest()==VHASH,'pinned frontier sha256')
    V={int(s.split()[0]) for s in raw.decode().splitlines() if s.strip() and not s.startswith('#')}
    check(len(V)==627,'frontier627');check(all(native_edge(s)[0] in V for s in V),'frontier forward closed')
    results=[];digests={}
    for N in (512,1024,4096):
        records=native_records(N,V)
        payload=json.dumps(records,sort_keys=True,separators=(',',':'))+'\n'
        path=HERE/f'native_sources_N{N}.json';path.write_text(payload)
        digests[path.name]=hashlib.sha256(payload.encode()).hexdigest()
        cutoffs=([1,3,10,max(r['d'] for r in records)] if N==4096 else [max(r['d'] for r in records)])
        for d in sorted(set(cutoffs)):
            results.append(native_ledger([r for r in records if r['d']<=d],records,N,d))
    first=results[0]
    check((first['N'],first['depth_cutoff'],first['S_full'],first['terminal_energy'],first['D_plus'],
           first['sum_positive_parent_bin_interactions'],first['sum_positive_labelled_pair_interactions'])
          ==(512,42,1244,1313,79,293,433),'authenticated 433 versus79 native instance')
    rows=json.loads((HERE/'native_sources_N512.json').read_text());by_s={r['s']:r for r in rows}
    source_799,source_1065=by_s[799],by_s[1065]
    marked_cross=2*(2*source_799['v'][0]-source_799['v'][1])*(2*source_1065['v'][0]-source_1065['v'][1])
    check(source_799['u']==source_1065['u']==799
          and beta(source_799['v'],source_1065['v'])==0 and marked_cross==6,
          'native zero-Q positive-marked pair')
    local_witnesses=[]
    for N,q,u,parent_residue in ((512,1,435,0),(4096,4,443,0)):
        rows=json.loads((HERE/f'native_sources_N{N}.json').read_text())
        children=[[z for z in rows if z['u']==u and z['m']%(4*q)==parent_residue+k*q] for k in range(4)]
        xs=[total(z['v'] for z in group) for group in children]
        P=Q(total(xs));C=energy(xs);ms=matching_energies(xs)
        check(P==C==3 and sorted(ms)==[2,2,5],'native three-source no-zero-debt witness')
        check(sum(map(len,children))==3 and sum(bool(g) for g in children)==3,'actual empty bin retained')
        local_witnesses.append({'N':N,'q':q,'terminal':u,'parent_residue':parent_residue,'child_records':children,
                                'child_sums':xs,'P':P,'C':C,'middle_energies':ms,
                                'local_debts':[debt([P,m,C]) for m in ms],'native_middle_index':1})
    whole=next(z for z in first['four_child_controls'] if z['q']==16)
    check((whole['P'],whole['C'],whole['global_same_pattern_middle_energies'])
          ==(1255,1252,[1251,1242,1266]),'native global correlated-three-pattern obstruction')
    result={'source_pin':PIN,'frontier_sha256':VHASH,'source_artifact_sha256':digests,'controls':results,
            'native_minimal_support_pairing_hostiles':local_witnesses,
            'verified_433_vs_79':'N512, complete finite first arrival, max depth42; historical provenance not identified.',
            'native_marked_cross_hostile':{'N':512,'records':[source_799,source_1065],
                                         'beta_Q':0,'marked_cross':marked_cross}}
    # Regression: validate the final serialized consequence object, after every
    # local-witness loop, so source records cannot be overwritten by residues.
    serialized_records=json.loads(json.dumps(result))['native_marked_cross_hostile']['records']
    expected_records=[
        {'s':799,'u':799,'p':0,'t':0,'d':0,'A':[-2,-1],'v':[-2,-1],'m':399},
        {'s':1065,'u':799,'p':1065,'t':2,'d':1,'A':[0,1],'v':[0,1],'m':532}]
    check(serialized_records==expected_records,'serialized marked-hostile source records regression')
    return result

if __name__=='__main__':
    out={'status':'PASS','algebra':algebra_controls(),'trees':tree_controls(),'native':native_controls()}
    print(json.dumps(out,indent=2,sort_keys=True))

"""Independent complete-section controller and native physical-row replay.

The supply checker authenticates complete local domains. This checker rebuilds
all chosen path words, resources, charge sections and arrival row currents.
It uses no assert, so -O has the same checking behavior.
"""
from collections import Counter
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import argparse,json,random,sys

HERE=Path(__file__).resolve().parent
RESULTS=HERE.parents[1]
sys.path.insert(0,str(RESULTS/'injective_native_blocks_20261023/search'))
import explore as e
COUNT=0
U,V=1965,2234
REV=[e.CHORDS.index((7-j,7-i)) for i,j in e.CHORDS]

def need(ok,reason):
    global COUNT
    COUNT+=1
    if not ok:raise RuntimeError(reason)

def plus(a,b):return ((((a>>4)+(b>>4))%128)<<4)|((a^b)&15)
def minus(a,b):return ((((a>>4)-(b>>4))%128)<<4)|((a^b)&15)
def label(s):return ((e.X[s]%128)<<4)|sum((e.Y[j][s]%2)<<(j-1) for j in range(1,5))
def charge(ss):
    out=0
    for s in ss:out=plus(out,label(s))
    return out
def reflect(m):return sum(((m>>k)&1)<<REV[k] for k in range(21))
def tile(m):return sorted({m,reflect(m),e.FULL^m,e.FULL^reflect(m)})
def matrix(s):
    A=[[0]*8 for _ in range(8)];c=int(e.NODES[s]['canon_a'])
    for k,(i,j) in enumerate(combinations(range(8),2)):
        A[i][j]=(c>>(27-k))&1;A[j][i]=1-A[i][j]
    return A

def compile_route(a,ha,b,hb,section):
    a,b,ha,hb=map(frozenset,(a,b,ha,hb))
    need(len(a)==len(b) and len(ha)==len(hb)==3,'controller exact separate margins')
    q=plus(charge(a),charge(ha))
    need(q==plus(charge(b),charge(hb)),'controller accepted full endpoints')
    route=[(a,ha)];current=set(a)
    exchanges=list(zip(sorted(a-b),sorted(b-a)))
    for index,(x,y) in enumerate(exchanges):
        current.remove(x);current.add(y)
        hs=hb if index+1==len(exchanges) else section[minus(q,charge(current))]
        route.append((frozenset(current),frozenset(hs)))
    if not exchanges and ha!=hb:route.append((a,hb))
    need(route[-1]==(b,hb),'prescribed helper endpoint retained')
    for core,helpers in route:
        need(len(core)==len(a) and len(helpers)==3,'every arrival margins')
        need(plus(charge(core),charge(helpers))==q,'every arrival fullQ')
    for (c,h),(cc,hh) in zip(route,route[1:]):
        need(len(c^cc)+len(h^hh)<=8,'uniform eight-source support')
    return route

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true');args=parser.parse_args()
    old_path=RESULTS/'native_packet_decomposition_20261023/native/domains.json'
    new_path=HERE.parent/'supply/selected_supply.json'
    old=json.loads(old_path.read_text())['sources'];new=json.loads(new_path.read_text())
    bank=new['banks'][new['main_bank_index']]
    cores=sorted(z['source'] for z in old);helpers=bank['source_ids']
    need(len(cores)==28 and len(helpers)==40 and not set(cores)&set(helpers),'forty genuine new slots')
    need(bank['V_occupancy']==3 and bank['atomic_support_bound']==8,'selected main bank')
    data={}
    for z in old:
        choices={r['target']:r for r in z['exact_all_shared_clocks_target_domain']}
        data[z['source']]={'records':[choices[U],choices[V]],'fixed':z['fixed_clocks']}
    for z in new['helpers']:
        if z['source'] in helpers:
            m=z['selected'];data[z['source']]={'records':[m['before'],m['after']],'fixed':m['common_clocks']}
    reservoir=set();dist={}
    for s,d in sorted(data.items()):
        A=matrix(s)
        for b,r in enumerate(d['records']):
            p=r['path'];need(sorted(p)==list(range(8)) and all(A[i][j] for i,j in zip(p,p[1:])),'literal Hamiltonian source path')
            m=sum(A[p[i]][p[j]]<<k for k,(i,j) in enumerate(e.CHORDS))
            need(m==r['mask'] and e.MERGED[e.OWNER[m]]==s,'literal path source frame')
            need(e.MERGED[e.OWNER[e.FULL^m]]==(U,V)[b],'literal target owner')
            need(all(p[t]==v for t,v in d['fixed']),'all saved common pins')
            need((e.X[s]-e.X[(U,V)[b]])%8==r['phase'],'chosen actual original phase')
            orbit=set(tile(m));need(len(orbit)==4 and not orbit&reservoir,'full alternative port reservoir injective')
            reservoir.update(orbit)
            reverse=list(reversed(p));op=[[A[j][i] for j in range(8)]for i in range(8)]
            need(all(op[i][j] for i,j in zip(reverse,reverse[1:])),'literal converse path')
            rm=sum(op[reverse[i]][reverse[j]]<<k for k,(i,j) in enumerate(e.CHORDS))
            need(rm==reflect(m),'reflected converse word')
            need(all(reverse[7-t]==v for t,v in d['fixed']),'literal converse shared clocks')
        dist[s]=sum(x!=y for x,y in zip(d['records'][0]['path'],d['records'][1]['path']))
    need(len(reservoir)==544,'complete core/helper reservation544closedwords')

    # Independent triple enumeration: direct integer sum and XOR of native
    # marked potentials, not the producer's recursive section routine.
    section={};multiplicity=Counter()
    for ids in combinations(helpers,3):
        parity=0
        for s in ids:parity^=sum((e.Y[j][s]%2)<<(j-1) for j in range(1,5))
        q=((sum(e.X[s] for s in ids)%128)<<4)|parity
        section.setdefault(q,frozenset(ids));multiplicity[q]+=1
    need(len(section)==2048 and sum(multiplicity.values())==9880,'complete three-active charge section')
    for q in range(2048):
        mask=bank['charge_section'][str(q)]
        need(mask>=0 and mask<(1<<40) and mask.bit_count()==3,'saved section uses exactly three actual slots')
        chosen=frozenset(helpers[i] for i in range(40) if mask>>i&1)
        need(charge(chosen)==q and chosen==section[q],'independent canonical section agreement')

    # Rebuild the packet Laplacian action from actual selected black edges.
    # No scalar-code surrogate is used in this physical-row computation.
    potentials=[e.X]+[e.Y[j] for j in range(5)]
    moduli=[256]+[4]*5
    owner_ids=sorted(data)+[U,V]
    def rows(core,hs):
        active=core|hs;out={s:[0]*6 for s in owner_ids}
        for s in data:
            t=V if s in active else U
            for j,p in enumerate(potentials):
                w=2*(p[s]-p[t]);out[s][j]+=w;out[t][j]-=w
        return out
    def replay(route,save=False):
        baseline=rows(*route[0]);states=[];events=[]
        for c,h in route:
            current=rows(c,h)
            need(all((current[s][j]-baseline[s][j])%moduli[j]==0 for s in owner_ids for j in range(6)),'all actual original row currents, including diagonals')
            selected=c|h
            if save:
                records=[];ports=set()
                for s,d in sorted(data.items()):
                    r=d['records'][int(s in selected)];orb=tile(r['mask'])
                    need(not ports&set(orb),'no concurrent port reuse at saved arrival');ports.update(orb)
                    records.append(dict(source=s,target=r['target'],path=r['path'],word=r['mask'],normalized_resource=min(orb),closed_words=orb,
                        fixed_clocks=d['fixed'],converse_path=list(reversed(r['path'])),converse_fixed_clocks=[[7-t,v]for t,v in d['fixed']]))
                need(len(ports)==272,'actual arrival occupies68resources')
                states.append(dict(core_V=sorted(c),helper_V=sorted(h),integer_X_sum=sum(e.X[s]for s in selected),
                    physical_delta_rows={str(s):[current[s][j]-baseline[s][j]for j in range(6)]for s in owner_ids if current[s]!=baseline[s]},records=records))
        for (c,h),(cc,hh) in zip(route,route[1:]):
            changed=sorted((c^cc)|(h^hh));cost=sum(dist[s]for s in changed)
            need(cost<=sum(sorted((dist[s]for s in cores),reverse=True)[:2])+sum(sorted((dist[s]for s in helpers),reverse=True)[:6]),'actual integer-clock uniform bound')
            events.append(dict(changed_sources=changed,support=len(changed),owner_clock_cost=cost,
                changed_clocks={str(s):[i for i,(x,y)in enumerate(zip(data[s]['records'][0]['path'],data[s]['records'][1]['path']))if x!=y]for s in changed}))
        return states,events

    # All2,048 total charges exercised on a fixed core pair, then independently
    # varied endpoints at EVERY core margin, with noncanonical helper endpoints.
    cases=0;events_total=0
    a=frozenset(cores[:14]);b=frozenset(cores[14:])
    for q in range(2048):
        route=compile_route(a,section[minus(q,charge(a))],b,section[minus(q,charge(b))],section)
        # Full-row replay is deliberately done for every charge, not just q=0.
        replay(route);cases+=1;events_total+=len(route)-1
    rng=random.Random(2026102340)
    for k in range(29):
        for repeat in range(12):
            a=frozenset(rng.sample(cores,k));b=frozenset(rng.sample(cores,k));ha=frozenset(rng.sample(helpers,3))
            q=plus(charge(a),charge(ha));hb=section[minus(q,charge(b))]
            route=compile_route(a,ha,b,hb,section);replay(route);cases+=1;events_total+=len(route)-1

    a=frozenset([176,1425,1549,2087]);b=frozenset([1394,1923,2837,3369]);h=frozenset(helpers[:3])
    need(charge(a)==charge(b) and sum(e.X[s]for s in a)==sum(e.X[s]for s in b)==527,'smallest native hostile inherited')
    route=compile_route(a,h,b,h,section);states,events=replay(route,True)
    receipt=dict(status='FINITE-EXACT',checks=COUNT,core_sources=cores,helpers=helpers,helper_occupancy=3,
        charge_values=len(section),helper_triples=9880,route_cases=cases,checked_events=events_total,
        reserved_closed_words=len(reservoir),arrival_closed_words=272,
        max_core_pair_clock_cost=sum(sorted((dist[s]for s in cores),reverse=True)[:2]),
        max_six_helper_clock_cost=sum(sorted((dist[s]for s in helpers),reverse=True)[:6]),
        minimum_hostile_route=dict(states=states,events=events,total_owner_clock_cost=sum(x['owner_clock_cost']for x in events)),
        source_hashes={p.relative_to(RESULTS).as_posix():sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest()for p in (old_path,new_path)},
        scope='Complete fixed-margin fullQ fibres of the declared independent28+40slot packet; not full atlas or universal per-port phases')
    out=HERE/'receipt.json'
    if args.write:out.write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8')
    else:need(json.loads(out.read_text())==receipt,'frozen complete receipt agrees')
    print(json.dumps({k:v for k,v in receipt.items()if k not in ('minimum_hostile_route','core_sources','helpers','source_hashes')},indent=2))
    print('PASS',COUNT)

if __name__=='__main__':main()

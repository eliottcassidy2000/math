"""Build a definition-first sequence registry from the complete bounded census.
Usage: python3 -B analyze.py <census JSONL>; writes only this result directory.
No third-party packages. No novelty claim is attached to a registry identifier.
"""
from collections import Counter, defaultdict
from itertools import combinations
from pathlib import Path
import csv, gzip, hashlib, json, sys

HERE=Path(__file__).resolve().parent
SOURCE=HERE.parent/'session67_metagraph_enumeration_information_20261026/session67-root/output'
def need(x,s):
    if not x: raise RuntimeError(s)
def dump(name,obj): (HERE/name).write_text(json.dumps(obj,indent=2,sort_keys=True)+'\n',encoding='utf-8')
def freeze(x):
    if isinstance(x,dict):return tuple(sorted((k,freeze(v)) for k,v in x.items()))
    if isinstance(x,(tuple,list)):return tuple(map(freeze,x))
    return x
def histogram(x):return sorted(Counter(x).items())
records=[json.loads(s) for s in Path(sys.argv[1]).read_text(encoding='utf-8-sig').splitlines()]
with (HERE/'records.jsonl.gz').open('wb') as f:
    with gzip.GzipFile(fileobj=f,mode='wb',mtime=0,filename='') as z:
        for r in records:z.write((json.dumps(r,separators=(',',':'))+'\n').encode())
registry={}; summaries=[]; profiles=[]; witnesses={}; histograms={}; roles={}; adjacency_witnesses=[]
def register(key,definition,universe,kind='tracked here'):
    if key not in registry:registry[key]={'definition':definition,'universe':universe,'kind':kind,'offset':1,'status':'FINITE-EXACT','terms':[]}
def emit(key,value,definition,universe='ordinary isomorphism classes',kind='tracked here'):
    register(key,definition,universe,kind);registry[key]['terms'].append(value)
def signature_counts(name,values,definition,universe='ordinary isomorphism classes'):
    c=Counter(map(freeze,values));N=sum(c.values())
    emit(name+'_distinct',len(c),'Number of distinct '+definition,universe)
    emit(name+'_ambiguous_pairs',sum(v*(v-1)//2 for v in c.values()),'Number of unordered distinct class pairs with equal '+definition,universe)
    emit(name+'_largest_fibre',max(c.values(),default=0),'Maximum number of classes sharing '+definition,universe)

for n in range(1,9):
    rs=[r for r in records if r['n']==n];byrep={r['rep']:r for r in rs};groups=defaultdict(list)
    sig=defaultdict(list);rows=[];line_coefficients=set();unit_orbits=0;currents={}
    for r in rs:
        H=r['H'];groups[H].append(r);deg=[v.bit_count() for v in r['rows']];ep=r['endpoint']
        starts=list(map(sum,ep));ends=[sum(ep[i][j] for i in range(n)) for j in range(n)]
        E=sorted(set(r['extension']));hist=histogram(r['extension'])
        orbithist=histogram(r['extension'][s] for s,o in enumerate(r['attachment_orbit']) if s==o)
        role=sorted(zip(deg,starts,ends,r['delete_H']))
        dual=sorted(zip([n-1-d for d in deg],ends,starts,r['delete_H']))
        flow=Counter()
        for i in range(n):
            for j in range(n):
                if ep[i][j]:flow[deg[i]-deg[j]-1]+=ep[i][j]
        current=sum(k*v for k,v in flow.items());need(current%r['aut']==0,'integral current')
        currents[r['rep']]=current//r['aut'];seedprof=histogram((c,s) for _,c,s in r['seeds'])
        for _,c,s in r['seeds']:line_coefficients.add(c)
        unit_orbits+=len(r['seeds'])//r['aut']
        p={'n':n,'rep':r['rep'],'H':H,'aut':r['aut'],'tc':r['tc'],'sc':r['sc'],
           'extension_width':len(E),'extension_collision_pairs':sum(t*(t-1)//2 for _,t in hist),
           'neutral_attachments':r['extension'].count(H),'attachment_orbits':sum(s==o for s,o in enumerate(r['attachment_orbit'])),
           'unit_selections':len(r['seeds']),'unit_orbits':len(r['seeds'])//r['aut'],
           'endpoint_support':sum(v>0 for row in ep for v in row),'zero_flux_paths':flow[0],
           'current':current//r['aut'],'endpoint_pairs7':sum(v==7 for v in starts+ends),
           'ternary2':int(3*H+2 in E),'ternary4':int(3*H+4 in E)}
        rows.append(p);profiles.append(p)
        sig['H'].append(H);sig['tc'].append(r['tc']);sig['H_aut'].append((H,r['aut']))
        sig['extension_support'].append((H,E));sig['extension_distribution'].append((H,hist))
        sig['rooted_extension_distribution'].append((H,orbithist));sig['deletion_H'].append((H,sorted(r['delete_H'])))
        sig['endpoint_roles'].append((H,role));sig['endpoint_roles_converse'].append((H,min(role,dual)))
        sig['flux'].append((H,sorted(flow.items())));sig['unit_profile'].append((H,seedprof))
        sig['joint_response'].append((H,hist,orbithist,sorted(r['delete_H']),min(role,dual),seedprof))
    emit('ordinary',len(rs),'Number of ordinary classes',kind='baseline')
    emit('merged',len(set(r['merged'] for r in rs)),'Number of converse-merged classes','converse-merged classes','baseline')
    emit('strong',sum(r['strong'] for r in rs),'Number of strongly connected classes; singleton is strong',kind='baseline')
    emit('self_converse',sum(r['sc'] for r in rs),'Number of self-converse ordinary classes',kind='baseline')
    emit('symmetric',sum(r['aut']>1 for r in rs),'Number of ordinary classes with a nontrivial automorphism',kind='baseline')
    emit('path_presentations',sum(r['tc'] for r in rs),'Sum H/Aut; fixed-Hamilton-frame presentations',kind='baseline')
    definitions={'H':'H','tc':'H/Aut','H_aut':'(H,Aut order)',
        'extension_support':'(H, set of one-vertex child H values)',
        'extension_distribution':'(H, multiplicity distribution over all named attachment subsets)',
        'rooted_extension_distribution':'(H, distribution over Aut(parent)-orbits of attachment subsets)',
        'deletion_H':'(H, multiset of H after each vertex deletion)',
        'endpoint_roles':'(H, multiset of (outdegree, start count, end count, deletion H))',
        'endpoint_roles_converse':'converse-normalized endpoint-role signature',
        'flux':'(H, endpoint-score flux polynomial)',
        'unit_profile':'(H, multiset of (odd coefficient, first ray exponent) over unit selections)',
        'joint_response':'joint extension, rooted-extension, deletion, converse-normalized endpoint and unit signature'}
    for key,definition in definitions.items():signature_counts(key,sig[key],definition)
    for key,description in [
        ('extension_width','distinct child counts from one source'),
        ('extension_collision_pairs','pairs of distinct named attachments having equal child H'),
        ('neutral_attachments','attachments preserving H'),
        ('attachment_orbits','Aut(parent)-orbits of attachment subsets'),
        ('unit_selections','selected subsets with W=H-1 and U>0'),
        ('unit_orbits','Aut(parent)-orbits of positive unit selections'),
        ('endpoint_support','ordered endpoint pairs supporting a Hamiltonian path'),
        ('zero_flux_paths','Hamiltonian paths with d(start)-d(end)-1=0')]:
        emit(key+'_sum',sum(p[key] for p in rows),'Sum over classes of '+description)
        emit(key+'_max',max(p[key] for p in rows),'Maximum over classes of '+description)
    emit('unit_sources',sum(bool(r['seeds']) for r in rs),'Classes admitting at least one positive unit selection')
    emit('unit_odd_coefficients',len(line_coefficients),'Distinct odd parts of U among positive unit selections')
    emit('endpoint7_sources',sum(p['endpoint_pairs7']>0 for p in rows),'Classes with a start or end path count equal to 7')
    emit('ternary_both_sources',sum(p['ternary2'] and p['ternary4'] for p in rows),'Actual sources supporting both 3H+2 and 3H+4 in separate one-vertex attachments')
    mixed=0;falsepairs=0;fibrewitness=[];union_both=common_both=0
    for h,g in sorted(groups.items()):
        supports=[set(r['extension']) for r in g];targets=sorted(set.union(*supports));index={v:i for i,v in enumerate(targets)}
        masks=[sum(1<<index[v] for v in s) for s in supports];owners=defaultdict(list)
        for i,s in enumerate(supports):
            for v in s:owners[v].append(i)
        bad=0;example=None;allmask=(1<<len(targets))-1
        for t in targets:
            companions=0
            for i in owners[t]:companions|=masks[i]
            missing=(allmask^companions)&~((1<<(index[t]+1))-1)
            bad+=missing.bit_count()
            if missing and example is None:example=(t,targets[(missing&-missing).bit_length()-1])
        mixed+=bad>0;falsepairs+=bad
        union_both+=3*h+2 in index and 3*h+4 in index
        common_both+=any(3*h+2 in s and 3*h+4 in s for s in supports)
        if example:
            a,b=example;ra=next(r for r in g if a in r['extension']);rb=next(r for r in g if b in r['extension'])
            fibrewitness.append({'H':h,'pair':example,'sources':[ra['rep'],rb['rep']],
                'attachments':[ra['extension'].index(a),rb['extension'].index(b)],'incompatible_pairs':bad,'source_count':len(g)})
    emit('incompatible_target_pairs',falsepairs,'Sum over H-fibres of unordered distinct target pairs individually attainable but never from one common source','same-order H-fibres')
    emit('incoherent_H_fibres',mixed,'H-fibres containing at least one incompatible target pair','same-order H-fibres')
    emit('ternary_union_fibres',union_both,'H-fibres whose union of source supports contains both 3H+2 and 3H+4','same-order H-fibres')
    emit('ternary_common_fibres',common_both,'H-fibres with one source supporting both ternary targets','same-order H-fibres')
    witnesses[str(n)]=fibrewitness
    histograms[str(n)]={'H':histogram(r['H'] for r in rs),'tc':histogram(r['tc'] for r in rs),
        'extension_width':histogram(p['extension_width'] for p in rows),'unit_coefficients':sorted(line_coefficients),
        'H_aut':histogram((r['H'],r['aut']) for r in rs)}
    roles[str(n)]={str(v):{'H_classes':sum(r['H']==v for r in rs),'tc_classes':sum(r['tc']==v for r in rs),
        'aut_classes':sum(r['aut']==v for r in rs),'endpoint_entries':sum(x==v for r in rs for row in r['endpoint'] for x in row),
        'unit_coefficient_selections':sum(c==v for r in rs for _,c,s in r['seeds'])} for v in [7,11,33]}
    # Exact inherited coloured line carrier. Loops count twice in degrees.
    mergedgroups=defaultdict(list)
    for r in rs:mergedgroups[r['merged']].append(r)
    mr={m:g[0] for m,g in mergedgroups.items()};adj={m:Counter() for m in mr};flows=Counter();colour=Counter();degree=Counter();loops=Counter();equalH=0;tcdescent=set();Hdescent=set();pairs=set()
    path=SOURCE/f'lines_{n}.tsv';op=open
    if n==8:path=path.with_suffix('.tsv.gz');op=gzip.open
    for row in csv.DictReader(op(path,'rt'),delimiter='\t'):
        u,v,b,m=map(int,row.values());degree[u]+=m;degree[v]+=m;colour[b]+=m
        adj[u][v]+=m;adj[v][u]+=m
        flows[u]+=m*(mr[v]['c3']-mr[u]['c3']);flows[v]+=m*(mr[u]['c3']-mr[v]['c3'])
        if u==v:loops[u]+=m
        else:pairs.add((u,v))
        if mr[u]['H']==mr[v]['H']:equalH+=m
        for x,y in [(u,v),(v,u)]:
            if mr[x]['H']>mr[y]['H']:Hdescent.add(x)
            if mr[x]['tc']>mr[y]['tc']:tcdescent.add(x)
        if mr[u]['H'] in (11,33) and mr[v]['H'] in (11,33):adjacency_witnesses.append({'n':n,'u':u,'v':v,'blue':b,'multiplicity':m,'H':[mr[u]['H'],mr[v]['H']]})
    if n>=3:
        need(sum(colour.values())==1<<((n-1)*(n-2)//2-1),'line mass')
        for m,g in mergedgroups.items():
            need(degree[m]==sum(r['tc'] for r in g),'half-edge degree')
            need(flows[m]==sum(currents[r['rep']] for r in g),'endpoint/Kirchhoff current')
    emit('line_count',sum(colour.values()),'Literal complement lines, retaining multiplicity','converse-merged line multigraph','baseline')
    emit('line_blue',colour[1],'Literal blue (mirror-fixed frame) lines','converse-merged line multigraph','baseline')
    emit('line_loops',sum(loops.values()),'Literal loop lines','converse-merged line multigraph','baseline')
    emit('line_loop_nodes',len(loops),'Merged nodes carrying a literal loop','converse-merged line multigraph')
    emit('line_support_edges',len(pairs),'Distinct nonloop endpoint pairs, colours and multiplicities forgotten','simple support of merged line graph')
    emit('line_equal_H',equalH,'Literal lines whose two endpoint H values agree, including loops','converse-merged line multigraph')
    emit('line_H_descent_failures',sum(mr[m]['H']>1 and m not in Hdescent for m in mr),'Nodes with H>1 and no neighbour of smaller H','converse-merged line multigraph')
    emit('line_tc_descent_failures',sum(mr[m]['tc']>1 and m not in tcdescent for m in mr),'Nodes with tc>1 and no neighbour of smaller tc','converse-merged line multigraph')
    emit('line_zero_current_nodes',sum(flows[m]==0 for m in mr),'Nodes with zero cyclic-triangle current','converse-merged line multigraph')
    emit('line_current_square_sum',sum(flows[m]**2 for m in mr),'Sum of squared integer cyclic-triangle currents','converse-merged line multigraph')
    partner=[];coarse=[]
    for m,r in mr.items():
        partner.append((r['H'],r['aut'],r['sc'],sorted(Counter({}).items())))
        c=Counter()
        for v,mult in adj[m].items():c[mr[v]['H']]+=mult
        partner[-1]=(r['H'],r['aut'],r['sc'],sorted(c.items()))
        coarse.append((r['H'],r['aut'],r['sc']))
    signature_counts('merged_H_aut_sc',coarse,'(H, Aut, self-converse flag)','converse-merged classes')
    signature_counts('partner_H_response',partner,'(H, Aut, self-converse flag, multiplicity-weighted neighbour H profile); loops contribute twice','converse-merged classes')
    summaries.append({'n':n,'classes':len(rs),'named_attachments':len(rs)*(1<<n),'unit_selections':sum(len(r['seeds']) for r in rs),'new_registry_terms':len(registry)})
    print(json.dumps(summaries[-1]),flush=True)

dump('registry.json',registry);dump('summary.json',summaries);dump('histograms.json',histograms)
dump('source_obstructions.json',witnesses);dump('number_roles.json',roles);dump('number_line_witnesses.json',adjacency_witnesses)
with (HERE/'profiles.tsv').open('w',newline='',encoding='utf-8') as f:
    w=csv.DictWriter(f,fieldnames=list(profiles[0]),delimiter='\t');w.writeheader();w.writerows(profiles)
with (HERE/'sequences.tsv').open('w',newline='',encoding='utf-8') as f:
    w=csv.writer(f,delimiter='\t');w.writerow(['id','kind','universe','definition']+[f'n={n}' for n in range(1,9)])
    for k,v in registry.items():w.writerow([k,v['kind'],v['universe'],v['definition']]+v['terms'])
dump('inputs.json',{'universe':'complete inherited ordinary classes and coloured merged lines, n=1..8',
    'source_directory':str(SOURCE.relative_to(HERE.parents[2])).replace('\\','/'),
    'hash_convention':'gzip bytes unchanged; text CRLF normalized to LF',
    'sha256':{p.name:hashlib.sha256(p.read_bytes() if p.suffix=='.gz' else p.read_bytes().replace(b'\r\n',b'\n')).hexdigest() for p in SOURCE.iterdir() if p.name.startswith(('nodes_','lines_'))}})
print('Registered',len(registry),'exact sequences; no novelty or all-order recurrence inferred.')

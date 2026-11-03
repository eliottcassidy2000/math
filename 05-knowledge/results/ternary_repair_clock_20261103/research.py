"""Reproduce exact bank saturation, phase-preserving repairs and native programs.

Requires Python's standard library and g++17. Generated atlases and binaries
stay in a task-specific temporary directory; only compact certificates are saved.
"""
import argparse,hashlib,json,math,subprocess,sys,tempfile
from collections import Counter
from itertools import product,combinations,permutations
from pathlib import Path
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
sys.path.insert(0,str(HERE.parent/'golden_inverse_defects_20261103'))
import inverse as old
import native
TAU=(2,-1,0)
FAIL=((287,122,225),(289,121,225),(293,119,225),(295,118,225))
SEEDS={31:('a','a','r5'),95:('a','a','q1','a','a'),97:('a','a','a','a','a','q0')}


def digest(path):return hashlib.sha256(Path(path).read_text(encoding='utf-8').encode()).hexdigest()


def chart(n,code):
    rows=[0]*n;bit=0
    for a in range(n):
        for b in range(a+1,n):
            yes=b==a+1 or code>>bit&1;bit+=b!=a+1
            rows[a if yes else b]|=1<<(b if yes else a)
    return tuple(rows)


def fringe():
    seen=set(FAIL);stack=list(FAIL)
    while stack:
        v=stack.pop()
        for det,A in old.INVERSES.values():
            q=old.step(v,A)
            if any(x%det for x in q):continue
            u=tuple(x//det for x in q)
            if u[0]<1 or not u[0]%2 or min(u[1:])<0 or min(old.defect(u))<0:continue
            if sum(u[1:])>=sum(v[1:]):continue
            if u not in seen:seen.add(u);stack.append(u)
    return sorted(seen)


def run(command):
    return subprocess.run(list(map(str,command)),check=True,capture_output=True,text=True).stdout.strip()


def bank_file(path,bank):
    with path.open('w',encoding='utf-8',newline='\n') as f:
        for c,P in sorted(bank.items()):f.write(c+'\t'+'\t'.join(str(x) for row in P for x in row)+'\n')


def read_states(path):
    rows=[];index={}
    for line in path.read_text().splitlines():
        h,x,y,parent,c=line.split();v=(int(h),int(x),int(y));index[v]=len(rows);rows.append((v,int(parent),c))
    def word(v):
        j=index[v];out=[]
        while j:
            _,parent,c=rows[j];out.append(c);j=parent
        return tuple(out)
    return rows,index,word


def analyse_bank(index):
    states=set(index);answer={}
    for H in (31,95,97):
        fibre=[v for v in states if v[0]==H]
        both=[v for v in fibre if all((3*v[0]+2*j,3*v[1]-j,3*v[2]) in states for j in(1,2))]
        answer[str(H)]={'state_fibre':len(fibre),'both_exact_offsets':sorted(both)}
    answer['original_failed_states_present']=[v for v in FAIL if v in states]
    answer['original_clock_repair_fibres']={str(v):sorted(z for z in states if z[0]==v[0] and all((a-b)%7==0 for a,b in zip(z,v))) for v in FAIL}
    return answer


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true');args=parser.parse_args()
    cache=Path(tempfile.gettempdir())/'math-ternary-repair-clock-20261103';cache.mkdir(exist_ok=True)
    transcript={}
    for name in ('atlas7','audit_atlas','saturate','check_closure','root_probe'):
        run(['g++','-O3','-std=c++17',HERE/(name+'.cpp'),'-o',cache/(name+'.exe')])
    atlaspath=cache/'atlas7.tsv';transcript['compile7']=run([cache/'atlas7.exe',atlaspath])
    boundary=fringe();fringepath=cache/'fringe.tsv'
    fringepath.write_text(''.join('\t'.join(map(str,v))+'\n' for v in boundary),newline='\n')
    transcript['compile8fringe']=run([cache/'atlas7.exe','--initial8',fringepath,cache/'eight_hits.tsv'])
    transcript['independent_atlas']=run([cache/'audit_atlas.exe',atlaspath,fringepath])
    full={c:tuple(map(tuple,a['matrix'])) for c,a in old.ATLAS.items()}
    sources={c:tuple(a['source_rows']) for c,a in old.ATLAS.items()}
    for line in atlaspath.read_text().splitlines():
        z=list(map(int,line.split()));c='t'+str(z[0]);sources[c]=chart(7,z[0])
        full[c]=tuple(tuple(z[1+3*j:4+3*j]) for j in range(3))
    old.need(len(full)==31733,'complete seven-chart alphabet')
    for c,P in full.items():
        old.need(all(P[r][s]>=old.BANK['a'][r][s] for r in range(3) for s in range(3)),('singleton domination',c))
        old.need(P[1][0]%2==P[2][0]%2==0,('integral correction parity',c))
    counts=[];strong=[]
    for n in range(1,9):
        total=2**((n-1)*(n-2)//2);s=total-sum(strong[k-1]*counts[n-k-1] for k in range(1,n))
        counts.append(total);strong.append(s)
    old.need(strong==[1,0,1,5,50,903,30773,2032504],'independent SCC word recurrence')
    tau_symbols={c:(1 if old.step(TAU,P)==TAU else -1) for c,P in full.items()
                 if old.step(TAU,P) in (TAU,tuple(-x for x in TAU))}
    old.need(tau_symbols=={'b':-1,'q2':1,'s78':1},'complete signed tau inventory through7')
    allpath=cache/'all7.tsv';bank_file(allpath,full)
    banks={'65':old.BANK,'960':{c:P for c,P in full.items() if not c.startswith('t')},
           '67':old.BANK|{c:full[c] for c in ('t18413','t12083')}}
    analyses={};smallstates={}
    for label,bank in banks.items():
        p=cache/(label+'.tsv');out=cache/(label+'-states.tsv');bound=2671 if label=='67' else 1000
        bank_file(p,bank);transcript[label]=run([cache/'saturate.exe',bound,p,out])
        transcript[label+'_closure']=run([cache/'check_closure.exe',p,out,bound])
        rows,index,w=read_states(out);analyses[label]=analyse_bank(index)
        smallstates[label]=set(index)
        if label=='65':old.need(set(index)==old.forward_states(1000),'independent65 saturation')
    allstates=cache/'all7-states.tsv'
    transcript['31733']=run([cache/'saturate.exe',2671,allpath,allstates])
    transcript['31733_closure']=run([cache/'check_closure.exe',allpath,allstates,2671])
    rows,index,witness=read_states(allstates);analyses['31733']=analyse_bank(index)
    old.need(not any(v in index for v in FAIL),'failed exact states absent through7 at every order')
    old.need(not analyses['31733']['95']['both_exact_offsets'] and not analyses['31733']['97']['both_exact_offsets'],'no exact clock repair by reencoding')
    cost=lambda w:sum(len(sources[c]) for c in w)
    def state(w):
        v=old.ROOT
        for c in reversed(w):v=old.step(v,full[c])
        return v
    def literal(w):
        r=()
        for c in w:
            atom=sources[c];n,d=len(r),len(atom)
            r=tuple(x|(((1<<d)-1)<<n) for x in r)+tuple(x<<n for x in atom)
        return native.native(r,tuple(range(len(r))))
    def check_word(w):
        v=state(w);T=literal(w);H,E,C=old.endpoint_stats(T)
        old.need(v==(H,(C-E)//2,(C+E)//2),('independent literal state',w))
        old.need(len(w)>=2 and T[-1]&1 and all(T[j]>>(j+1)&1 for j in range(len(T)-1)),('strong Hamilton cycle',w))
        return T
    short={287:('t18413','a'),289:('r21','a','a','b'),293:('r30','r2'),295:('t12083','a')}
    compact=[]
    for H in (95,97):
        source=SEEDS[H];v=state(source)
        for j in (1,2):
            out=short[3*H+2*j];u=state(out);ideal=tuple(3*a+j*b for a,b in zip(v,TAU));delta=tuple(a-b for a,b in zip(u,ideal))
            old.need(delta[0]==0 and all(x%7==0 for x in delta),'compact residue repair')
            compact.append({'source_word':source,'target_word':out,'source_state':v,'target_state':u,
                            'digit':j,'carry':tuple(x//7 for x in delta),
                            'program':native.compile_repair(check_word(source),check_word(out))})
    # Actual same-order, same-count, same-F49/charge twins with different next H.
    twinwords=(('a','t11547'),('t16510','a'))
    twins=[state(w) for w in twinwords]
    old.need(twins==[(293,98,218),(293,147,204)],'literal clock twins')
    old.need(all(cost(w)==8 for w in twinwords),'same vertex clock')
    for w in twinwords:check_word(w)
    twin_next=[state(('q2',)+w)[0] for w in twinwords]
    old.need(twin_next==[3557,3865],'clock twin continuation')
    for w in twinwords:check_word(('q2',)+w)
    # Every Q4 phase acts on the complete carry, without identifying charts.
    q4checks=[]
    for c in ('q0','q1','q2','q4','q5'):
        a=compact[0];v=old.step(a['source_state'],full[c]);u=old.step(a['target_state'],full[c])
        correction=old.step(TAU,full[c]);carry=old.step(a['carry'],full[c])
        old.need(u==tuple(3*x+y+7*z for x,y,z in zip(v,correction,carry)),('Q4 transport',c))
        check_word((c,)+tuple(a['source_word']));check_word((c,)+tuple(a['target_word']))
        q4checks.append({'phase':c,'correction':correction,'carry':carry,'target_state':u})
    front=[(31,14%7,25%7)];keys=[];levels=[]
    for depth in range(1,5):
        front=[(3*h+2*j,(3*x-j)%7,3*y%7) for h,x,y in front for j in(1,2)]
        keys+=front;levels.append(front[:])
    wanted=set(keys);best={}
    for v in index:
        key=(v[0],v[1]%7,v[2]%7)
        if key in wanted:
            w=witness(v);candidate=(cost(w),w,v)
            if key not in best or candidate<best[key]:best[key]=candidate
    old.need(len(best)==30 and max(t[0] for t in best.values())==12,'four complete branching levels')
    chosen={H:(state(w),w) for H,w in SEEDS.items()}
    for key,(_,w,v) in best.items():
        if key[0] not in chosen:chosen[key[0]]=(v,w)
    graph={H:check_word(w) for H,(v,w) in chosen.items()}
    tree=[]
    for depth,layer in enumerate(levels,1):
        for H,x,y in layer:
            j=1 if H%3==2 else 2;parent=(H-2*j)//3;v,pw=chosen[parent];u,w=chosen[H]
            delta=tuple(a-3*b-j*c for a,b,c in zip(u,v,TAU))
            old.need(delta[0]==0 and all(t%7==0 for t in delta),'branch carry')
            field=lambda z:(4*(2*z[1]+z[0]-z[2])%7,z[2]%7)
            old.need(field(u)==tuple(3*t%7 for t in field(v)),'literal F49 scalar transport')
            old.need((u[0]-u[2]-3*(v[0]-v[2])-2*j)%7==0,'literal charge transport')
            tree.append({'depth':depth,'parent_H':parent,'H':H,'digit':j,'word':w,'state':u,
                         'carry':tuple(t//7 for t in delta),
                         'program':native.compile_repair(graph[parent],graph[H])})
    smallkeys={(h,x%7,y%7) for h,x,y in smallstates['67']}
    missing67=[[v for v in layer if v not in smallkeys] for layer in levels]
    q=full['q2'];q2=old.mm(q,q)
    O=tuple((int(j==0),q[j][0],q2[j][0]) for j in range(3))
    old.need(old.adjugate(O)[0]==96,'three q2 probes are faithful over Q and F7')
    def smith(P):
        d1=math.gcd(*(x for row in P for x in row))
        d2=math.gcd(*(P[a][c]*P[b][d]-P[a][d]*P[b][c] for a,b in combinations(range(3),2) for c,d in combinations(range(3),2)))
        return (d1,d2//d1,abs(old.adjugate(P)[0])//d2)
    closedq=((14,70,1022),(0,56,1008),(0,42,840))
    allphase=((14,14,14),(0,14,14),(0,0,14))
    old.need(smith(closedq)==(14,14,336) and smith(allphase)==(14,14,14),'closed carry observation Smith factors')
    kappas=[tuple(t['carry'][1:]) for t in compact]
    old.need(math.gcd(*(a*d-b*c for (a,b),(c,d) in combinations(kappas,2)))==1,'realized carries span integer plane')
    s=full['s78'];J=((0,1,1),(0,2,0),(1,-1,-1));R=((77,152,34),(28,55,12),(0,0,1))
    old.need(old.mm(s,J)==old.mm(J,R),'new tau-preserving split')
    sv=state(('s78',)+SEEDS[31]);old.need(sv==(4141,2203,2877),'zero phase literal source')
    check_word(('s78',)+SEEDS[31])
    phase=[];charges={6};alpha=(3,4)
    for d in range(13):
        phase.append({'depth':d,'field':alpha,'charges':sorted(charges)})
        alpha=tuple(3*x%7 for x in alpha);charges={(3*g+2*j)%7 for g in charges for j in(1,2)}
    old.need(all(len(r['charges'])==7 for r in phase[3:]),'all seven charge phases from depth3')
    # Incoming unequal-port mechanism: cheap local probe, with a literal control.
    control=[0]*5
    for bit,(a,b) in enumerate(combinations(range(5),2)):
        u,v=(a,b) if 41>>bit&1 else(b,a);control[u]|=1<<v
    h,e,c=old.endpoint_stats(control);old.need(h==13,'incoming H13 control')
    rootinput=cache/'root-input.tsv'
    with rootinput.open('w') as f:
        f.write(' '.join(map(str,[5,h,(c-e)//2,(c+e)//2]+control))+'\n')
        for H in(95,97):
            T=literal(SEEDS[H]);f.write(' '.join(map(str,(len(T),)+state(SEEDS[H])+T))+'\n')
    transcript['one_exception_root_probe']=run([cache/'root_probe.exe',rootinput,cache/'root-output.tsv'])
    expected=['H 13 tested 180 Hhits 2 0 residue 0 0 exact 0 0',
              'H 95 tested 1232 Hhits 0 0 residue 0 0 exact 0 0',
              'H 97 tested 1992 Hhits 0 0 residue 0 0 exact 0 0']
    old.need(transcript['one_exception_root_probe'].splitlines()==expected,'bounded unequal-port probe')
    independent=Counter();tested=0
    for root in range(1,5):
        kept=[u for u in range(5) if u!=root];base=sum(1<<j for j,u in enumerate(kept) if control[root]>>u&1)
        options=[base]+[base^(1<<j) for j in range(4)]
        for ports in product(options,repeat=3):
            if ports>ports[1:]+ports[:1] or ports>ports[2:]+ports[:2]:continue
            tested+=1;T=[sum(1<<j for j,v in enumerate(kept) if control[u]>>v&1) for u in kept]+[0]*3
            for k,mask in enumerate(ports):
                T[4+k]|=1<<(4+(k+1)%3)
                for j in range(4):
                    if mask>>j&1:T[4+k]|=1<<j
                    else:T[j]|=1<<(4+k)
            H=sum(all(T[a]>>b&1 for a,b in zip(p,p[1:])) for p in permutations(range(7)))
            # Strongness by reachability; no reliance on the source-word cycle guard.
            strong=True
            for u in range(7):
                seen=1<<u
                while True:
                    more=seen
                    for j,r in enumerate(T):
                        if seen>>j&1:more|=r
                    if more==seen:break
                    seen=more
                strong&=seen==127
            if strong and H in(41,43):independent[H]+=1
    old.need(tested==180 and independent==Counter({41:2}),'literal permutation control for root probe')
    sourcefiles=['golden_inverse_defects_20261103/THEORY.md','golden_inverse_defects_20261103/TRANSPORT.md',
                 'golden_inverse_defects_20261103/inverse.py','golden_seam_recurrence_20261103/atom_atlas.json',
                 'golden_hamilton_transfer_20261103/RECURRENCE.md','adaptive_metagraph_reformulation_20261103/METAGRAPH.md',
                 'inverse65_zeckendorf_20261103/INCOMING.md',
                 'session154_gadgets_20261103/closure/NOTE.md',
                 'session154_gadgets_20261103/interfaces/README.md']
    data={'status':'PROVED scoped / FINITE-EXACT','transcript':transcript,'bank_analysis':analyses,
          'strong_chart_counts':strong,'fringe':boundary,'tau_symbols':tau_symbols,
          'compact_repairs':compact,'four_level_tree':tree,'67_missing_by_depth':missing67,
          'clock_twins':{'words':twinwords,'states':twins,'next_H':twin_next},'Q4_transport':q4checks,
          'carry_Smith':{'q2_only':smith(closedq),'all_phases':smith(allphase)},
          'q2_observation_matrix':O,'s78_split':R,'s78_zero_phase_state':sv,'ternary_clock':phase,
          'generated_sha256_lf':{'atlas7':digest(atlaspath),'closed_set':digest(allstates)},
          'source_sha256_lf':{s:digest(HERE.parent/s) for s in sourcefiles}}
    path=HERE/'receipt.json'
    if args.write:path.write_text(json.dumps(data,sort_keys=True,separators=(',',':'))+'\n',encoding='utf-8',newline='\n')
    else:old.need(json.loads(path.read_text())==json.loads(json.dumps(data)),'complete receipt reproduction')
    print(json.dumps({'status':'PASS','atlas7':30773,'eight_fringe':2032504,'closed_states':len(index),
                      'branch_edges':len(tree),'maximum_endpoint_order':max(map(len,graph.values())),
                      'native_moves':sum(t['program']['native_moves'] for t in tree),
                      'maximum_peak_order':max(t['program']['peak_order'] for t in tree)}))


if __name__=='__main__':main()

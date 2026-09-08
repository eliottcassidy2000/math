"""Independent native-grid referee for last-B and the full7200 composition."""
from collections import Counter
from fractions import Fraction as Q
from functools import lru_cache
from itertools import combinations,combinations_with_replacement,product
from math import gcd,lcm
from pathlib import Path
import hashlib,json,sys
sys.stdout.reconfigure(newline='\n')
HERE=Path(__file__).resolve().parent
T=7200
G=0
def gate(ok,label):
    global G
    G+=1
    if not ok: raise RuntimeError(label)
def allowed_sum(n):
    for p in range(2,n+1):
        if n%p: continue
        exponent=0
        while n%p==0: exponent+=1;n//=p
        if p%3!=2 or exponent>2:return False
    return n==1
atlas=[(p,total-p) for total in range(3,357) if allowed_sum(total)
       for p in range(1,(total+1)//2) if gcd(p,total-p)==1]
gate(len(atlas)==5855,'complete strict coprime atlas')
@lru_cache(None)
def spatial(p,q):
    L=14*p*q; cuts={0,L}
    for v in (p,q):
        unit=L//(14*v)
        for k in range(v):
            cuts.update(((14*k-1)*unit%L,(14*k+1)*unit%L))
    cut=sorted(cuts); intervals=[]; den=2*L
    for lo,hi in zip(cut,cut[1:]):
        num=lo+hi
        if all(14*min(v*num%den,(-v*num)%den)<den for v in (p,q)):
            intervals.append((lo,hi))
    gate(intervals[0][0]==0 and intervals[-1][1]==L,'literal open origin component')
    return L,[(intervals[-1][0]-L,intervals[0][1])]+intervals[1:-1]
def literal(n,phase,p,q):
    den=n*phase.denominator; num=phase.numerator; out=0
    for j in range(n):
        out+=all(14*min(v*num%den,(-v*num)%den)<den for v in (p,q))
        num+=phase.denominator
    return out
@lru_cache(None)
def profile(n,p,q):
    p,q=sorted((p,q)); L,I=spatial(p,q)
    enters=Counter(n*a%L for a,b in I); exits=Counter(n*b%L for a,b in I)
    walls=sorted(enters.keys()|exits.keys())
    phase=Q(walls[-1]+walls[0]+L,2*L)
    initial=cur=literal(n,phase,p,q)
    minimum=maximum=cur; minphase=maxphase=phase
    for i,w in enumerate(walls):
        wall=cur-exits[w]; after=wall+enters[w]
        nextw=walls[i+1] if i+1<len(walls) else walls[0]+L
        for value,phase0 in ((wall,Q(w,L)),(after,Q(w+nextw,2*L))):
            if value<minimum:minimum=value;minphase=phase0
            if value>maximum:maximum=value;maxphase=phase0
        cur=after
    gate(cur==initial,'native circular event balance')
    gate(literal(n,minphase,p,q)==minimum and literal(n,maxphase,p,q)==maximum,
         'independent literal extremum owners')
    return minimum,maximum,minphase,maxphase,len(walls)
def arms(left,right):
    e=gcd(left,right); n=T//e; result=[]
    for p,q in atlas:
        for a,b in ((p,q),(q,p)):
            if e*gcd(n,a)==left and e*gcd(n,b)==right:
                credit=e*profile(n,p,q)[0]
                result.append((a,b,credit))
    return sorted(result)

def locate(name):
    for candidate in (HERE/name,HERE.parent/'05-knowledge/results'/name,
                      Path('C:/w/s0905/04-computation')/name,
                      Path('C:/w/s0905/05-knowledge/results')/name):
        if candidate.exists():return candidate
    raise FileNotFoundError(name)
def pinned(name,expected):
    raw=locate(name).read_bytes()
    gate(hashlib.sha256(raw).hexdigest()==expected,'exact frozen dependency '+name)
    return raw
stem='continuing10_20260907_lrc_last_b'
pinned(stem+'.py','3cb1c77dadcdccf2778e24df58dfad9a0bc8a24395698ee973c85afc0a0ecef0')
pinned(stem+'.out','571884a84d72e7a55c89c7a89835bf86577066b375774935101b0e9f5e1e63fd')
cert=json.loads(pinned(stem+'_certificate.json','01381d94aaeaf92ca7dfb37b902c29db16da0aa36c5587b602dd836b5616cabd'))
mst=json.loads(pinned('continuing8_20260906_lrc_minimum_tree_certificate.json',
    '580a7c930103aab3bea867ad463a90b0e0208323a90ee95a685ff811a761582d'))
old=next(row for row in mst['clocks'] if row['t']==T)
weights={tuple(pair):row[0] for pair,row in old['weights']}
word=(5,8,9,30,32,36,48)
E=sum(d*((T+7*d-1)//(7*d)) for d in word)-T
gate(E==103 and cert['word']==list(word) and cert['E']==E,'literal full-word excess')
full=[(weights[tuple(sorted((word[i],word[j])))],i,j) for i,j in combinations(range(7),2)]
cheap=[row for row in full if row[0]<=E]
gate(sorted(full)==sorted(map(tuple,cert['all_edges'])),'all21 inherited possible edge credits')
gate(sorted(cheap)==sorted(map(tuple,cert['cheap_graph'])) and len(cheap)==7,'entire affordable graph')
def connected(edge_list):
    found={0}
    while True:
        enlarged=found|{j for _,i,j in edge_list if i in found}|{i for _,i,j in edge_list if j in found}
        if enlarged==found:return len(found)==7
        found=enlarged
connected_graphs=[];trees=[]
for mask in range(1<<7):
    present=[cheap[i] for i in range(7) if mask>>i&1]
    if not connected(present):continue
    connected_graphs.append(present)
    actual={tuple(sorted((word[i],word[j]))) for _,i,j in present}
    gate((9,30) in actual and (36,48) in actual,'two actual mandatory bridges')
    gate((8,48) in actual or (32,48) in actual,'actual doubled zero-arm path forced')
    gate(all({9,30}.isdisjoint({a,48,36}) for a in (8,32)),
         '102 edge disjoint from either forced path')
    if len(present)==6:trees.append(tuple(sorted(present)))
gate(len(connected_graphs)==4 and len(trees)==3,'complete connected subgraph and tree universe')
gate(sorted(trees)==sorted(tuple(sorted(map(tuple,tree))) for tree in cert['all_cheap_spanning_trees']),
     'every cheap tree compared independently by edge subsets')
gate([c for c,i,j in cheap if {word[i],word[j]}=={9,30}]==[102],
     'mandatory native102 credit')
domains={a:arms(a,48) for a in (8,32,36)}
zeros={a:[(p,q) for p,q,c in rows if c==0] for a,rows in domains.items()}
gate([len(domains[a]) for a in (8,32,36)]==[316,633,231],'all1180 freshly computed directed arms')
gate([len(zeros[a]) for a in (8,32,36)]==[10,51,18],'all fresh zero arms')
for a,rows in domains.items():
    gate(all(c==0 or (c%gcd(a,48)==0 and c+102>103) for _,_,c in rows),
         'every positive arm closes with disjoint102 edge')
    saved=next(row for row in cert['arm_domains'] if row['margins']==[a,48])
    gate(sorted(map(tuple,saved['arms']))==rows,'entire arm list and all exact native credits')

results={};minimums=[]
for a in (8,32):
    rows=[]
    for u,v in zeros[a]:
        for w,z in zeros[36]:
            L=lcm(v,z);R=[u*(L//v),L,w*(L//z)];G0=gcd(*R);R=tuple(x//G0 for x in R)
            gate(len(set(R))==3 and tuple(gcd(T,4*x) for x in R)==(a,48,36),
                 'every native scale realization passes full joint valuation tree')
            D=gcd(R[0],R[2]);p,q=sorted((R[0]//D,R[2]//D));e=gcd(T,4*D)
            gate(e==4,'every endpoint really has quotient1800')
            c=e*profile(T//e,p,q)[0]
            gate(c>103,'actual nonatlas endpoint alone closes the word')
            rows.append(((u,v),(w,z),R,(e,p,q,c)))
    saved=next(row for row in cert['wedges'] if row['margins']==[a,48,36])
    expected=sorted((tuple(r['arms'][0]),tuple(r['arms'][1]),tuple(r['primitive']),tuple(r['pair'])) for r in saved['rows'])
    gate(sorted(rows)==expected,'every joint product and endpoint credit compared')
    minimum=min(row[-1][-1] for row in rows);minimums.append(minimum)
    gate(minimum==saved['minimum'],'sharp fresh family endpoint minimum')
    results[a]=rows
gate([len(results[a]) for a in (8,32)]==[180,918] and minimums==[112,140],
     'complete1098 joint universe and exact sharp minima')
gate(2*profile(3600,578,801)[0]==114 and 4*profile(1800,578,801)[0]==112,
     'hostile: doubled margins do not preserve old credit114')
print('DOUBLED_WEDGES arm_counts',[len(domains[a]) for a in (8,32,36)],
      'zero_counts',[len(zeros[a]) for a in (8,32,36)],'products',[len(results[a]) for a in (8,32)],
      'minimum_credits',minimums,flush=True)

gate(len(cert['profiles'])==2277,'complete frozen profile cardinality')
keys=set()
for row in cert['profiles']:
    key=(row['n'],row['p'],row['q']);keys.add(key)
    actual=profile(*key)
    gate((actual[0],actual[1],actual[4])==(row['minimum'],row['maximum'],row['walls']),
         'every native profile extrema and wall count')
    L,I=spatial(row['p'],row['q'])
    gate((L,len(I))==(row['L'],row['components']),'all geometric ruler and component counts')
    gate(literal(row['n'],Q(*row['minimizer']),row['p'],row['q'])==row['minimum'] and
         literal(row['n'],Q(*row['maximizer']),row['p'],row['q'])==row['maximum'],
         'literal native grids at all frozen extremum owners')
gate(len(keys)==2277,'all saved profile keys distinct')

control=cert['positive_control'];A=tuple(T*k for k in range(1,7));B=word
gate(tuple(control['A'])==A and tuple(control['B'])==B and len(set(A+B))==13 and gcd(*(A+B))==1,
     'actual primitive distinct unitless thirteen-speed positive control')
actual_edges=[]
for i,j in combinations(range(7),2):
    d=gcd(B[i],B[j]);p,q=sorted((B[i]//d,B[j]//d))
    if (p,q) in set(atlas):actual_edges.append([i,j])
gate(actual_edges==control['edges'] and connected([(0,i,j) for i,j in actual_edges]),
     'actual native strict graph connected, not inferred from margin graph')
safe=[]
for j in range(T):
    den=7*T;num=7*j+1
    if all(14*min(v*num%den,(-v*num)%den)>=den for v in A+B):safe.append(j)
gate(len(safe)==control['safe_count']==2449,'literal full-row weak-safe lifts')
gate(hashlib.sha256(json.dumps(safe,separators=(',',':')).encode()).hexdigest()==control['safe_indices_sha256'],
     'entire literal safe-index set')
print('TOPOLOGY cheap_edges',len(cheap),'connected_subgraphs',len(connected_graphs),'spanning_trees',len(trees))
print('NATIVE_PROFILES',len(keys),'all strict extrema and literal owners verified')
print('SCALE_HOSTILE ratio578:801 old_sheet2_clock3600_credit114 new_sheet4_clock1800_credit112')
print('FULL_ROW_SAFE_LIFTS',len(safe),'at alpha1/7')

# Composition uses the separately proved/audited A theorem, not a second A
# census. Match every indexed word before subtracting the single clock.
third=json.loads(pinned('continuing10_20260907_lrc_third_wedge_certificate.json',
    '60a6de57a623dbe64a23a77113f140b312f5496247c2253fbe695872740039e9'))
last_a=json.loads(pinned('continuing10_20260907_lrc_last_a_certificate.json',
    'bce97b849036ba89ca8e8d1593b2ab03ff7aa9ecc647bdc72d55bd4962383274'))
gate(old['word_count']==76814 and old['survivor_count']==len(old['survivors'])==15,
     'complete inherited7200 necessary word universe')
gate([(r['word'],r['E']) for r in third['topology']]==[(w,E) for w,E,_ in old['survivors']],
     'entire topology scope exactly equals inherited residual universe')
still=[r['word'] for r in third['topology'] if not r['closed']]
gate(sum(r['closed'] for r in third['topology'])==13 and
     still==third['remaining_words']==[last_a['word'],list(word)],
     'thirteen deletions plus disjoint A and B exhaust all fifteen words')
gate(last_a['clock']==7200 and last_a['E']==116 and last_a['inherited_pin']==cert['inherited_sha256'],
     'A and B use the same actual clock and inherited quantitative supplier')
scales=mst['new_scales']
canonical=lambda x:hashlib.sha256(json.dumps(x,separators=(',',':')).encode()).hexdigest()
gate(len(scales)==len(set(scales))==7646 and scales.count(7200)==1 and
     canonical(scales)=='8ffc6d14b3883cf7e02c3ab02ddca5339d909a8051411def9096dee83b0aaed7',
     'exact inherited necessary clock array')
new=[t for t in scales if t!=7200]
gate(len(new)==7645 and max(new)==max(scales)==11995 and
     canonical(new)=='f8c42793c4a5081d40a3937c9a79d9cd307e6c08c8b9ea3e94bac5388e47d16c',
     'single-clock subtraction preserves every other inherited candidate')
print('COMPOSITION inherited_words76814 residual15 third_deletes13 A_and_B_close_remaining2')
print('NEW_NECESSARY_CLOCKS',len(new),'maximum',max(new),'semantic_sha256',canonical(new))
print('PASS',G,'always-active exact gates')

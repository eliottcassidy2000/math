"""Independent raw-wall/native-grid referee for third wedge and topology."""
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
left=arms(24,18);right=arms(30,18)
gate(len(left)==231 and len(right)==141,'complete directed arm universes')
low=[];inconsistent=[];consistent=[];high=0
for a,b,ac in left:
    for c,d,bc in right:
        if ac+bc>107:high+=1;continue
        L=lcm(b,d); R=[a*(L//b),L,c*(L//d)]; common=gcd(*R); R=tuple(v//common for v in R)
        gate(len(set(R))==3,'all low normalized triples have distinct speeds')
        row=(a,b,ac,c,d,bc,R)
        low.append(row)
        if tuple(gcd(T,6*v) for v in R)!=(24,18,30):
            inconsistent.append(row);continue
        D=gcd(R[0],R[2]);p,q=R[0]//D,R[2]//D;e=gcd(T,6*D)
        endpoint=e*profile(T//e,min(p,q),max(p,q))[0]
        forest=sum(sorted((ac,bc,endpoint))[-2:])
        gate(forest>107,'exact actual two-edge forest pays whole-word bound')
        consistent.append(row+(endpoint,forest))
gate((len(low),high,len(inconsistent),len(consistent))==(30,32541,26,4),'all32571 path products and exact split')
gate(min(row[-1] for row in consistent)==162,'sharp low-path forest minimum')
print('THIRD_WEDGE arms',len(left),len(right),'all products',len(left)*len(right),'high',high,'low',len(low),'inconsistent',len(inconsistent),'valid',len(consistent))
print('VALID_PATHS',consistent)
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
stem='continuing10_20260907_lrc_third_wedge'
pinned(stem+'.py','a768f216da83a8882cd2ef2d52abe08f344a971035e0618aa74242d0b623e6fe')
pinned(stem+'.out','f8ebfb3f7d302bdf24274934c80d06d8c663df3bcc52c0feae6e400cfca9e044')
cert=json.loads(pinned(stem+'_certificate.json','60a6de57a623dbe64a23a77113f140b312f5496247c2253fbe695872740039e9'))
wedges=json.loads(pinned('continuing10_20260907_lrc_composite_wedges_certificate.json',
    '890d00d44f0195765d62fe1d40b59ad102311f34a84a42a5fddc229d037209e9'))
supplier=json.loads(pinned('overnight12_20260906_lrc_decoder_descent_inherited_profiles.json',
    '935f3f687b6d7c89cc099e536f536238fd753bcc4c1747906d213cef387ca93f'))
mst=json.loads(pinned('continuing8_20260906_lrc_minimum_tree_certificate.json',
    '580a7c930103aab3bea867ad463a90b0e0208323a90ee95a685ff811a761582d'))
gate(sorted(map(tuple,cert['arm_domains'][0]['arms']))==left,'all first-arm rows')
gate(sorted(map(tuple,cert['arm_domains'][1]['arms']))==right,'all second-arm rows')
gate(cert['arm_sum_closed']==high,'all high products')
observed_bad={tuple(row['primitive']):(tuple(map(tuple,row['arms'])),tuple(row['actual_margins']))
              for row in cert['valuation_rejections']}
expected_bad={row[6]:(((row[0],row[1],row[2]),(row[3],row[4],row[5])),
                     tuple(gcd(T,6*v) for v in row[6])) for row in inconsistent}
gate(observed_bad==expected_bad and len(observed_bad)==26,'all rejected joint rows')
for R,(_,actual) in expected_bad.items():
    gate(all(gcd(a,modulus)==gcd(b,modulus) for a,b in zip(actual,(24,18,30))
             for modulus in (32,25)) and actual!=(24,18,30),
         'each rejected row fails only prime3 depth')
gate(sorted((tuple(c['primitive']),c['pair'][-1],c['forest']) for c in cert['accepted'])==
     sorted((row[6],row[-2],row[-1]) for row in consistent),'all accepted triples and exact credits')
gate(len(cert['profiles'])==376,'complete saved native profiles')
keys=set()
for row in cert['profiles']:
    key=(row['n'],row['p'],row['q']);keys.add(key)
    actual=profile(*key)
    gate((actual[0],actual[1],actual[4])==(row['minimum'],row['maximum'],row['walls']),
         'all native extrema and exact wall cardinalities')
    L,I=spatial(row['p'],row['q'])
    gate((L,len(I))==(row['L'],row['components']),'literal geometric component normalization')
    gate(literal(*[row['n'],Q(*row['minimizer']),row['p'],row['q']])==row['minimum'] and
         literal(*[row['n'],Q(*row['maximizer']),row['p'],row['q']])==row['maximum'],
         'producer extremal phases verified by independent literal grids')
gate(len(keys)==376,'saved native profile uniqueness')
# An exact weak-endpoint hostile: the two endpoints are outside strict danger.
gate(literal(7,Q(1,2),1,1)==0,'weak endpoint control retained')
closed=sum(14*min(Q(2*j+1,14),1-Q(2*j+1,14))<=1 for j in range(7))
gate(closed==2,'closing strict danger would create two false credits')

# Independent full-profile evaluation by literal subsets, not the producer's
# dynamic gcd array. The four free slots are an exhaustive multiset domain.
tables={int(k):{(c,tuple(word)) for c,word in row['profiles']}
        for k,row in supplier['levels'].items()}
domain=[d for d in supplier['levels']['6']['gcds'] if T%d==0]
gate(len(domain)==26 and len(set(domain))==26,'complete 7200 divisor alphabet')
def word_valid(word):
    if gcd(*word)!=1:return False
    for size in range(1,7):
        for indices in combinations(range(7),size):
            c=gcd(*(word[i] for i in indices))
            rest=tuple(sorted(gcd(c,word[j]) for j in range(7) if j not in indices))
            if (c,rest) not in tables[7-size]:return False
    return True
def cost(word):return sum(d*((T+7*d-1)//(7*d)) for d in word)-T
conditioned=[];attempted=0
for free in combinations_with_replacement(domain,4):
    attempted+=1;word=tuple(sorted((24,18,30)+free))
    if word_valid(word):conditioned.append((word,cost(word)))
gate(attempted==23751 and len(conditioned)==162,'entire conditional full-profile universe')
gate(max(E for word,E in conditioned)==107,'complete sharp conditional ceiling')
gate(((9,16,18,24,25,30,32),107) in conditioned,'exact attaining full word')
saved=next(row for row in wedges['conditional_words'] if sorted(row['forced'])==[18,24,30])
gate(sorted(conditioned)==sorted((tuple(w),E) for w,E in saved['words']),
     'every conditioned word compared with frozen primary')

# An independent exhaustive tree path: all 7^5 labeled Pruefer words, without
# a greedy spanning-tree algorithm or the producer's cut-condition verifier.
edges=list(combinations(range(7),2));edge_index={edge:i for i,edge in enumerate(edges)}
trees=[]
for sequence in product(range(7),repeat=5):
    degree=[1]*7
    for i in sequence:degree[i]+=1
    selected=[]
    for i in sequence:
        leaf=next(j for j in range(7) if degree[j]==1)
        selected.append(edge_index[tuple(sorted((leaf,i)))])
        degree[leaf]-=1;degree[i]-=1
    selected.append(edge_index[tuple(j for j in range(7) if degree[j]==1)])
    trees.append((sum(1<<j for j in selected),tuple(selected)))
gate(len(trees)==len({mask for mask,_ in trees})==16807,'all labeled spanning trees exactly once')
old=next(row for row in mst['clocks'] if row['t']==T)
weights={tuple(pair):record[0] for pair,record in old['weights']}
gate(len(old['survivors'])==15,'entire inherited necessary-word universe')
all_branches=0;deleted=[];remaining=[];branch_report=[]
for wi,(word,E,old_floor) in enumerate(old['survivors']):
    gate(word_valid(word) and cost(word)==E,'inherited word native full-profile and cost')
    # Every old 24 policy is retained, even if its second target class is absent.
    # The new 18 policy is only non-vacuous when both 24 and 30 occur.
    policies=[(i,{18},{4,16}) for i,v in enumerate(word) if v==24]
    policies += [(i,{24},{30}) for i,v in enumerate(word) if v==18 and 24 in word and 30 in word]
    row=cert['topology'][wi]
    gate((row['word'],row['E'],row['word_index'])==(word,E,wi),'topology positions including repeated margins')
    gate(row['policies']==[[i,sorted(a),sorted(b)] for i,a,b in policies],'complete independent branch policies')
    base={edge_index[i,j]:weights[tuple(sorted((word[i],word[j]))) ]
          for i,j in edges if tuple(sorted((word[i],word[j]))) in weights}
    tree_costs=[(mask,sum(base.get(j,10**9) for j in js)) for mask,js in trees]
    costs=[]
    for choices in product((0,1),repeat=len(policies)):
        all_branches+=1;removed=set()
        for choice,(i,a,b) in zip(choices,policies):
            removed.update(tuple(sorted((i,j))) for j in range(7) if j!=i and word[j] in (a if choice==0 else b))
        allowed={index:weight for index,weight in base.items() if edges[index] not in removed}
        mask=sum(1<<index for index in allowed)
        options=[v for tree_mask,v in tree_costs if tree_mask&mask==tree_mask]
        tree_min=min(options) if options else None
        branch=next(b for b in row['branches'] if b['choices']==list(choices))
        gate(branch['removed']==[list(e) for e in sorted(removed)],'all forbidden positional edges')
        gate(sorted(branch['allowed_edges'])==sorted([weight,*edges[j]] for j,weight in allowed.items()),
             'entire generous graph compared without using producer graph construction')
        gate(tree_min==branch['cost'],'exhaustive Pruefer minimum equals certified tree cost')
        chosen=[edge_index[tuple(sorted((i,j)))] for w,i,j in branch['tree']]
        gate(len(chosen)==6 and sum(1<<j for j in chosen) in {m for m,_ in trees} and
             all(j in allowed for j in chosen) and sum(allowed[j] for j in chosen)==tree_min,
             'saved native attaining tree is a literal labeled spanning tree')
        gate(branch['closed']==(tree_min is None or tree_min>E),'strict sufficient credit gate')
        costs.append(tree_min)
    gate(len(costs)==len(row['branches']),'all branches with no omissions or duplicate producer entries')
    closed=all(x is None or x>E for x in costs)
    gate(row['closed']==closed,'complete role-disjunction quantifier')
    if closed:deleted.append(wi)
    else:remaining.append(word)
    branch_report.append((wi,E,costs,closed))
gate(all_branches==37 and len(deleted)==13,'all37 branches and thirteen word deletions')
gate(deleted==[1,2,3,4,5,6,7,8,9,10,11,13,14],'exact deletion indices')
gate(remaining==cert['remaining_words']==[[1,9,16,18,24,32,60],[5,8,9,30,32,36,48]],
     'two necessary-condition survivors; no unsafe-row inference')
print('CONDITIONED_WORDS attempted',attempted,'valid',len(conditioned),'maximum',max(E for _,E in conditioned))
print('NATIVE_PROFILES',len(keys),'strict walls and literal extrema verified')
print('TOPOLOGY labeled_trees',len(trees),'branches',all_branches,'deleted',deleted)
for row in branch_report:print('WORD',row)
print('REMAINING',remaining)
print('PASS',G,'always-active exact gates')

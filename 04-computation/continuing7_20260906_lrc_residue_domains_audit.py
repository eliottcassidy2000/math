"""Independent full residue-domain referee; no producer implementation imports."""
from bisect import bisect_left, bisect_right
from collections import Counter
from fractions import Fraction as F
from itertools import combinations_with_replacement
from math import comb, gcd, lcm
from pathlib import Path
import hashlib
import json
import sys
sys.stdout.reconfigure(newline='\n')
HERE=Path(__file__).resolve().parent
STEM='continuing7_20260906_lrc_residue_domains'
G=0
def gate(ok, name):
    global G
    G+=1
    if not ok: raise RuntimeError(name)
def find(name):
    for path in (HERE/name,HERE.parent/'05-knowledge/results'/name,
                 Path('C:/w/s0905/04-computation')/name,
                 Path('C:/w/s0905/05-knowledge/results')/name):
        if path.is_file(): return path
    raise FileNotFoundError(name)
def digest(obj):
    return hashlib.sha256(json.dumps(obj,separators=(',',':')).encode()).hexdigest()
for ext,pin in (
    ('.py','1afc8c3438c9db2ad6bad53b771d84a9b6e46d4f6464722aa4ce240d73084e1b'),
    ('.out','f7554bdb450c2f6a042c0caef9270eba4832934682d72b1c8cfaaff5db223858'),
    ('_certificate.json','4764d3def5d24702465f702180bad715b14402e9c45ddcd72ae1968e676a8086')):
    raw=find(STEM+ext).read_bytes()
    gate(hashlib.sha256(raw).hexdigest()==pin,'frozen producer '+ext)
    if ext=='.out': gate(b'\r' not in raw,'producer literal LF')
cert=json.loads(find(STEM+'_certificate.json').read_bytes())
raw=find('overnight12_20260906_lrc_decoder_descent_inherited_profiles.json').read_bytes()
gate(hashlib.sha256(raw).hexdigest()=='935f3f687b6d7c89cc099e536f536238fd753bcc4c1747906d213cef387ca93f','full profile input')
data=json.loads(raw)
profiles={int(k):{(c,tuple(w)) for c,w in v['profiles']} for k,v in data['levels'].items()}
allowed=data['levels']['6']['gcds']
gate(len(allowed)==42 and lcm(*allowed)==5388292800 and all(gcd(d,7)==1 for d in allowed),'exact inherited state domain')
global_caps=[]
for line in find('third_20260906_grid_full_words.out').read_text().splitlines():
    if line.startswith('GLOBAL '):
        row=json.loads(line[7:]); global_caps.append([row[0],row[2]])
gate(global_caps==[[1,210],[2,270],[3,192],[4,239],[5,197],[6,224]]==cert['inherited_global_caps'],'proved global cap supplier')

def full_word(w):
    """Subset bitmask gcd recurrence; different from producer subset tuples."""
    gs=[0]*128
    for mask in range(1,128):
        bit=mask & -mask
        gs[mask]=gcd(gs[mask^bit],w[bit.bit_length()-1])
    if gs[127]!=1: return False
    for mask in range(1,127):
        c=gs[mask]
        rem=tuple(sorted(gcd(c,w[j]) for j in range(7) if not(mask>>j & 1)))
        if (c,rem) not in profiles[7-mask.bit_count()]: return False
    return True

def frequency_words(domain, size):
    """Enumerate multiplicity vectors, not positional combinations."""
    word=[]
    def visit(i,left):
        if i==len(domain)-1:
            yield tuple(word+[domain[i]]*left)
            return
        d=domain[i]
        for count in range(left+1):
            word.extend([d]*count)
            yield from visit(i+1,left-count)
            if count: del word[-count:]
    yield from visit(0,size)

domain_results={}
for rec in cert['clocks']:
    T=rec['clock']; domain=tuple(d for d in allowed if T%d==0)
    gate(list(domain)==rec['allowed_divisors'],'complete clock divisor domain')
    gate(gcd(T,lcm(*allowed))==rec['state_gcd'],'exact arithmetic domain key')
    L=lcm(*domain)
    # Actual ceiling costs on seven chosen multiples of L. No inverse-mod-seven formula.
    clocks=[next(L*k for k in range(1,8) if L*k%7==tau) for tau in range(7)]
    costs=[{d:d*((t+7*d-1)//(7*d)) for d in domain} for t in clocks]
    tables=[{key:None for key in combinations_with_replacement(domain,2)} for _ in range(7)]
    words=[]; total=0
    for w in frequency_words(domain,7):
        total+=1
        if not full_word(w): continue
        native=sum(d*((T+7*d-1)//(7*d)) for d in w)-T
        words.append([list(w),native])
        Es=[sum(cost[d] for d in w)-t for cost,t in zip(costs,clocks)]
        gate(Es[T%7]==native,'native equivalent-residue ceiling transport')
        have=Counter(w)
        pairs=[key for key in tables[0] if all(have[d]>=count for d,count in Counter(key).items())]
        for tau,E in enumerate(Es):
            for key in pairs:
                old=tables[tau][key]
                if old is None or E>old[0] or (E==old[0] and w<old[1]): tables[tau][key]=(E,w)
    words.sort()
    gate(total==comb(len(domain)+6,7)==rec['word_universe'],'complete frequency-vector universe')
    gate(words==rec['valid_words_and_costs'] and len(words)==rec['valid_word_count'],'all full words/costs, no omitted types')
    gate(digest(words)==rec['valid_word_sha256'],'full-word semantic pin')
    maxima=[]
    for tau,(table,produced) in enumerate(zip(tables,rec['residue_tables'])):
        actual=[[a,b,None if value is None else value[0],None if value is None else list(value[1])]
                for (a,b),value in table.items()]
        gate(actual==produced['conditional_maxima'],'whole residue conditional table including nulls and lex owners')
        M=max(value[0] for value in table.values() if value is not None)
        maxima.append(M)
        gate(M==produced['maximum'] and produced['residue']==tau,'residue maximum and correct index')
        for a,b,E,owner in actual:
            if owner is None: continue
            gate(full_word(owner),'every attaining owner independently satisfies all profiles')
            gate(Counter(owner)[a]>=1+(a==b) and Counter(owner)[b]>=1,'forced pair positional multiplicity')
    gate(maxima[T%7]==rec['global_maximum'],'actual global maximum')
    gate(sum(v is not None for v in tables[0].values())==rec['realized_margin_count'],'complete feasible margin count')
    domain_results[T]=tables[T%7]
    print('DOMAIN',T,'universe',total,'valid',len(words),'residue maxima',maxima)

# Independent atlas generated by admissible prime-product sums, then reduced splits.
prime_list=[p for p in range(2,357) if all(p%d for d in range(2,int(p**.5)+1))]
def good_sum(s):
    for p in prime_list:
        e=0
        while s%p==0: e+=1; s//=p
        if e and (p%3!=2 or e>2): return False
    return s==1
atlas=sorted((p,s-p) for s in range(3,357) if good_sum(s)
             for p in range(1,(s+1)//2) if p<s-p and gcd(p,s-p)==1)
gate(len(atlas)==5855 and [list(x) for x in atlas]==cert['atlas'],'complete exact atlas from sum-first enumeration')

def spatial(p,q):
    """Raw boundary arrangement and literal midpoint predicates; no arc intersections."""
    L=14*p*q; cuts={0,L}
    for v in (p,q):
        width=L//(14*v)
        for k in range(v):
            for sign in (-1,1): cuts.add((14*k+sign)*width%L)
    cuts=sorted(cuts); intervals=[]; den=2*L
    for lo,hi in zip(cuts,cuts[1:]):
        num=lo+hi
        if all(14*min(v*num%den,(-v*num)%den)<den for v in (p,q)):
            intervals.append((lo,hi))
    gate(intervals[0][0]==0 and intervals[-1][1]==L,'native circle origin joins')
    return L,[(intervals[-1][0]-L,intervals[0][1])]+intervals[1:-1]
geometry={key:spatial(*key) for key in atlas}

def grid(n,alpha,p,q):
    """Binary-search actual grid in the independent geometric intervals, including endpoints."""
    L,I=geometry[p,q]
    a,d=alpha.numerator,alpha.denominator
    points=[(a+j*d)*L for j in range(-n,n)]
    return sum(bisect_left(points,hi*n*d)-bisect_right(points,lo*n*d) for lo,hi in I)

def event_profile(n,p,q):
    L,I=geometry[p,q]
    ent,ext=Counter(),Counter()
    for a,b in I:
        ent[n*a%L]+=1; ext[n*b%L]+=1
    walls=sorted(ent.keys()|ext.keys())
    alpha=F(walls[-1]+walls[0]+L,2*L)%1
    current=grid(n,alpha,p,q); initial=current
    values={}
    for i,w in enumerate(walls):
        wall_value=current-ext[w]
        current=wall_value+ent[w]
        a=F(w,L); b=F(w+walls[(i+1)%len(walls)]+(L if i==len(walls)-1 else 0),2*L)%1
        values[a]=wall_value; values[b]=current
        gate(grid(n,a,p,q)==wall_value and grid(n,b,p,q)==current,'independent binary-search wall and open-cell counts')
    gate(initial==current,'full native event sweep returns to initial cell')
    return values,len(walls)

profiles_by_n={}
for n in range(2478,2486):
    values,walls=event_profile(n,1,355)
    lo,hi=min(values.values()),max(values.values())
    gate((lo,hi)==((50,51) if n<2485 else (0,51)),'all selected-edge near/resonance controls')
    profiles_by_n[n]=(values,walls)
    print('NATIVE',n,'minimum',lo,'maximum',hi,'walls',walls)
    if n==2485: gate(values[F(1,2)]==0,'exact resonant half phase is actually empty')
for q,s,expect in ((43,6,6),(43,7,5),(15,7,1)):
    if (1,q) not in geometry: geometry[1,q]=spatial(1,q)
    values,walls=event_profile(7*q-s,1,q)
    gate(min(values.values())==expect and max(values.values())==2*(q//14)+1,'strict/equality inherited hostile with endpoints')

total_rows=0
for rec in cert['clocks']:
    T=rec['clock']; table=domain_results[T]; rows=[]; counts=Counter()
    for e in range(1,7):
        if T%e: continue
        n=T//e
        for p,q in atlas:
            L,I=geometry[p,q]
            sep=e*sum((n*(b-a)+L-1)//L-1 for a,b in I)
            dp,dq=e*gcd(n,p),e*gcd(n,q)
            value=table.get(tuple(sorted((dp,dq))))
            bound=None if value is None else value[0]
            mode='infeasible' if bound is None else 'separate' if sep>bound else 'located'
            if mode=='located':
                gate((e,p,q)==(6,1,355),'exact sole location consumer')
                gate(e*min(profiles_by_n[n][0].values())>bound,'actual located credit beats same owner cost')
            counts[mode]+=1
            rows.append([e,p,q,dp,dq,sep,bound,mode])
    gate(rows==rec['all_candidates'],'every edge/sheet input compared entry by entry')
    gate(len(rows)==rec['base_count'] and dict(counts)==rec['conditional_first_counts'],'complete accounting')
    gate([row[:3]+[row[5]] for row in rows if row[5]<=rec['global_maximum']]==[[6,1,355,0]],'every noncritical edge closes even by global cost')
    gate(min((row[5],*row[:3]) for row in rows if row[:3]!=[6,1,355])==tuple(rec['least_noncritical_credit']),'sharp least noncritical credit')
    stored=rec['located_profile']; values,walls=profiles_by_n[T//6]
    gate(walls==stored['wall_count']==102 and len(stored['walls_and_cells'])==2*walls,'entire wall/cell phase universe')
    gate({F(a):c for a,c in stored['walls_and_cells']}==values,'all native phase counts compared')
    gate(min(values.values())==stored['minimum'] and max(values.values())==stored['maximum'],'exact located extrema')
    total_rows+=len(rows)
    print('CLOCK',T,'all candidates',len(rows),'decisions',dict(counts),'credit',6*min(values.values()))
gate(total_rows==81970,'all three complete candidate universes')

# General selected-edge proof parameters. This is not a further clock scan.
gate(7*351<2485 and 6*50>max(E for tau,E in global_caps),'entire strict s<=7 range')
gate(6*gcd(2485,355)==2130 and 2130 not in allowed,'exact resonance violates full-profile state domain')
L,I=geometry[1,355]
gate(len(I)==51 and all(F(b-a,L)==F(1,2485) for a,b in I),'all full teeth have exact common length')
gate(6*sum((2486*(b-a)+L-1)//L-1 for a,b in I)==306>270,'monotone far-range lower bound')
for T,h,beta,raw in cert['multiplicity_controls']:
    n=T//6; alpha=(h*F(beta))%1
    counts=Counter(h*j%n for j in range(T))
    gate(gcd(h,n)==1 and set(counts.values())=={6} and len(counts)==n,'literal residue map has exact sixfold multiplicity')
    gate(6*grid(n,alpha,1,355)==raw>=300,'independent geometric physical lift count')

old=json.loads(find('continuing6_20260906_lrc_near_resonance_certificate.json').read_bytes())['new_scales']
gate(digest(old)=='c146d73c149c448e744e138aa3bb6f8c286748be81353b8a6605adc47a129117','full old necessary set pin')
new=[t for t in old if t not in (14886,14880,14874)]
gate(len(old)==8198 and len(new)==8195 and max(new)==14868,'actual count and maximum after exact deletions')
gate(new==cert['new_scales'] and set(old)-set(new)=={14886,14880,14874},'full exact inherited-set comparison')
gate(digest(new)=='c67f5e98eff3fe406a4416057c6063095290330a50f039e731ced0fc2ca4657a','complete new necessary set pin')
print('NEW SET',len(new),'maximum',max(new),'SHA256',digest(new))
print('PASS',G,'always-active independent exact gates; raw LF')
print('Scope: three complete clocks; uniform selected edge only; inherited proper-six input; weak safety, LRC14 open.')

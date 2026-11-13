"""Selective ternary factorization forests and literal one-shot source controls."""
from collections import Counter
from functools import cache
from itertools import permutations
from pathlib import Path
import hashlib
import json
import sys

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]
SEEDS=(1,3,5,23,49,63,67,343)
LIMIT=100001

def need(ok,message):
    if not ok: raise RuntimeError(message)

def factor(n):
    result=[]
    d=3
    while d*d<=n:
        while n%d==0: result.append(d);n//=d
        d+=2
    if n>1: result.append(n)
    return result

def atoms(n):
    need(n>0 and n%2 and n not in (7,21),'allowed arithmetic target')
    if n==1:return []
    a=0
    while n%7==0:a+=1;n//=7
    other=factor(n)
    if a==0:return other
    if a%2==0:return [49]*(a//2)+other
    if a>=3:return [343]+[49]*((a-3)//2)+other
    chosen=next((p for p in other if p>=5),None)
    if chosen is not None:
        other.remove(chosen)
        return [7*chosen]+other
    need(len(other)>=2,'63 carry boundary')
    return [63]+other[2:]

@cache
def forest_stats(n):
    if n in SEEDS:return (0,0,((n,1),))
    pieces=atoms(n)
    if len(pieces)>1:
        gates=depth=0;leaves=Counter()
        for value in pieces:
            g,d,s=forest_stats(value)
            gates+=g;depth=max(depth,d);leaves.update(dict(s))
        return gates,depth,tuple(sorted(leaves.items()))
    need(pieces==[n] and n%3,'unseeded relative atom')
    digit=2 if n%3==2 else 4
    parent=(n-digit)//3
    need(parent>=3 and parent not in (7,21) and parent<n,'guarded smaller source')
    gates,depth,leaves=forest_stats(parent)
    return gates+1,depth+1,leaves

def certificate(n):
    if n in SEEDS:return {'kind':'seed','value':n}
    pieces=atoms(n)
    if len(pieces)>1:return {'kind':'product','value':n,'children':[certificate(a) for a in pieces]}
    digit=2 if n%3==2 else 4
    parent=(n-digit)//3
    return {'kind':'one_shot_gate','value':n,'digit':digit,'source_count':parent,'source':certificate(parent)}

def divisors(n):
    values=[1]
    for p,e in Counter(factor(n)).items():
        values=[d*p**j for d in values for j in range(e+1)]
    return values

@cache
def optimal(n):
    if n in SEEDS:return 0,('seed',)
    if atoms(n)==[n]:
        digit=2 if n%3==2 else 4
        parent=(n-digit)//3
        return optimal(parent)[0]+1,('gate',parent,digit)
    options=[]
    for d in divisors(n):
        q=n//d
        if 1<d<=q and d not in (7,21) and q not in (7,21):
            options.append((optimal(d)[0]+optimal(q)[0],('product',d,q)))
    need(options,'nonatom factor split')
    return min(options)

def optimal_certificate(n):
    cost,choice=optimal(n)
    if choice[0]=='seed':return {'kind':'seed','value':n,'gate_cost':cost}
    if choice[0]=='gate':
        _,parent,digit=choice
        return {'kind':'one_shot_gate','value':n,'gate_cost':cost,'digit':digit,'source_count':parent,
                'source':optimal_certificate(parent)}
    _,d,q=choice
    return {'kind':'product','value':n,'gate_cost':cost,'children':[optimal_certificate(d),optimal_certificate(q)]}

def valid(rows):
    n=len(rows)
    return all(not (rows[i]>>i)&1 and rows[i]<(1<<n)
               and all(((rows[i]>>j)&1)+((rows[j]>>i)&1)==1 for j in range(i))
               for i in range(n))

def hp(rows):
    need(valid(rows),'actual tournament')
    n=len(rows);dp=[[0]*n for _ in range(1<<n)]
    for mask in range(1,1<<n):
        for v in range(n):
            if not (mask>>v)&1:continue
            previous=mask^(1<<v)
            dp[mask][v]=1 if not previous else sum(dp[previous][u] for u in range(n) if (rows[u]>>v)&1)
    return sum(dp[-1])

def strong(rows):
    for start in range(len(rows)):
        seen=1<<start
        while True:
            more=seen
            for v,row in enumerate(rows):
                if seen>>v&1:more|=row
            if more==seen:break
            seen=more
        if seen!=(1<<len(rows))-1:return False
    return True

def apex(word):
    n=len(word)
    return tuple(sum(1<<j for j in range(i+1,n))|((bit=='0')<<n) for i,bit in enumerate(word))+(sum((bit=='1')<<i for i,bit in enumerate(word)),)

def two_apex(word):
    n=len(word)+2;rows=[0]*n;rows[0]|=2
    for i,letter in enumerate(word,2):
        c=int(letter)
        for a in (0,1):
            if (c>>a)&1:rows[a]|=1<<i
            else:rows[i]|=1<<a
        for j in range(i+1,n):rows[i]|=1<<j
    return tuple(rows)

def transitive(n):return tuple(sum(1<<j for j in range(i+1,n)) for i in range(n))

def inflate(frame,root,module):
    slots=[];cursor=0
    for v in range(len(frame)):
        size=len(module) if v==root else 1
        slots.append(tuple(range(cursor,cursor+size)));cursor+=size
    rows=[0]*cursor
    for v in range(len(frame)):
        for w in range(len(frame)):
            if (frame[v]>>w)&1:
                for a in slots[v]:
                    for b in slots[w]:rows[a]|=1<<b
    for i,row in enumerate(module):
        for j in range(len(module)):
            if (row>>j)&1:rows[slots[root][i]]|=1<<slots[root][j]
    return tuple(rows)

def join(a,b):
    n=len(a);later=((1<<len(b))-1)<<n
    return tuple(row|later for row in a)+tuple(row<<n for row in b)

def permutation_hp(rows):
    return sum(all((rows[a]>>b)&1 for a,b in zip(p,p[1:])) for p in permutations(range(len(rows))))

def main():
    maxima={'gates':(-1,None),'depth':(-1,None),'leaves':(-1,None)}
    checked=0;critical=0;improved=0;largest_saving=(0,None)
    for n in range(1,LIMIT+1,2):
        if n in (7,21):continue
        a=atoms(n);product=1
        for value in a:product*=value
        need(product==n,'relative-atom product')
        gates,depth,leaves=forest_stats(n)
        seed_product=1
        for value,multiplicity in leaves:seed_product*=value**multiplicity
        need(3**gates*seed_product<=n,'multiplicative gate budget')
        need(gates==0 or 3**gates*seed_product<n,'strict positive-digit budget')
        need(3**depth<=n,'depth budget')
        optimum=optimal(n)[0]
        need(optimum<=gates and 3**optimum<=n,'optimal forest cost')
        if optimum<gates:
            improved+=1
            if gates-optimum>largest_saving[0]:largest_saving=(gates-optimum,n)
        metrics={'gates':gates,'depth':depth,'leaves':sum(c for value,c in leaves if value!=1)}
        for label,value in metrics.items():
            if value>maxima[label][0]:maxima[label]=(value,n)
        checked+=1
        critical+=int(len(a)==1 and n not in (49,63,343))
    words={3:('one','10'),5:('one','100'),23:('one','11000'),49:('two','10212'),
           63:('one','1001100'),67:('two','12301'),343:('two','1012000')}
    seeds=[{'H':1,'rows':[0],'strong':True}]
    independent=0
    for target,(kind,word) in words.items():
        rows=apex(word) if kind=='one' else two_apex(word)
        need(hp(rows)==target,'inherited literal seed')
        if len(rows)<=7:
            need(permutation_hp(rows)==target,'independent seed permutations');independent+=1
        seeds.append({'H':target,'rows':rows,'word':word,'kind':kind,'strong':strong(rows)})
    frames=[((2,20,9,35,45,7),0),((2,12,9,49,135,215,159,15),4),(apex('10010110'),4)]
    physical=[]
    for m in (1,2,3,7):
        rows=[inflate(frame,root,transitive(m)) for frame,root in frames]
        values=[hp(t) for t in rows]
        h=15*2**m-1
        need(values==[h,3*h+2,3*h+4] and all(strong(t) for t in rows),'incoming actual fork')
        physical.append({'module_order':m,'values':values,'orders':[len(t) for t in rows],
                         'rows':rows if m==7 else None})
    c3=(2,4,1)
    nonunit_rows=[inflate(frame,root,c3) for frame,root in frames]
    nonunit=[hp(t) for t in nonunit_rows]
    need(nonunit==[177,537,543],'scaled-digit hostile')
    h19=apex('10100')
    h101=(250,508,505,240,224,192,128,256,121)
    need(hp(h19)==19 and hp(h101)==101,'actual factor sources')
    canonical=join(h19,h101)
    need(hp(canonical)==1919 and not strong(canonical),'canonical source is reducible')
    eligible=physical[-1]['rows'][0]
    need(strong(eligible),'alternate source is strong')
    need(factor(823)==[823] and physical[-1]['values'][2]==7*823,'critical seven-prime child')
    example=certificate(5761)
    need(forest_stats(5761)==(4,3,((3,2),(5,1))),'selected forest control')
    need(factor(1367)==[1367] and optimal(1367)[0]==3 and forest_stats(1367)[0]==4,'greedy factorization hostile')
    need(optimal(455)==(2,('product',5,91)),'alternative seven-factor allocation')
    h91=inflate(*frames[2],transitive(1))
    h455=join(apex('100'),h91)
    need(hp(h91)==91 and hp(h455)==455,'actual alternative source455')
    sources=[
        '05-knowledge/results/session193_triangulated_proof_20261105/prime_induction/PROOF.md',
        '05-knowledge/results/prime_ternary_lift_20261113/CORE_REPLACEMENTS.md',
        '05-knowledge/results/sixty_three_bridge_20261113/semigroup/PROOF.md',
        '05-knowledge/results/chanpak_five_core_bridges_20261110/history/NEAR_SQUARE_AND_ADJACENT.md']
    pins=[{'path':p,'sha256_lf':hashlib.sha256((ROOT/p).read_bytes().replace(b'\r\n',b'\n')).hexdigest()} for p in sources]
    report={'status':'PASS','scope':'Formal arithmetic forests; physical gates remain conditional except the separately verified finite controls',
            'arithmetic_cap':LIMIT,'allowed_targets':checked,'critical_targets':critical,
            'maximum_symbolic_costs':maxima,'literal_seeds':seeds,'independent_permutation_seed_checks':independent,
            'optimal_factorization':{'improved_targets':improved,'largest_gate_saving':largest_saving,
                                     'certificate_1367':optimal_certificate(1367),
                                     'actual_source455_rows':h455,'H91_rows':h91},
            'incoming_frame_controls':physical,'nonunit_module_hostile':{'module':c3,'values':nonunit,'actual_digits':[6,12]},
            'same_count_source_hostile':{'H':1919,'factor_counts':[19,101],'factor_rows':[h19,h101],
                                         'canonical_rows':canonical,'canonical_strong':False,
                                         'eligible_strong':True,'eligible_rows':eligible,'critical_child':5761},
            'certificate_5761':example,'certificate_5761_statistics':forest_stats(5761),'source_pins':pins}
    suffix='.optimized' if sys.flags.optimize else ''
    (HERE/('receipt'+suffix+'.json')).write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({key:value for key,value in report.items() if key in
                     ('status','arithmetic_cap','allowed_targets','critical_targets','maximum_symbolic_costs','certificate_5761_statistics')}))
    print('Actual source1919: reducible canonical join, strong alternative; strong child5761=7*823.')


if __name__=='__main__':main()

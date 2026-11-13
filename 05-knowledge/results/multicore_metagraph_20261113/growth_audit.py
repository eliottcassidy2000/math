"""Independent rooted deletion-card formula and all literal lift profiles through6."""
from rooted_growth import *
from cycle_index import census

def run():
    d=json.loads(gzip.decompress((HERE/'rooted_growth.json.gz').read_bytes()))
    records={(r['n'],tuple(r['rows'])):r for r in d['records']}
    observed=defaultdict(Counter)
    for t in d['transitions']:observed[t['n'],tuple(t['parent'])][tuple(t['child'])]=t['multiplicity']
    expected=defaultdict(Counter);cards=0
    for rec in d['records']:
        U=tuple(rec['rows']);n=len(U)-1
        if not 1<=n<=6:continue
        T=relabel(U,tuple(range(n)));root=n
        counts=Counter(rooted(T,v)[0] for v in range(n) if U[v]>>root&1 and endpoints(T)[v])
        for parent,I in counts.items():
            a=records[n,parent]['aut_root'];b=rec['aut_root']
            need(I*a%b==0,'card quotient integrality')
            expected[n,parent][U]=I*a//b
        cards+=1
    for key,val in expected.items():need(val==observed[key],'complete rooted deletion-card transfer')
    literal=0;unmarked_hostile=None
    for n in range(1,7):
        for bits in range(1<<((n-1)*(n-2)//2)):
            T=tiling(n,bits);key=rooted(T,n-1)[0]
            profile=Counter(rooted(tuple(cc.child(T,mask)),n)[0] for mask in range(1<<(n-1)))
            need(profile==observed[n,key],'individual tiling lift profile equals rooted transfer')
            literal+=1
    for n in range(1,7):
        grouped=defaultdict(list)
        for (size,T),profile in observed.items():
            if size!=n:continue
            histogram=Counter()
            for U,m in profile.items():histogram[sum(endpoints(U))]+=m
            grouped[ordinary(T)].append((T,sorted(histogram.items())))
        for code,values in grouped.items():
            if any(v[1]!=values[0][1] for v in values[1:]):
                unmarked_hostile=dict(n=n,ordinary_code=code,rooted_rows_and_next_H_histograms=values);break
        if unmarked_hostile:break
    burnside=census(8)
    for r,b in zip(d['levels'],burnside):
        need(r['end_rooted_types']==b['end_rooted'] and r['strong_rooted_types']==b['strong_rooted'],'independent Burnside and SCC census')
    save('growth_audit',dict(scope='Complete rooted-card formula through child7; individual tiling forward profiles through parent6; independent Burnside counts through8',child_cards=cards,literal_presentations=literal,unmarked_hostile=unmarked_hostile))
    print('PASS',cards,'child cards',literal,'individual profiles; unmarked hostile',unmarked_hostile)
if __name__=='__main__':run()

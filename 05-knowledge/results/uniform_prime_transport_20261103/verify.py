"""Exact replay; standard library, checks active under -O."""
from collections import Counter
from math import gcd,comb
from pathlib import Path
import json
from transport import *

HERE=Path(__file__).resolve().parent


def main():
    bank=profiles();rec={};atlas={}
    for name,rows in bank.items():
        c=covers(rows);d=[];g=0
        require(c==permutation_covers(rows),'independent covers '+name)
        for k in c[1:]:g=gcd(g,k);d.append(g)
        atlas[name]=native_atlas(rows)
        rec[name]={'rows':rows,'C':c[1:],'gcd_filtration':d,
                   'slot_polynomial':slot_polynomial(c),
                   'native_frame_counts':atlas[name]['frame_counts'],
                   'one_arm_coefficients':[c[1]+c[2],c[1]-c[2]],
                   'two_arm_coefficients':[2*(c[2]+c[3]),3*c[1]+c[2]-4*c[3],-c[1]-5*c[2]+2*c[3]]}
    require(rec['P7']['gcd_filtration'][-1]==63,'Paley content')
    require(rec['R7']['gcd_filtration'][-1]==7,'circular content')
    require(rec['B13']['gcd_filtration'][:3]==[13,13,1],'third run needed')
    require(rec['B15']['C'][1:]==rec['B11']['C'][1:],'mixed H15 response line')
    require(rec['R5']['C'][1:]!=rec['B15']['C'][1:],'H15 hostile: scalar count is not a module type')
    for name,v in rec.items():
        c=[0]+v['C'];poly=v['slot_polynomial'];content=v['gcd_filtration'][-1]
        for r in range(1,14):
            require(sum(a*r**j for j,a in enumerate(poly))==
                    sum(c[j]*comb(r,j) for j in range(1,min(r,len(c)-1)+1)),
                    'Newton/monomial identity '+name)
        for p in [2,3,5,7,11,13]:
            zero=all(sum(a*r**j for j,a in enumerate(poly))%p==0 for r in range(p))
            require(zero==(content%p==0),'fixed-divisor field criterion '+name)
    # Complete labelled five-vertex score fibre, independently canonicalized.
    fibre=Counter();edges=list(combinations(range(5),2));b13_four_external=Counter()
    quotient_words={j:[] for j in (1,2,3)}
    for j in quotient_words:
        for positions in combinations(range(4+j),j):
            if any(b==a+1 for a,b in zip(positions,positions[1:])):continue
            for exterior in permutations(range(1,5)):
                it=iter(exterior)
                quotient_words[j].append(tuple(0 if i in positions else next(it) for i in range(4+j)))
    for mask in range(1<<10):
        rows=[0]*5
        for j,(a,b) in enumerate(edges):
            if mask>>j&1:rows[a]|=1<<b
            else:rows[b]|=1<<a
        if sorted(r.bit_count() for r in rows)==[1,2,2,2,3]:
            fibre[canonical(rows)]+=1
        if strong(rows):
            expanded=substitute_root(rows,bank['B13']);value=H(expanded)
            bj=[sum(all(rows[a]>>b&1 for a,b in zip(word,word[1:])) for word in quotient_words[j])
                for j in (1,2,3)]
            require(value==sum(a*b for a,b in zip(bj,rec['B13']['C'])),'independent quotient-word contraction')
            require(not prime(value),'no prime with exactly four exterior B13 vertices')
            b13_four_external[value]+=1
    require(set(fibre)=={canonical(bank[x]) for x in ['B11','B13','B15']},'complete mixed-five fibre')
    rewrites={}
    for name,target in [('B11',13),('B11',15),('P7',175),('R7',171),('U7',189)]:
        step=atlas[name]['examples'][target];source=bank[name]
        require(covers(step['rows'])[2:]==rec[name]['C'][1:],'entire higher-run response fixed')
        checks=[]
        for b in (1,2):
            for m in range(1,5):
                out=family_rewrite(source,step,m,b)
                old=cyclic(source,transitive(m),transitive(b))
                old_h,new_h=H(old),H(out['rows'])
                require(old_h==family_value([0]+rec[name]['C'],m,b),'uniform family identity')
                quotient=cyclic((0,),transitive(m),transitive(b))
                require(new_h-old_h==(target-rec[name]['C'][0])*H(quotient),'native affine exchange')
                checks.append({'m':m,'b':b,'H_before':old_h,'H_after':new_h,'quotient_H':H(quotient)})
        rewrites[name+'_'+str(target)]={'certificate':step,'contexts':checks}
    # An inherited fixed-attachment mean carrier cannot transport legal guards.
    p=rewrites['B11_13']['certificate']['path_before']
    i,j=rewrites['B11_13']['certificate']['relocation']
    parity=[{'path_valid':0,'macro_valid':0,'H_sum':0,'ear_sums':[0]*32} for _ in range(2)]
    guard_rows=[]
    for mask in range(8):
        rows=list(bank['B11'])
        for bit,(a,b) in enumerate([(0,1),(0,2),(1,2)]):
            if mask>>bit&1:rows[a]^=1<<b;rows[b]^=1<<a
        group=parity[mask.bit_count()%2];valid=path_ok(rows,p)
        legal=valid and rewrite(rows,p,i,j) is not None
        group['path_valid']+=int(valid);group['macro_valid']+=int(legal)
        group['H_sum']+=H(rows)
        for cut in range(32):group['ear_sums'][cut]+=H(ear(rows,cut))
        guard_rows.append({'mask':mask,'rows':rows,'path_valid':valid,'macro_valid':legal})
    require(parity[0]['H_sum']==parity[1]['H_sum'],'triangle parity ordinary mean')
    require(parity[0]['ear_sums']==parity[1]['ear_sums'],'triangle parity every one-ear mean')
    require([v['path_valid'] for v in parity]==[1,1],'same initial path guard mass')
    require([v['macro_valid'] for v in parity]==[1,0],'selected macro defeats mean carrier')
    ears={}
    for name in ['U7','R7','P7']:
        r=bank[name];values=[H(ear(r,s)) for s in range(128)]
        primes=[(h,s) for s,h in enumerate(values) if prime(h)]
        ears[name]={'gcd':gcd(*values),'histogram':dict(sorted(Counter(values).items())),
                    'prime_witnesses':sorted(primes)}
        if name=='P7':
            require(all(h==189+39*s.bit_count()*(7-s.bit_count()) for s,h in enumerate(values)),'Paley ear law')
    step=atlas['P7']['examples'][175]
    repaired=ear(step['rows'],1)
    require(H(repaired)==409 and prime(409) and strong(repaired),'native Paley-to-prime construction')
    # Minimal number of added vertices when every Paley edge must remain.
    first=ear(bank['P7'],3);positive=[]
    for cut in range(256):
        h=H(ear(first,cut))
        if prime(h):positive.append((h,cut))
    require(bool(positive),'two-ear positive prime witness')
    smallest=min(positive);two=ear(first,smallest[1])
    require(strong(two),'strong two-ear prime witness')
    # Literal paths independent of the subset recurrence, one 9-vertex witness.
    brute=sum(all(two[a]>>b&1 for a,b in zip(p,p[1:])) for p in permutations(range(9)))
    require(brute==smallest[0],'independent prime witness')
    prime_families={}
    for name in ['B11','B13','R7']:
        c=[0]+rec[name]['C'];divisor=7 if name=='R7' else 1
        prime_families[name]=[{'m':m,'H':family_value(c,m),'prime':family_value(c,m)//divisor}
                              for m in range(1,11) if prime(family_value(c,m)//divisor)]
    large_direct=[]
    for name,m in [('B13',3),('B11',10),('R7',8)]:
        rows=cyclic(bank[name],transitive(m),transitive(2));h=H(rows)
        require(h==family_value([0]+rec[name]['C'],m),'direct order17 family check')
        large_direct.append({'name':name,'m':m,'n':len(rows),'H':h})
    # The literal hostile: equal complete cover responses, different group orders inherited.
    tri=(2,4,1);a=ordinal(tri,(0,));b=ordinal((0,),tri)
    host_a=cyclic(a,a,a);host_b=cyclic(a,a,b)
    ca,cb=covers(host_a),covers(host_b)
    require(ca==cb and ca[1]==462915,'81/27 hostile full covers')
    result={'status':'FINITE-EXACT controls; all-size statements proved in THEORY.md',
            'profiles':rec,'rewrites':rewrites,'ears':ears,
            'guard_parity':{'sums_over_four_sources':parity,'sources':guard_rows},
            'complete_five_fibre':[{'rows':r,'H':H(r),'labelled_mass':v} for r,v in sorted(fibre.items())],
            'B13_four_external':{'strong_quotients':sum(b13_four_external.values()),
                                 'H_histogram':dict(sorted(b13_four_external.items()))},
            'large_direct':large_direct,
            'prime409':{'rows':repaired,'base_rewrite':'P7_175','cut':1},
            'two_ear_prime':{'first_cut':3,'second_cut':smallest[1],'H':smallest[0],'rows':two,
                             'scope':'least count among the256 second ears of this fixed first cut only'},
            'two_arm_prime_examples':prime_families,
            'hostile_81_27':{'rows':[host_a,host_b],'H':ca[1],'C':ca[1:],
                             'content':gcd(*ca[1:]),'Aut_orders_inherited':[81,27]}}
    (HERE/'verification.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'profiles':{k:(v['C'][0],v['gcd_filtration']) for k,v in rec.items()},
                      'prime_families':prime_families,'two_ear_prime':smallest,
                      'hostile_content':result['hostile_81_27']['content']},indent=2))


if __name__=='__main__':main()

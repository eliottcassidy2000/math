"""Independent complete word/bridge, native-bank and arbitrary-product referee."""
from pathlib import Path
from fractions import Fraction
from itertools import combinations
from collections import Counter,defaultdict
from math import gcd,lcm,comb
import argparse,hashlib,json,os,shutil,subprocess,sys,tempfile

sys.stdout.reconfigure(encoding='utf-8',newline='\n')
HERE=Path(__file__).resolve().parent
STEM=Path(__file__).stem
PRIMARY='continuing15_20260909_lrc_mandatory_wedges'
EXPECTED_PRIMARY={
    '.md':'0fb26e30697f69c4ce7346a0e21da8bf1ff15cd6aea4b7a706ad532e9c1bd7a8',
    '.py':'695d8ae6f73236db83a0c75500d15ddae28d2819043ed10d8338d38be7c59a00',
    '.cpp':'15e9bba0e4331beff34fa1ffb7a4edbbb087e3838503294e0e6fcd1d443cd96e',
    '.out':'629b81839af568a50433c783c878fdeb4bb25a35a9eea9bd55ea50a94dc2961b',
    '_certificate.json':'6927c5d3c64d5983acd149b675ad4df334e38e90c8952e793be5772b45bda277',
}
GATES=Counter()

def check(ok,label):
    GATES[label]+=1
    if not ok:raise ArithmeticError('independent gate failed: '+label)

def canonical(value):return json.dumps(value,sort_keys=True,separators=(',',':')).encode()
def semantic(value):return hashlib.sha256(canonical(value)).hexdigest()
def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def rows(path):return [list(map(int,l.split())) for l in path.read_text().splitlines()]
def normalized(t,a,b,c):
    g=gcd(gcd(a,b),c)
    return(t//g,a//g,b//g,c//g)

def literal_bridges(word,t,table):
    edges=[]
    for i,j in combinations(range(7),2):
        g=gcd(word[i],word[j]);a,b=sorted((word[i]//g,word[j]//g))
        if table.get((t//g,a,b),0):edges.append((i,j))
    def reach(excluded):
        reached={0}
        while True:
            enlarged=reached|{v for e in edges if e!=excluded for u,v in [e,e[::-1]] if u in reached}
            if enlarged==reached:return reached
            reached=enlarged
    check(len(reach(None))==7,'eligible graph really connected')
    return [list(e) for e in edges if len(reach(e))!=7]

def main():
    filed=HERE.name=='04-computation'
    parser=argparse.ArgumentParser()
    parser.add_argument('--root',type=Path,default=HERE.parent if filed else Path('C:/w/s0905'))
    parser.add_argument('--producer',type=Path,default=HERE.parent/'05-knowledge/results' if filed else Path('C:/w/continuing15_20260909_lrc'))
    args=parser.parse_args();root=args.root.resolve();producer=args.producer.resolve()
    primary=json.loads((producer/(PRIMARY+'_certificate.json')).read_bytes())
    primary_pins={suffix:digest((root/'04-computation' if filed and suffix in ['.py','.cpp'] else producer)/(PRIMARY+suffix)) for suffix in EXPECTED_PRIMARY}
    check(primary_pins==EXPECTED_PRIMARY,'all frozen primary target bytes pinned')
    old_path=root/'05-knowledge/results/continuing14_20260908_lrc_zero_classification_certificate.json'
    check(digest(old_path)=='863d33accf19a5515b3ecd2048aac35d77f30f4e004371c73f5760d345bbb240','complete inherited classification pin')
    old=json.loads(old_path.read_bytes());baseline=old['new_scales']
    check(len(baseline)==6889 and max(baseline)==11934 and semantic(baseline)=='b603ca0186c1ebaa2318782348a7446a16f529543147eb53062f1ba360e8e57b','whole inherited necessary array')
    inherited_audit=root/'05-knowledge/results/continuing14_20260908_lrc_independent_audit_certificate.json'
    check(digest(inherited_audit)=='0d2572099bc7d314ee8fc50ffac32ba1789f4361586a9ed174726430a7631f4f','prior complete independent native acceptance pin')
    raw=root/'04-computation/overnight12_20260906_lrc_decoder_descent_inherited_profiles.json'
    check(digest(raw)=='935f3f687b6d7c89cc099e536f536238fd753bcc4c1747906d213cef387ca93f','all inherited profiles pinned')
    profile=json.loads(raw.read_bytes())['levels'];alphabet=profile['6']['gcds']
    clocks=[r for r in old['clocks'] if r[5]];ids=sorted({r[1] for r in clocks})
    check([r[0] for r in clocks]==[t for t in baseline if t%7==0]==primary['declared_zero_clocks'],'whole 356-clock input, not outcome-selected')
    check(len(clocks)==356 and len(ids)==121,'complete clock and domain counts')
    table={(n,a,b):opened for n,a,b,bits,closed,opened in old['quotient_zero_pairs']}
    for t,d,*_ in clocks:
        check(old['domains'][d]['domain']==[a for a in alphabet if t%a==0],'complete allowed divisor alphabet')
    pr=[(int(k),c,w) for k,lev in profile.items() for c,w in lev['profiles']]
    inp=[str(len(pr))]+[' '.join(map(str,[k,c,*w])) for k,c,w in pr]
    inp += [str(len(ids))]+[' '.join(map(str,[d,len(old['domains'][d]['domain']),*old['domains'][d]['domain']])) for d in ids]
    inp += [str(len(clocks))]+[f'{r[0]} {r[1]}' for r in clocks]
    inp += [str(len(table))]+[' '.join(map(str,[*key,count])) for key,count in sorted(table.items())]
    compiler=shutil.which('g++') or 'C:/Users/Eliott/scoop/apps/gcc/current/bin/g++.exe'
    env=os.environ.copy();env['PATH']=str(Path(compiler).resolve().parent)+os.pathsep+env.get('PATH','')
    native_transcripts={};domain_records=[];independent_products=[]
    with tempfile.TemporaryDirectory(prefix=STEM+'_') as directory:
        work=Path(directory).resolve()
        check(work.parent==Path(tempfile.gettempdir()).resolve(),'private outputs within explicit temporary parent')
        exe=work/('audit.exe' if os.name=='nt' else 'audit')
        opts=['-std=c++17','-O3','-DNDEBUG'] if sys.flags.optimize else ['-std=c++17','-O2']
        build=subprocess.run([compiler,*opts,str(HERE/(STEM+'.cpp')),'-o',str(exe)],capture_output=True)
        check(build.returncode==0,'independent native build '+build.stderr.decode(errors='replace'))
        def run(mode,data):
            path=work/(mode+'_input.txt');path.write_bytes(('\n'.join(data)+'\n').encode())
            result=subprocess.run([str(exe),mode,str(path),str(work)],capture_output=True,env=env)
            check(result.returncode==0,'independent '+mode+' execution '+result.stderr.decode(errors='replace'))
            check(b'\r' not in result.stdout,'raw native LF output '+mode)
            native_transcripts[mode]=result.stdout.decode()
        run('words',inp)
        total_raw=total_accepted=0
        for d,raw_count,count in rows(work/'domains.txt'):
            data=(work/f'words_{d}.json').read_bytes();expected=old['domains'][d]
            check(raw_count==comb(len(expected['domain'])+6,7)==expected['unpruned_multisets'],'every unpruned multiset independently visited')
            check(count==expected['word_count'] and hashlib.sha256(data).hexdigest()==expected['words_sha256'],'entire accepted word bank independently regenerated')
            total_raw+=raw_count;total_accepted+=count
            domain_records.append({'index':d,'raw_multisets':raw_count,'word_count':count,'sha256':hashlib.sha256(data).hexdigest()})
        check((total_raw,total_accepted)==(7105596,291811),'full independent profile totals')
        residuals=defaultdict(list);word_keys={};singles=Counter();bounded=Counter()
        for values in rows(work/'residuals.txt'):
            t=values[0];word=values[1:8];j=8;n=values[j];j+=1
            si=[tuple(values[k:k+3]) for k in range(j,j+3*n,3)];j+=3*n;n=values[j];j+=1
            bo=[tuple(values[k:k+3]) for k in range(j,j+3*n,3)]
            check(j+3*n==len(values),'whole native bridge record parsed')
            residuals[t].append(word);word_keys[t,tuple(word)]=(si,bo)
            for key in si:singles[t,*key]+=1
            for key in bo:bounded[t,*key]+=1
        check(sum(map(len,residuals.values()))==47568,'all complete inherited residual words')
        census=rows(work/'clocks.txt')
        check(sum(r[2] for r in census)==353571,'all independent word-clock evaluations')
        for t,d,total,rc,se,ce in census:
            prior=next(r for r in clocks if r[0]==t)
            check((d,total,rc)==(prior[1],prior[2],prior[5]),'complete per-clock inherited predicates')
        check(sum(r[4] for r in census)==1533 and sum(r[5] for r in census)==5424,'complete eligibility counts')
        single_keys=sorted({normalized(*k) for k in singles});bounded_keys=sorted({normalized(*k) for k in bounded})
        check((len(singles),len(single_keys),len(bounded),len(bounded_keys))==(94,53,561,304),'both entire mandatory-bridge catalogues')
        check([[list(k),v] for k,v in sorted(singles.items())]==primary['singleton_signatures'],'all singleton clock-labelled multiplicities')
        check([[list(k),v] for k,v in sorted(bounded.items())]==primary['bounded_signatures'],'all bounded clock-labelled multiplicities')
        check([list(k) for k in single_keys]==primary['singleton_normalized_keys'] and [list(k) for k in bounded_keys]==primary['bounded_normalized_keys'],'all primitive normalization keys')
        check(len(primary['full_clock_topology_survey'])==356 and {r['t'] for r in primary['full_clock_topology_survey']}==set(residuals),'all original topology rows present once')
        for row in primary['full_clock_topology_survey']:
            t=row['t'];check(row['residual_sha256']==semantic(residuals[t]),'whole clock residual digest')
            actual_eligible={tuple(w):(si,bo) for w in residuals[t] for si,bo in [word_keys[t,tuple(w)]] if si or bo}
            check(set(actual_eligible)=={tuple(w['word']) for w in row['eligible']},'complete eligible original words')
            for w in row['eligible']:
                word=w['word'];si,bo=actual_eligible[tuple(word)]
                check(si==[tuple(k) for k,p in w['singleton']] and bo==[tuple(k) for k,p in w['bounded']],'literal edge deletion reconstructs all wedge signatures')
                bridge=literal_bridges(word,t,table)
                check(bridge==w['bridges'],'every claimed bridge is necessary by actual edge deletion')
                for name in ['singleton','bounded']:
                    for key,positions in w[name]:
                        i,j,k=positions
                        check(len(set(positions))==3 and [word[i],word[j],word[k]]==key,'retained original wedge owner labels')
                        check(sorted((i,j)) in bridge and sorted((j,k)) in bridge,'both owner edges are mandatory bridges')

        # Replay the declared arithmetic union, then independently verify that
        # the adaptive policy indeed closes6804 and first stops at5376.
        arithmetic_keys=sorted(set(single_keys)|{normalized(*k) for k in bounded if k[0] in [6804,5376]})
        check(len(arithmetic_keys)==64,'exact declared arithmetic union')
        arms=sorted({(t,a,b) for t,a,b,c in arithmetic_keys}|{(t,c,b) for t,a,b,c in arithmetic_keys})
        run('banks',[str(len(arms))]+[' '.join(map(str,k)) for k in arms])
        banks=defaultdict(set)
        for T,A,B,p,q in rows(work/'banks.txt'):
            check(gcd(p,q)==1,'every native ratio primitive')
            banks[T,A,B].add(Fraction(p,q))
        expected_banks={tuple(r['key']):r for r in primary['native_banks']}
        check(len(primary['native_banks'])==91 and set(expected_banks)==set(arms),'complete native bank requests, no extra or missing arm')
        typed=0
        for T,A,B,native_count,zero_count in rows(work/'bank_counts.txt'):
            key=(T,A,B);e=gcd(A,B);expected=expected_banks[key]
            check(expected['quotient_clock']==T//e,'actual native arm quotient')
            check(expected['typed_native_pairs']==native_count,'all compatible atlas ratios independently evaluated without skips')
            check(len(banks[key])==zero_count==table.get((T//e,min(A,B)//e,max(A,B)//e),0),'full oriented bank cardinality and inherited comparison')
            check(list(map(str,sorted(banks[key])))==expected['ratios'],'every oriented zero ratio, not only aggregate count')
            typed+=native_count
        check(typed==57438==primary['typed_native_scans'] and len(arms)==91,'entire independent native arm-bank universe')
        pending=[];endpoint_keys=set()
        for T,A,B,C in arithmetic_keys:
            products=[]
            for left in sorted(banks[T,A,B]):
                for right in sorted(banks[T,C,B]):
                    denominator=lcm(left.denominator,right.denominator)
                    raw=[int(denominator*left),denominator,int(denominator*right)];g=gcd(gcd(raw[0],raw[1]),raw[2]);triple=[a//g for a in raw]
                    actual=[gcd(T,a) for a in triple];endpoint=left/right;capkey=None
                    if len(set(triple))<3:kind='duplicate'
                    elif actual!=[A,B,C]:kind='depth'
                    else:
                        p,q=sorted((endpoint.numerator,endpoint.denominator));capkey=(T//gcd(A,C),p,q);endpoint_keys.add(capkey);kind='pending'
                    products.append({'r':str(left),'s':str(right),'primitive_triple':triple,'actual_margins':actual,'endpoint_ratio':str(endpoint),'endpoint_quotient':capkey[0] if capkey else None,'kind':kind,'capacity_key':capkey})
            pending.append({'key':[T,A,B,C],'banks':[list(map(str,sorted(banks[T,A,B]))),list(map(str,sorted(banks[T,C,B])))],'products':products})
        check(len(endpoint_keys)==2408 and sum(len(r['products']) for r in pending)==3055,'every unpruned ratio Cartesian product and unique endpoint capacity')
        check(all(28*n*p*q<2**63 for n,p,q in endpoint_keys),'all signed interval products fit the declared int64 range; speed products use int128')
        run('capacities',[str(len(endpoint_keys))]+[' '.join(map(str,k)) for k in sorted(endpoint_keys)])
        capacities={tuple(r[:3]):r[3:5] for r in rows(work/'capacities.txt')}
        allstats=Counter()
        for row in pending:
            stats=Counter()
            for product in row['products']:
                capkey=product.pop('capacity_key');cap=capacities[capkey] if capkey else None
                product['capacity']=cap
                if cap is not None:
                    check(bool(cap[0])==bool(cap[1]),'no chamber-only composite positive in declared bank')
                    product['kind']='positive' if cap[0]>0 else 'residual'
                stats[product['kind']]+=1
                triple=product['primitive_triple']
                check(gcd(gcd(*triple[:2]),triple[2])==1 and Fraction(triple[0],triple[1])==Fraction(product['r']) and Fraction(triple[2],triple[1])==Fraction(product['s']),'exact primitive triple and both ratio orientations')
            row['counts']=dict(sorted(stats.items()));row['closed']=not stats['residual'];allstats.update(stats)
            independent_products.append(row)
        check(independent_products==primary['evaluated_signatures'],'all independent products, clipped depths, duplicates and BOTH endpoint minima')
        check(dict(allstats)=={'positive':2158,'depth':630,'duplicate':9,'residual':258},'entire four-way product classification')
        proved={tuple(r['key']) for r in independent_products if r['closed']}
        check(len(proved)==37 and [list(k) for k in sorted(proved)]==primary['closed_primitive_signatures'],'all closed arithmetic signatures')
        initial_proved=proved&set(single_keys)
        def surviving(t,closed,which):
            return [w for w in residuals[t] if not any(normalized(t,*k) in closed for k in word_keys[t,tuple(w)][which])]
        singleton_closed=[t for t in sorted(residuals) if not surviving(t,initial_proved,0)]
        check(singleton_closed==[5880,7056]==primary['singleton_clock_closures'],'complete singleton stage clock closures')
        check(len(initial_proved)==33 and sum(len(r['products']) for r in independent_products if tuple(r['key']) in single_keys)==1980,'complete singleton arithmetic universe')
        candidates=sorted([t for t in residuals if all(word_keys[t,tuple(w)][1] for w in residuals[t]) and t not in singleton_closed],reverse=True)
        check(candidates==[6804,5376,5124]==primary['boundary_declared_order'],'predeclared descending newly wholly eligible candidates')
        visited=set(single_keys);boundary=[]
        for t in candidates:
            local=sorted({normalized(t,*k) for w in residuals[t] for k in word_keys[t,tuple(w)][1]})
            check(set(local)<=set(arithmetic_keys),'only visited boundary arithmetic requested')
            visited.update(local);remaining=surviving(t,proved&visited,1)
            boundary.append({'t':t,'complete_keys':[list(k) for k in local],'remaining_words':remaining})
            if remaining:break
        check(boundary==primary['boundary_results'] and [r['t'] for r in boundary]==[6804,5376],'whole boundary results and first-residual stopping policy')
        check(boundary[-1]['remaining_words']==[[3,3,6,8,12,16,16]],'exact stopping original word')
        untested={normalized(5124,*k) for w in residuals[5124] for k in word_keys[5124,tuple(w)][1]}
        check(not untested&set(arithmetic_keys) and primary['arithmetic_untested']==[5124],'5124 really remains arithmetically untested')
        changes=[];clock_results=[];removed_set=set()
        for t in sorted(residuals):
            remaining=surviving(t,proved,1);removed={tuple(w) for w in residuals[t]}-{tuple(w) for w in remaining}
            removed_set.update((t,w) for w in removed)
            if not remaining:changes.append(t)
            clock_results.append({'t':t,'old_count':len(residuals[t]),'removed':len(removed),'new_count':len(remaining),'new_residual_sha256':semantic(remaining),'first_remaining':remaining[0] if remaining else None})
        check(clock_results==primary['clock_results'],'all complete surviving clock banks and first owners')
        check(removed_set=={(r['t'],tuple(r['word'])) for r in primary['removed_word_witnesses']},'all removed original words, not just grand totals')
        check(len(removed_set)==1119 and changes==[5880,6804,7056]==primary['new_clock_closures'],'only complete new clock exclusions')
        for witness in primary['removed_word_witnesses']:
            t=witness['t'];word=witness['word'];i,j,k=witness['positions']
            check(normalized(t,word[i],word[j],word[k])==tuple(witness['key']) in proved,'every removal has its correctly oriented proved wedge')
            bridge=literal_bridges(word,t,table)
            check(sorted((i,j)) in bridge and sorted((j,k)) in bridge,'every printed removal owner has two actual mandatory bridges')
        new=[t for t in baseline if t not in changes];zero=[t for t in new if t%7==0]
        check(new==primary['new_scales'] and semantic(new)==primary['new_array_sha256'],'entire new array and semantic hash')
        check((len(new),max(new),len(zero),max(zero))==(6886,11934,353,6552),'exact final necessary frontiers')

    # Literal full seven-speed original-grid hostile, independent of all pair
    # minimizers and of the native floating-free sweep.
    physical=primary['physical_stopping_hostile'];t=physical['clock'];speeds=physical['speeds'];pn,pd=physical['phase'];den=t*pd
    check(t==5376 and len(set(speeds))==7 and gcd(*speeds)==1,'literal physical row is primitive and distinct')
    check([gcd(t,u) for u in speeds]==[3,3,6,8,12,16,16]==physical['margins'],'all actual physical margins')
    dangers=[{j for j in range(t) if 14*min((u*(pd*j+pn))%den,den-(u*(pd*j+pn))%den)<den} for u in speeds]
    def strict_sum(n):
        p=2
        while p*p<=n:
            e=0
            while n%p==0:n//=p;e+=1
            if e and (p%3!=2 or e>2):return False
            p+=1
        return n==1 or n%3==2
    edges=[];allpairs=[]
    for i,j in combinations(range(7),2):
        g=gcd(speeds[i],speeds[j]);p,q=sorted((speeds[i]//g,speeds[j]//g));strict=p+q<=356 and strict_sum(p+q);overlap=len(dangers[i]&dangers[j])
        allpairs.append([i,j,p,q,strict,overlap])
        if strict:edges.append([i,j,p,q]);check(overlap==0,'every actual strict edge zero at one COMMON original phase')
    check(edges==physical['actual_strict_edges'] and allpairs==physical['all_actual_pairs'],'full actual atlas graph and all omitted-pair overlaps')
    reached={0}
    while True:
        enlarged=reached|{v for i,j,p,q in edges for u,v in [(i,j),(j,i)] if u in reached}
        if enlarged==reached:break
        reached=enlarged
    check(len(edges)==6 and reached==set(range(7)),'physical full strict graph is a connected tree')
    check([len(z) for z in dangers]==[768]*7==physical['danger_sizes'],'every physical marginal attained')
    check(t-len(set.union(*dangers))==1600==physical['safe_points'],'common strict-edge zeros still leave exactly1600 safe points')
    check([semantic(sorted(z)) for z in dangers]==physical['danger_set_sha256'],'all seven literal danger sets, not only sizes')
    check(any(not strict and count>0 for i,j,p,q,strict,count in allpairs),'non-atlas actual overlap is retained')
    native_gates=sum(int(next(l.split()[-1] for l in text.splitlines() if l.startswith('NATIVE_ALWAYS_ACTIVE_GATES'))) for text in native_transcripts.values())
    certificate={'status':'INDEPENDENT ANALYTIC ACCEPTANCE + COMPLETE FINITE-EXACT REPLAY',
        'primary_pins':primary_pins,'native_source_sha256':digest(HERE/(STEM+'.cpp')),
        'native_transcripts':native_transcripts,'native_gates':native_gates,'python_gate_counts':dict(sorted(GATES.items())),
        'python_gates':sum(GATES.values()),'total_gates':native_gates+sum(GATES.values()),
        'domain_records':domain_records,'unpruned_multisets':total_raw,'accepted_domain_words':total_accepted,
        'word_clock_evaluations':353571,'whole_residual_words':47568,
        'singleton_clock_signatures':94,'singleton_primitive_signatures':53,
        'bounded_clock_signatures':561,'bounded_primitive_signatures':304,
        'evaluated_signatures':independent_products,'product_counts':dict(sorted(allstats.items())),
        'boundary_declared_order':candidates,'boundary_results':boundary,'arithmetic_untested':[5124],
        'removed_word_count':len(removed_set),'new_clock_closures':changes,'new_scales':new,'new_array_sha256':semantic(new),'remaining_E0':zero,
        'physical_stopping_hostile':physical,'scope':'Primitive thirteen distinct positive speeds; selected-six gcd t; actual connected strict complementary seven. Remaining signatures do not assert a common unsafe phase.'}
    dest=root/'05-knowledge/results' if filed else HERE
    target=dest/(STEM+'_certificate.json');target.write_bytes(canonical(certificate)+b'\n')
    print('PASS: independent complete mandatory-wedge LRC audit')
    for key in ['words','banks','capacities']:print(native_transcripts[key].strip())
    print('PRODUCTS3055: positive2158 depth630 duplicate9 residual258; closed signatures37')
    print('CLOSED_CLOCKS5880,6804,7056; 6886 remaining/max11934; E0 353/max6552')
    print('PHYSICAL5376: every actual strict tree edge zero at one phase; safe points1600')
    print('PYTHON_ALWAYS_ACTIVE_GATES',sum(GATES.values()),'TOTAL_GATES',native_gates+sum(GATES.values()))
    print('CERTIFICATE_SHA256',digest(target))

if __name__=='__main__':main()

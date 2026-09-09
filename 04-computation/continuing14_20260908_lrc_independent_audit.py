"""Independent full-domain and fully unpruned native-ratio LRC referee.

No primary producer import/execution. The C++ companion reconstructs all
geometry from raw strict walls and signed events, and exhausts every word.
"""
from pathlib import Path
from hashlib import sha256
from itertools import combinations
from collections import Counter,defaultdict
from math import gcd,comb
import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile

sys.stdout.reconfigure(encoding='utf-8',newline='\n')
HERE=Path(__file__).resolve().parent
STEM=Path(__file__).stem
PRIMARY='continuing14_20260908_lrc_zero_classification'
GATES=Counter()

def check(ok,label):
    GATES[label]+=1
    if not ok:raise ArithmeticError('always-active audit gate failed: '+label)

def canonical(value):return json.dumps(value,sort_keys=True,separators=(',',':')).encode()
def semantic(value):return sha256(canonical(value)).hexdigest()
def read_rows(path):return [list(map(int,line.split())) for line in path.read_text().splitlines()]

def connected_parts(word,table,t,component):
    parent=list(range(7))
    def root(i):
        while parent[i]!=i:i=parent[i]
        return i
    for i,j in combinations(range(7),2):
        e=gcd(word[i],word[j]);a,b=sorted((word[i]//e,word[j]//e))
        if table.get((t//e,a,b),(0,0,0))[component]:parent[root(i)]=root(j)
    groups=defaultdict(list)
    for i in range(7):groups[root(i)].append(i)
    return sorted(groups.values())

def ratio_count(table,t,a,b):
    e=gcd(a,b);aa,bb=sorted((a//e,b//e))
    return table.get((t//e,aa,bb),(0,0,0))[2]

def first_slot_failure(word,table,t):
    mult=Counter(word)
    for a in sorted(mult):
        if ratio_count(table,t,a,a):continue
        terms=[[b,mult[b],ratio_count(table,t,a,b)] for b in sorted(mult) if b!=a and ratio_count(table,t,a,b)]
        bound=sum(number*slots for b,number,slots in terms)
        if mult[a]>bound:return {'margin':a,'multiplicity':mult[a],'slot_bound':bound,'terms':terms}
    return None

def main():
    filed=HERE.name=='04-computation'
    parser=argparse.ArgumentParser()
    parser.add_argument('--root',type=Path,default=HERE.parent if filed else Path('C:/w/s0905'))
    parser.add_argument('--producer',type=Path,default=HERE.parent/'05-knowledge/results' if filed else Path('C:/w/continuing14_20260908_lrc'))
    args=parser.parse_args();root=args.root.resolve();producer=args.producer.resolve()
    expected_pins={'.py':'77b7edaa50e1a067a554ac330ee41d27b2aa1201797abd37b0091c8c636596e5',
        '.cpp':'dc7b5aec55e727963523cba0f976b9f9caaaa8c8d162ee4e4ad09b83b6a49aaf',
        '.out':'70b45b15b37f545f1d5af826abd6a1a2d77f4e65d9f8b8582a723a610b010b4c',
        '_certificate.json':'863d33accf19a5515b3ecd2048aac35d77f30f4e004371c73f5760d345bbb240'}
    for suffix,digest in expected_pins.items():
        directory=root/'04-computation' if filed and suffix in ('.py','.cpp') else producer
        check(sha256((directory/(PRIMARY+suffix)).read_bytes()).hexdigest()==digest,'frozen primary artifact '+suffix)
    primary=json.loads((producer/(PRIMARY+'_certificate.json')).read_bytes())
    raw=(root/'04-computation/overnight12_20260906_lrc_decoder_descent_inherited_profiles.json').read_bytes()
    check(sha256(raw).hexdigest()=='935f3f687b6d7c89cc099e536f536238fd753bcc4c1747906d213cef387ca93f','complete inherited profile pin')
    profiles=json.loads(raw)['levels'];alphabet=profiles['6']['gcds']
    check(alphabet==sorted({r[0] for r in profiles['6']['profiles']}) and len(alphabet)==42,'full original sheet alphabet')
    check(all(gcd(a,7)==1 for a in alphabet),'every sheet prime to seven')
    baseline_raw=(root/'05-knowledge/results/continuing13_20260908_lrc_zero_cuts_certificate.json').read_bytes()
    check(sha256(baseline_raw).hexdigest()=='23ae81ba6cb39bc14aecf2d93fc3fcf6d4e3c587b5240db35ea1d8f3992551f1','full predecessor certificate')
    baseline=json.loads(baseline_raw)['new_scales']
    check(len(baseline)==7618 and max(baseline)==11935 and semantic(baseline)=='89a3e7545cc77467c5f85fbe0bdaab1f71071cff65184bbc1212eb86c9f07177','complete predecessor necessary array')
    clocks=[t for t in baseline if t%7==0]
    check(len(clocks)==1085 and clocks==primary['declared_clocks'],'whole E0 survey selected before outcomes')
    domains=sorted({tuple(a for a in alphabet if t%a==0) for t in clocks});domain_index={d:i for i,d in enumerate(domains)}
    check(len(domains)==203 and [list(d) for d in domains]==[r['domain'] for r in primary['domains']],'every full divisor alphabet')
    quotients={t//gcd(a,b) for t in clocks for a in alphabet if t%a==0 for b in alphabet if t%b==0}
    check(quotients==set(clocks),'all actual pair quotient clocks covered exactly')
    profile_rows=[(int(k),c,word) for k,level in profiles.items() for c,word in level['profiles']]
    rows=[str(len(profile_rows))]+[' '.join(map(str,[k,c,*word])) for k,c,word in profile_rows]
    rows += [str(len(domains))]+[' '.join(map(str,[len(d),*d])) for d in domains]
    rows += [str(len(clocks))]+[f'{t} {domain_index[tuple(a for a in alphabet if t%a==0)]}' for t in clocks]
    temp_parent=Path(tempfile.gettempdir()).resolve()
    work=Path(tempfile.mkdtemp(prefix=STEM+'_',dir=temp_parent)).resolve()
    check(work.parent==temp_parent and work.name.startswith(STEM+'_'),'private native outputs confined to declared temporary parent')
    input_path=work/'input.txt';input_path.write_bytes(('\n'.join(rows)+'\n').encode())
    compiler=shutil.which('g++') or 'C:/Users/Eliott/scoop/apps/gcc/current/bin/g++.exe'
    executable=work/('audit.exe' if os.name=='nt' else 'audit')
    options=['-std=c++17','-O3','-DNDEBUG'] if sys.flags.optimize else ['-std=c++17','-O2']
    build=subprocess.run([compiler,*options,str(HERE/(STEM+'.cpp')),'-o',str(executable)],capture_output=True)
    check(build.returncode==0,'independent native compilation '+build.stderr.decode(errors='replace'))
    env=os.environ.copy();env['PATH']=str(Path(compiler).resolve().parent)+os.pathsep+env.get('PATH','')
    run=subprocess.run([str(executable),str(input_path),str(work)],capture_output=True,env=env)
    check(run.returncode==0,'independent native survey '+run.stderr.decode(errors='replace'))
    native=run.stdout.replace(b'\r\n',b'\n').decode()
    table_rows=read_rows(work/'zero_pairs.txt')
    check(table_rows==primary['quotient_zero_pairs'] and len(table_rows)==728,'complete unpruned native oriented ratio counts')
    table={tuple(r[:3]):tuple(r[3:]) for r in table_rows}
    raw_zero=read_rows(work/'all_zero_ratios.txt');ratio_sets=defaultdict(lambda:[set(),set()])
    for n,p,q,a,b,bits in raw_zero:
        check(gcd(p,q)==1 and p<q and gcd(n,p)==a and gcd(n,q)==b,'every individual native ratio typed')
        for left,right,numerator,denominator in ((a,b,p,q),(b,a,q,p)):
            ratio_sets[n,left,right][0].add((numerator,denominator))
            if bits&2:ratio_sets[n,left,right][1].add((numerator,denominator))
    for n,a,b,bits,c,o in table_rows:
        check(len(ratio_sets[n,a,b][0])==c and len(ratio_sets[n,a,b][1])==o,'distinct oriented sets reconstruct every slot count')
        check(ratio_sets[n,a,b][1]<=ratio_sets[n,a,b][0],'complete chamber ratio subset')
    critical=read_rows(work/'critical_banks.txt')
    check(critical==primary['critical_native_banks'],'every retained critical native bank from unpruned atlas')
    domain_rows=read_rows(work/'domains.txt');wordbanks=[];raw_count=0;valid_count=0
    for i,count,raw_multisets in domain_rows:
        data=(work/f'words_{i}.json').read_bytes();wordbanks.append(json.loads(data))
        expected=primary['domains'][i]
        check(len(wordbanks[-1])==count==expected['word_count'],'entire independently accepted word-bank length')
        check(raw_multisets==comb(len(domains[i])+6,7)==expected['unpruned_multisets'],'unpruned complete multiset universe')
        check(sha256(data).hexdigest()==expected['words_sha256'],'whole independently regenerated word-bank digest')
        raw_count+=raw_multisets;valid_count+=count
    check(raw_count==9978787 and valid_count==415506,'full exact domain census totals')
    records=read_rows(work/'survey.txt')
    check(records==primary['clocks'],'all 1085 independent clock records, stages and first owners')
    check(sum(r[2] for r in records)==630818,'all complete word-clock evaluations')
    stages=[[r[0] for r in records if r[3]==0],
            [r[0] for r in records if r[3]>0 and r[4]==0],
            [r[0] for r in records if r[4]>0 and r[5]==0]]
    check(list(map(len,stages))==[713,8,8],'three disjoint deletion stages')
    deleted=set();intermediate=[]
    for actual,expected in zip(stages,primary['stages']):
        check(actual==expected['removed_clocks'] and not deleted.intersection(actual),'stage has only its independently proved new removals')
        deleted.update(actual);remaining=[t for t in baseline if t not in deleted]
        check(remaining==expected['new_scales'] and semantic(remaining)==expected['new_scales_sha256'],'complete intermediate array and semantic pin')
        intermediate.append(len(remaining))
    check(intermediate==[6905,6897,6889],'all intermediate counts')
    check(max(remaining)==11934 and semantic(remaining)=='b603ca0186c1ebaa2318782348a7446a16f529543147eb53062f1ba360e8e57b','complete final array and maximum')
    check(remaining==primary['new_scales'],'full final necessary array matches certificate')
    e0=[t for t in remaining if t%7==0]
    check(e0==primary['remaining_E0'] and len(e0)==356 and max(e0)==7056,'whole residual E0 slice')
    refinement_pairs=0
    for j,(T,d,*_) in enumerate(records):
        for t,small_d,*_ in records[:j]:
            if small_d!=d or T%t:continue
            refinement_pairs+=1
            for a in domains[d]:
                for b in domains[d]:
                    e=gcd(a,b)
                    small=ratio_sets[t//e,a//e,b//e];large=ratio_sets[T//e,a//e,b//e]
                    check(large[0]<=small[0] and large[1]<=small[1],'full oriented-set inclusion for every lawful divisor refinement')
    check(refinement_pairs==primary['lawful_divisible_refinement_pairs']==1055,'complete lawful refinement universe')
    refined=[]
    for expected in primary['refined_complete_residuals']:
        t=expected['clock'];d=domain_index[tuple(a for a in alphabet if t%a==0)]
        closed=[];opened=[];slots=[];proofs=[]
        for word in wordbanks[d]:
            parts0=connected_parts(word,table,t,1);parts1=connected_parts(word,table,t,2)
            if len(parts0)==1:closed.append(word)
            if len(parts1)>1:
                if len(parts0)==1:proofs.append({'word':word,'chamber_components':parts1})
                continue
            opened.append(word);failure=first_slot_failure(word,table,t)
            if failure:proofs.append({'word':word,'slot_failure':failure})
            else:slots.append(word)
        actual={'clock':t,'closed_residuals':closed,'chamber_residuals':opened,'slot_residuals':slots,'proofs':proofs}
        check(actual==expected,'entire refined residual bank and every explicit cut or slot witness')
        refined.append(actual)
    stop=read_rows(work/'stop.txt')
    check(stop==[[4,8,8,9,9,18,24],[4,8,9,9,16,18,24],[4,9,9,16,16,18,24]],'complete 7056 stopping words')
    check(len([w for r in refined if r['clock'] in stages[1] for w in r['closed_residuals']])==18,'all 18 ordinary-only refined words')
    check(len([w for r in refined if r['clock'] in stages[2] for w in r['chamber_residuals']])==31,'all 31 slot-only refined words')
    native_gates=int(next(line.split()[-1] for line in native.splitlines() if line.startswith('NATIVE_ALWAYS_ACTIVE_GATES ')))
    certificate={'status':'INDEPENDENTLY ACCEPTED complete scoped E0 classification',
        'reviewed_primary_report_sha256':'0ae9aa1b4d52989b69da202b8981ac41f39075cda8da90d557668512cbd0b1d3',
        'primary_executable_pins':expected_pins,'native_audit_source_sha256':sha256((HERE/(STEM+'.cpp')).read_bytes()).hexdigest(),
        'native_transcript':native,'native_gates':native_gates,'python_gate_counts':dict(sorted(GATES.items())),
        'python_gates':sum(GATES.values()),'total_gates':native_gates+sum(GATES.values()),
        'complete_domains':203,'unpruned_multisets':raw_count,'distinct_valid_words':valid_count,
        'word_clock_evaluations':630818,'all_zero_ratio_rows':len(raw_zero),
        'all_zero_ratios_sha256':sha256((work/'all_zero_ratios.txt').read_bytes()).hexdigest(),
        'quotient_zero_pairs':table_rows,'separate_removed_clocks':stages,'intermediate_counts':intermediate,
        'new_scales':remaining,'new_scales_sha256':semantic(remaining),'remaining_E0':e0,
        'lawful_full_ratio_set_refinements':refinement_pairs,'stop_words':stop,
        'scope':'Distinct positive speeds, selected-six gcd clock, actual connected strict complement. Remaining graphs are not joint ratio/phase realizations.'}
    destination=root/'05-knowledge/results' if filed else HERE
    (destination/(STEM+'_certificate.json')).write_bytes(canonical(certificate)+b'\n')
    print('Independent complete zero-budget clock classification: PASS')
    print(native.rstrip())
    print('All 203 word banks, 728 complete oriented zero-pair rows, 1085 clock records, and 1055 full ratio-set refinements agree')
    print('SEPARATE_REMOVALS 713+8+8; NECESSARY_ARRAY 6889 MAX 11934; E0_ARRAY 356 MAX 7056')
    print('FINAL_ARRAY_SHA256',semantic(remaining))
    print('PYTHON_ALWAYS_ACTIVE_GATES',sum(GATES.values()))
    print('TOTAL_ALWAYS_ACTIVE_GATES',native_gates+sum(GATES.values()))

if __name__=='__main__':main()

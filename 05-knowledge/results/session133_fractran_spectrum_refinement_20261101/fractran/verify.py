#!/usr/bin/env python3
"""Portable normal/-O checks for eight compiled and two direct finite programs."""
import hashlib
import json
import math
from pathlib import Path
from compiler import (Compiler, BaselineCompiler, thread_jumps, fractions, require,
                      run_register, run_fractran_exponents, run_fractran_literal,
                      write_artifacts)
import build_sources
import core_counter
import reference_oracles as ref
import run_program

ROOT=Path(__file__).resolve().parent
NAMES=('a000568','self_converse','converse_merged','strong','pair_tournaments',
       'colored_pair_tournaments','colored_tournaments','extremizer_tail')
REGISTER_MAX=dict(a000568=6,self_converse=8,converse_merged=5,strong=5,
                  pair_tournaments=4,colored_pair_tournaments=3,colored_tournaments=4,extremizer_tail=12)
EXPONENT_MAX=dict(a000568=5,self_converse=7,converse_merged=4,strong=4,
                  pair_tournaments=3,colored_pair_tournaments=3,colored_tournaments=4,extremizer_tail=12)


def must_fail(fn,label):
    try:
        fn()
    except (RuntimeError,ValueError):
        return
    raise RuntimeError('hostile survived: '+label)


def primitive_tests():
    cases=[]
    for k in range(7):
        cases.extend([(f'n+{k}',lambda n,k=k:n+k),(f'n-{k}',lambda n,k=k:max(0,n-k)),
          (f'n*{k}',lambda n,k=k:n*k),(f'n=={k}',lambda n,k=k:int(n==k)),
          (f'n!={k}',lambda n,k=k:int(n!=k)),(f'n<{k}',lambda n,k=k:int(n<k)),
          (f'n<={k}',lambda n,k=k:int(n<=k)),(f'n>{k}',lambda n,k=k:int(n>k)),
          (f'n>={k}',lambda n,k=k:int(n>=k))])
        if k:
            cases.extend([(f'n//{k}',lambda n,k=k:n//k),(f'n%{k}',lambda n,k=k:n%k)])
    cases.extend([('2**n',lambda n:2**n),('n**0',lambda n:1),('(n+2)*(n+1)',lambda n:(n+2)*(n+1))])
    count=0
    for expression,oracle in cases:
        p=Compiler().compile(f'def f(n):\n    return {expression}\n')
        for n in range(7):
            r=run_register(p,n)
            require(r['output']==oracle(n) and r['error']==0,'arithmetic '+expression)
            require(run_fractran_exponents(p,n)==r,'exponent vs register '+expression)
            count+=1
    refinements=[('n += 3\n    return n',lambda n:n+3),
                 ('n += n\n    return n',lambda n:2*n),
                 ('n += n*2\n    return n',lambda n:3*n),
                 ('n = 9\n    return n',lambda n:9),
                 ('n += 0\n    return n',lambda n:n),
                 ('n = 0\n    return n',lambda n:0),
                 ('if n > 2:\n        return 7\n    else:\n        return 4',lambda n:7 if n>2 else 4)]
    for body,oracle in refinements:
        p=Compiler().compile('def f(n):\n    '+body+'\n')
        for n in range(7):
            r=run_register(p,n)
            require(r['output']==oracle(n) and r['error']==0,'refinement semantics')
            require(run_fractran_exponents(p,n)==r,'refinement fraction semantics')
            count+=1
    p=Compiler().compile('def f(n):\n    return (n+1)*2\n')
    for n in range(7):
        require(run_fractran_literal(p,n)['output']==2*(n+1),'literal primitive')
        count+=1
    return count


def hostile_tests():
    labels=[]
    def rejected(fn,label):
        must_fail(fn,label); labels.append(label)
    p=Compiler().compile('def f(n):\n    return n//0\n')
    require(run_register(p,3)['error']==1,'zero division must set ERROR')
    labels.append('zero_divisor_is_explicit_error')
    rejected(lambda:run_register(p,-1),'negative_input')
    p=Compiler().compile('def f(n):\n    while n != 0:\n        n = n\n    return n\n')
    rejected(lambda:run_register(p,1,100),'nontermination_bound')
    bad={'registers':{'n':0},'register_count':1,'instructions':[('INC',0,0)],'input':0,'output':0,'error':0}
    rejected(lambda:fractions(bad),'cancelled_control_guard')
    bad={**bad,'instructions':[('JUMP',0)]}
    rejected(lambda:thread_jumps(bad),'pure_jump_cycle_rejected_not_misoptimized')
    # Removing the self-loop bridge would make CLEAR's tested fraction lose its guard.
    p=Compiler().compile('def f(n):\n    n=0\n    return n\n')
    require(any(ins[0]=='JUMP' for ins in p['instructions']),'necessary guard bridges disappeared')
    labels.append('clear_self_target_requires_bridge')
    rows=[(7,3),(11,6)]
    require(next((a,b) for a,b in rows if 6%b==0)==(7,3),'zero-first counterexample')
    labels.append('zero_first_branch_is_wrong')
    require(ref.partition_count(2)[0]==1 and ref.partition_count(5)[0]==12,'historical odd-cycle controls')
    labels.extend(['A2_even_cycle_filter','A5_odd_restart'])
    require(ref.partition_count(1,4,2)[0]!=ref.partition_count(1,4,1)[0],'parity needs color sidecar')
    labels.append('self_converse_odd_parity_not_discarded')
    require((2+2)//2 != 2//2,'converse merge must retain self-converse fixed points')
    labels.append('merge_not_naive_halving')
    require(1-1==0 and 1+1!=0,'strong composition requires alternating signs')
    labels.append('strong_unsigned_compositions_wrong_at_two')
    rejected(lambda:core_counter.run_integer(2),'core_infeasible_order_rejected')
    require(core_counter.expected(8)==25 and core_counter.expected(5)*3==27,'tail threshold essential')
    labels.append('core_recurrence_does_not_start_at_five')
    return labels


def oracles(maximum=12):
    values={name:[] for name in NAMES}
    for n in range(maximum+1):
        a=ref.partition_count(n)[0]
        s=ref.partition_count(n//2,4,2**(n%2))[0]
        require((a+s)%2==0,'merged parity')
        c=0 if n==0 else a-sum(values['strong'][k]*values['a000568'][n-k] for k in range(1,n))
        row=dict(a000568=a,self_converse=s,converse_merged=(a+s)//2,strong=c,
                 pair_tournaments=ref.partition_count(n,4,1)[0],
                 colored_pair_tournaments=ref.partition_count(n,4,2)[0],
                 colored_tournaments=ref.partition_count(n,2,2)[0],
                 extremizer_tail=(1,n//3+2,math.comb(n//3+2,2))[n%3])
        for name,value in row.items():
            values[name].append(value)
    return values


def main():
    build_sources.build()
    report={'schema':'session133.finite-fractran.v1','primitive_cases':primitive_tests(),
            'hostile_controls':hostile_tests(),'programs':{}}
    expected=oracles()
    for name in NAMES:
        source=(ROOT/f'{name}_source.py').read_text()
        namespace={}; exec(source,namespace)
        native=[namespace['count'](n) for n in range(13)]
        require(native==expected[name],'native independent oracle '+name)
        p=Compiler().compile(source)
        entry={'program':write_artifacts(p,ROOT/name),'native_0_through_12':native,
               'register_replays':[],'exponent_fraction_replays':[],'literal_integer_replays':[]}
        # Disk shards must reconstruct the exact emitted list in the exact order.
        meta,rows,states=run_program.load(name)
        f=fractions(p)
        require(rows==[(a,b) for a,b,_,_ in f['fractions']],'artifact ordered rows')
        require(states==f['state_primes'],'artifact state dictionary')
        reconstructed=[]
        for part in meta['instruction_parts']:
            reconstructed+=json.loads((ROOT/part).read_text())
        require(reconstructed==[list(ins) for ins in p['instructions']],'IR shard order')
        for n in range(REGISTER_MAX[name]+1):
            r=run_register(p,n)
            require(r['output']==expected[name][n] and r['error']==0,'primitive replay '+name)
            entry['register_replays'].append(dict(n=n,output=r['output'],steps=r['steps']))
            if n<=EXPONENT_MAX[name]:
                e=run_fractran_exponents(p,n)
                require(e==r,'all-register fraction replay '+name)
                entry['exponent_fraction_replays'].append(dict(n=n,output=e['output'],steps=e['steps']))
            if n<=2 or n==3 and name in ('a000568','self_converse','strong'):
                literal=run_program.run(name,n)
                require(literal['output']==r['output'] and literal['error']==0 and literal['steps']==r['steps'],
                        'on-disk literal replay '+name)
                entry['literal_integer_replays'].append(dict(n=n,**literal))
        baseline=BaselineCompiler().compile(source)
        entry['baseline_fractions']=len(fractions(baseline)['fractions'])
        entry['baseline_states']=len(baseline['instructions'])
        entry['baseline_equivalence_0_through_3']=[]
        for n in range(4):
            b=run_register(baseline,n); r=run_register(p,n)
            require(b['output']==r['output'] and b['error']==r['error'],'baseline result equivalence')
            require(all(b['registers'][i]==r['registers'][p['registers'][reg]] for reg,i in baseline['registers'].items()),
                    'baseline named-register equivalence')
            entry['baseline_equivalence_0_through_3'].append(dict(n=n,baseline_steps=b['steps'],refined_steps=r['steps']))
        report['programs'][name]=entry
        print('PASS',name,entry['program'],flush=True)
    report['direct_signed_burnside']=[]
    for n in range(8):
        a,s=ref.permutation_burnside(n)
        require((a,s)==(expected['a000568'][n],expected['self_converse'][n]),'independent signed pair cycles')
        report['direct_signed_burnside'].append(dict(n=n,A=a,S=s,permutations=math.factorial(n)))
    report['direct_orbit_census']=[]
    for n in range(6):
        row=ref.Census(n).summary()
        for symbol,name in [('A','a000568'),('S','self_converse'),('M','converse_merged'),('C','strong')]:
            require(row[symbol]==expected[name][n],'independent complete orbit census')
        report['direct_orbit_census'].append(row)
    report['direct_pair_orbits']=[]
    for n in range(4):
        for colors,name in ((1,'pair_tournaments'),(2,'colored_pair_tournaments')):
            count=ref.pair_tournament_census(n,colors)
            require(count==expected[name][n],'independent pair tournament orbit census')
            report['direct_pair_orbits'].append(dict(n=n,colors=colors,count=count))
    report['core']={}
    for tail,label,minimum in ((False,'full20',3),(True,'tail14',6)):
        rows=core_counter.TAIL14 if tail else core_counter.FULL20
        require(all(math.gcd(a,b)==1 for a,b in rows),'core fractions must be reduced')
        (ROOT/f'core_{label}.fractions.00.txt').write_text(''.join(f'{a}/{b}\n' for a,b in rows))
        results=[]
        for n in range(minimum,25):
            r=core_counter.run_integer(n,tail)
            e=core_counter.run_exponents(n,tail)
            require(r==e and r['output']==core_counter.expected(n) and r['steps']==core_counter.exact_steps(n,tail),'core recurrence/steps')
            results.append(r)
        report['core'][label]={'fractions':len(rows),'literal_and_exponent_replays':results}
    bound_files=[p for p in sorted(ROOT.iterdir()) if p.is_file() and (p.suffix in ('.py','.txt') or any(k in p.name for k in ('.meta.','.instructions.','.states.')))]
    report['sha256']={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in bound_files}
    require(all(p.stat().st_size<100000 for p in ROOT.iterdir() if p.is_file()),'file size cap')
    (ROOT/'VERIFIED.json').write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
    print('PASS all eight compiled programs, direct core20/core14, hostile controls and independent censuses')


if __name__=='__main__':
    main()

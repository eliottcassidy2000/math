#!/usr/bin/env python3
"""Own stdlib checks of actual local CA against independent Fraction PRIMEGAME."""
from fractions import Fraction
from hashlib import sha256
from pathlib import Path
from random import Random
from itertools import combinations
import json
from primegame_ca import *

ROOT=Path(__file__).resolve().parent
RACK=tuple(map(Fraction,('17/91','78/85','19/51','23/38','29/33','77/29','95/23','77/19','1/17','11/13','13/11','15/2','1/7','55/1')))
CHECKS=0
MICRO_UPDATES=0

def require(x,why):
    global CHECKS
    CHECKS+=1
    if not x:raise RuntimeError(why)

def integer_next(n):
    for i,f in enumerate(RACK):
        y=n*f
        if y.denominator==1:return y.numerator,i
    raise RuntimeError('unconditional fraction absent')

def canonical(configuration):return ';'.join(f'{x}:{v}' for x,v in sorted(configuration.items()))

def shifted(configuration,k):return {x+k:v for x,v in configuration.items()}

def checked_macro(configuration,compare_full=False,check_micro=True,hash_object=None):
    global MICRO_UPDATES
    initial=observe(configuration);expected,rule=integer_next(initial['integer'])
    touched=DENOMINATOR_SUPPORTS[rule]+NUMERATOR_SUPPORTS[rule]
    expected_time=1+len(touched)+2*sum(initial['exponents'][i] for i in touched)
    initial_spectator={x:v&SPECTATOR for x,v in configuration.items() if v!=ERROR and v&SPECTATOR}
    initial_origin=initial['origin'];limit=10*max(initial['exponents'])+6
    for elapsed in range(1,limit+1):
        previous=configuration;configuration=ca_step(previous);MICRO_UPDATES+=1
        if compare_full:require(configuration==ca_step_full(previous),'sparse/full CA mismatch')
        if check_micro:
            s=observe(configuration,False)
            require(s['origin']==initial_origin and s['head_position']>=initial_origin,'marker/head boundary failure')
            require({x:v&SPECTATOR for x,v in configuration.items() if v!=ERROR and v&SPECTATOR}==initial_spectator,'spectator changed')
        if hash_object is not None:hash_object.update((canonical(configuration)+'\n').encode())
        if ready(configuration):
            result=observe(configuration)
            require(result['integer']==expected,'CA/Fraction macro mismatch')
            require(elapsed==expected_time and elapsed<=limit,'exact macro time/bound mismatch')
            require(configuration==encode(expected,initial_origin),'noncanonical READY encoding')
            return configuration,elapsed,rule,result
    raise RuntimeError('bounded macro failed to return READY')

def run_to_prime_outputs(seed,count):
    cfg=encode(seed);n=seed;macro=0;micro=0;outputs=[];h=sha256();records=[]
    while len(outputs)<count:
        cfg,elapsed,rule,info=checked_macro(cfg,hash_object=h)
        macro+=1;micro+=elapsed;n=info['integer']
        records.append({'macro':macro,'micro':micro,'instruction':rule+1,'integer':n,'duration':elapsed})
        if n>1 and n&(n-1)==0:
            exponent=n.bit_length()-1
            if exponent>=2 and all(exponent%d for d in range(2,1+int(exponent**0.5))):
                outputs.append({'prime_exponent':exponent,'integer':n,'FRACTRAN_updates':macro,'CA_updates':micro})
    return {'seed':seed,'outputs':outputs,'FRACTRAN_updates':macro,'CA_updates':micro,'microconfiguration_sha256':h.hexdigest()},records

# Fully local CA replay through the first five prime outputs.
main,records=run_to_prime_outputs(2,5)
require([r['FRACTRAN_updates'] for r in main['outputs']]==[19,69,281,710,2375],'known prime milestone mismatch')
with (ROOT/'ready_trace.tsv').open('w') as f:
    f.write('FRACTRAN_updates\tCA_updates\tinstruction_1based\tinteger\tmacro_duration\n')
    for r in records:f.write(f"{r['macro']}\t{r['micro']}\t{r['instruction']}\t{r['integer']}\t{r['duration']}\n")

# Exact cycles/joins for requested seeds, and first three prime outputs for16.
seed_results=[]
for seed in (7,17,55,25,15,23):
    cfg=encode(seed);seen={canonical(cfg):(0,0)};sequence=[seed];micro=0;durations=[]
    for macro in range(1,100):
        cfg,t,rule,info=checked_macro(cfg,compare_full=True)
        micro+=t;sequence.append(info['integer']);durations.append(t)
        key=canonical(cfg)
        if key in seen:
            a,b=seen[key]
            seed_results.append({'seed':seed,'sequence':sequence,'outcome':'exact CA cycle','FRACTRAN_period':macro-a,'CA_period':micro-b,'FRACTRAN_transient':a,'CA_transient':b});break
        seen[key]=(macro,micro)
        if seed in (15,23) and info['integer']==2 or seed==15 and info['integer']==4:
            seed_results.append({'seed':seed,'sequence':sequence,'outcome':'exact standard-orbit join','FRACTRAN_updates':macro,'CA_updates':micro});break
    else:raise RuntimeError('seed join/cycle bound exceeded')
seed16,_=run_to_prime_outputs(16,3);seed_results.append(seed16)

# Every single instruction is selected by its own denominator, including B's5 ops.
rule_durations=[]
for i,(_,den) in enumerate(FRACTIONS):
    require(integer_next(den)[1]==i,'denominator not a priority witness')
    cfg,t,j,_=checked_macro(encode(den),compare_full=True)
    require(j==i,'instruction coverage failure');rule_durations.append(t)

# All380 balanced support shapes over the ten registers, including multi-control.
balanced=set(PRIMES)|{p**3 for p in PRIMES}
for p in PRIMES:
    for q,r in combinations((q for q in PRIMES if q!=p),2):balanced.add(p*p*q*r)
require(len(balanced)==380,'balanced-family input count')
multi_control=0
for n in sorted(balanced):
    cfg=encode(n);initial=observe(cfg)
    if sum(initial['exponents'][4:])>1:multi_control+=1
    checked_macro(cfg,compare_full=True)

# Random general exponents, spectators, shifted origins; every macro local.
rng=Random(126)
random_cases=240
for j in range(random_cases):
    exponents=[rng.randrange(5) for _ in PRIMES]
    spectator=rng.choice((1,31,37,31*37,41**2,43*47))
    n=spectator*prod(p**e for p,e in zip(PRIMES,exponents));origin=rng.randrange(-9,10)
    cfg=encode(n,origin)
    require(observe(cfg)['integer']==n,'all-positive-integer codec')
    for _ in range(3):cfg,_,_,_=checked_macro(cfg,compare_full=True)

# Translation-equivariance and literal/full scheduler agreement on microsteps.
cfg=encode(2)
for j in range(600):
    next_cfg=ca_step(cfg)
    require(next_cfg==ca_step_full(cfg),'micro full evaluator mismatch')
    require(ca_step(shifted(cfg,7))==shifted(next_cfg,7),'translation mismatch')
    cfg=next_cfg

# Exhaust every finite head/tape pair. This is not all A^3 triples.
head_cases=0;error_head_cases=0
for head in range(1,HEAD_COUNT+1):
    for data in range(TAPE_SIZE):
        delta=head_delta(data,head);head_cases+=1
        if delta is None:error_head_cases+=1;continue
        nd,nh,direction=delta
        require(0<=nd<TAPE_SIZE and 1<=nh<=HEAD_COUNT and direction in (-1,0,1),'finite head codomain')
        require(nd&(ORIGIN|SPECTATOR)==data&(ORIGIN|SPECTATOR),'persistent tape markers')

# Sparse-map representation normalizes explicit background entries.
for cfg in ({0:0},{-2:0,0:ORIGIN,3:0},{-2:0,0:ERROR,3:0}):
    require(ca_step(cfg)==ca_step_full(cfg),'explicit background normalization')

# Every tape-only center is fixed, independently of its nonhead neighbors.
for data in range(TAPE_SIZE):
    require(local_rule(rng.randrange(TAPE_SIZE),data,rng.randrange(TAPE_SIZE))==data,'head-free triple changed')

# Invalid-local configurations: explicit errors; globally malformed encodings
# still have a total rule, but need not be recognized as malformed locally.
hostiles=[]
def error_case(name,cfg):
    out=ca_step_full(cfg)
    require(ERROR in out.values(),name+' failed to flag')
    require(ca_step(cfg)==out,name+' sparse/full mismatch')
    hostiles.append(name)
error_case('READY without origin',{0:symbol(0,READY)})
dec_seek=HEAD_CODE[('SEEK',0,0)]
error_case('decrement zero at origin',{0:symbol(ORIGIN,dec_seek)})
dec_left=HEAD_CODE[('DECLEFT',0,0)]
error_case('DECLEFT without required1',{0:symbol(ORIGIN,dec_left)})
error_case('adjacent heads',{0:symbol(ORIGIN,READY),1:symbol(0,READY)})
error_case('two heads distance2',{0:symbol(ORIGIN,READY),2:symbol(ORIGIN,READY)})
error_case('absorbing error neighborhood',{0:ERROR})
for invalid in (-1,ALPHABET_SIZE,'x'):
    try:local_rule(0,invalid,0)
    except ValueError:pass
    else:raise RuntimeError('out-of-alphabet input accepted')
    CHECKS+=1
# Random alphabet triples and arbitrary sparse configurations cover error and
# multiple-head branches without materializing alphabet_size**3 triples.
for _ in range(10000):
    triple=tuple(rng.randrange(ALPHABET_SIZE) for _ in range(3))
    out=local_rule(*triple);require(0<=out<ALPHABET_SIZE,'rule outside alphabet')
for _ in range(500):
    cfg={x:rng.randrange(ALPHABET_SIZE) for x in range(-3,4)}
    a=ca_step(cfg);require(a==ca_step_full(cfg),'arbitrary configuration scheduler mismatch')
    require(ca_step(shifted(cfg,-5))==shifted(a,-5),'arbitrary configuration translation mismatch')

report={'status':'PASS','alphabet':{'tape_symbols':TAPE_SIZE,'head_states':HEAD_COUNT,'plain_or_head_options':HEAD_COUNT+1,'error_symbols':1,'total_symbols':ALPHABET_SIZE,'radius':1,'all_triples_materialized':False},'operation_counts':[len(x) for x in OPS],'exact_macro_time':'1+k+2*sum(input exponent over touched coordinates)','uniform_macro_bound':'10*M+6','seed2_prime_outputs':main,'requested_seed_results':seed_results,'instruction_witness_macro_durations':rule_durations,'balanced_inputs_checked':len(balanced),'balanced_multi_control_inputs':multi_control,'random_general_inputs':random_cases,'random_general_macro_updates':3*random_cases,'all_head_tape_pairs_checked':head_cases,'locally_invalid_head_tape_pairs':error_head_cases,'random_alphabet_triples':10000,'arbitrary_sparse_configurations':500,'hostile_tests':hostiles,'micro_updates_executed':MICRO_UPDATES,'explicit_require_checks':CHECKS,'ready_trace_sha256':sha256((ROOT/'ready_trace.tsv').read_bytes()).hexdigest(),'scope':'Actual finite-alphabet homogeneous radius1 rule. Exact integer dynamics only at READY times; not one-step conjugacy, a ten-cell finite CA, or a tournament metagraph map.'}
(ROOT/'ca_results.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))

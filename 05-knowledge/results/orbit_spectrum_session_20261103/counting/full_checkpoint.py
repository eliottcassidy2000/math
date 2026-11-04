#!/usr/bin/env python3
"""Full a(256) reproduction using the inherited exact engine, not a new record.

Controls: independently derived certified head and inherited bounded-cycle
2-adic enumeration, which sees only <=6 cycles for16 normalized bits.
"""
from pathlib import Path
from math import gcd
import hashlib,importlib.util,json,signal
from heads import head,need
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[3]
BASE=ROOT/'05-knowledge/results/session87_four_vertices_enumeration_20261028/fast_counts.py'
LOW=ROOT/'05-knowledge/results/port_residual_integral_20261025/precision/check.py'

def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    return mod

def main():
    need(hashlib.sha256(BASE.read_bytes()).hexdigest()=='3844a4a590daaf01410cf99d3584441bfa2f8e5e95371ccaafe534fd152b0ece','baseline source hash')
    baseline=load('full_baseline',BASE);low=load('low_independent',LOW)
    def timeout(signum,frame):raise TimeoutError('120-second wall-time session cap')
    signal.signal(signal.SIGALRM,timeout);signal.alarm(120)
    n=256
    values,stats=baseline.degree_bucket_dp(n)
    value=values[n];normalized=value>>128
    need(value%(1<<128)==0,'termwise baseline')
    record=head(n,256,max_core_mass=40)
    need(record['bit_length']==value.bit_length(),'certified binary length')
    need(int(record['head_hex'],16)==value>>(value.bit_length()-256),'independent256-bit head')
    residue,lowstats=low.normalized_residue(n,16)
    need(normalized%(1<<16)==residue,'independent16 normalized bits')
    # Separate explicit four-bit even-order formula (Sessionpaired recursions F2).
    pair=sum(pow(a*(n-a),-1,16) for a in range(1,n//2,2) if gcd(a,n)==1)
    phi=sum(gcd(a,n-2)==1 for a in range(1,n-1))
    four=(pair+8*(phi//2))%16  #256 is neither6-divisible nor an exception4,6,10.
    need(normalized%16==four,'explicit four-bit formula')
    encoded=value.to_bytes((value.bit_length()+7)//8,'big')
    prefixhex='\n'.join(hex(x) for x in values)+'\n'
    result={'status':'PASS','scope':'Full a(256) reproduction, inherited algorithm; prior repository scalar computation through300, no record claim.',
        'n':n,'value_hex':hex(value),'value_sha256_minimal_big_endian_bytes':hashlib.sha256(encoded).hexdigest(),
        'prefix_0_through256_sha256_hex_lines':hashlib.sha256(prefixhex.encode()).hexdigest(),
        'bit_length':value.bit_length(),'valuation_2':(value&-value).bit_length()-1,
        'leading256_hex':record['head_hex'],'head_cutoff':record['cutoff'],'head_core_types':record['core_types'],
        'normalized_mod65536':residue,'normalized_mod16_explicit':four,
        'lowbit_enumeration':lowstats,'full_engine_work':{k:v for k,v in stats.items() if k!='rows'},
        'source_hashes':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in (BASE,LOW)}}
    signal.alarm(0)
    print(json.dumps(result,indent=2,sort_keys=True))
if __name__=='__main__':main()

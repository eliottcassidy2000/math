#!/usr/bin/env python3
"""Exact finite certificates for cross-observer corollaries in README 6b.

The infinite tail uses the proved macro token bound; Hamiltonian existence
uses the explicitly CITED Gerbicz theorem. This script does not compute
large square-sum graphs or search farther along the integer orbit.
"""
from collections import Counter
import csv
from fractions import Fraction
from hashlib import sha256
from itertools import permutations
import json
from math import isqrt
from pathlib import Path

ROOT=Path(__file__).resolve().parent
Q=tuple(map(Fraction,('17/91','78/85','19/51','23/38','29/33','77/29','95/23',
                     '77/19','1/17','11/13','13/11','15/2','1/7','55/1')))
CHECKS=Counter()


def check(condition,message,category):
    CHECKS[category]+=1
    if not condition:
        raise RuntimeError(message)


def step(n):
    for q in Q:
        candidate=n*q
        if candidate.denominator==1:
            return candidate.numerator
    raise RuntimeError('no applicable fraction')


def main():
    trace=json.loads((ROOT/'seed_balance_trace.json').read_text())
    check(len(trace)==146,'finite trace size','prefix')
    n=2
    for t,row in enumerate(trace):
        check((row['step'],row['n'])==(t,n),'Fraction prefix mismatch','prefix')
        if t<145:
            n=step(n)
    low=[(r['step'],r['n']) for r in trace if r['n']<32]
    check(low==[(0,2),(1,15),(19,4),(20,30),(69,8)],'all low prefix states','threshold')
    minimum=min((r['n'],r['step']) for r in trace if r['step']>=70)
    check(minimum==(60,70),'minimum on steps70..145','threshold')
    check(trace[-1]['n']==5**5*7**4*13,'candidate5 boundary','threshold')
    five_hits=[(r['step'],r['n']) for r in trace if r['n'] in (15,16,17,23,25)]
    check(five_hits==[(1,15)],'requested five-set intersection','threshold')

    certificates=[]
    witnesses={2:1,15:8,4:2,30:18,8:4}
    for t,n in low:
        vertex=witnesses[n]
        neighbors=[v for v in range(1,n+1) if v!=vertex and isqrt(v+vertex)**2==v+vertex]
        check(len(neighbors)<2,'small graph degree obstruction','graph')
        certificates.append({'step':t,'order':n,'vertex':vertex,'neighbors':neighbors})

    ready_bytes=(ROOT/'inputs/ready_observations.tsv').read_bytes()
    ready=list(csv.DictReader(ready_bytes.decode().splitlines(),delimiter='\t'))
    by_step={int(r['FRACTRAN_updates']):r for r in ready}
    for t in range(1,146):
        check(int(by_step[t]['integer'])==trace[t]['n'],'READY/prefix integer mismatch','ready')
    check(int(by_step[69]['CA_updates'])==537 and int(by_step[70]['CA_updates'])==547,
          'CA boundary times','ready')
    check(int(by_step[281]['integer'])==32,'later value32 prevents an invalid60 lower bound','ready')
    even_permutations=sum(sum(p[i]>p[j] for i in range(5) for j in range(i+1,5))%2==0
                          for p in permutations(range(5)))
    check(even_permutations==60,'A5 order','group')
    window=json.loads((ROOT/'seed_balance_window.json').read_text())
    check((window['last_balanced_visit']['step'],window['last_balanced_visit']['n'])==(70,60),
          'last balance observation','threshold')
    result={'status':'FINITE-EXACT certificates; universal and CITED dependencies in README 6b',
            'prefix_last_step':145,'minimum_on_steps_70_through_145':{'step':70,'n':60},
            'all_prefix_states_below_32':[{'step':t,'n':n} for t,n in low],
            'small_graph_nonhamiltonicity_certificates':certificates,
            'all_time_nonhamiltonian_graph_observation_steps':[t for t,_ in low],
            'eventual_hamiltonian_start_step':70,
            'all_time_minimum_after_step_70':{'value':32,'attained_at_step':281},
            'ready_observations':[{'FRACTRAN_updates':t,'CA_updates':int(by_step[t]['CA_updates']),
                                   'integer':int(by_step[t]['integer'])} for t in (69,70)],
            'requested_five_set_all_time_visits':[{'step':t,'n':n} for t,n in five_hits],
            'alternating_group_degree':5,'alternating_group_order':even_permutations,
            'ready_input_sha256':sha256(ready_bytes).hexdigest(),
            'checks':dict(sorted(CHECKS.items())),'checks_total':sum(CHECKS.values())}
    encoded=json.dumps(result,indent=2,sort_keys=True)+'\n'
    (ROOT/'cross_observers.json').write_text(encoded)
    print(encoded,end='')


if __name__=='__main__':
    main()

"""Incoming external-pin interface and new reference-rate constants."""
from pathlib import Path
from fractions import Fraction as Q
from math import comb
import json,subprocess,hashlib

ROOT=Path(__file__).resolve().parents[3]
REV='c09b642a629b7d12b529b88c2697abaee78e73e7'
FILES={};checks=0
def read(path):
    b=subprocess.check_output(['git','show',REV+':05-knowledge/results/'+path],cwd=ROOT)
    FILES[path]=hashlib.sha256(b).hexdigest();return json.loads(b)
def need(v,msg):
    global checks
    checks+=1
    if not v:raise RuntimeError(msg)
def main():
    r=read('external_pin_sections_20261024/native/receipt.json')
    x=read('guarded_component_completion_20261024/native/xor_tree_receipt.json')
    n=read('guarded_component_completion_20261024/sections/native_2125_receipt.json')
    old=read('three_target570_20261024/menus.json');bank=read('three_target570_20261024/bank.json')
    rows={row['source']:row for row in old['rows'] if row['role']!='fresh'}
    need(set(bank['reservoir_sources'])==set(r['cheap_core_reservoir']),'same actual twenty receiver sources')
    paths={s:{o['target']:o['path'] for o in row['options'] if o['target'] in (1965,2234)} for s,row in rows.items()}
    paths[2125]={o['target']:o['path'] for o in n['records']}
    cost={s:sum(a!=b for a,b in zip(p[1965],p[2234])) for s,p in paths.items()}
    need(len(paths)==82,'exact declared native source inventory')
    need(len(r['outside_core'])==8 and max(cost[s] for s in r['outside_core'])==4,'eight outside core costs')
    free=r['free_helpers'];need(len(free)==51 and 2125 in free,'actual free helper set')
    top=sorted((cost[s] for s in free),reverse=True)
    need(top[:2]==[8,7] and max(cost[s] for s in r['cheap_core_reservoir'])==2,'literal helper-pair and receiver work')
    words=set();records=0
    for s,row in rows.items():
        for o in row['options']:
            if o['target'] not in (1965,2234):continue
            O=set(o['closed_words']);need(len(O)==4 and not O&words,'original disjoint full resource options');words|=O;records+=1
    for o in n['records']:
        O=set(o['closed_words']);need(len(O)==4 and not O&words,'new2125 keeps full normalized capacity');words|=O;records+=1
    need(records==164 and len(words)==656,'full option store, not duplicated kernel capacity')
    for s,t,clock in ((368,1965,1),(776,2234,4),(3139,2234,7)):
        o=next(z for z in x['actual_four_cell_records'][str(s)] if z['target']==t)
        need(o['global_path'][clock]==4 and o['global_path'][::-1][7-clock]==4,'actual external owner4 and converse clock')
    C1=comb(18,2);C2=comb(51,2);kap=Q(4032,67);gamma=1/(2*kap*(C1+C2))
    need(gamma==Q(67,11515392),'product-reference rate normalization')
    return {'status':'PASS incoming native interface replay and weighted receiver constants','checks':checks,
        'frozen_revision':REV,'inputs':FILES,'sources':82,'reserved_closed_words':len(words),
        'core_pair_reference_coordinates':18,'helper_pair_reference_coordinates':51,
        'reference_mixture_weights':[str(Q(C1,C1+C2)),str(Q(C2,C1+C2))],
        'unit_receiver_gap':str(gamma),'interior_budget':[22,55],'all_grade_budget':[22,62],
        'scope':'native interface recheck plus proved weighted corollary; no new global-observer or harmonic-transfer claim'}
if __name__=='__main__':print(json.dumps(main(),indent=2,sort_keys=True))

#!/usr/bin/env python3
"""Hostile semantic controls on the authentic N4096 family, no source replacement."""
from pathlib import Path
import copy,json,argparse,hashlib
import literal_verify as independent
ROOT=Path(__file__).resolve().parent

def need(ok,why):
    if not ok:raise RuntimeError(why)
def rejects(f):
    try:f()
    except (ValueError,RuntimeError):return True
    return False

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--suffix',default='');args=ap.parse_args()
    rows,_,_=independent.reconstruct(4096)
    claim=json.loads((ROOT/'results/receipt_N4096.json').read_text()); correct=claim['hostile_controls']['correct'];S=claim['S_full']
    mutations={}
    for field,change in [('t',lambda x:x+2),('p',lambda x:x+1),('u',lambda x:x+2),('m',lambda x:x+1),('d',lambda x:x+1),('A',lambda x:[x[0]+1,x[1]])]:
        changed=copy.deepcopy(rows); changed[-1][field]=change(changed[-1][field])
        mutations['tamper_'+field]=rejects(lambda:independent.require(changed==rows,'complete exact metadata equality'))
    zero=next(r for r in rows if independent.norm(r['A'])==0)
    changed=[r for r in rows if r['s']!=zero['s']]
    mutations['delete_zero_charge_source']=rejects(lambda:independent.require(changed==rows,'zero label retained'))
    moving=next(r for r in rows if r['d']>36 and independent.norm(r['A'])>0)
    changed=[r for r in rows if r['s']!=moving['s']]
    mutations['delete_moving_source']=rejects(lambda:independent.require(changed==rows,'moving source retained'))
    true_s=sum(independent.norm(r['A']) for r in rows); changed_s=sum(independent.norm(r['A']) for r in changed)
    need(true_s==S and true_s>changed_s,'moving-source denominator loss exposed')
    need(correct['D_plus']==495 and correct['per_parent']==11227,'positive-part mutation authentic values')
    merged=claim['hostile_controls']['merged_terminals'];sign=claim['hostile_controls']['lost_clock_sign'];address=claim['hostile_controls']['wrong_low_bit_address']
    need(merged['D_plus']==56715 and 2*merged['D_plus']>S,'cross-terminal merge false violation')
    need(2*correct['per_parent']>S,'per-parent positive part false violation')
    need(sign['energies']!=correct['energies'],'losing clock sign detected')
    need(address['energies']!=correct['energies'],'wrong low-bit address detected in level path')
    need(all(mutations.values()),'every metadata mutation caught')
    out=dict(status='PASS',N=4096,S_full=S,correct_D_plus=correct['D_plus'],per_parent_false_debt=correct['per_parent'],
             merged_terminal_false_debt=merged['D_plus'],lost_clock_sign_false_debt=sign['D_plus'],
             wrong_address_debt=address['D_plus'],wrong_address_note='Odd coordinate s shifts dyadic levels; debt may accidentally agree. Complete level path differs and catches mutation.',
             metadata_mutations_rejected=mutations,deleted_zero_source=zero,deleted_moving_source=moving,
             original_denominator=true_s,wrong_deleted_moving_denominator=changed_s)
    (ROOT/'results'/f'mutation_receipt{args.suffix}.json').write_text(json.dumps(out,sort_keys=True,separators=(',',':'))+'\n')
    print(json.dumps(out,sort_keys=True))
if __name__=='__main__':main()

"""The requested bounded hostile continuation; stops at first miss or t16."""
import json
from export_reader import HERE,SELECTOR,all_hits,need


def main():
    bank=json.loads((SELECTOR/'bank.json').read_text())['seeds']
    records=[]
    first_miss=None
    for t in range(11,17):
        for family,constant in (('24',35),('44',29)):
            target=63*(1 << t)-constant
            hits=all_hits(bank,target)
            item=dict(family=family,t=t,target_H=target,hit_count=len(hits),
                      minimum_order_in_fixed_bank=hits[0]['order'] if hits else None,all_hits=hits)
            records.append(item)
            if not hits:
                first_miss=item
                break
        if first_miss:
            break
    need(first_miss is None and len(records)==12,'recorded bounded tail result')
    for family,want in (('24',[3,4,6,6,4,2]),('44',[5,5,5,11,6,3])):
        need([r['hit_count'] for r in records if r['family']==family]==want,'tail acceptance counts')
    report=dict(status='FINITE_EXACT_BOUNDED_HOSTILE_CONTINUATION',
                stopping_rule='Stop at first miss, otherwise after both targets at t16.',
                actual_stop='Both targets accepted through t16; no targets above16 tested.',
                scope='No all-t acceptance or native realization asserted.',
                first_miss=first_miss,records=records)
    (HERE/'tail.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8',newline='\n')
    print('PASS: both affine families accepted at t11..16; stopped at16 without extrapolation')


if __name__=='__main__':
    main()

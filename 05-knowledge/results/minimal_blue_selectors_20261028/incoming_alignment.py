"""Replay incoming wreath-family alignment; compare two distinct targets."""
from pathlib import Path
import sys,json
import blue_minimum as m
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
INCOMING=ROOT/'05-knowledge/results/session90_even_graphs_odd_descent_20261029/carriers'


def main():
    sys.path.insert(0,str(INCOMING))
    incoming=m.r.old.source('incoming_wreath_family',str((INCOMING/'check_wreath_family.py').relative_to(ROOT)))
    receipt=incoming.run()
    rows=[]
    for item in receipt['family']:
        n=item['vertices'];c=item['centralizer_order']
        if n>=5:
            rows.append(dict(order=n,blue_states=4*n-8,blue_fixed_length_bits=(4*n-9).bit_length(),
                             reversal_marks_at_named_root=c,exact_mark_fixed_length_bits=(c-1).bit_length()))
    m.need(rows[1]['order']==27 and rows[1]['reversal_marks_at_named_root']==243,'named-root memory control')
    result=dict(status='PASS',scope='Fresh replay of incoming algorithm, not an independent implementation.',
                inherited=receipt,distinct_observer_comparison=rows)
    (HERE/'incoming-receipt.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(status=result['status'],active_checks=receipt['active_checks'],comparison=rows)))
    return result


if __name__=='__main__':main()

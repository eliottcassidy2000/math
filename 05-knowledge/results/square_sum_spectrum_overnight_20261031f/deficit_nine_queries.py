#!/usr/bin/env python3
"""Actual changing-word witnesses beyond the fixed-affine/periodic obstructions."""
import json,sys
from pathlib import Path
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'square_sum_spectrum_overnight_20261031d'))
from bank_query import query
sys.path.insert(0,str(HERE.parent/'square_sum_spectrum_overnight_20261031c'))
from audit_port_words import adjacency,count,word_ports

def main():
    records=[]
    for k in range(9,37):
        target=(1<<k)-9;r=query(target)
        if r['x'] is None:raise RuntimeError(('unrepresented probe',k,r))
        rows=adjacency(r['word']);actual=count(rows,word_ports(r['word']))
        if actual!=target:raise RuntimeError('actual adjacency disagreement')
        r.update({'exponent':k,'vertices':len(rows),'actual_adjacency_count':actual});records.append(r)
    print(json.dumps({'status':'passed','complete_exponent_range':[9,36],'records':records,
                      'scope':'Exact finite changing-word witnesses; no eventual or universal acceptance theorem.'},indent=2,sort_keys=True))

if __name__=='__main__':main()

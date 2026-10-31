#!/usr/bin/env python3
"""Count actual adjacency for the point-query witnesses, without their inverse."""
from pathlib import Path
import json
import sys
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'square_sum_spectrum_overnight_20261031c'))
from audit_port_words import adjacency,count,word_ports
records=[]
for line in (HERE/'large_point_queries.jsonl').read_text().splitlines():
    r=json.loads(line);rows=adjacency(r['word']);actual=count(rows,word_ports(r['word']))
    if actual!=r['target']:raise RuntimeError(r['target'])
    records.append({'target':r['target'],'order':len(rows),'actual':actual})
print(json.dumps({'status':'passed','records':records,'scope':'Actual adjacency and port partitions, independent of arithmetic inversion.'},indent=2))

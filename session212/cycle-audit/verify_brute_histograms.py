#!/usr/bin/env python3
"""Compare the independent permutation oracle with the complete class census."""
import csv,json
from pathlib import Path
from audit_packet import LEVELS,ck
HERE=Path(__file__).resolve().parent
records=[]
for line in (HERE/'brute_census.out').read_text().splitlines():
    fields=line.split();n,seen=map(int,fields[:2])
    actual=dict(tuple(map(int,x.split(':'))) for x in fields[2:])
    expected={int(r['cycles']):int(r['isomorphism_classes']) for r in csv.DictReader((LEVELS/f'n{n}_hist.tsv').open(),delimiter='\t')}
    ck(actual==expected,('brute class histogram mismatch',n))
    records.append({'order':n,'classes':seen,'distinct_cycle_counts':len(actual)})
ck([r['order'] for r in records]==list(range(3,9)),'order coverage')
print(json.dumps({'status':'PASS','independent_permutation_histograms':records,'total_classes':sum(r['classes'] for r in records)},indent=2))

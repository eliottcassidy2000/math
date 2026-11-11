"""Read exact sequence definitions, source profiles, or typed numerical roles."""
from pathlib import Path
import argparse,csv,gzip,json
p=argparse.ArgumentParser();g=p.add_mutually_exclusive_group(required=True)
g.add_argument('--sequence');g.add_argument('--number',type=int,choices=[7,11,33]);g.add_argument('--class-key',nargs=2,type=int,metavar=('ORDER','FRAME_CODE'));g.add_argument('--coefficient',type=int)
a=p.parse_args();here=Path(__file__).resolve().parent
if a.sequence:
    reg=json.loads((here/'registry.json').read_text());result=reg[a.sequence]
elif a.number:
    roles=json.loads((here/'number_roles.json').read_text());result={n:r[str(a.number)] for n,r in roles.items()}
elif a.class_key:
    result=next(json.loads(s) for s in gzip.open(here/'records.jsonl.gz','rt') if (lambda r:[r['n'],r['rep']]==a.class_key)(json.loads(s)))
else:
    result=next((dict(r) for r in csv.DictReader((here/'seed_witnesses9.tsv').open(),delimiter='\t') if int(r['c'])==a.coefficient),None)
    if result is None:result={'status':'not in complete nine-base unit-seed support','coefficient':a.coefficient,'larger_orders':'not decided by this query'}
print(json.dumps(result,indent=2))

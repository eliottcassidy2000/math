#!/usr/bin/env python3
"""Compile a supplied signed rectangle menu without splitting shared physical rows.

This checks electrical algebra and explicit resource-key capacity. It does NOT
turn unauthenticated keys or paths into native ownership certificates.
"""
import argparse,json
from pathlib import Path
from check import compile_rectangles,label_defects,carry_certificate

def key(x):return json.dumps(x,sort_keys=True,separators=(',',':'))
def unique(items):return len(items)==len({key(x)for x in items})
def main():
 p=argparse.ArgumentParser();p.add_argument('--labels',type=Path,default=Path(__file__).with_name('node_labels.json'));p.add_argument('--menu',type=Path,required=True);p.add_argument('--out',type=Path);a=p.parse_args()
 data=json.loads(a.labels.read_text());menu=json.loads(a.menu.read_text());labels=data['labels'];moduli=data['moduli']
 rects=menu['rectangles'];F=compile_rectangles(rects);defects=label_defects(F,labels,moduli)
 X=[r[-1]for r in labels];carry=carry_certificate(F,X,7)
 result={'status':'EXACT algebraic certificate for supplied rectangle amplitudes','all_dyadic_original_bank_inclusion':not any(any(row)for row in defects.values()),'physical_integer_rows':[[i,sorted(row.items())]for i,row in sorted(F.items())],'owned_row_defects':defects,'higher_carry':carry,'resource_authentication':'NOT PERFORMED: verify paths, owners, shared coordinate meaning, converse orbits and unchanged complement separately.'}
 if 'resource_ledger'in menu:
  ledger=menu['resource_ledger'];before=ledger['before'];after=ledger['after'];fixed=ledger.get('unchanged',[])
  result['explicit_resource_key_capacity']={'before_injective':unique(before+fixed),'after_injective':unique(after+fixed),'old_new_union_injective':unique(before+after+fixed),'unchanged_count':len(fixed)}
 if 'shared_records'in menu:
  records=menu['shared_records'];conflicts=[];maps={}
  for state in ['before','after']:
   values={}
   for r in records.get(state,[]):
    k=key(r['key']);v=key(r['value'])
    if k in values and values[k]!=v:conflicts.append({'state':state,'key':r['key']})
    values[k]=v
   maps[state]=values
  for r in records.get('fixed',[]):
   k=key(r['key']);v=key(r['value'])
   for state in ['before','after']:
    if maps[state].get(k)!=v:conflicts.append({'state':state,'fixed_key':r['key']})
  result['explicit_shared_record_consistency']={'passes':not conflicts,'conflicts':conflicts}
 raw=json.dumps(result,indent=2,sort_keys=True)+'\n'
 if a.out:a.out.write_text(raw)
 else:print(raw,end='')
if __name__=='__main__':main()

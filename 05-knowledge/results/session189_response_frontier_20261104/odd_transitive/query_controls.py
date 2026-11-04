#!/usr/bin/env python3
import json,itertools,sys,hashlib
from pathlib import Path
HERE=Path(__file__).resolve().parent

def need(x,msg):
 if not x:raise RuntimeError(msg)
def bank():return sum([json.loads(p.read_text())for p in sorted(HERE.glob('bank.part*.json'))],[])
def weights(m):
 import math
 return {m-1:1,**{r:(-1)**r*math.comb(m-2,r-1)for r in range(2,(m-1)//2+1)}}
def value(B,a,c):
 masks=B['masks'];co=B['coefficients'];L=len(masks)
 return sum(co[i*L+j]for i,x in enumerate(masks)if x&a==x for j,y in enumerate(masks)if y&c==y)
def prepare():
 controls=json.loads((HERE/'controls.json').read_text());banks=bank();records=[];queries=[]
 for case in controls['cases']:
  if'witness_groups'not in case:continue
  m=case['m'];src=banks[case['graph_ids'][0]]['rows'];ext=case.get('ambient_names',[]);fixed=sum(1<<v for v in ext if src[v]&1)
  for group in case['witness_groups']:
   for r,weight in weights(m).items():
    rec={'case':case['id'],'r':r,'weight':weight,'group':group,'outside_fixed_one_mask':fixed,'queries':[]}
    for flags in range(16):
     a=fixed|(int(bool(flags&1))<<0)|(int(bool(flags&2))<<1)
     c=fixed|((1<<r)if flags&4 else 0)|(sum(1<<i for i in group)if flags&8 else 0)
     q={'flags':flags,'a_mask':a,'c_mask':c,'finite_difference_sign':(-1)**(4-flags.bit_count()),'graph_records':[]}
     for gid in case['graph_ids']:
      B=banks[gid];entry={'id':len(queries),'graph_id':gid,'kind':controls['graphs'][gid]['kind'],'expected_full_polynomial_count':value(B,a,c)}
      queries.append({'rows':B['rows'],'a':a,'c':c,'expected':entry['expected_full_polynomial_count']});q['graph_records'].append(entry)
     rec['queries'].append(q)
    records.append(rec)
 (HERE/'physical_queries.json').write_text(json.dumps(records,indent=2)+'\n')
 with(HERE/'query_input.txt').open('w')as f:
  f.write(str(len(queries))+'\n')
  for q in queries:f.write(str(len(q['rows']))+' '+' '.join(map(str,q['rows']))+' '+str(q['a'])+' '+str(q['c'])+'\n')
 print(json.dumps({'physical_query_records':len(queries),'quartic_extractions':len(records)}))
def check():
 records=json.loads((HERE/'physical_queries.json').read_text());counts=list(map(int,(HERE/'query_output.txt').read_text().split()));seen=0;combo={}
 for rec in records:
  vals={}
  for q in rec['queries']:
   for g in q['graph_records']:
    need(counts[g['id']]==g['expected_full_polynomial_count'],'actual high-weight query differs from full polynomial');seen+=1
    vals[g['kind']]=vals.get(g['kind'],0)+q['finite_difference_sign']*counts[g['id']]
  key=(rec['case'],tuple(rec['group']));z=combo.setdefault(key,{'local':0,'quotient':0})
  z['local']+=rec['weight']*(vals['source']-vals['local_converse']);z['quotient']+=rec['weight']*(vals['source']-vals['quotient_only'])
  need(vals['source']-vals['local_converse']==vals['source']-vals['quotient_only'],'literal grouped restoration identity')
 physical=json.loads((HERE/'physical_receipt.json').read_text())['cases'];out=[]
 for(case,group),z in combo.items():
  ref=next(q for q in physical if q['case']==case);expected=next(q['local_contrast']for q in ref['witness']['grouped']if q['group']==list(group))
  need(z['local']==z['quotient']==expected,'binomial actual query extraction');out.append({'case':case,'group':group,**z})
 need(seen==len(counts),'extra or missing actual counts')
 print(json.dumps({'status':'PASS','actual_hamilton_queries':seen,'quartic_extractions':len(records),'grouped_binomial_results':out,'source_sha256':{f:hashlib.sha256((HERE/f).read_bytes()).hexdigest()for f in['physical_queries.json','query_input.txt','query_output.txt']}},indent=2,sort_keys=True))
if __name__=='__main__':prepare()if len(sys.argv)>1 and sys.argv[1]=='prepare'else check()

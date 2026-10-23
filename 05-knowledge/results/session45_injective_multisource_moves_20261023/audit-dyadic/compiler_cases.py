#!/usr/bin/env python3
"""CLI black-box audit of producer compiler; no producer imports."""
import argparse,json,subprocess,sys,tempfile
from pathlib import Path

p=argparse.ArgumentParser();p.add_argument('--packet',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();packet=a.packet.resolve();checks=0

def need(ok,label):
 global checks
 checks+=1
 if not ok:raise RuntimeError(label)

with tempfile.TemporaryDirectory(prefix='session45-dyadic-compiler-')as temp:
 temp=Path(temp);reports={}
 def run(name,menu,labels=None):
  path=temp/(name+'.json');path.write_text(json.dumps(menu));args=['-B',str(packet/'compile_menu.py'),'--menu',str(path)]
  if labels is not None:
   lab=temp/(name+'.labels.json');lab.write_text(json.dumps(labels));args+=['--labels',str(lab)]
  outputs=[subprocess.check_output([sys.executable]+flag+args)for flag in[[],['-O']]]
  need(outputs[0]==outputs[1],name+' normal/-O equality');out=json.loads(outputs[0]);reports[name]=out;return out
 menu=json.loads((packet/'algebraic_circuit_menu.json').read_text());out=run('algebraic_shared_circuit',menu)
 need(out['all_dyadic_original_bank_inclusion'],'shared circuit passes')
 for i,r in enumerate(menu['rectangles']):need(not run('unsafe_component_'+str(i),{'rectangles':[r]})['all_dyadic_original_bank_inclusion'],'component individually fails')
 hostile={'rectangles':[{'sources':[0,28],'targets':[398,528]}]};out=run('higher_carry_hostile',hostile)
 need(not out['all_dyadic_original_bank_inclusion']and out['higher_carry']['first_failed_bit']==1,'higher defect precisely bit1')
 overlap={'rectangles':[{'sources':[0,1],'targets':[0,1]}]};labels={'labels':[[0,0,0,0,0],[0,0,0,0,64]],'moduli':[2,2,2,2,128]};out=run('overlap_doubled_rows',overlap,labels)
 need(out['all_dyadic_original_bank_inclusion'],'overlap row contributions cancel')
 need(out['physical_integer_rows']==[[0,[[0,2],[1,-2]]],[1,[[0,-2],[1,2]]]],'overlap diagonal and factor2 retained')
 resources={'rectangles':[],'resource_ledger':{'before':['a'],'after':['b'],'unchanged':['a']},'shared_records':{'before':[{'key':'p','value':0},{'key':'p','value':1}],'after':[{'key':'p','value':0}],'fixed':[{'key':'p','value':0}]}}
 out=run('resource_and_shared_hostile',resources)
 need(out['explicit_resource_key_capacity']=={'before_injective':False,'after_injective':True,'old_new_union_injective':False,'unchanged_count':1},'unchanged complement reserves capacity')
 need(not out['explicit_shared_record_consistency']['passes'],'physical shared-key conflict and pin rejected')
 valid={'rectangles':[],'resource_ledger':{'before':['a'],'after':['b'],'unchanged':['c']},'shared_records':{'before':[{'key':'p','value':0},{'key':'p','value':0}],'after':[{'key':'p','value':0}],'fixed':[{'key':'p','value':0}]}}
 out=run('consistent_explicit_ledger',valid)
 need(out['explicit_shared_record_consistency']['passes']and all(out['explicit_resource_key_capacity'][k]for k in['before_injective','after_injective','old_new_union_injective']),'consistent ledger passes')
 need(all('NOT PERFORMED' in r['resource_authentication']for r in reports.values()),'compiler does not overclaim native authentication')
 result={'status':'PASS','checks':checks,'cases':reports,'scope':'Only algebraic currents and consistency of supplied keys; native paths, coefficient-to-resource matching, and semantic ownership are not authenticated by compiler.'}
 a.out.write_text(json.dumps(result,sort_keys=True,indent=2)+'\n');print(json.dumps({'status':'PASS','checks':checks}))

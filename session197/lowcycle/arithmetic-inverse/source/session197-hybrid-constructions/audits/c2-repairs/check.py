import json,sys
from pathlib import Path
P=Path(__file__).resolve().parent
sys.path.insert(0,str(P.parents[1]/'pair-entry-bridge'))
from core import *
sys.path.insert(0,str(P.parents[1]/'audits/lowcycle-two-port'))
from compiler import compile_two
src=P.parents[2]/'session197-lowcycle/arithmetic-c2-witnesses.json'
recipes=json.loads(src.read_text());wanted={13,41,43,55,77,83,163,173,231,233,247,251,277,299,421,449,503,533,569,583,637,821,1067,1339,1861,2647}
require({x['target'] for x in recipes}==wanted and len(recipes)==26,'all26 exact target set')
records=[]
for rec in recipes:
 T=rec['source_rows'];a,b=rec['pair'];p=inventories(T);require(cycles(T)==1 and p['s'][a]==p['s'][b]==1,'C1 pair source')
 response=pair_response(T,a,b);require(list(response)==rec['response'],'actual pair response')
 T,r=glue(T,{a,b},gadget(rec['pair_L']));state=pair_transfer(response,rec['pair_L']);require(root_profile(T,r)==state,'pair seed transfer')
 for L in rec['ordinary_Ls']:
  T,r=glue(T,{r},gadget(L));state=single_transfer(state,L)
 p=inventories(T);require((p['H'],p['e'][r],p['R'][r])==state and p['H']==rec['target'] and p['s'][r]==2 and cycles(T)==2,'final actual profile')
 word=compile_two(T,r);require(r in word['ports'],'retained root')
 records.append(dict(target=rec['target'],rows=T,root=r,profile=list(state),ports=word['ports'],birth_order=word['birth_order'],moves=[[x['mask'],x['drop']] for x in word['steps'][2:]]))
(P/'physical-words.json').write_text(json.dumps(records,indent=2,sort_keys=True)+'\n')
(P/'walk-input.txt').write_text(str(len(records))+'\n'+'\n'.join(' '.join(map(str,[len(x['rows']),x['target']]+x['rows'])) for x in records)+'\n')
print(json.dumps({'status':'PASS','targets':sorted(wanted),'physical_witnesses':len(records),'max_order':max(len(x['rows']) for x in records),'same_root_two_port_compilations':len(records)},indent=2,sort_keys=True))

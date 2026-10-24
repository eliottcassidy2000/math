#!/usr/bin/env python3
"""Independent complete helper-pin census and target570 local hostile."""
import argparse,csv,gzip,json,sys
from array import array
from collections import Counter
from itertools import combinations
from pathlib import Path
H=Path(__file__).resolve().parent;O=H/'inputs/original';checks=0
p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);a=p.parse_args()
def need(x,m):
 global checks
 checks+=1
 if not x:raise RuntimeError(m)
N=list(csv.DictReader((O/'metagraph_order8_20261023/merged_nodes.tsv').read_text().splitlines(),delimiter='\t'))
own=array('H');own.frombytes(gzip.decompress((O/'metagraph_order8_20261023/owner_map.u16.gz').read_bytes()))
if sys.byteorder!='little':own.byteswap()
merge={int(r[k]):int(r['merged_id']) for r in N for k in ['class_a','class_b']};pair=list(combinations(range(8),2));chords=[(i,j) for i in range(8) for j in range(i+2,8)];full=2**21-1
E=json.loads((O/'metagraph_order8_electrical_20261023/receipt.json').read_text())['order8'];X=E['higher_layers'][-1]['lifts'][0];Y=E['radical_lifts_mod4']
def label(s):return (sum((Y[j][s]%2)<<j for j in range(4)),X[s]%128)
def plus(a,b):return (a[0]^b[0],(a[1]+b[1])%128)
B=json.loads((O/'compatible_charged_helpers_20261023/supply/selected_supply.json').read_text());main=set(B['banks'][B['main_bank_index']]['source_ids']);out=[];extra=[]
for h in B['helpers']:
 s=h['source'];M=[[0]*8 for _ in range(8)];c=int(N[s]['canon_a'])
 for k,(i,j) in enumerate(pair):M[i][j]=(c>>(27-k))&1;M[j][i]=1-M[i][j]
 F=dict(h['selected']['common_clocks']);stack=[((),255)];raw=[]
 while stack:
  p,left=stack.pop();k=len(p)
  if not left:
   w=sum(M[p[i]][p[j]]<<b for b,(i,j) in enumerate(chords));raw.append((p,w,merge[own[full^w]]));continue
  for v in range(8):
   if left>>v&1 and (k not in F or F[k]==v) and (not p or M[p[-1]][v]):stack.append((p+(v,),left^(1<<v)))
 count=Counter(t for p,w,t in raw);need(dict(count)=={int(k):v for k,v in h['selected']['pins_target_count'].items()},'full raw helper target census')
 need(len(raw)==h['selected']['complete_pinned_count'],'full raw helper path count')
 need(count[1965]==count[2234]==1,'selected U/V records are complete in this pin context')
 hits=[{'path':p,'word':w,'target':t} for p,w,t in raw if t==570]
 if hits:extra.append({'source':s,'pins':sorted(F.items()),'label':label(s),'in_main40':s in main,'paths':hits})
 out.append({'source':s,'raw_count':len(raw),'raw570':count[570],'raw1965':count[1965],'raw2234':count[2234]})
old=json.loads((O/'session45_injective_multisource_moves_20261023/multisource-trades/target_pair_inventory.json').read_text())['records']
old_census=[]
for r in old:
 need(r['before']['path'][0]==r['after']['path'][0] and r['before']['path'][-1]==r['after']['path'][-1],'old source retains endpoint owner pins')
 s=r['source'];p0=r['before']['path'];p1=r['after']['path'];F={i:v for i,v in enumerate(p0) if v==p1[i]};M=[[0]*8 for _ in range(8)];c=int(N[s]['canon_a'])
 for k,(i,j) in enumerate(pair):M[i][j]=(c>>(27-k))&1;M[j][i]=1-M[i][j]
 stack=[((),255)];raw=[]
 while stack:
  p,left=stack.pop();k=len(p)
  if not left:
   w=sum(M[p[i]][p[j]]<<b for b,(i,j) in enumerate(chords));raw.append((p,w,merge[own[full^w]]));continue
  for v in range(8):
   if left>>v&1 and (k not in F or F[k]==v) and (not p or M[p[-1]][v]):stack.append((p+(v,),left^(1<<v)))
 count=Counter(t for p,w,t in raw)
 need(count==Counter({1965:1,2234:1,1559:1} if s==1923 else {1965:1,2234:1}),'complete old all-common-pin raw target census')
 need({p for p,w,t in raw if t in [1965,2234]}=={tuple(p0),tuple(p1)},'complete old selected record domain')
 old_census.append({'source':s,'pins':sorted(F.items()),'raw_count':len(raw),'targets':dict(count)})
need(sum(r['raw_count'] for r in old_census)==57,'all old28 raw paths accounted for')
need(len(extra)==1 and extra[0]['source']==2393 and len(extra[0]['paths'])==1,'only incoming2393 can choose570 in complete pins')
S=[1130,1673,1754,2678];wanted=plus(label(1130),label(1673));need(wanted==(11,126) and label(2393)==(11,0),'actual local source charge target')
cases=[]
for s in S:
 total=plus(label(s),label(2393));need(total!=wanted,'one fresh plus2393 cannot pay original570 source sum')
 # Current change at570 is -2 times the changed sum of source chi, because margin is fixed.
 delta=-2*(X[s]+X[2393]-X[1130]-X[1673]);need(delta%256!=0,'actual original X alone separates2393 case')
 cases.append({'remaining_fresh_source':s,'paired_sum':total,'target570_X_current_residue_mod256':delta%256})
result={'status':'PASS: complete53-helper raw pin census; only2393 has570 alternative, excluded by the original570 row in the margin2 baseline','checks':checks,'helper_raw_census':out,'old28_retains_endpoint_pins':True,'old28_complete_common_pin_census':old_census,'extra570_records':extra,'required_target570_source_label':wanted,'2393_plus_one_fresh_rejections':cases,'scope':'Original fresh four, supplied old28 and any subset of supplied53, exact common-clock pins, fixed global target570 margin2, zero margins outside{570,1965,2234}, baseline helper/old570 populationzero. Does not establish arbitrary mixed-pin fibres or allow new supplies/target margins.'}
a.out.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n');print(json.dumps({'status':'PASS','checks':checks,'extra570_sources':[r['source'] for r in extra]}))

#!/usr/bin/env python3
"""Independent exact extension census at the unique strong H=13 core."""
import hashlib,itertools,json
from pathlib import Path
from verify import graph,dp,strong_rows,require,value,strong_word
CORE=graph((1,2,0))

def brute(rows):
 return sum(all(rows[u]>>v&1 for u,v in zip(w,w[1:])) for w in itertools.permutations(range(len(rows))))
def from_mask(n,mask):
 rows=[0]*n;k=0
 for u in range(n):
  for v in range(u+1,n):
   if mask>>k&1:rows[u]|=1<<v
   else:rows[v]|=1<<u
   k+=1
 return rows

def extend(rows,mask):
 n=len(rows);out=rows+[mask]
 for u in range(n):
  if not(mask>>u&1):out[u]|=1<<n
 return out

def relabel(rows,p):
 return tuple(sum(1<<j for j in range(len(rows)) if rows[p[i]]>>p[j]&1) for i in range(len(rows)))

def permutation_completion():
 # Count every 7-name path into ALL compatible choices of the 11 new arcs.
 # This independently avoids Hamilton subset DP and the two-apex transfer.
 counts=[0]*2048;accepted_orders=0
 for w in itertools.permutations(range(7)):
  ones=0;zeros=0;valid=True
  for u,v in zip(w,w[1:]):
   hi=max(u,v);lo=min(u,v)
   if hi<5:
    if not(CORE[u]>>v&1):valid=False;break
   else:
    bit=lo if hi==5 else 5+lo
    if u==hi:ones|=1<<bit
    else:zeros|=1<<bit
  if not valid:continue
  require(not(ones&zeros),'contradictory repeated edge');accepted_orders+=1
  free=2047^(ones|zeros);s=free
  while True:
   counts[ones|s]+=1
   if s==0:break
   s=(s-1)&free
 return counts,accepted_orders

def main():
 require(5**5>13**3 and 5**7>39**3,'all-order strong cutoffs')
 parents=[]
 for mask in range(1024):
  rows=from_mask(5,mask);h=dp(rows)
  require(h==brute(rows),'all order-five path counts')
  if h==13:parents.append(tuple(rows));require(strong_rows(rows),'13 strong')
 orbit={relabel(CORE,p) for p in itertools.permutations(range(5))}
 require(len(parents)==120 and set(parents)==orbit,'unique H13 isomorphism class')
 one=[]
 for a in range(32):
  rows=extend(CORE,a);h=dp(rows);require(h==brute(rows),'one-birth independent count')
  one.append({'mask':a,'H':h,'strong':strong_rows(rows)})
 expected=[13,23,25,29,31,33,37,43,45]
 require(sorted({r['H'] for r in one})==expected,'complete one-birth support')
 independent,orders=permutation_completion();two=[];hist={}
 for code in range(2048):
  a=code&31;b=code>>5;rows=extend(extend(CORE,a),b);h=dp(rows)
  require(h==independent[code],'complete two-birth independent count')
  hist[h]=hist.get(h,0)+1
  if h in (39,41,43):two.append({'first_mask':a,'second_mask':b,'H':h,'strong':strong_rows(rows)})
 require(39 not in hist,'no two-birth H39')
 for w,h in [((1,2,0),13),((1,2,0,1),23),((1,2,0,1,1),41),((1,0,2,2,2),39)]:
  require(dp(graph(w))==brute(graph(w))==h and strong_word(w),'strong positive witness')
 require([graph((1,2,0,1,1))[i]&31 for i in range(5)]==CORE,'literal two-birth containment')
 end_fan=[((1,2,0),13),((1,2,0,1,1),41),((3,1,2,0),43),((3,1,2,0,1,1),179)]
 for w,h in end_fan:require(dp(graph(w))==brute(graph(w))==h and strong_word(w),'asymmetric end fan')
 require((3,)+((1,2,0)+(1,1))==((3,)+(1,2,0))+(1,1),'literal end commutation')
 require(3*41+4==127 and 3*43+2==131,'incompatible mixed targets')
 result={'status':'PASS','core_rows':CORE,'all_order_five_tournaments_checked':1024,
 'H13_labelled_parents':120,'H13_isomorphism_classes':1,'one_birth_support':expected,'one_birth_records':one,
 'two_birth_extensions_checked':2048,'independent_permutation_orders_surviving_old_arcs':orders,
 'two_birth_histogram':hist,'two_birth_target_records':two,
 'strong_chain':{'words':['120','1201','12011'],'H':[13,23,41]},
 'strong_H39_witness':{'word':'10222','rows':graph((1,0,2,2,2))},
 'strong_H39_never_contains_induced_H13':True,'asymmetric_end_fan':{'parent':13,'append11':41,'prepend3':43,'mixed_word':'312011','mixed_H':179,'incompatible_targets':[127,131]}}
 raw=(json.dumps(result,indent=2,sort_keys=True)+'\n').encode();Path(__file__).with_name('extension_receipt.json').write_bytes(raw)
 print(json.dumps({'status':'PASS','H13_isomorphism_classes':1,'one_birth_support':expected,'two_birth_extensions_checked':2048,
 'receipt_sha256':hashlib.sha256(raw).hexdigest()},sort_keys=True))
if __name__=='__main__':main()

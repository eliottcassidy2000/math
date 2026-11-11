"""Exact paired-count restoration, not an unrestricted spectrum search."""
from collections import defaultdict
from pathlib import Path
import importlib.util,json
HERE=Path(__file__).resolve().parent
sp=importlib.util.spec_from_file_location('paired_probe',HERE/'probe.py')
p=importlib.util.module_from_spec(sp);sp.loader.exec_module(p)


def main():
 cap=32767
 old=defaultdict(list)
 for C in range(1,cap+1):old[13*p.c.g(C)+334*C+198].append(C)
 matches=[]
 for D in range(1,cap+1):
  h=25*p.c.g(D)+244*D+186
  for C in old[h]:
   a=p.word_state((211,13,347,13),C)
   b=p.word_state((211,25,269,25),D)
   if a[0]!=b[0] or a[0]!=h:raise RuntimeError(('synchronization',C,D))
   matches.append({'H_pair':[h,h+2],'old_address':C,'new_address':D,
      'old_state':a,'new_state':b,'U_saving':a[2]-b[2],
      'old_original_tail':4+8*C,'new_original_tail':1+2*D})
 matches.sort(key=lambda x:(x['H_pair'][0],x['old_address'],x['new_address']))
 first=matches[0]
 N=first['H_pair'][0]
 if (N-211)//334>cap or (N-211)//244>cap:raise RuntimeError('leastness cap insufficient')
 # Two fixed-digit infinite-family attempts; only primitive odd t matters.
 p.c.g.cache_clear()
 ones=[]
 for t in range(1,200001,2):
  if 13*p.c.g(122*t-1)-25*p.c.g(167*t-1)==78:ones.append(t)
 p.c.g.cache_clear()
 def f(n):return n+1-(n&1)-p.c.g(n)
 zeros=[]
 for t in range(1,500001,2):
  if 13*f(269*t)==25*f(347*t):zeros.append(t)
 report={'status':'FINITE-EXACT actual paired resynchronizations; no all-height family established',
  'address_box':[1,cap],'match_count':len(matches),
  'new_U_no_larger':sum(x['U_saving']>=0 for x in matches),
  'new_U_larger':sum(x['U_saving']<0 for x in matches),
  'least_nonempty_continuation_sync':first,
  'leastness_target_bounds':{'old':(N-211)//334,'new':(N-211)//244},
  'all_one_extension_attempt':{'odd_t_bound':200000,'equation':'13g(122t-1)-25g(167t-1)=78','solutions':ones},
  'all_zero_extension_attempt':{'odd_t_bound':500000,'equation':'13f(269t)=25f(347t)','solutions':zeros},
  'matches':matches}
 (HERE/'synchronization.json').write_text(json.dumps(report,indent=2)+'\n')
 print(json.dumps({k:v for k,v in report.items() if k!='matches'},indent=2))


if __name__=='__main__':main()

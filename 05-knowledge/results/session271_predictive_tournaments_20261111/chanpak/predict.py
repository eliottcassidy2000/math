#!/usr/bin/env python3
"""Bounded predictive probes; no exhaustive large-target census. Python + NumPy."""
from pathlib import Path
import importlib.util,itertools,json,math
import numpy as np
from decimal import Decimal,localcontext,ROUND_CEILING
HERE=Path(__file__).resolve().parent
OLD=HERE.parents[1]/'session270_chanpak_enumeration_20261111'/'chanpak'/'spine.py'
spec=importlib.util.spec_from_file_location('spine',OLD);spine=importlib.util.module_from_spec(spec);spec.loader.exec_module(spine)

def growth_rank(core,support):
 k=len(core);best=0;witness=None
 for p in itertools.permutations(range(k)):
  allowed=[any(not(c>>p[0]&1) for c in support)]
  allowed += [any(c>>p[j-1]&1 for c in support) and any(not(c>>p[j]&1) for c in support) for j in range(1,k)]
  allowed += [any(c>>p[-1]&1 for c in support)]
  for active in itertools.product((0,1),repeat=k+1):
   if any(x and not y for x,y in zip(active,allowed)):continue
   if any(not active[j] and not(core[p[j-1]]>>p[j]&1) for j in range(1,k)):continue
   if sum(active)>best:best=sum(active);witness={'permutation':p,'active_slots':active}
 return best,witness

def main():
 rng=np.random.default_rng(27120261111)
 v,a,edges=spine.machine((2,0));a=np.array(a)
 matrices=[]
 for es in edges:
  M=np.zeros((36,36))
  for i,j in es:M[j,i]+=1
  matrices.append(M)
 samples=1200;lengths=[16,32,64,128,256,512]
 cases={'uniform':[.25]*4,'skewed':[.7,.1,.1,.1],'binary03':[.5,0,0,.5],'binary01':[.5,.5,0,0],'constant2':[0,0,1,0]}
 stats=[]
 for name,prob in cases.items():
  support=[i for i,p in enumerate(prob) if p];r,w=growth_rank((2,0),support)
  V=np.tile(v,(samples,1)).astype(float)
  for m in range(1,max(lengths)+1):
   letters=rng.choice(4,size=samples,p=prob)
   for c in support:
    mask=letters==c
    V[mask]=V[mask]@matrices[c].T/r
   if m in lengths:
    scaled=V@a;d=np.log(scaled)
    stats.append({'law':name,'m':m,'samples':samples,'rank':r,'mean_H_div_r_pow_m':float(np.mean(scaled)),'sd_log_H':float(np.std(d,ddof=1)),'quantiles_H_div_r_pow_m':dict(zip(['p05','p50','p95'],map(float,np.quantile(scaled,[.05,.5,.95])))),'mean_logH_per_m':float(math.log(r)+np.mean(d)/m)})
 support_table=[]
 for mask in range(1,16):
  S=[i for i in range(4) if mask>>i&1];r,w=growth_rank((2,0),S)
  word=S*20;H=spine.value(word)
  support_table.append({'support':S,'rank':r,'witness':w,'periodic_probe_length':len(word),'periodic_logH_per_m':math.log(H)/len(word)})
 # Small random finite support coverage diagnostics, intentionally not a large atlas.
 hist=[]
 for m in [8,12,16,24,48]:
  vals=[spine.value(rng.integers(0,4,size=m).tolist()) for _ in range(2000)]
  hist.append({'m':m,'samples':len(vals),'distinct':len(set(vals)),'median_H_floor':sum(sorted(vals)[999:1001])//2,'residues':{str(q):sorted({h%q for h in vals}) for q in [8,9,15,25]},'max_H':max(vals)})
 # Rigorous full-support growth bracket and mean-scale proxy, not witness claims.
 targets=[]
 for digits in [6,12,30,100,1000,1000000]:
  with localcontext() as ctx:
   ctx.prec=60
   L=Decimal(digits)*Decimal(10).ln()
   lower=2+int(((L-Decimal(2).ln())/Decimal(3).ln()).to_integral_value(rounding=ROUND_CEILING))
   proxy=2+int(((L+Decimal(8).ln())/Decimal(3).ln()).to_integral_value(rounding=ROUND_CEILING))
  targets.append({'target':'10^'+str(digits),'universal_family_order_lower_bound':lower,'mean_scale_order_proxy':proxy})
 result={'seed':27120261111,'status':'finite stochastic predictions; theorem-derived rank; no exact coverage inference','growth_samples':stats,'support_classification':support_table,'bounded_coverage_samples':hist,'large_target_order_scales':targets}
 (HERE/'predictions.json').write_text(json.dumps(result,indent=2)+'\n')
 rows=['# Bounded probes and large-target scales','','Seed: 27120261111. Each growth row uses 1,200 iid words; standard deviations are empirical, not confidence guarantees. No atlas was regenerated.','','| Law | m | r | mean H/r^m | median H/r^m | sd(log H) |','|---|---:|---:|---:|---:|---:|']
 for x in stats:
  rows.append(f"| {x['law']} | {x['m']} | {x['rank']} | {x['mean_H_div_r_pow_m']:.6g} | {x['quantiles_H_div_r_pow_m']['p50']:.6g} | {x['sd_log_H']:.6g} |")
 rows+=['','## Target order scales, not existence certificates','','For H=10^d (even: deliberately only a numerical scale), use a nearby eligible odd target if seeking a witness. The lower bound holds for every word. The proxy is the order where the ensemble mean reaches H; it is neither a minimum nor an upper bound.','','| d | rigorous language-specific order lower bound | ensemble mean-scale order proxy |','|---:|---:|---:|']
 for x in targets: rows.append(f"| {x['target'][3:]} | {x['universal_family_order_lower_bound']} | {x['mean_scale_order_proxy']} |")
 rows+=['','## Complete support classification for the two-core iid carrier','','| Positive-probability letters | exponential base r |','|---|---:|']
 for x in support_table:rows.append(f"| {x['support']} | {x['rank']} |")
 (HERE/'TABLES.md').write_text('\n'.join(rows)+'\n')
 print(json.dumps({'growth_rows':len(stats),'support_cases':len(support_table),'coverage_sample_words':10000,'targets':targets},indent=2))
if __name__=='__main__':main()

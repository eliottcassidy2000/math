from pathlib import Path
DATA_DIR = Path(__file__).resolve().parent
from itertools import product

def req(x,m):
 if not x:raise RuntimeError(m)
def step(z,b):
 F,q,g,h=z
 return {'0':(F+q+g,2*q,0,h),'1':(F,2*q+h,h,h),
 'C0':(3*F+9*q+3*g,12*q,0,3*h),'C3':(3*F,12*q+9*h,3*h,3*h),
 'C1':(3*F+6*q+2*g+2*h,12*q+3*h,h,3*h),'C2':(3*F+3*q+g+2*h,12*q+6*h,2*h,3*h)}[b]
zero=lambda z:all(a%3==0 for a in z)
for z in product(range(3),repeat=4):
 for b in ('C0','C3'):req(zero(step(z,b)),('constant',z,b))
 for a,b,c in product(('C1','C2'),repeat=3):req(zero(step(step(step(z,a),b),c)),('triple',z,a,b,c))
 for a,c in product(('C1','C2'),repeat=2):
  for b in ('0','1'):req(zero(step(step(step(z,a),b),c)),('separated',z,a,b,c))
print('mod3_all_81_input_states_constant_triple_separated_PASS')
miss=[34869,85885,103061,116587,153451,245923,535405,594793,595543,609873,618723,767135,774073,782103,824119,826373,868367,949777]
factors={34869:[3,59,197],85885:[5,89,193],116587:[23,37,137],153451:[173,887],245923:[31,7933],535405:[5,13,8237],595543:[13,61,751],609873:[3,11,18481],618723:[49,9,23,61],767135:[5,153427],782103:[3,260701],824119:[149,5531],826373:[353,2341],868367:[157,5531]}
import math,json
for k,fs in factors.items():
 req(math.prod(fs)==k,('factorization',k,fs))
 req(all(0<v<=1000000 and v%2 and v not in miss+[7,21] for v in fs),('factor_not_direct',k,fs))
repairs={103061,594793,774073,949777}
req(set(factors)|repairs==set(miss),'miss_partition')
rows=[json.loads(x)for x in open(DATA_DIR / 'block_target_witnesses.jsonl')]
req({x['target']for x in rows}==repairs,'repair_target_set')
print('all_18_misses_partition_14_direct_factorizations_4_explicit_repairs_PASS')
print('CONDITIONAL_ON_seeded_words_and_check_block_targets_receipts_allowed_prefix_999999')

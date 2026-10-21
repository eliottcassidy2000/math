#!/usr/bin/env python3
import json,pathlib,hashlib,datetime,time
p=pathlib.Path(__file__).resolve();t=time.monotonic();start=datetime.datetime.now(datetime.timezone.utc).isoformat()
orders=set(range(3,48,2));frames=[]
for q in range(2,9):
 for m in (1,3,5,7):
  n=q*m
  if n==2:continue
  if q==2 and m>=5:continue
  if q in (6,8) and m>=5:continue
  if q in (5,7) and m>=5 and n>=49:continue
  if n>47:continue
  orders.add(n);frames.append([q,m,n])
even=sorted(n for n in orders if n%2==0)
if even!=[4,6,8,12,18,20,24,28]:raise RuntimeError('even support mismatch')
unknown=sorted(n for n in orders if n>9)
expected=sorted(set(range(11,48,2))|{12,18,20,24,28})
if unknown!=expected or len(unknown)!=24:raise RuntimeError('unknown coefficient support mismatch')
R={'status':'parameter arithmetic only; low-order coefficients inherited','start_utc':start,'end_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'wall_seconds':time.monotonic()-t,'source_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'surviving_q_gt_1_parameter_rows':frames,'possible_even_orders':even,'unknown_coefficient_orders':unknown,'unknown_coefficient_count':len(unknown)}
(p.parent/'support_receipt.json').write_text(json.dumps(R,indent=2)+'\n');print(json.dumps(R,indent=2))

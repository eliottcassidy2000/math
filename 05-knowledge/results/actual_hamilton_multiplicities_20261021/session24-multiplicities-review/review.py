import csv, json, math, pathlib, hashlib
root=pathlib.Path(__file__).resolve().parent.parent
p=root/'session24-multiplicities'
rows=[{k:int(v) for k,v in r.items()} for r in csv.DictReader((p/'n8_classes.csv').open())]
checks=0
def need(ok,label):
 global checks
 checks+=1
 if not ok: raise RuntimeError(label)
need(len(rows)==6880,'class count')
need(len({r['code'] for r in rows})==6880,'distinct canonical codes')
for r in rows:
 a=[(r['code']>>(8*i))&255 for i in range(8)]
 need(all(not(a[i]>>i&1) for i in range(8)),'loopless')
 need(all(((a[i]>>j)&1)+((a[j]>>i)&1)==1 for i in range(8) for j in range(i+1,8)),'tournament')
 need(r['aut']%2==1 and math.factorial(8)%r['aut']==0,'group order')
 need(r['H']%r['aut']==0,'free Hamilton action')
need(sum(math.factorial(8)//r['aut'] for r in rows)==2**28,'full labelled universe mass')
need(sum(r['H']//r['aut'] for r in rows)==2**21,'fixed Hamilton presentation mass')
need(sum(r['self_converse'] for r in rows)==176,'SC count')
need(sum(r['strong'] for r in rows)==6008,'strong count')
need(5**8>73**3,'global low-H cutoff')
replay=(root/'session24-multiplicities-review/n8-independent-replay.csv').read_bytes()
need(replay==(p/'n8_classes.csv').read_bytes(),'fresh compiled augmentation matches')
print(json.dumps({'status':'PASS','checks':checks,'n8_sha256':hashlib.sha256(replay).hexdigest(),'scope':'Validity, freeness, full labelled mass and freshly rebuilt canonical augmentation; not an all-order spectrum-completeness claim.'},indent=2,sort_keys=True))

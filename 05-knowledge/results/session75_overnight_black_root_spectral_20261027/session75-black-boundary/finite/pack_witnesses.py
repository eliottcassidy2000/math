#!/usr/bin/env python3
"""One-time receipt compressor; args are raw complete producer JSONL files."""
import json,gzip,sys,collections,hashlib
counts=collections.Counter();nindex=collections.Counter();digest=hashlib.sha256()
with gzip.GzipFile(filename='',mode='wb',fileobj=sys.stdout.buffer,mtime=0,compresslevel=9)as out:
 for path in sys.argv[1:]:
  for line in open(path):
   r=json.loads(line);n=len(r['rows']);z={'n':n,'i':nindex[n],'c':r['code']};nindex[n]+=1
   if r.get('verdict')=='NOT_STRONG':z['t']='N'
   elif r.get('verdict')=='PURE_BLUE':z['t']='P'
   else:
    w=next(c for c in r['B']if c['verdict']=='WITNESS');z.update(t='W',d=w['dual'],v=w['port'],p=w['route'])
   s=json.dumps(z,separators=(',',':'),sort_keys=True)+'\n';out.write(s.encode())

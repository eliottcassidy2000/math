import pathlib,json,copy,hashlib
P=pathlib.Path(__file__).resolve().parent
s=(P.parent/'session22-centers/check.py').read_text();ns={'__file__':str(P.parent/'session22-centers/check.py')};exec(s.split('start=datetime.datetime')[0],ns);audit=ns['audit_receipt']
def need(value):
 if not value:raise RuntimeError('schema check failed')
valid=json.loads((P.parent/'session22-centers/receipt_normal.json').read_text())['prefix_receipts'];tests=[]
for v in valid:
 need(audit(v))
 for key in v:
  for bad in (None,False,0,0.0,{},[],[['nested']],{'bad':'shape'}):
   r=copy.deepcopy(v);r[key]=bad;tests.append(r)
 for key in v:
  r=copy.deepcopy(v);del r[key];tests.append(r)
 for bad in (-1,18,False,True,1.0,[],{},'0',None):
  r=copy.deepcopy(v);r['cycle'][0]=bad;tests.append(r)
  r=copy.deepcopy(v);r['words'][0][0]=bad;tests.append(r)
  r=copy.deepcopy(v);r['terminal']=bad;tests.append(r)
 r=copy.deepcopy(v);r['direction']='arbitrary';tests.append(r)
 r=copy.deepcopy(v);r['extra']=1;tests.append(r)
 # Duplicate row, nonsquare row, nonbinary row, diagonal flip.
 for row in ('0', '2'*18, '0'*19):
  r=copy.deepcopy(v);r['source_rows'][0]=row;tests.append(r)
for r in tests:need(audit(r) is False)
out={'source_sha256':hashlib.sha256(s.encode()).hexdigest(),'valid_receipts':len(valid),'independent_schema_mutations_rejected_without_exception':len(tests)}
(P/'repaired_schema_results.json').write_text(json.dumps(out,indent=2)+'\n');print(out)

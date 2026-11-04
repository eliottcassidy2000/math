import pathlib,json,re,hashlib
D=pathlib.Path(__file__).resolve().parent
R=D.parent
s=(D/'NOTE.md').read_text()
def need(x,msg):
 if not x:raise RuntimeError(msg)
links=re.findall(r'\]\(([^)]+)\)',s)
for link in links:need((D/link).is_file(),'Broken relative link '+link)
c=json.load(open(R/'missing44/candidate44.json'));a=json.load(open(R/'missing44_audit/native_audit44.json'));p=json.load(open(R/'ports/candidate.json'));r=next(q for q in p['ports']if q['name']==-1)
need((c['h'],c['F'],c['J'],c['target_H'],c['root_E2'])==(40085,633897,499420,4611841,5111261),'Main fields')
need(a['source_names']==c['names']==a['target_names'] and a['target']==c['rows'],'Physical lineage')
need(a['counts']==[1537279,2014907,1495091,1559067,2144877,2910413,4610117,5056377,5311829,4216865,4611841],'Native trace')
need(a['native_complements']==20 and a['temporary_helpers']==10 and a['deleted_old_names']==[],'Native inventory')
need([q['helper']for q in a['trace']]==list(range(100,110)),'Helpers')
need(p['in_ports']==[2,5,6,10,12,13,15,17,20,22],'Alternate cut')
need((r['E1'],r['S1'],r['deletion_H'],r['E2'],r['S2'],r['defect'])==(1,1,1,5420133,7137627,808292),'Alternate root fields')
need(r['golden']==[4611841,3568813,3568814],'Golden')
need(r['charts']['F49_charge']==[6,4,6] and r['charts']['F49_frobenius']==[1,3,6],'F49')
need(r['charts']['F343_vector_codec']==[[3,3,4],[4,6,4],[3,5,6]],'F343')
q=json.load(open(R/'ports_audit/receipt.json'));need(q['distinct_named_adjacencies']==4 and q['min_family_order']==24 and len(q['all_six'])==6,'Representation scope')
receipt={'status':'PASS','relative_links':links,'synthesis_sha256':hashlib.sha256((D/'NOTE.md').read_bytes()).hexdigest(),'scope':'Synthesis consistency and link checks; symbolic review in AUDIT.md.'}
(D/'audit_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,sort_keys=True))

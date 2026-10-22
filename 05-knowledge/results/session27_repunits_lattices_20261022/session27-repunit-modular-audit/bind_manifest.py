import pathlib,hashlib,json
O=pathlib.Path(__file__).parent
P=pathlib.Path('/workspace/shared/session27-repunits');M=pathlib.Path('/workspace/shared/session27-modular');R=pathlib.Path('/workspace/scratch/07d5532f2fce/math')
files=[P/n for n in ['REPORT.md','verify_weighted.py','verify_weighted.json','catalogue.json','catalogue_initial.json','catalogue_breakdown.json','discover.py','discovery_summary.txt','fetch_manifest.json','corrections.patch']]+[M/n for n in ['report.md','check.py','check.json','check_camion.py','check_camion.json','historical_1729_correction_proposal.md']]+list((P/'fetched').iterdir())
files += [R/p for p in ['05-knowledge/results/actual_hamilton_multiplicities_20261021/session24-multiplicities/PROOF.md','05-knowledge/results/actual_hamilton_multiplicities_20261021/session24-multiplicities/n8_classes.csv','01-canon/theorems/THM-227-knacci-mersenne.md','01-canon/theorems/THM-868-e8-bridge-score-lattice.md','01-canon/theorems/THM-489-eta-discriminant-bridge.md','07-reflections/the-1729-resonance-is-isolated-the-tournament-ratio-has-no-modular-structure.md']]
files += list((R/'01-canon/theorems').glob('THM-4746-*'))+list((R/'01-canon/theorems').glob('THM-4749-*'))
blob_checks=[]
for e in json.loads((P/'fetch_manifest.json').read_text()):
 b=(P/e['local_copy']).read_bytes();actual=hashlib.sha1(b'blob '+str(len(b)).encode()+b'\0'+b).hexdigest();blob_checks.append({'file':e['local_copy'],'expected_remote_blob_sha':e['sha'],'actual_local_git_blob_sha':actual,'equal':actual==e['sha']})
assert all(e['equal'] for e in blob_checks)
result={'status':'PASS','freeze':'Producer freeze confirmed; original packet only. Remote-search and E8-holes addenda excluded.','source_files':[{'path':str(p),'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in files],'connector_blob_checks':blob_checks,'review_files':[{'path':p.name,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in O.iterdir() if p.is_file() and p.name!='manifest.json']}
(O/'manifest.json').write_text(json.dumps(result,indent=2)+'\n')
print('PASS; bound',len(files),'source files')

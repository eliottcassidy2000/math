from pathlib import Path, PurePosixPath
import json,hashlib,re,subprocess
root=Path(__file__).resolve().parent.parent
sha=lambda b:hashlib.sha256(b).hexdigest()
def need(c,msg):
 if not c: raise RuntimeError(msg)
m=json.loads((root/'session53-root/first/PUBLICATION_MAP.json').read_text()); dests={}
for f in m['files']:
 b=Path(f['source']).read_bytes(); need(len(b)==f['bytes'] and sha(b)==f['sha256'],'map bytes '+f['source']); need(f['destination'] not in dests,'duplicate'); dests[f['destination']]=f
counts={}
for name,pin in m['packet_manifests'].items():
 fs=[f for d,f in dests.items() if str(PurePosixPath(d).parent)==m['destination']+'/'+name and PurePosixPath(d).name.lower()=='manifest.json'];need(len(fs)==1,'manifest unique')
 f=fs[0];need(f['sha256']==pin,'pin mismatch'); md=json.loads(Path(f['source']).read_text()); entries=md.get('files',md)
 if isinstance(entries,list): entries={x['path']:x['sha256'] for x in entries}
 actual={str(PurePosixPath(d).relative_to(m['destination']+'/'+name)):v['sha256'] for d,v in dests.items() if d.startswith(m['destination']+'/'+name+'/') and d!=f['destination']}
 need(entries==actual,'closure '+name);counts[name]=len(entries)
summary=(root/'session53-root/first/README.md').read_bytes()
links=re.findall(r'\]\(([^)]+)\)',summary.decode())
for l in links: need(m['destination']+'/'+l in dests,'summary link '+l)
r=json.loads((root/'session53-publication/navigation-compaction-receipt.json').read_text()); src=subprocess.check_output(['git','-C',str(root/'math'),'show','HEAD:'+r['source_path']]); compact=(root/'session53-publication/navigation-option'/r['source_path']).read_bytes()
need(sha(src)==r['source_sha256'],'navigation source');need(sha(compact)==r['result_sha256'],'navigation result')
s=compact.decode();expanded=s; occurrences={}
for x in r['replacements']:
 label,target=x['label'],x['target'];definition=f'[{label}]: {target}\n';need(expanded.count(definition)==1,'definition unique');expanded=expanded.replace(definition,'');pat=r'\[([^\]\n]+)\]\['+label+r'\]';n=len(re.findall(pat,expanded));need(n==x['occurrences'],'occurrences '+label); occurrences[label]=n;expanded=re.sub(pat,lambda mt:'['+mt[1]+']('+target+')',expanded)
need(expanded.endswith('\n\n'),'appended separator');expanded=expanded[:-1];need(expanded.encode()==src,'exact expansion');need(len(src)-len(compact)==r['bytes_saved'],'savings');need(len(src.splitlines())==r['lines_before'] and len(compact.splitlines())==r['lines_after'],'lines')
result={'status':'PASS','summary_sha256':sha(summary),'publication_map_sha256':sha((root/'session53-root/first/PUBLICATION_MAP.json').read_bytes()),'mapped_files':len(dests),'packet_payload_counts':counts,'summary_links':links,'navigation_source_git_head':subprocess.check_output(['git','-C',str(root/'math'),'rev-parse','HEAD']).decode().strip(),'navigation_source_sha256':sha(src),'navigation_result_sha256':sha(compact),'navigation_expansion_byte_exact':True,'navigation_reference_uses':sum(occurrences.values()),'navigation_bytes_saved':len(src)-len(compact)}
print(json.dumps(result,indent=2))

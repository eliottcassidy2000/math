"""Read only physically present tracked UTF-8 files. Never fetch lazy git blobs.
Usage: python3 recover_inventory.py /path/to/math > inventory.json
"""
import sys,subprocess,pathlib,re,json,hashlib
root=pathlib.Path(sys.argv[1]);manifest=subprocess.check_output(['git','-C',str(root),'ls-files','-z']);paths=[p.decode() for p in manifest.split(b'\0') if p]
patterns={'fraction':r'143\s*/\s*900|21\s*/\s*90|4\s*/\s*9','hessian':'Hessian','representation':'representation theor','triple':r'\{2,\s*3,\s*11\}','repunit':'repunit','1729':'1729'}
patterns={k:re.compile(v,re.I) for k,v in patterns.items()};matches={k:[] for k in patterns};stats=dict(tracked_paths=len(paths),read_files=0,missing=0,oversized=0,unreadable=0)
for p in paths:
 f=root/p
 if not f.is_file(): stats['missing']+=1;continue
 if f.stat().st_size>5000000: stats['oversized']+=1;continue
 try:t=f.read_text(encoding='utf-8')
 except (OSError,UnicodeError):stats['unreadable']+=1;continue
 stats['read_files']+=1
 for k,rx in patterns.items():
  if not rx.search(t):continue
  hits=[{'line':i,'text':s[:1000]} for i,s in enumerate(t.splitlines(),1) if rx.search(s)]
  matches[k].append({'path':p,'hits':hits})
print(json.dumps({'stats':stats,'tracked_manifest_sha256':hashlib.sha256(manifest).hexdigest(),'queries':{k:r.pattern for k,r in patterns.items()},'matches':matches},indent=2))

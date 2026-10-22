"""Search present tracked text only; never fetch Git blobs. Usage: python discover.py REPO OUT."""
import pathlib,subprocess,re,json,sys,collections,hashlib
root=pathlib.Path(sys.argv[1]); dest=pathlib.Path(sys.argv[2]); git=lambda *a:subprocess.check_output(['git','-C',str(root),*a]).decode()
paths=[x for x in git('ls-files','-z').split('\0') if x]
pat=re.compile(r'repunit|mersenne|wagstaff|fibonacci[ -_]+fan|binary[ -_]+ones|all[ -_]+ones|trailing[ -_]+ones|2\s*(?:\^|\*\*)\s*\{?[knm]\}?\s*[-−]\s*1',re.I)
entries=[]; scanned=0;missing=0;large=[];nontext=[]
for p in paths:
 f=root/p;fn=bool(pat.search(p));hits=[];sha=None
 if f.is_file():
  if f.stat().st_size>10_000_000:large.append(p);continue
  try:
   raw=f.read_bytes();s=raw.decode();scanned+=1
   hits=[{'line':i,'text':l[:700]}for i,l in enumerate(s.splitlines(),1)if pat.search(l)]
   if hits:sha=hashlib.sha256(raw).hexdigest()
  except(UnicodeError,OSError):nontext.append(p)
 else:missing+=1
 if fn or hits:entries.append({'path':p,'present':f.is_file(),'filename_match':fn,'sha256':sha,'matches':hits})
obj={'head':git('rev-parse','HEAD').strip(),'tracked_paths':len(paths),'text_files_scanned':scanned,'missing_worktree_files':missing,'large_files_skipped':large,'nontext_or_read_error_count':len(nontext),'regex':pat.pattern,'entries':entries}
dest.write_text(json.dumps(obj,indent=2));print(json.dumps({k:v for k,v in obj.items()if k!='entries'},indent=2));print('matching entries',len(entries),'matching lines',sum(len(e['matches'])for e in entries))

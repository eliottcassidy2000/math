#!/usr/bin/env python3
"""Repository-wide text discovery, not a claim to have audited every hit."""
import subprocess,json,re,collections,hashlib
from pathlib import Path
root=Path(__file__).resolve().parent
base='66cf347030dc'
pattern=r'square[- _–]sum|squareful|A090461|A071983|A071984|Gerbicz|ninja pair|nice pair|sums (are|to be) (perfect )?squares'
raw=subprocess.check_output(['git','grep','-n','-i','-E',pattern,base,'--','*.md'],text=True)
rows={}
for line in raw.splitlines():
 p,ln,txt=line[len(base)+1:].split(':',2)
 row=rows.setdefault(p,{'path':p,'matches':[]})
 row['matches'].append({'line':int(ln),'text':txt[:1800]})
terms={
 'interval_paths':r'hamilton|Gerbicz|A07198|ninja|nice pair|Q_[nN0-9]|G_[nN0-9]',
 'ports_and_gluing':r'port|seam|glue|endpoint|matching|forest|boundary',
 'elliptic_and_magic':r'elliptic|magic|biclique|K.?3.?3|quartic|genus',
 'dynamics_and_clocks':r'collatz|Pell|clock|Berggren|FRACTRAN|cat.map|Rule.?90',
 'group_and_spectrum':r'monster|supersingular|tournament|metagraph|Fano|Hamming',
 'labelling_and_arithmetic':r'graceful|aliquot|sigma|Fermat|Zeckendorf|triangular',
}
cache={}
for p,row in rows.items():
 blob=subprocess.check_output(['git','rev-parse',base+':'+p],text=True).strip();row['blob']=blob
 if blob not in cache:
  text=subprocess.check_output(['git','show',base+':'+p],text=True)
  cache[blob]=text
 text=cache[blob];row['bytes']=len(text.encode());row['title']=next((x[2:] for x in text.splitlines() if x.startswith('# ')),p)
 row['themes']=[k for k,v in terms.items() if re.search(v,text,re.I)]
 row['copy_or_audit']=bool(re.search(r'/(sources?|inputs?|audit[^/]*|review[^/]*)/',p,re.I))
 row['session_hint']=(re.findall(r'^session:\s*(.+)',text,re.M) or [''])[0]
 row['date_hint']=(re.findall(r'^date:\s*(.+)',text,re.M) or re.findall(r'20(?:25|26)[-_]?\d{2}[-_]?\d{2}',p) or [''])[0]
 row['direct_subject']=bool(re.search(r'square[- _]sum|Gerbicz|ninja|nice pair|A07198',row['title'],re.I) or re.search(r'square[- _]sum|Gerbicz|ninja',p,re.I))
 row['headings']=[x.lstrip('# ') for x in text.splitlines() if x.startswith('##')][:25]
 # Large logs retain only located hit excerpts, not the whole text.
report={'baseline':base,'pattern':pattern,'scope':'All tracked Markdown at the pinned revision, including audit copies and incidental matches. Discovery/indexing is not proof review. Untracked files, other branches, deleted-only files and unmentioned synonyms are not exhausted.','files':list(rows.values()),'counts':{'files':len(rows),'distinct_blobs':len(cache),'direct_subject_files':sum(r['direct_subject'] for r in rows.values()),'direct_noncopy_files':sum(r['direct_subject'] and not r['copy_or_audit'] for r in rows.values())}}
(root/'CATALOG.json').write_text(json.dumps(report,indent=2,ensure_ascii=False)+'\n')
lines=['# Square-sum discovery index','',report['scope'],'',str(report['counts']),'','The reviewed mechanism map is in RECOVERY.md. The list below is discovery, not acceptance.','']
for p,row in rows.items():
 if row['direct_subject'] and not row['copy_or_audit']:
  lines.append('- `'+p+'` — '+row['title'][:200]+'. Themes: '+', '.join(row['themes']))
(root/'CATALOG.md').write_text('\n'.join(lines)+'\n')
print(json.dumps(report['counts']))

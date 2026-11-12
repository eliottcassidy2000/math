#!/usr/bin/env python3
"""Separate shallow reachable commit-subject search, not content history."""
import subprocess,re,json,pathlib,collections
root=pathlib.Path(subprocess.check_output(['git','rev-parse','--show-toplevel'],text=True).strip());out=pathlib.Path(__file__).resolve().parent
ref='7510e38f9c777a61a1b3210334077e0ebe0df2ad'
rules={'five':r'\bfive\b|fivefold','named':r'platonic|pentagon|Fermat|Ogg|moonshine|golden|icosahed','central':r'isothermal|stress traject|Nem.nyi','meta_search':r'(five|Platonic|Fermat).*(web|history|atlas|search|bridge)|(history|atlas|search|bridge).*(five|Platonic|Fermat)'}
lines=subprocess.check_output(['git','log','--format=%H%x09%as%x09%s',ref],cwd=root,text=True).splitlines();rows=[];counts=collections.Counter()
for line in lines:
 sha,date,subject=line.split('\t',2);tags=[k for k,v in rules.items() if re.search(v,subject,re.I)]
 if tags: rows.append({'sha':sha,'date':date,'subject':subject,'tags':tags});counts.update(tags)
(out/'commit_history.json').write_text(json.dumps(rows,indent=2)+'\n')
(out/'commit_scope.json').write_text(json.dumps({'ref':ref,'shallow':subprocess.check_output(['git','rev-parse','--is-shallow-repository'],cwd=root,text=True).strip(),'reachable_commit_subjects':len(lines),'matched_subjects':len(rows),'tag_counts':dict(counts),'rules':rules,'oldest_reachable':lines[-1],'newest_reachable':lines[0],'limitation':'Only reachable commit subjects; no exhaustive older/deleted content scan. Old dated bodies in current tree are separately included.'},indent=2)+'\n')

"""Validate the completed frozen-ref scan, then index retrieval leads.

Superseded partial passes were stopped without publishing their results.
Only history_hits.jsonl.gz from the completed frozen pass is authoritative.
The scan covers complete blob tag presence and bounded early excerpts.
"""
from pathlib import Path
import json,gzip,re,hashlib,collections
import history_scan as hs
HERE=Path(__file__).resolve().parent
def need(x,msg):
    if not x:raise RuntimeError(msg)
def run():
    stats=json.loads((HERE/'history_stats.json').read_text(encoding='utf8'))
    raw=hs.frozen_git('rev-list','--objects').decode('utf8','replace').splitlines()
    candidates=[q for q in raw if ' ' in q and q.split(' ',1)[1].lower().endswith(('.md','.txt','.rst'))]
    need(len(raw)==stats['reachable_object_count'] and len(candidates)==stats['candidate_text_objects'],'frozen initial universe mismatch')
    manifest=json.loads(gzip.decompress((HERE/'history_universe_manifest.json.gz').read_bytes()))
    need([(q['oid'],q['path']) for q in manifest]==[tuple(q.split(' ',1)) for q in candidates],'exact frozen manifest mismatch')
    need(sum(q['bytes'] for q in manifest if q['type']=='blob')==stats['counts']['bytes_scanned'],'manifest byte scope mismatch')
    baseline=stats['baseline'];current={}
    for line in hs.git('ls-tree','-r',baseline).decode('utf8','replace').splitlines():
        meta,path=line.split('\t',1);current[path]=meta.split()[2]
    with gzip.open(HERE/'history_hits.jsonl.gz','rt',encoding='utf8') as f:lines=f.read().splitlines()
    hits=[json.loads(q) for q in lines]
    for h in hits:
        h.pop('current_at_discovery_path',None)
        h['baseline_at_discovery_path']=current.get(h['path'])==h['oid']
    packed=('\n'.join(json.dumps(h,ensure_ascii=False) for h in hits)+'\n').encode('utf8')
    (HERE/'history_hits.jsonl.gz').write_bytes(gzip.compress(packed,mtime=0))
    commits=hs.frozen_git('log','--format=%H%x09%aI%x09%s').decode('utf8','replace').splitlines()
    matches=[q for q in commits if any(rx.search(q.split('\t',2)[-1]) for rx in hs.RX.values())]
    (HERE/'history_commits.txt').write_text('\n'.join(matches)+'\n',encoding='utf8')
    stats.pop('current_discovery_hits',None)
    stats.update({'frozen_refs_file':'history_refs.txt','baseline_discovery_hits':sum(h['baseline_at_discovery_path'] for h in hits),'reachable_commits':len(commits),'matching_commit_subjects':len(matches),'is_shallow_repository':hs.git('rev-parse','--is-shallow-repository').decode().strip()=='true','raw_jsonl_sha256':hashlib.sha256(packed).hexdigest(),'compressed_sha256':hashlib.sha256((HERE/'history_hits.jsonl.gz').read_bytes()).hexdigest(),'universe_check':'Completed pass uses saved initial refs and baseline tree; independently regenerated frozen object-universe counts agree.'})
    stats['limits'][0]='Only the saved initial local Git refs; shallow boundaries can omit earlier history.'
    (HERE/'history_stats.json').write_text(json.dumps(stats,indent=2,ensure_ascii=False)+'\n',encoding='utf8')
    # Every candidate excerpt is indexed, not just the high-priority subset.
    grouped={};markers=re.compile(r"(?:owner.?s?|user.?s?|human.?s?)\s+(?:directive|prompt|request|instruction)|prompt\s*\(?user|verbatim",re.I)
    for h in hits:
        for p in h['prompt_candidates']:
            normalized=' '.join(p['text'].split());key=hashlib.sha256(normalized.encode()).hexdigest()
            if key not in grouped:grouped[key]={'key':key,'text':normalized,'explicit_prompt_marker':bool(markers.search(normalized)),'provenance':[],'topics':set()}
            grouped[key]['provenance'].append({'path':h['path'],'oid':h['oid'],'line':p['line'],'baseline_at_discovery_path':h['baseline_at_discovery_path']})
            grouped[key]['topics'].update(h['topics'])
    prompts=[]
    for p in grouped.values():p['topics']=sorted(p['topics']);prompts.append(p)
    prompts.sort(key=lambda p:(not p['explicit_prompt_marker'],p['provenance'][0]['path'],p['key']))
    data=(json.dumps(prompts,ensure_ascii=False,separators=(',',':'))+'\n').encode()
    (HERE/'history_prompt_candidates.json.gz').write_bytes(gzip.compress(data,mtime=0))
    path_topics=collections.defaultdict(set)
    for h in hits:path_topics[h['path']].update(h['topics'])
    topics=collections.Counter(t for ts in path_topics.values() for t in ts)
    summary={'unique_context_candidates':len(prompts),'with_explicit_prompt_marker':sum(p['explicit_prompt_marker'] for p in prompts),'candidate_discovery_paths':len({q['path'] for p in prompts for q in p['provenance']}),'topic_hit_paths':dict(topics),'limits':'Normalized contextual excerpts, not verified quotations or counts of distinct sessions. A prompt can mention five outside its quoted text. All source variants retained by blob.'}
    (HERE/'history_prompt_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    # Bounded readable queue for manual inspection; the complete catalog is above.
    eligible=[p for p in prompts if p['explicit_prompt_marker']]
    queue={'scope':'First120 explicit-marker contexts in the deterministic sort; the complete version-aware collection is history_prompt_candidates.json.gz. Retrieval leads, not authenticated user messages.','total_eligible':len(eligible),'records':eligible[:120]}
    (HERE/'history_prompt_review_queue.json').write_text(json.dumps(queue,indent=2,ensure_ascii=False)+'\n',encoding='utf8')
    print(json.dumps({'history':stats,'prompts':summary}),flush=True)
if __name__=='__main__':run()

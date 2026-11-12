"""Full-file prompt-attribution harvest on the largest revision of each path.

This supplements (does not replace) the all-revision bounded-excerpt catalog.
Largest revision is a retrieval choice, not a statement of current truth.
"""
from pathlib import Path
import json,gzip,subprocess,re,hashlib,collections,time
import history_scan as h
MARKER=re.compile(r"\b(?:owner(?:['’]s)?\s+(?:directive|request|prompt|seed|asked|asks)|user(?:['’]s)?\s+(?:prompt|request|directive)|original\s+request|(?:directive|prompt)\s*\(verbatim\))",re.I)
def run():
    start=time.time();items=json.loads(gzip.decompress((h.HERE/'history_universe_manifest.json.gz').read_bytes()));by={}
    for q in items:
        if q['path'] not in by or (q['bytes'],q['oid'])>(by[q['path']]['bytes'],by[q['path']]['oid']):by[q['path']]=q
    chosen=sorted(by.values(),key=lambda q:q['path']);catalog={};source_count=0;markers=0
    p=subprocess.Popen(['git','cat-file','--batch'],cwd=h.ROOT,stdin=subprocess.PIPE,stdout=subprocess.PIPE)
    for idx,q in enumerate(chosen):
        p.stdin.write((q['oid']+'\n').encode());p.stdin.flush();header=p.stdout.readline().decode().split()
        if len(header)!=3 or header[1]!='blob':raise RuntimeError(header)
        data=p.stdout.read(int(header[2]));need=p.stdout.read(1)
        if need!=b'\n':raise RuntimeError('batch framing')
        s=data.decode('utf8','replace');low=s.lower().replace('\u0131','i').replace('\u017f','s')
        if not any(w in low for w in ['owner','user','original','verbatim']):continue
        source_tags=[k for k,rx in h.RX.items() if (len(s)!=len(low) or any(w in low for w in h.CHEAP[k])) and h.fast_search(rx,h.STARTS[k],s,low)]
        if not source_tags:continue
        ms=list(MARKER.finditer(s))
        if not ms:continue
        source_count+=1
        for m in ms:
            markers+=1;left=s.rfind('\n',0,m.start())+1;right=left
            for _ in range(12):
                right=s.find('\n',right+1)
                if right<0:right=len(s);break
            text=s[left:right][:2600];normalized=' '.join(text.split());key=hashlib.sha256(normalized.encode()).hexdigest()
            direct=[k for k,rx in h.RX.items() if rx.search(text)]
            if key not in catalog:catalog[key]={'key':key,'text':text,'five_in_excerpt':bool(direct),'excerpt_tags':direct,'provenance':[]}
            catalog[key]['provenance'].append({'oid':q['oid'],'path':q['path'],'line':s.count('\n',0,left)+1,'source_tags':source_tags})
        if idx%5000==0:print(json.dumps({'at':idx,'source_paths':source_count,'contexts':len(catalog)}),flush=True)
    p.stdin.close();p.wait()
    contexts=sorted(catalog.values(),key=lambda q:(not q['five_in_excerpt'],q['provenance'][0]['path'],q['key']))
    raw=(json.dumps(contexts,ensure_ascii=False,separators=(',',':'))+'\n').encode()
    (h.HERE/'history_full_prompt_attributions.json.gz').write_bytes(gzip.compress(raw,mtime=0))
    scope={'selection':'Largest byte-size revision per Git discovery path, ties by object ID, within frozen manifest. Complete marker search inside each chosen text; at most twelve lines/2600 characters per context.','selected_paths':len(chosen),'selected_bytes':sum(q['bytes'] for q in chosen),'source_paths_with_five_and_prompt_marker':source_count,'marker_occurrences':markers,'unique_normalized_contexts':len(contexts),'contexts_with_five_in_excerpt':sum(q['five_in_excerpt'] for q in contexts),'seconds':round(time.time()-start,2),'limits':['Attribution markers are not authentication of original chats.','Five may occur only elsewhere in a source; such contexts are separately flagged.','Largest revision is not necessarily the current or corrected one.','All-revision catalog remains separate; this pass does not search every line of every revision for prompt markers.'],'selected_objects':[(q['oid'],q['path']) for q in chosen]}
    (h.HERE/'history_prompt_harvest_scope.json.gz').write_bytes(gzip.compress((json.dumps(scope,separators=(',',':'))+'\n').encode(),mtime=0))
    summary={k:v for k,v in scope.items() if k!='selected_objects'}
    (h.HERE/'history_prompt_harvest_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary),flush=True)
if __name__=='__main__':run()

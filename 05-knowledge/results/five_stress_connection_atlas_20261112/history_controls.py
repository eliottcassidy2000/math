"""Ordinary-regex control of localized searching and bounded excerpts."""
from pathlib import Path
import json,gzip,subprocess,random,time
import history_scan as h
def need(x,msg):
    if not x:raise RuntimeError(msg)
def reference(s):
    lines=s.splitlines();out=[];prompts=[]
    for ln,line in enumerate(lines,1):
        tags=[k for k,rx in h.RX.items() if rx.search(line)]
        if not tags:continue
        out.append({'line':ln,'tags':tags,'text':line[:700]})
        context=' '.join(lines[max(0,ln-3):min(len(lines),ln+1)])
        if h.PROMPT.search(context) and len(prompts)<6:prompts.append({'line':ln,'text':context[:1100]})
        if len(out)==8:break
    return out,prompts
def run():
    t=time.time();items=json.loads(gzip.decompress((h.HERE/'history_universe_manifest.json.gz').read_bytes()))
    samples=random.Random(512).sample([q for q in items if q['bytes']<100000],120)+sorted(items,key=lambda q:q['bytes'],reverse=True)[:1]+[next(q for q in items if q['path']=='05-knowledge/results/INDEX.md')]
    synthetic=['Hecke checked fifteen fivefold five\n','family\nof\n5\nuser request five\n', 'İSOTHERMAL ſtress trajectories ıSOTHERMAL\nfive\n','\n'.join('owner request five '+str(j) for j in range(20)), 'user request\n\nfive\n\nnext','5-vertex 5 element 5-set\npentagon\n','a family of     5\na system of 5\n', 'five\nfamily of\n5\nrequest\nfive']
    for q in samples:
        s=subprocess.check_output(['git','cat-file','blob',q['oid']]).decode('utf8','replace')
        check(s,q['oid'],False)
    for i,s in enumerate(synthetic):check(s,'synthetic'+str(i),True)
    out={'scope':'All eighteen ordinary-regex and localized first-match positions agree on122 deterministic Git blobs and8 hostile strings. First-eight-line excerpts agree exactly, and contextual prompts agree modulo whitespace with a separate linewise implementation. LF line convention.','sample_rule':'Random(512),120 blobs under100000 bytes; largest blob; first results INDEX blob.','sample_oids':[q['oid'] for q in samples],'synthetic_cases':synthetic,'seconds':round(time.time()-t,3)}
    (h.HERE/'history_search_controls.json').write_text(json.dumps(out,indent=2,ensure_ascii=False)+'\n',encoding='utf8')
    print(json.dumps({'blobs':len(samples),'hostiles':len(synthetic),'seconds':out['seconds']}),flush=True)
def check(s,label,synthetic):
    low=s.lower().replace('\u0131','i').replace('\u017f','s')
    tags=[]
    for k,rx in h.RX.items():
        a=rx.search(s);b=h.fast_search(rx,h.STARTS[k],s,low)
        need((a.start() if a else None)==(b.start() if b else None),label+' '+k)
        guarded=(len(s)!=len(low) or any(t in low for t in h.CHEAP[k])) and b is not None
        need(bool(a)==guarded,label+' cheap '+k)
        if a:tags.append(k)
    for k,rx in h.TRX.items():
        a=rx.search(s);b=h.fast_search(rx,h.TCHEAP[k],s,low)
        need((a.start() if a else None)==(b.start() if b else None),label+' '+k)
    a,b=h.early_excerpts(s,low,tags),reference(s)
    need(a[0]==b[0],label+' excerpt/control')
    normalize=lambda rows:[(q['line'],' '.join(q['text'].split())) for q in rows]
    need(normalize(a[1])==normalize(b[1]),label+' contextual prompt/control')
if __name__=='__main__':run()

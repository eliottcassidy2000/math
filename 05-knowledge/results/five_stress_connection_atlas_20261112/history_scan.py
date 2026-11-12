"""Bounded-excerpt, word-boundary census of all reachable Markdown/text blobs.

No mathematical claim is validated by a lexical hit. One blob may have many
paths/revisions, and a hit is not automatically a distinct user session.
Run from this worktree: python3 -X utf8 .../history_scan.py
"""
from pathlib import Path
import subprocess, re, json, collections, hashlib, time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
PATTERNS = {
    'five_word': r'\bfive\b',
    'five_object': r'(?<!\w)5[- ](?:element|vertex|point|parameter|constant|dimensional|fold|core|group|set|tuple|cycle|adic)\b',
    'set_of_5': r'\b(?:set|group|collection|family|list|system)s?\s+of\s+5\b',
    'five_geometry': r'\b(?:pentagon(?:al)?|pentagram|pentatope|quintic|quintuplet|icosahedr(?:on|al)|dodecahedr(?:on|al))\b',
    'stress': r'\b(?:nem[eé]nyi|isothermal|isometric\s+stress|stress\s+trajector(?:y|ies))\b',
}
RX = {k: re.compile(v, re.I) for k,v in PATTERNS.items()}
PROMPT = re.compile(r'\b(?:user|prompt|owner|request|directive|instructions?|asked|asks|search|recover|find|session)\b',re.I)
TOPICS = {
    'tournaments': r'\b(?:tournaments?|chan[– -]?pak|hamilton(?:ian)?|metagraph)\b',
    'primes': r'\b(?:primes?|primality|sieve|lehmer|carmichael)\b',
    'collatz': r'\b(?:collatz|syracuse|hailstone)\b',
    'fibonacci': r'\b(?:fibonacci|zeckendorf|golden|lucas|tribonacci)\b',
    'lrc': r'\b(?:lrc|lonely\s+runner|runners?)\b',
    'jacobian': r'\b(?:jacobian|keller|gmc|gaussian\s+moments?)\b',
    'pentagon_coherence': r'\b(?:associahedr\w*|pentagon|coherence|catalan)\b',
    'groups': r'\b(?:groups?|representations?|characters?|involutions?|symmetries)\b',
    'complex_geometry': r'\b(?:hermitian|m[oö]bius|lorentz|conformal|riemann|isothermal)\b',
    'stress': r'\b(?:stress|elasticity|airy|biharmonic|nem[eé]nyi)\b',
    'graphs': r'\b(?:graphs?|flows?|incidence|spectral|laplacian)\b',
    'computation': r'\b(?:automata|turing|cellular|finite.state|recursion)\b',
    'number_fields': r'\b(?:elliptic|bsd|chowla|eisenstein|cyclotomic|modular)\b',
}
TRX={k:re.compile(v,re.I) for k,v in TOPICS.items()}

def git(*args):
    return subprocess.check_output(['git',*args],cwd=ROOT)

def scan():
    t=time.time()
    head=git('rev-parse','HEAD').decode().strip()
    # rev-list yields each reachable object once and one observed path. Paths
    # are discovery names, not an exhaustive rename history.
    raw=git('rev-list','--all','--objects').decode('utf8','replace').splitlines()
    objects=[]
    for line in raw:
        parts=line.split(' ',1)
        if len(parts)==2 and parts[1].lower().endswith(('.md','.txt','.rst')):
            objects.append(parts)
    print(json.dumps({'phase':'inventory','reachable_objects':len(raw),'candidate_text_objects':len(objects)}),flush=True)
    p=subprocess.Popen(['git','cat-file','--batch'],cwd=ROOT,stdin=subprocess.PIPE,stdout=subprocess.PIPE)
    hits=[]; counts=collections.Counter(); skip=[]; by_topic=collections.Counter()
    for idx,(oid,path) in enumerate(objects):
        p.stdin.write((oid+'\n').encode());p.stdin.flush()
        header=p.stdout.readline().decode().split()
        if len(header)!=3: raise RuntimeError(header)
        size=int(header[2]);data=p.stdout.read(size);trail=p.stdout.read(1)
        if trail!=b'\n':raise RuntimeError('batch framing')
        if header[1]!='blob':continue
        counts['text_blobs_scanned']+=1;counts['bytes_scanned']+=size
        s=data.decode('utf8','replace')
        tags=[k for k,rx in RX.items() if rx.search(s)]
        if not tags:continue
        lines=s.splitlines();excerpts=[];prompt=[]
        for ln,line in enumerate(lines,1):
            lt=[k for k,rx in RX.items() if rx.search(line)]
            if lt:
                e={'line':ln,'tags':lt,'text':line[:700]}
                if len(excerpts)<8:excerpts.append(e)
                # Two-line context helps multiline prompts without asserting
                # that every contextual hit is an actual human instruction.
                context=' '.join(lines[max(0,ln-3):min(len(lines),ln+1)])
                if PROMPT.search(context) and len(prompt)<6:
                    prompt.append({'line':ln,'text':context[:1100]})
        topics=[k for k,rx in TRX.items() if rx.search(s)]
        by_topic.update(topics);counts.update(tags)
        title=next((x.lstrip('# ').strip()[:200] for x in lines[:70] if x.startswith('# ')),path)
        hits.append({'oid':oid,'path':path,'bytes':size,'title':title,'tags':tags,'topics':topics,'excerpts':excerpts,'prompt_candidates':prompt})
        if idx%10000==0:print(json.dumps({'phase':'scan','at':idx,'hits':len(hits),'seconds':round(time.time()-t,1)}),flush=True)
    p.stdin.close();p.wait()
    commits=git('log','--all','--format=%H%x09%aI%x09%s').decode('utf8','replace').splitlines()
    cm=[]
    for c in commits:
        if any(rx.search(c.split('\t',2)[-1]) for rx in RX.values()):cm.append(c)
    # Current paths and exact object identities provide stable local retrieval.
    current={}
    for line in git('ls-files','-s').decode('utf8','replace').splitlines():
        meta,path=line.split('\t',1);current[path]=meta.split()[1]
    for h in hits:h['current_at_discovery_path']=current.get(h['path'])==h['oid']
    with (HERE/'history_hits.jsonl').open('w',encoding='utf8',newline='\n') as f:
        for h in hits:f.write(json.dumps(h,ensure_ascii=False)+'\n')
    (HERE/'history_commits.txt').write_text('\n'.join(cm)+'\n',encoding='utf8')
    stats={'baseline':head,'reachable_object_count':len(raw),'candidate_text_objects':len(objects),'counts':dict(counts),'hit_blobs':len(hits),'hit_discovery_paths':len({h['path'] for h in hits}),'current_discovery_hits':sum(h['current_at_discovery_path'] for h in hits),'topic_cooccurrence_blobs':dict(by_topic),'reachable_commits':len(commits),'matching_commit_subjects':len(cm),'seconds':round(time.time()-t,2),'patterns':PATTERNS,'topic_patterns':TOPICS,'limits':['Reachable local Git refs only; fetched origin/main at startup.','Each object has one Git discovery path; rename aliases not exhaustively reconstructed.','Lexical hits and prompt candidates are not validated claims or distinct sessions.','Binary attachments and inaccessible chats are outside the text universe.','Topic co-occurrence is evidence for retrieval, not a mathematical map.']}
    (HERE/'history_stats.json').write_text(json.dumps(stats,indent=2,ensure_ascii=False)+'\n',encoding='utf8')
    print(json.dumps(stats,ensure_ascii=False),flush=True)

if __name__=='__main__':scan()

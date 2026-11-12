"""Bounded-excerpt, word-boundary census of all reachable Markdown/text blobs.

No mathematical claim is validated by a lexical hit. One blob may have many
paths/revisions, and a hit is not automatically a distinct user session.
Run from this worktree: python3 -X utf8 .../history_scan.py
"""
from pathlib import Path
import subprocess, re, json, collections, hashlib, time, gzip

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
CHEAP={
 'five_word':['five'],
 # Necessary suffix words reject huge numeric b-files before a regex pass.
 'five_object':['element','vertex','point','parameter','constant','dimensional','fold','core','group','set','tuple','cycle','adic'],
 'set_of_5':['set','group','collection','family','list','system'],
 'five_geometry':['pentagon','pentagram','pentatope','quintic','quintuplet','icosahedr','dodecahedr'],
 'stress':['nem','isothermal','isometric','stress'],
}
TCHEAP={
 'tournaments':['tournament','chan','hamilton','metagraph'],
 'primes':['prime','primality','sieve','lehmer','carmichael'],
 'collatz':['collatz','syracuse','hailstone'],
 'fibonacci':['fibonacci','zeckendorf','golden','lucas','tribonacci'],
 'lrc':['lrc','lonely','runner'], 'jacobian':['jacobian','keller','gmc','gaussian'],
 'pentagon_coherence':['associahedr','pentagon','coherence','catalan'],
 'groups':['group','representation','character','involution','symmetries'],
 'complex_geometry':['hermitian','mobius','möbius','lorentz','conformal','riemann','isothermal'],
 'stress':['stress','elasticity','airy','biharmonic','nem'],
 'graphs':['graph','flow','incidence','spectral','laplacian'],
 'computation':['automata','turing','cellular','finite','recursion'],
 'number_fields':['elliptic','bsd','chowla','eisenstein','cyclotomic','modular'],
}
STARTS={**CHEAP,'five_object':['5'],'set_of_5':['set','group','collection','family','list','system']}

def fast_search(rx,starts,s,lower,start=0):
    """Exact regex search localized to necessary literal initial prefixes.

    Each regex alternative starts with one of these literals. Preserve the
    original regex for boundaries, suffixes and arbitrary whitespace. Unicode
    lowercase length changes use the ordinary regex as a safe fallback.
    """
    if len(s)!=len(lower):return rx.search(s,start)
    best=None
    for word in starts:
        p=lower.find(word,start)
        while p>=0 and (best is None or p<best.start()):
            m=rx.match(s,p)
            if m is not None:best=m;break
            p=lower.find(word,p+1)
    return best

def early_excerpts(s,lower,tags):
    # The first eight matching lines of the union can be found among the
    # first eight matching lines of each constituent pattern.
    positions=set()
    for k in tags:
        start=0;found=0
        while found<8:
            m=fast_search(RX[k],STARTS[k],s,lower,start)
            if m is None:break
            left=s.rfind('\n',0,m.start())+1
            right=s.find('\n',m.start())
            if right<0 or m.end()<=right:positions.add(left);found+=1
            if right<0:break
            start=right+1
    excerpts=[];prompt=[]
    for left in sorted(positions)[:8]:
        right=s.find('\n',left)
        if right<0:right=len(s)
        line=s[left:right];ln=s.count('\n',0,left)+1
        lt=[k for k in tags if RX[k].search(line)]
        excerpts.append({'line':ln,'tags':lt,'text':line[:700]})
        before=left
        for _ in range(2):before=s.rfind('\n',0,max(0,before-1))+1
        after=s.find('\n',right+1)
        if after<0:after=len(s)
        context=' '.join(s[before:after].splitlines())
        if PROMPT.search(context) and len(prompt)<6:prompt.append({'line':ln,'text':context[:1100]})
    return excerpts,prompt

def git(*args):
    return subprocess.check_output(['git',*args],cwd=ROOT)

def frozen_git(*args):
    refs=(HERE/'history_refs.txt').read_text(encoding='utf-8-sig').splitlines()
    tips=sorted({line.split()[0] for line in refs if line.strip()})
    return subprocess.check_output(['git',*args,'--stdin'],input=('\n'.join(tips)+'\n').encode(),cwd=ROOT)

def scan():
    t=time.time()
    head=git('rev-parse','7510e38f9c77').decode().strip()
    # rev-list yields each reachable object once and one observed path. Paths
    # are discovery names, not an exhaustive rename history.
    raw=frozen_git('rev-list','--objects').decode('utf8','replace').splitlines()
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
        lower=s.lower().replace('\u0131','i').replace('\u017f','s')
        tags=[k for k,rx in RX.items() if (len(s)!=len(lower) or any(t in lower for t in CHEAP[k])) and fast_search(rx,STARTS[k],s,lower)]
        if idx%10000==0:print(json.dumps({'phase':'scan','at':idx,'hits':len(hits),'seconds':round(time.time()-t,1)}),flush=True)
        if not tags:continue
        excerpts,prompt=early_excerpts(s,lower,tags)
        topics=[k for k,rx in TRX.items() if fast_search(rx,TCHEAP[k],s,lower)]
        by_topic.update(topics);counts.update(tags)
        title=next((x.lstrip('# ').strip()[:200] for x in s[:60000].splitlines()[:70] if x.startswith('# ')),path)
        hits.append({'oid':oid,'path':path,'bytes':size,'title':title,'tags':tags,'topics':topics,'excerpts':excerpts,'prompt_candidates':prompt})
    p.stdin.close();p.wait()
    commits=frozen_git('log','--format=%H%x09%aI%x09%s').decode('utf8','replace').splitlines()
    cm=[]
    for c in commits:
        if any(rx.search(c.split('\t',2)[-1]) for rx in RX.values()):cm.append(c)
    # Current paths and exact object identities provide stable local retrieval.
    current={}
    for line in git('ls-tree','-r',head).decode('utf8','replace').splitlines():
        meta,path=line.split('\t',1);current[path]=meta.split()[2]
    for h in hits:h['baseline_at_discovery_path']=current.get(h['path'])==h['oid']
    with (HERE/'history_hits.jsonl.gz').open('wb') as rawfile:
        with gzip.GzipFile(filename='',fileobj=rawfile,mode='wb',mtime=0) as f:
            for h in hits:f.write((json.dumps(h,ensure_ascii=False)+'\n').encode('utf8'))
    (HERE/'history_commits.txt').write_text('\n'.join(cm)+'\n',encoding='utf8')
    stats={'baseline':head,'frozen_refs_file':'history_refs.txt','reachable_object_count':len(raw),'candidate_text_objects':len(objects),'counts':dict(counts),'hit_blobs':len(hits),'hit_discovery_paths':len({h['path'] for h in hits}),'baseline_discovery_hits':sum(h['baseline_at_discovery_path'] for h in hits),'topic_cooccurrence_blobs':dict(by_topic),'reachable_commits':len(commits),'matching_commit_subjects':len(cm),'seconds':round(time.time()-t,2),'patterns':PATTERNS,'topic_patterns':TOPICS,'limits':['Reachable local Git refs at frozen history_refs.txt only; shallow boundaries may omit older commits.','Each object has one Git discovery path; rename aliases not exhaustively reconstructed.','Whole-blob tag presence is scanned, but excerpts are the first eight matching lines only; contextual prompt candidates come from these lines. Later prompts can be absent from the excerpt index.','Lexical hits and prompt candidates are not validated claims or distinct sessions.','Binary attachments and inaccessible chats are outside the text universe.','Topic co-occurrence is evidence for retrieval, not a mathematical map.']}
    (HERE/'history_stats.json').write_text(json.dumps(stats,indent=2,ensure_ascii=False)+'\n',encoding='utf8')
    print(json.dumps(stats,ensure_ascii=False),flush=True)

if __name__=='__main__':scan()

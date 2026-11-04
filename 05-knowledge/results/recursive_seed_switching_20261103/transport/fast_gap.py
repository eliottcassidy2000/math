"""Exact compiled gap counter, n<=30 and bins<=3; binary lives only in /tmp.

API count matches the inherited gap.count. Engine provides persistent or batch
queries. This module does not change callers or assert target acceptance.
"""
from functools import lru_cache
from hashlib import sha256
from pathlib import Path
import atexit
import os
import subprocess
import tempfile

SOURCE = Path(__file__).with_suffix('.cpp')


def executable():
    digest = sha256(SOURCE.read_bytes()).hexdigest()[:20]
    binary = Path(tempfile.gettempdir()) / ('math-fast-gap-' + digest)
    if not binary.exists():
        temporary = str(binary) + '.' + str(os.getpid())
        subprocess.run(['c++','-O3','-std=c++17',str(SOURCE),'-o',temporary],check=True)
        os.replace(temporary,binary)
    return binary


def encode(rows,feedback,bins=1,starts=(),ends=()):
    values = [len(rows),bins,len(feedback),len(starts),len(ends),*rows,*feedback]
    for pair in (*starts,*ends):
        values.extend(pair)
    return ' '.join(map(str,values))+'\n'


class Engine:
    def __init__(self,sparse=False):
        self.process = subprocess.Popen([str(executable())]+(['--sparse'] if sparse else []),stdin=subprocess.PIPE,
                                        stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
        atexit.register(self.close)

    def count(self,rows,feedback,bins=1,starts=(),ends=()):
        self.process.stdin.write(encode(rows,feedback,bins,starts,ends))
        self.process.stdin.flush()
        result = self.process.stdout.readline()
        if not result:
            error = self.process.stderr.read()
            raise RuntimeError('Compiled gap counter failed: '+error)
        return int(result)

    def batch(self,queries):
        return [self.count(*query) for query in queries]

    def close(self):
        if self.process.poll() is None:
            self.process.stdin.close()
            self.process.wait()


_engine = None


@lru_cache(None)
def count(rows,feedback,bins=1,starts=(),ends=()):
    global _engine
    if _engine is None:
        _engine = Engine()
    return _engine.count(rows,feedback,bins,starts,ends)


def self_test():
    """Literal small controls plus a few inherited gap-counter comparisons."""
    import importlib.util
    from itertools import combinations,permutations
    import json
    import time
    here=SOURCE.parent
    old=here.parents[1]/'zero_ready_component_returns_20261103/interface/gap.py'
    spec=importlib.util.spec_from_file_location('old_gap_for_audit',old)
    gap=importlib.util.module_from_spec(spec);spec.loader.exec_module(gap)
    def need(ok,message):
        if not ok:raise RuntimeError(message)
    def literal(rows,k,starts,ends):
        answer=0
        for word in permutations(range(len(rows))):
            for sizes in gap.compositions(len(rows),k):
                words=[];pos=0
                for size in sizes:words.append(word[pos:pos+size]);pos+=size
                if any(not words[b] or words[b][0]!=v for b,v in starts):continue
                if any(not words[b] or words[b][-1]!=v for b,v in ends):continue
                if all(all(rows[a]>>b&1 for a,b in zip(w,w[1:])) for w in words):answer+=1
        return answer
    def feedback(rows):
        n=len(rows)
        for size in range(n+1):
            for fb in combinations(range(n),size):
                rest=[v for v in range(n) if v not in fb]
                order=sorted(rest,key=lambda v:-sum(rows[v]>>u&1 for u in rest))
                if all(rows[a]>>b&1 for i,a in enumerate(order) for b in order[i+1:]):return fb
        raise RuntimeError('Missing feedback')
    engine=Engine();sparse=Engine(sparse=True);tested=0;graphs=[]
    for n in range(6):
        pairs=list(combinations(range(n),2));limit=1<<len(pairs)
        codes=range(limit) if n<=3 else sorted({0,1,limit//3,limit//2,limit-1})
        for code in codes:
            rows=[0]*n
            for j,(a,b) in enumerate(pairs):
                u,v=(a,b) if code>>j&1 else (b,a);rows[u]|=1<<v
            rows=tuple(rows);fb=feedback(rows);graphs.append((n,code,fb))
            queries=[(k,(),()) for k in (1,2,3)]
            if n:
                queries += [(1,((0,0),),()),(2,(),((0,0),)),(3,((0,0),),((0,0),)),(2,((0,0),),((1,0),))]
            if n>=2:queries += [(2,((1,1),),((0,0),)),(3,((0,0),(1,1)),())]
            for k,starts,ends in queries:
                got=engine.count(rows,fb,k,starts,ends)
                want=literal(rows,k,starts,ends)
                need(got==want,'Literal ordered-bin mismatch')
                need(got==sparse.count(rows,fb,k,starts,ends),'Sparse/dense agreement')
                need(got==gap.count(rows,fb,k,starts,ends),'Old gap mismatch')
                tested+=1
    sparse.close()
    rejected=0
    for bad in ('31 1 0 0 0\n','0 4 0 0 0\n',encode((2,4,1),()),encode((2,4,1),(0,),1,((0,0),(0,1)),())):
        result=subprocess.run([str(executable())],input=bad,text=True,capture_output=True)
        need(result.returncode==2,'Required input guard rejection')
        rejected+=1
    additional=[]
    base=(2,12,9,1)
    for tail in ('010101','111000'):
        rows=base
        for bit in tail:
            out=1 if bit=='0' else 0;n=len(rows)
            rows=tuple(x|(0 if out>>i&1 else 1<<n) for i,x in enumerate(rows))+(out,)
        for query in ((1,(),()),(2,(),()),(3,((0,0),),((1,3),))):
            k,starts,ends=query
            expected=gap.count(rows,(0,1,2,3),k,starts,ends)
            got=engine.count(rows,(0,1,2,3),k,starts,ends)
            need(got==expected,'Feedback-four comparison')
            additional.append(dict(tail=tail,bins=k,starts=starts,ends=ends,count=got))
    # Actual first repaired target, using its independently declared feedback.
    inp=here.parent/'selector/named_repair.json';data=json.loads(inp.read_text())
    target=data['target'];rows=tuple(target['rows']);names=target['names']
    fb=tuple(names.index(v) for v in (4,7,14,39))
    t0=time.monotonic();H=engine.count(rows,fb);P2=engine.count(rows,fb,2)
    elapsed=time.monotonic()-t0
    need(H==target['H']==13835447,'Actual repaired target H')
    # Existing Python implementation supplies an independent path on this H.
    need(H==gap.count(rows,fb),'Actual repaired target versus Python gap')
    engine.close()
    return dict(status='FINITE-EXACT implementation controls; no new coverage claim',
        universe='All labelled tournaments orders0..3 plus five fixed codes at each order4,5; all bins1..3 and listed endpoint patterns.',
        small_graphs=len(graphs),literal_queries=tested,feedback_four_cases=additional,
        sparse_dense_agreement_queries=tested,rejected_invalid_inputs=rejected,
        actual_repair=dict(order=len(rows),H=H,P2=P2,feedback_names=[4,7,14,39]),
        arithmetic_bound='Each intermediate and final count <=30!*496<2^117<2^128.',
        source_sha256={str(p.relative_to(here.parents[1])):sha256(p.read_bytes()).hexdigest() for p in (SOURCE,Path(__file__),old,inp)})


if __name__ == '__main__':
    import sys
    # JSON-lines query format: [rows,feedback,bins,starts,ends].
    import json
    if sys.argv[1:]==['--self-test']:
        print(json.dumps(self_test(),sort_keys=True,indent=2))
    else:
        engine=Engine()
        for line in sys.stdin:
            print(engine.count(*json.loads(line)),flush=True)

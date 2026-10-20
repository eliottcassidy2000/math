"""Bounded endpoint-two coverage; exact labelled certificates, no class census.

All source presentations have standard blue closed P and endpoint scores
(n-3,2), n=6..9. Construct the universe from all blue middle presentations
and every marked exceptional port. First try the proved local ledger bank;
complete three-K search is used ONLY on remaining small presentations.
"""
from collections import Counter
from hashlib import sha256
from pathlib import Path
import json
import probe as P
import probe_search as SEARCH


def candidates(n, exceptional):
    a,b,x,c,z = 0,1,2,n-2,n-1
    if exceptional == n-3:
        d = tuple(range(3,n-2))
        yield 'B_open_penultimate', ((a,b,c,z,x)+d,
                                    (b,c,a,z,x)+d,
                                    (b,x)+d+(c,a,z))
    if exceptional == c:
        interior = tuple(range(3,n-2))
        s = (a,b,x,c)+interior[::-1]+(z,)
        yield 'L2_last_port', ((a,x)+interior+(c,z,b),
                              s[::-1],(b,z)+interior+(c,x,a))
        for cut in range(4,n-1):
            d, u = tuple(range(3,cut)),tuple(range(cut,n-1))
            yield 'A_last_port_cut_'+str(cut), (
                (a,x)+d+u+(z,b), d+(a,x)+u+(z,b),d+(a,b,x)+u+(z,))
    if exceptional == x and n>=7:
        interior,y = tuple(range(3,n-3)),n-3
        yield 'C_closed_second_port', (
            (a,b,x)+interior+(z,y,c),
            (a,x)+interior+(z,y,b,c),
            (a,x)+interior+(y,b,c,z))
    for name,words in P.cyclic_words(n):
        yield name,words
        yield name+'_reflected',P.reflect_words(n,words)


def independent_replay(rows, words):
    """Set-of-directed-arcs replay, independent of P's bitmask K routine."""
    n = len(rows)
    initial = {(u,v) for u in range(n) for v in range(n) if rows[u]>>v&1}
    now = initial
    witnesses = []
    for word in words:
        P.need(set(word)==set(range(n)) and len(word)==n,'permutation')
        directed_path = set(zip(word,word[1:]))
        P.need(directed_path <= now,'independent path directedness')
        reversal = dict(zip(word,word[::-1]))
        failures = [(u,v,reversal[v],reversal[u]) for u,v in sorted(now)
                    if (reversal[v],reversal[u]) not in now]
        P.need(bool(failures),'independent blackness')
        witnesses.append(failures[0])
        reverse = {(v,u) for u,v in now}
        now = (reverse - {(v,u) for u,v in directed_path}) | directed_path
    original_path = {(i,i+1) for i in range(n-1)}
    target = ({(v,u) for u,v in initial} - {(v,u) for u,v in original_path}) | original_path
    P.need(now==target,'independent literal target')
    return witnesses


def main():
    summaries, residual_certificates, unresolved = [],[],[]
    source_digest, certificate_digest = sha256(),sha256()
    for n in range(6,10):
        stats, names = Counter(),Counter()
        seen = set()
        for t,u,v in P.endpoint_two(n):
            P.need(t not in seen,'duplicate source')
            seen.add(t)
            p = tuple(range(n))
            P.need(P.path(t,p) and P.blue(t,p) and P.arc(t,n-1,0),'blue closed source')
            P.need(t[0].bit_count()==n-3 and t[-1].bit_count()==2,'scores')
            P.need(v==n-1-u,'marked port reflection')
            stats['universe'] += 1
            source_digest.update(json.dumps([n,t,u,v],separators=(',',':')).encode()+b'\n')
            partner = P.complement(t,p)
            if not P.strong(t) or not P.strong(partner):
                stats['excluded_nonstrong_endpoint'] += 1
                continue
            if not any(not P.blue(t,q) for q in P.paths(t)) or not any(not P.blue(partner,q) for q in P.paths(partner)):
                stats['excluded_pure_endpoint'] += 1
                continue
            stats['strong_mixed'] += 1
            selected = None
            for name,words in candidates(n,u):
                result = P.evaluate(t,words)
                if result['success']:
                    selected = (name,words)
                    names[name] += 1
                    stats['template_covered'] += 1
                    break
            if selected is None:
                print('template residual',n,u,list(t),flush=True)
                stats['template_residual'] += 1
                answer = SEARCH.search(t)
                if not answer['best']:
                    unresolved.append({'n':n,'rows':t,'u':u,'v':v,'search':answer})
                    continue
                cert = answer['best'][0]
                words = tuple(tuple(q) for q in cert['words'])
                selected = ('finite_three_K',words)
                residual_certificates.append({
                    'n':n,'rows':t,'partner_rows':partner,'u':u,'v':v,
                    'words':words,'independent_anti_failures':independent_replay(t,words),
                    'endpoint_excesses':[s['endpoint_excess'] for s in cert['evaluation']['steps']],
                    'complete_search_counts':{k:answer[k] for k in
                        ['first_black_paths','last_black_paths','pairs_examined','three_K_certificates']}})
                stats['finite_certified'] += 1
            name,words = selected
            independent_replay(t,words)
            certificate_digest.update(json.dumps([n,t,u,v,name,words],separators=(',',':')).encode()+b'\n')
        P.need(stats['universe']=={6:12,7:64,8:320,9:3072}[n],
               'the complete endpoint-two universe has the expected size')
        P.need(stats['universe']==stats['strong_mixed'],
               'the finite theorem covers every source, without filtering any out')
        summaries.append({'n':n,**dict(stats),'selected_templates':dict(sorted(names.items()))})
        print(json.dumps(summaries[-1]),flush=True)
    P.need(not unresolved,'a strong mixed endpoint-two presentation has no certificate')
    receipt = {
        'status':'FINITE-EXACT: every retained labelled presentation has three black K moves to its literal blue partner',
        'universe':'All standard closed blue endpoint-score (n-3,2) presentations at n=6,7,8,9; every blue middle and every allowed marked incoming port u=2,...,n-2; no isomorphism quotient',
        'inherited_filters':'Both literal endpoints strong and mixed; these are tested, not assumed',
        'checks':'Every selected ledger is independently replayed with directed arc sets, actual reversal failure, and literal target equality. Complete labelled three-K search is used only after template failure.',
        'source_hashes':{Path(f).name:sha256(Path(f).read_bytes()).hexdigest() for f in [__file__,P.__file__,SEARCH.__file__]},
        'universe_sha256':source_digest.hexdigest(),'certificate_stream_sha256':certificate_digest.hexdigest(),
        'summaries':summaries,'finite_residual_certificates':residual_certificates,'unresolved':unresolved}
    Path(__file__).with_suffix('.json').write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8')
    print('PASS',sum(s.get('strong_mixed',0) for s in summaries),'presentations;',len(residual_certificates),'finite exceptions',flush=True)


if __name__ == '__main__':
    main()

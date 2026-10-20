"""Exact endpoint-(n-4,3) marked presentations and inherited ledger diagnostics.

No isomorphism quotient. This imports the frozen endpoint-two bit operations,
not its theorem assumptions. All blackness predicates are retested literally.
"""
from collections import Counter, defaultdict
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import importlib.util
import json

OLD_PATH = Path(__file__).parents[1]/'endpoint_two_repairs_20261020g'/'probe.py'
spec = importlib.util.spec_from_file_location('endpoint_two_kernel',OLD_PATH)
OLD = importlib.util.module_from_spec(spec)
spec.loader.exec_module(OLD)
need,arc,edge,edges,path,blue = OLD.need,OLD.arc,OLD.edge,OLD.edges,OLD.path,OLD.blue
anti_witness,complement,strong,paths = OLD.anti_witness,OLD.complement,OLD.strong,OLD.paths
compile_ledger,evaluate = OLD.compile_ledger,OLD.evaluate


def endpoint_three(n):
    for middle in OLD.blue_presentations(n-2):
        for incoming in combinations(range(2,n-1),2):
            outgoing = tuple(sorted(n-1-v for v in incoming))
            rows = [sum(1<<v for v in range(1,n-1))]
            rows += [(row<<1)|(1<<(n-1)) for row in middle]
            rows += [1]
            for u in incoming:
                rows[0] ^= 1<<u
                rows[u] ^= 1
            for v in outgoing:
                rows[v] ^= 1<<(n-1)
                rows[-1] ^= 1<<v
            yield tuple(rows),incoming,outgoing


def bank(n):
    a,b,x,c,z = 0,1,2,n-2,n-1
    raw = list(OLD.cyclic_words(n))
    interior = tuple(range(3,n-2))
    s = (a,b,x,c)+interior[::-1]+(z,)
    raw.append(('L2',((a,x)+interior+(c,z,b),s[::-1],(b,z)+interior+(c,x,a))))
    for cut in range(4,n-1):
        d,r = tuple(range(3,cut)),tuple(range(cut,n-1))
        raw.append(('A_cut_'+str(cut),((a,x)+d+r+(z,b),d+(a,x)+r+(z,b),d+(a,b,x)+r+(z,))))
        raw.append(('B_cut_'+str(cut),((a,b)+r+(z,x)+d,(b,)+r+(a,z,x)+d,(b,x)+d+r+(a,z))))
    if n>=7:
        interior,y = tuple(range(3,n-3)),n-3
        raw.append(('C',((a,b,x)+interior+(z,y,c),
                         (a,x)+interior+(z,y,b,c),
                         (a,x)+interior+(y,b,c,z))))
    seen = set()
    for name,words in raw:
        for suffix,w in [('',words),('_reflected',OLD.reflect_words(n,words))]:
            if w in seen:
                continue
            seen.add(w)
            compiled = compile_ledger(n,w)
            need(not compiled['guard_conflicts'] and compiled['path_XOR_is_standard'], 'structural ledger consistency')
            yield name+suffix,w


def profile(t,u,v):
    n = len(t)
    return ' '.join([
        'open' if arc(t,1,n-2) else 'closed',
        'overlap='+str(len(set(u)&set(v))),
        'Unear='+''.join('1' if w in u else '0' for w in [2,3,n-3,n-2])])


def pair_kind(n,pair):
    u,v = pair
    if 0 in pair:
        return 'first-endpoint'
    if n-1 in pair:
        return 'last-endpoint'
    if abs(u-v)==1:
        return 'middle-path'
    return 'middle-chord'


def source_checks(t,u,v):
    n = len(t)
    p = tuple(range(n))
    need(path(t,p) and blue(t,p) and arc(t,n-1,0),'closed standard blue source')
    need(t[0].bit_count()==n-4 and t[-1].bit_count()==3,'endpoint scores')
    need(tuple(w for w in range(1,n-1) if arc(t,w,0))==u,'incoming mark exact')
    need(tuple(w for w in range(1,n-1) if arc(t,n-1,w))==v,'outgoing mark exact')
    need(strong(t),'closed spanning word gives a strong source')


def main():
    summaries, failures, examples, minimal_residuals, excluded = [],{}, {},[],[]
    all_digest = sha256()
    for n in range(6,10):
        stats, profiles, selected = Counter(),defaultdict(Counter),Counter()
        candidates = list(bank(n))
        candidate_stats = defaultdict(Counter)
        seen = set()
        for t,u,v in endpoint_three(n):
            source_checks(t,u,v)
            need(t not in seen,'duplicate marked source')
            seen.add(t)
            all_digest.update(json.dumps([n,t,u,v],separators=(',',':')).encode()+b'\n')
            stats['universe'] += 1
            other = complement(t,tuple(range(n)))
            if not strong(other):
                stats['nonstrong_partner'] += 1
                excluded.append({'n':n,'rows':t,'partner_rows':other,'incoming':u,'outgoing':v,
                                 'reason':'nonstrong partner','partner_HP':sum(1 for q in paths(other)),
                                 'partner_black_HP':sum(not blue(other,q) for q in paths(other))})
                continue
            if not any(not blue(t,q) for q in paths(t)) or not any(not blue(other,q) for q in paths(other)):
                stats['pure_endpoint'] += 1
                continue
            stats['strong_mixed'] += 1
            tag = profile(t,u,v)
            profiles[tag]['total'] += 1
            chosen = None
            for name,words in candidates:
                result = evaluate(t,words)
                if result['success']:
                    candidate_stats[name]['all_black'] += 1
                    if chosen is None:
                        chosen = (name,words)
                elif 'invalid_step' in result:
                    step = result['invalid_step']
                    kinds = '+'.join(sorted({pair_kind(n,pair) for pair in result['missing_arcs']}))
                    candidate_stats[name]['illegal_'+str(step)+'_'+kinds] += 1
                else:
                    candidate_stats[name]['legal_blue_'+''.join(map(str,result['blue_steps']))] += 1
            if chosen is not None:
                stats['bank_covered'] += 1
                profiles[tag]['covered'] += 1
                selected[chosen[0]] += 1
            else:
                stats['bank_residual'] += 1
                profiles[tag]['residual'] += 1
                key = str(n)+' '+tag
                item = {'n':n,'rows':t,'partner_rows':other,'incoming':u,'outgoing':v,'profile':tag}
                if n==6:
                    minimal_residuals.append(item)
                if key not in examples:
                    examples[key] = item
        need(stats['universe']=={6:12,7:96,8:640,9:7680}[n],'expected marked universe size')
        summaries.append({'n':n,**dict(stats),'selected':dict(sorted(selected.items())),
                          'profiles':{k:dict(v) for k,v in sorted(profiles.items())}})
        failures[n] = {k:dict(v) for k,v in sorted(candidate_stats.items())}
        print(json.dumps({'n':n,**dict(stats)}),flush=True)
    receipt = {'status':'FINITE-EXACT diagnostic, not an all-order detour or bridge theorem',
        'universe':'Every closed standard-blue endpoint-score (n-4,3) presentation at n6..9; arbitrary blue middle plus every two-element U subset of positions2..n-2; V is the reflected set; overlaps retained',
        'filters':'For ledger diagnostics, both literal endpoints are required strong and mixed; exclusions counted explicitly',
        'source_hashes':{'current_probe':sha256(Path(__file__).read_bytes()).hexdigest(),
                         'inherited_endpoint_two_kernel':sha256(OLD_PATH.read_bytes()).hexdigest()},
        'universe_sha256':all_digest.hexdigest(),'summaries':summaries,
        'ledger_failure_types':failures,'one_residual_per_profile':list(examples.values()),
        'minimal_residuals':minimal_residuals,'excluded_presentations':excluded}
    Path(__file__).with_suffix('.json').write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8')


if __name__=='__main__':
    main()

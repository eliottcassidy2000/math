"""Exact support-cover test for the 21 marked Fano junction probes.

Enumerates all 3^7 labelled words, all105 junction targets and the21-row
bank, without graph-realizability assumptions. Explicit checks survive -O.
"""
from itertools import combinations,product
from collections import Counter
from hashlib import sha256
from pathlib import Path
import json

HERE=Path(__file__).resolve().parent


def need(x,msg):
    if not x:
        raise RuntimeError(msg)


def canonical_word(w):
    seen={}
    return tuple(seen.setdefault(c,len(seen)) for c in w)


def label(triple,pairs):
    return tuple(sorted(triple)),tuple(sorted(tuple(sorted(p)) for p in pairs))


def show(j):
    t,ps=j
    return ''.join(map(str,t))+';'+','.join(''.join(map(str,p)) for p in ps)


def junctions():
    for triple in combinations(range(1,8),3):
        rest=[v for v in range(1,8) if v not in triple]
        for mate in rest[1:]:
            yield label(triple,[(rest[0],mate),tuple(v for v in rest if v not in (rest[0],mate))])


def accepts(j,w):
    t,ps=j
    return len({w[v-1] for v in t})==3 and all(w[a-1]==w[b-1] for a,b in ps)


def mapped(j,mapping):
    t,ps=j
    return label([mapping[v] for v in t],[[mapping[v] for v in p] for p in ps])


def main():
    all_words=list(product(range(3),repeat=7))
    admissible=[w for w in all_words if all(w.count(c)%2 for c in range(3))]
    orbits=sorted({canonical_word(w) for w in admissible})
    targets=list(junctions())
    need((len(all_words),len(admissible),len(orbits),len(targets))==(2187,546,91,105),'universe sizes')
    supports={j:{i for i,w in enumerate(orbits) if accepts(j,w)} for j in targets}
    need(all(len(s)==9 for s in supports.values()),'junction orbit support size')
    bank={j for j in targets if j[0][0]^j[0][1]^j[0][2]==0}
    # Reconstruct the same bank from seven parallel matching directions.
    flags=set()
    for direction in range(1,8):
        pairs=sorted({tuple(sorted((v,v^direction))) for v in range(1,8) if v!=direction})
        for selected in pairs:
            flags.add(label((direction,)+selected,[p for p in pairs if p!=selected]))
    need(bank==flags and len(bank)==21,'Fano bank differs from parallel-flag bank')
    bank=sorted(bank)
    records=[]
    distribution=Counter()
    fingerprints={}
    examples={}
    for target in targets:
        support=supports[target]
        zero=[j for j in bank if not(supports[j]&support)]
        covered=set().union(*(supports[j] for j in zero))
        survivors=set(range(91))-covered
        outside=survivors-support
        need(support<=survivors,'zero probes excluded the target itself')
        # Independent full-word check, rather than only normalized orbits.
        full_survivors=[w for w in admissible if not any(accepts(j,w) for j in zero)]
        need({canonical_word(w) for w in full_survivors}=={orbits[i] for i in survivors},'word/orbit survivor mismatch')
        need(len(full_survivors)==6*len(survivors),'global-colour normalization mismatch')
        key=(target in bank,len(zero),len(survivors),len(outside))
        distribution[key]+=1
        witness=min(orbits[i] for i in outside)
        need(not accepts(target,witness) and all(not accepts(j,witness) for j in zero),'surviving witness incorrect')
        fingerprint=tuple(len(support&supports[j]) for j in bank)
        need(fingerprint not in fingerprints,'two target junctions have equal full-bank readouts')
        fingerprints[fingerprint]=target
        row=dict(target=show(target),in_bank=target in bank,
                 zero_probes=[show(j) for j in zero],surviving_orbits=len(survivors),
                 outside_target_orbits=len(outside),surviving_labelled_words=len(full_survivors),
                 least_outside_word=''.join(map(str,witness)),
                 outside_orbit_words=[''.join(map(str,orbits[i])) for i in sorted(outside)],
                 full_bank_orbit_measurements=list(fingerprint))
        records.append(row)
        examples.setdefault(str(target in bank),row)
    need(distribution==Counter({(True,10,33,24):21,(False,8,39,30):84}),'support-cover distribution changed')

    # Independently replay the incoming complete14-probe certificate.
    raw=['034;12,56','034;15,26','034;16,25','056;12,34',
         '134;02,56','134;05,26','134;06,25','156;02,34',
         '234;01,56','256;01,34','256;03,14','256;04,13',
         '345;01,26','346;01,25']
    cover=[]
    for text in raw:
        t,ps=text.split(';')
        cover.append(label([int(v)+1 for v in t],[[int(v)+1 for v in p] for p in ps.split(',')]))
    incoming_verified=0
    for target in targets:
        t,ps=target
        image=tuple(t)+tuple(ps[0])+tuple(ps[1])
        mapping={i+1:v for i,v in enumerate(image)}
        probes=[mapped(j,mapping) for j in cover]
        need(all(not(supports[j]&supports[target]) for j in probes),'incoming probe overlaps target')
        union=set().union(*(supports[j] for j in probes))
        need(union==set(range(91))-supports[target],'incoming cover is incomplete')
        incoming_verified+=1

    # GL(3,2) explains the two uniform target strata, with the address mark fixed.
    actions=[]
    for a in range(1,8):
        for b in range(1,8):
            if b in (a,):continue
            for c in range(1,8):
                if c in (a,b,a^b):continue
                mapping={v:((a if v&1 else 0)^(b if v&2 else 0)^(c if v&4 else 0)) for v in range(1,8)}
                need(len(set(mapping.values()))==7,'linear address map is singular')
                actions.append(mapping)
    need(len(actions)==168,'GL3 order')
    base_in=next(j for j in targets if j in bank)
    base_out=next(j for j in targets if j not in bank)
    orbit_in={mapped(base_in,g) for g in actions}
    orbit_out={mapped(base_out,g) for g in actions}
    need(orbit_in==set(bank) and orbit_out==set(targets)-set(bank),'GL3 target orbit decomposition')
    report=dict(scope='all105 labelled junction targets; all2187 words; all91 parity-colour orbits',
                port_labels=list(range(1,8)),bank=[show(j) for j in bank],
                cover_distribution=[dict(in_bank=k[0],zero_rows=k[1],surviving_orbits=k[2],outside_target_orbits=k[3],targets=v)
                                    for k,v in sorted(distribution.items(),reverse=True)],
                examples=examples,all_targets=records,
                full_bank_distinct_target_fingerprints=len(fingerprints),
                incoming_14_probe_cover_targets_verified=incoming_verified,
                GL3_order=len(actions),GL3_target_orbit_sizes=[len(orbit_in),len(orbit_out)],
                source_sha256_lf=sha256(Path(__file__).read_text().encode()).hexdigest())
    output=json.dumps(report,indent=2)+'\n'
    (HERE/'arithmetic_junction_certificate.json').write_text(output,encoding='utf-8')
    print('Golden21 zero-support distributions:',report['cover_distribution'])
    print('Distinct complete21 target fingerprints:',len(fingerprints))
    print('Incoming14-probe covers independently replayed:',incoming_verified)
    print('GL3 target orbits:',report['GL3_target_orbit_sizes'])
    print('PASS; this is not a physical counterexample or a full21 tensor-certification theorem.')


if __name__=='__main__':
    main()

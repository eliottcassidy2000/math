"""Entrywise reconciliation with two independently published order-eight banks."""
from pathlib import Path
from collections import Counter
import argparse, csv, gzip, hashlib, io, json, struct, subprocess, tempfile
import check as core

HERE = Path(__file__).resolve().parent
REF = '64912fc832'
FIRST = '05-knowledge/results/metagraph_order8_20261023/'
SECOND = '05-knowledge/results/session38_incidence_ports_electrical_20261023/wiring/'

def saved(path):
    return subprocess.check_output(['git', 'show', REF+':'+path], cwd=core.ROOT)

def table(data): return list(csv.DictReader(io.StringIO(data.decode()), delimiter='\t'))

def main():
    p = argparse.ArgumentParser(); p.add_argument('--prefix', type=Path); p.add_argument('--write', action='store_true'); args = p.parse_args()
    temporary = None
    if args.prefix is None:
        temporary = tempfile.TemporaryDirectory(prefix='incidence-cross-'); base = Path(temporary.name)
        fields = ('class', 'canon', 'H', 'Aut', 'tc', 'bm', 'sigmaClass')
        records = list(csv.DictReader(core.CENSUS.open(), delimiter='\t'))
        (base/'input.txt').write_text(''.join(' '.join(r[k] for k in fields)+'\n' for r in records))
        subprocess.run(['clang++', '-O2', '-std=c++17', str(HERE/'wiring.cpp'), '-o', str(base/'producer')], check=True)
        args.prefix = base/'native'
        subprocess.check_output([str(base/'producer'), str(base/'input.txt'), str(args.prefix)])
    own = args.prefix.with_suffix('.owners').read_bytes()
    rawports = args.prefix.with_suffix('.ports').read_bytes()
    ports = Counter()
    for k, v in struct.iter_unpack('<QI', rawports): ports[k&8191, k>>13&8191, k>>26&7, k>>29&7, k>>32&1] += v
    edges = Counter()
    for row in table(gzip.decompress((HERE/'n8_merged_wiring.tsv.gz').read_bytes())):
        edges[int(row['source_class_min']), int(row['target_class_min']), int(row['blue'])] += int(row['literal_lines'])
    native = json.loads((HERE/'receipt.json').read_text())
    vertices = {u for u, v, b in edges} | {v for u, v, b in edges}
    for scale in (1, 2):
        fake = edges.copy()
        for u, v, b, d in native['ports']['reflected_counterfeit']['merged_edits']: fake[u, v, b] += scale*d
        core.need(min(fake.values()) >= 0 and len(core.components(vertices, [(u, v, w) for (u, v, b), w in fake.items()])) == 1, 'counterfeit lost connectivity or nonnegativity')
    incoming = gzip.decompress(saved(FIRST+'owner_map.u16.gz'))
    core.need(own == incoming, 'first incoming owner map differs')
    incoming_ports = Counter()
    for k, v in struct.iter_unpack('<QI', gzip.decompress(saved(FIRST+'ports.bin.gz'))):
        k, blue = divmod(k, 2); k, b = divmod(k, 8); k, a = divmod(k, 8); c, t = divmod(k, 6880)
        incoming_ports[c, t, a, b, blue] += v
    core.need(incoming_ports == ports, 'first incoming full port bank differs')
    nodes = {int(r['merged_id']): min(int(r['class_a']), int(r['class_b'])) for r in table(saved(FIRST+'merged_nodes.tsv'))}
    ee = Counter()
    for r in table(gzip.decompress(saved(FIRST+'merged_edges.tsv.gz'))):
        u, v = sorted((nodes[int(r['source'])], nodes[int(r['target'])]))
        for col, blue in (('blue', 1), ('black', 0)):
            if int(r[col]): ee[u, v, blue] += int(r[col])
    core.need(ee == edges, 'first incoming merged wiring differs')
    # Session38 uses canonical-code ranks and a descending path. Convert both.
    classes = table(gzip.decompress(saved(SECOND+'classes.tsv.gz')))
    old = {int(r['class_id']): int(r['input_class_id']) for r in classes}
    raw = gzip.decompress(saved(SECOND+'class_of_mask.u16le.gz'))
    theirs = [old[v[0]] for v in struct.iter_unpack('<H', raw)]
    ours = [v[0] for v in struct.iter_unpack('<H', own)]
    cells = [(i, j) for i in range(8) for j in range(i+2, 8)]
    rev = [cells.index((7-j, 7-i)) for i, j in cells]
    chunks = [[sum(((m >> b)&1) << rev[7*k+b] for b in range(7)) for m in range(128)] for k in range(3)]
    F = (1 << 21)-1
    core.need(all(ours[m] == theirs[F ^ chunks[0][m&127] ^ chunks[1][m>>7&127] ^ chunks[2][m>>14]] for m in range(F+1)), 'second incoming mask gauge conversion differs')
    pp = Counter()
    for r in table(gzip.decompress(saved(SECOND+'ports.tsv.gz'))):
        pp[old[int(r['source_class'])], old[int(r['target_class'])], int(r['start_score']), int(r['end_score']), int(r['blue'])] += int(r['multiplicity'])
    core.need(pp == ports, 'second incoming full port bank differs')
    out = {'status': 'FINITE-EXACT entrywise agreement', 'incoming_commit': REF,
           'literal_owners_each': len(ours), 'port_entries_each': len(ports), 'merged_colored_entries': len(edges),
           'first_owner_sha256': hashlib.sha256(own).hexdigest(), 'second_owner_sha256': hashlib.sha256(raw).hexdigest(),
           'second_gauge': 'descending_mask = all_ones XOR reflection(ascending_mask); canonical rank -> input_class_id',
           'native_counterfeit_scales_connected_and_nonnegative': [1, 2],
           'checks': core.CHECKS}
    # The separate electrical certificate is preserved as a pinned source; its
    # report is compared mathematically in AUDIT.md, not guessed JSON keys.
    out['incoming_electrical_receipt_sha256'] = hashlib.sha256(saved('05-knowledge/results/metagraph_order8_electrical_20261023/receipt.json')).hexdigest()
    dest = HERE/'cross_check.json'
    if args.write: dest.write_text(json.dumps(out, indent=2)+'\n')
    else: core.need(json.loads(dest.read_text()) == out, 'cross-check receipt differs')
    print(json.dumps(out, indent=2))
    if temporary is not None: temporary.cleanup()

if __name__ == '__main__': main()

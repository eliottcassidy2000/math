"""
lichess_tools.py -- independent cross-checks against the Syzygy tablebases served by tablebase.lichess.ovh.
Every response is cached in a JSON file (argv) so that the archived file is the reproducibility artifact.

  python3 lichess_tools.py validate <scratch_dir> <cache.json> [n_per_class]
      sample positions from the locally solved classes and compare WDL (category) for the side to move.
  python3 lichess_tools.py tempo <cache.json> <class> <symmetry>
      census of symmetric boards (fixed by sigma = board isometry o colour swap) in a 4-man class:
      class in {KPvKP, KRvKR, KQvKQ, KBvKB, KNvKN}; symmetry in {rank, rot, file, diag, rot90}.
      Prints the tempo value distribution (mover wins / draw / mover loses = mutual zugzwang) and lists the
      mutual zugzwangs.
Rate: a 0.12 s pause after each request (about one request per second including latency), retry on HTTP 429.
"""
import sys, json, time, os, random, itertools, pickle, urllib.request, urllib.parse
import chess

API = 'https://tablebase.lichess.ovh/standard?fen='

def load_cache(path):
    if os.path.exists(path):
        with open(path) as f:
            return json.load(f)
    return {}

def save_cache(path, cache):
    tmp = path + '.tmp'
    with open(tmp, 'w') as f:
        json.dump(cache, f)
    os.replace(tmp, path)

def query(fen, cache, path, delay=0.12):
    if fen in cache:
        return cache[fen]
    for attempt in range(8):
        try:
            r = urllib.request.urlopen(API + urllib.parse.quote(fen), timeout=30).read().decode()
            d = json.loads(r)
            rec = {k: d.get(k) for k in ('category', 'dtz', 'dtm', 'checkmate', 'stalemate', 'insufficient_material')}
            cache[fen] = rec
            if len(cache) % 50 == 0:
                save_cache(path, cache)
            time.sleep(delay)
            return rec
        except urllib.error.HTTPError as e:
            if e.code == 429:
                time.sleep(60); continue
            raise
        except Exception:
            time.sleep(5)
    raise RuntimeError('lichess query failed: ' + fen)

def wdl_of(rec):
    c = rec['category']
    if c in ('win', 'cursed-win', 'maybe-win'):
        return 1
    if c in ('loss', 'blessed-loss', 'maybe-loss'):
        return -1
    return 0

def validate(scr, path, n=120):
    import numpy as np
    import chesslib as cl
    cache = load_cache(path)
    random.seed(12)
    total_bad = 0
    for name in ('KQvK', 'KRvK', 'KBvK', 'KNvK', 'KPvK'):
        with open(f'{scr}/table_{name}.pkl', 'rb') as f:
            T = pickle.load(f)
        U = cl.Universe(T['pieces']); keys, val = T['keys'], T['val']
        idx = random.sample(range(len(keys)), n)
        bad = 0
        for i in idx:
            fen = U.fen(int(keys[i]))
            rec = query(fen, cache, path)
            if wdl_of(rec) != int(val[i]):
                bad += 1; print('  MISMATCH', name, fen, 'local', int(val[i]), 'syzygy', rec)
        print(f'[{name}] {n} random positions vs Syzygy: mismatches {bad}', flush=True)
        total_bad += bad
    # the assumption-based KBPvK sub-class: sample positions (both sides) and compare
    p = f'{scr}/kbpvk_sub.pkl'
    if os.path.exists(p):
        with open(p, 'rb') as f:
            S = pickle.load(f)
        U = cl.Universe(['K', 'B', 'P', 'k']); keys, val = S['keys'], S['val']
        idx = random.sample(range(len(keys)), n)
        bad = 0
        for i in idx:
            fen = U.fen(int(keys[i])); rec = query(fen, cache, path)
            if wdl_of(rec) != int(val[i]):
                bad += 1; print('  MISMATCH KBPvK', fen, 'local', int(val[i]), 'syzygy', rec)
        print(f'[KBPvK dark-B a-pawn sub-class, promotion assumption] {n} random positions vs Syzygy: mismatches {bad}', flush=True)
        # stress the assumption: positions with the pawn on a7 and White to move
        pa7 = [i for i in range(len(keys)) if (keys[i] & 1) and U.unkey(int(keys[i]))[0][2] == chess.A7]
        idx = random.sample(pa7, min(n, len(pa7))); bad = 0
        for i in idx:
            fen = U.fen(int(keys[i])); rec = query(fen, cache, path)
            if wdl_of(rec) != int(val[i]):
                bad += 1; print('  MISMATCH KBPvK a7', fen, 'local', int(val[i]), 'syzygy', rec)
        print(f'[KBPvK pawn on a7, White to move] {len(idx)} positions vs Syzygy: mismatches {bad}', flush=True)
        total_bad += bad
    save_cache(path, cache)
    print('TOTAL mismatches', total_bad)

# ----- symmetric boards -----
def iso(sym, sq):
    f, r = chess.square_file(sq), chess.square_rank(sq)
    if sym == 'rank':  f2, r2 = f, 7 - r
    elif sym == 'file': f2, r2 = 7 - f, r
    elif sym == 'rot':  f2, r2 = 7 - f, 7 - r
    elif sym == 'diag': f2, r2 = r, f
    elif sym == 'rot90': f2, r2 = 7 - r, f
    else: raise ValueError(sym)
    return chess.square(f2, r2)

def tempo(path, cls, sym):
    cache = load_cache(path)
    white_piece = {'KPvKP': 'P', 'KRvKR': 'R', 'KQvKQ': 'Q', 'KBvKB': 'B', 'KNvKN': 'N'}[cls]
    if white_piece == 'P' and sym not in ('rank', 'rot'):
        raise ValueError('pawn classes only admit rank reflection or rotation')
    rows = []
    squares = range(64)
    psq = [s for s in squares if 1 <= chess.square_rank(s) <= 6] if white_piece == 'P' else list(squares)
    for wk in squares:
        bk = iso(sym, wk)
        if bk == wk or chess.square_distance(wk, bk) <= 1:
            continue
        for wp in psq:
            bp = iso(sym, wp)
            if len({wk, bk, wp, bp}) < 4:
                continue
            b = chess.Board(None)
            b.set_piece_at(wk, chess.Piece.from_symbol('K')); b.set_piece_at(bk, chess.Piece.from_symbol('k'))
            b.set_piece_at(wp, chess.Piece.from_symbol(white_piece)); b.set_piece_at(bp, chess.Piece.from_symbol(white_piece.lower()))
            b.turn = chess.WHITE
            if not b.is_valid():
                continue
            fen = b.fen()
            rec = query(fen, cache, path)
            rows.append((fen, wdl_of(rec), rec))
    save_cache(path, cache)
    from collections import Counter
    c = Counter(v for _, v, _ in rows)
    print(f'[{cls} sigma={sym}] symmetric boards: {len(rows)}; mover wins {c[1]}, draw {c[0]}, mover loses (mutual zugzwang) {c[-1]}')
    mz = [(fen, rec) for fen, v, rec in rows if v == -1]
    for fen, rec in mz[:40]:
        print('   mutual zugzwang:', fen, rec)
    out = path.replace('.json', f'_{cls}_{sym}.txt')
    with open(out, 'w') as f:
        for fen, v, rec in rows:
            f.write(f'{fen}\t{v}\t{json.dumps(rec)}\n')
    return rows

if __name__ == '__main__':
    cmd = sys.argv[1]
    if cmd == 'validate':
        validate(sys.argv[2], sys.argv[3], int(sys.argv[4]) if len(sys.argv) > 4 else 120)
    elif cmd == 'tempo':
        tempo(sys.argv[2], sys.argv[3], sys.argv[4])

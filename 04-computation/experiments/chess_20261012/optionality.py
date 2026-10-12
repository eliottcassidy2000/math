"""
optionality.py -- "outcome first, optionality second" tested as a move-selection rule.

Winning side (White) in KRvK and KQvK: among value-preserving moves (those keeping the forced win), choose by a
secondary optionality score, with a deterministic tie-break; Black defends DTM-optimally (maximizing DTM).
Scores:
  O1 geometric: number of legal White moves at White's next turn, minimized over Black's replies;
  O2 strategic: number of win-preserving White moves at White's next turn, minimized over Black's replies.
Also: O2 restricted to DTM-non-increasing moves (weak progress) and to DTM-optimal moves (tie-break only).
We count, over all White-to-move winning positions, how many are mated within 100 moves (200 plies),
how many cycle (position repeats), and the longest successful line.
Defending side (Black) in KPvK drawn positions: any drawing move keeps the draw (trap property); we verify
the greedy-optionality defender never loses against a White policy that minimizes Black's safe options.
Run: python3 optionality.py <scratch_dir_with_tables>
"""
import sys, pickle, collections, time
import numpy as np
import chess
import chesslib as cl

scr = sys.argv[1]

def load(name):
    with open(f'{scr}/table_{name}.pkl', 'rb') as f:
        T = pickle.load(f)
    U = cl.Universe(T['pieces'])
    return U, T

def run_winner(name, max_plies=200, only=None):
    U, T = load(name)
    keys, val, dtm, indptr, indices, nmoves = T['keys'], T['val'], T['dtm'], T['indptr'], T['indices'], T['nmoves']
    N = len(keys)
    tw = (keys & 1).astype(bool)
    succ = lambda i: indices[indptr[i]:indptr[i+1]]
    # precompute per position: number of win-preserving moves (White to move, winning): successors with val == -1
    nwin = np.zeros(N, dtype=np.int32)
    for i in np.nonzero(tw & (val == 1))[0]:
        s = succ(i)
        nwin[i] = int((val[s] == -1).sum())
    def black_reply(j):
        # DTM-optimal defence: among Black's moves (all lead to White wins when val[j]==-1) choose max dtm; tie: first
        s = succ(j)
        if len(s) == 0:
            return None
        d = dtm[s]
        return int(s[int(np.argmax(d))])
    def score(j, kind):
        # j: Black-to-move position after White's candidate move (val[j] == -1). Worst case over Black replies.
        s = succ(j)
        if len(s) == 0:
            return 10**6  # mate delivered
        if kind == 'O1':
            return int(nmoves[s].min())
        return int(nwin[s].min())
    results = {}
    for kind, restrict in (('O1', None), ('O2', None), ('O2', 'nonincreasing'), ('O2', 'optimal')):
        if only is not None and restrict != only:
            continue
        t0 = time.time()
        wins = 0; cycles = 0; timeouts = 0; longest = 0
        example_cycle = None
        starts = np.nonzero(tw & (val == 1))[0]
        for i0 in starts:
            i = int(i0); seen = {i}; plies = 0; ok = None
            while True:
                s = succ(i)
                cand = [int(j) for j in s if val[j] == -1]
                if restrict == 'nonincreasing':
                    cand = [j for j in cand if dtm[j] - 1 <= dtm[i]]   # DTM after Black's best reply does not exceed the current DTM
                elif restrict == 'optimal':
                    cand = [j for j in cand if dtm[j] + 1 == dtm[i]]
                best = max(cand, key=lambda j: (score(j, kind), -j))
                plies += 1
                if len(succ(best)) == 0:
                    ok = True; break
                r = black_reply(best); plies += 1
                if r in seen:
                    ok = False; example_cycle = example_cycle or (int(i0), r); break
                seen.add(r); i = r
                if plies > max_plies:
                    ok = None; break
            if ok is True:
                wins += 1; longest = max(longest, plies)
            elif ok is False:
                cycles += 1
            else:
                timeouts += 1
        lab = f'{kind}' + (f' among {restrict} moves' if restrict else ' among win-preserving moves')
        results[lab] = (wins, cycles, timeouts, longest)
        ex = f' first cycling start: {U.fen(int(keys[example_cycle[0]]))}' if example_cycle else ''
        print(f'[{name}] {lab}: starts {len(starts)}, mated {wins}, cycled {cycles}, >{max_plies} plies {timeouts}, longest win {longest} plies ({time.time()-t0:.0f}s){ex}', flush=True)
    return results

def run_defender_kpvk(max_plies=200):
    U, T = load('KPvK')
    keys, val, indptr, indices = T['keys'], T['val'], T['indptr'], T['indices']
    bnd_best = T['bnd_best']
    N = len(keys); tw = (keys & 1).astype(bool)
    succ = lambda i: indices[indptr[i]:indptr[i+1]]
    # Black greedy: among drawing moves (val[j] != 1 for the White-to-move successor), maximize the number of
    # drawing replies available at Black's next turn (min over White replies).
    ndraw = np.zeros(N, dtype=np.int32)
    for i in np.nonzero(~tw & (val == 0))[0]:
        s = succ(i); ndraw[i] = int((val[s] != 1).sum()) + (1 if bnd_best[i] >= 0 else 0)
    losses = 0; draws = 0; starts = np.nonzero(~tw & (val == 0))[0]
    for i0 in starts:
        i = int(i0); plies = 0; seen = set()
        while plies < max_plies:
            s = succ(i)
            cand = [int(j) for j in s if val[j] != 1]
            if not cand:
                break   # stalemate or capture-draw; stays a draw by the table
            def sc(j):
                ss = succ(j)
                return int(ndraw[ss].min()) if len(ss) else 10**6
            j = max(cand, key=lambda j: (sc(j), -j)); plies += 1
            ss = succ(j)
            if len(ss) == 0:
                break
            # White minimizes Black's safe options (adversarial within the drawn region)
            w = min((int(x) for x in ss), key=lambda x: (ndraw[x] if val[x] == 0 else 10**6, x)); plies += 1
            if val[w] == -1:
                losses += 1; break
            if w in seen:
                break
            seen.add(w); i = w
        else:
            pass
        draws += 1
    print(f'[KPvK] greedy-optionality defender from all {len(starts)} drawn Black-to-move positions: losses {losses}')

if __name__ == '__main__':
    only = sys.argv[2] if len(sys.argv) > 2 else None
    for name in ('KRvK', 'KQvK'):
        run_winner(name, only=only)
    if only is None:
        run_defender_kpvk()

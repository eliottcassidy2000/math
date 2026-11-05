// spine_grammar_sublanguages_20261105.cpp -- opus-2026-11-05-S81
// Copy of spine_grammar_census_20261104.cpp (THM-4878 engine) with a letter filter selected by the environment variable
// SPINE_MODE: 1 = no simultaneous double events (b <= 1 and |T| <= 1); 2 = mode 1 and first-in-first-out closing (with
// two pending slots only the older one, slot 0, may terminate: crossing chains); 3 = mode 1 and last-in-first-out closing
// (only the newer slot may terminate: nested chains).  Everything else (census, word mode, dump) is unchanged.
// spine_grammar_census_20261104.cpp -- opus-2026-11-04-S79
// Census of Hamiltonian-path counts over the graded spine grammars G(k, w): a transitive chain v_1 < ... < v_m whose
// backward arcs (j -> i, i < j) have cutwidth <= w (at most w cross any cut between consecutive chain vertices), plus
// k <= 2 apexes x_1 (< x_2, with x_1 -> x_2) placed after the chain, each apex beating an arbitrary subset of the chain
// (apex arcs do not count toward w).  Special cases already in the repository:
//   G(1,0) = the one-apex words (feedback_word_hamilton_spectrum_20261020), G(2,0) = the two-apex words
//   (square_sum_spectrum_overnight_20261030), G(0,1) = Moon's interval chains (THM-4744), G(0,2) contains fa <= 2 (THM-4748).
// Letters.  At each chain vertex: T = the set of pending slots (backward arcs whose head is an earlier chain vertex and
// whose tail is this vertex) that terminate here, b = the number of new slots opened here (their head is this vertex;
// their tails come later), A = the set of apexes beating this vertex.  Slots are labelled in creation order (a physical
// labelling, the same in every state of one word); slots opened at the same vertex form a head group, and a letter may
// terminate at most one slot of each group (two would be a double arc).
// Cut states.  The fragments of the (unknown) Hamiltonian path restricted to the prefix, each with a start type (global
// start / entered later through a pending slot / entered later from an apex) and an end type (global end / open, with the
// set of pending slots and apexes whose head is that end vertex: the future vertices forbidden as its successor).
// The depth-first search over chain letters carries the integer count vector over cut states.  Soundness devices:
//   * pruning by the insertion lemma: the prefix tournament with the apexes placed is an induced subtournament of every
//     completion, so its count is a lower bound; prefixes above X are abandoned;
//   * saturating arithmetic at CAP = 2^62: an entry that contributes to the final count never exceeds that count, so a
//     recorded value (<= X < CAP) is exact and a saturated final sum is correctly rejected;
//   * sink runs: r consecutive slot-free chain vertices beaten by every apex are padding if the word ends there (each is
//     a sink), and otherwise the next non-sink vertex u (u -> x for some apex x; for k = 0 the tail of the spanning
//     pending slot) closes 2^(r-1) odd cycles x -> S -> u -> x (S an odd subset of the run), so H >= 2^r + 1 by the
//     odd-cycle-collection formula (THM-4876).  Runs with 2^r + 1 > X are cut, and no value is recorded after a sink;
//   * leading source letters (no slot event, A = empty, at depth 0) are padding and skipped;
//   * identical count vectors have identical futures and are visited once.
// For k = 0 the chain must be connected (every internal cut crossed by a slot): an empty cut would be an ordered join, and
// products are added by the closure afterwards; for k >= 1 no connectivity is imposed.
// Usage: grammar k w X [maxm] [nodedupe-flag] [witness-limit] [word]
//   maxm = -1 dumps the automaton (states, letters, transitions) for spine_grammar_ranks_20261104.py;
//   word = comma-separated letter indices: evaluates that word only (indices refer to lettersFor(c) at each step).
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <vector>
#include <map>
#include <unordered_map>
#include <algorithm>
#include <string>
using namespace std;
typedef unsigned long long u64;
static const u64 CAP = 1ULL << 62;
static inline u64 satAdd(u64 a, u64 b) { u64 s = a + b; return (s >= CAP || s < a) ? CAP : s; }
static inline u64 satMul(u64 a, u64 m) { if (a == 0 || m == 0) return 0; if (a >= CAP / m) return CAP; return a * m; }
// fragment: s = -1 global start, 0..w-1 pending slot, 100+a apex a; e = GEND or mask: bits 0..7 slots, bits 8+a apexes
static const int GS = -1; static const u64 GEND = ~0ULL;
struct Frag { int s; u64 e; bool operator<(const Frag& o) const { return s != o.s ? s < o.s : e < o.e; } bool operator==(const Frag& o) const { return s == o.s && e == o.e; } };
struct CutState { int c; vector<int> grp; vector<Frag> fr;
    bool operator<(const CutState& o) const { if (c != o.c) return c < o.c; if (grp != o.grp) return grp < o.grp; return fr < o.fr; } };
static int K, W; static u64 X; static int MAXM; static bool DEDUPE = true;
static map<CutState, int> stateId; static vector<CutState> states;
struct Letter { int T; int b; int A; };
static int getId(const CutState& st) { auto it = stateId.find(st); if (it != stateId.end()) return it->second; int id = (int)states.size(); stateId[st] = id; states.push_back(st); return id; }
static bool isComplete(const CutState& st) { return st.fr.size() == 1 && st.fr[0].s == GS && st.fr[0].e == GEND; }
static vector<pair<int, u64>> applyLetter(const CutState& st, const Letter& L) {
    vector<pair<int, u64>> res; int c = st.c; int m = (int)st.fr.size();
    // a letter may terminate at most one slot per head group (otherwise two arcs between the same pair)
    { vector<int> cnt(c + 1, 0); for (int i = 0; i < c; i++) if (L.T >> i & 1) { if (++cnt[st.grp[i]] >= 2) return res; } }
    map<CutState, u64> out;
    u64 tailMask = (u64)L.T;
    u64 headMask = 0; for (int j = 0; j < L.b; j++) headMask |= 1ULL << (c + j);
    u64 apexMask = 0; for (int a = 0; a < K; a++) if (L.A >> a & 1) apexMask |= 1ULL << (8 + a);
    bool hasGS = false, hasGE = false; vector<bool> apexUsed(K, false);
    for (auto& f : st.fr) { if (f.s == GS) hasGS = true; if (f.e == GEND) hasGE = true; if (f.s >= 100) apexUsed[f.s - 100] = true; }
    vector<int> pred; if (!hasGS) pred.push_back(-1);
    for (int p = 0; p < m; p++) if (st.fr[p].e != GEND && (st.fr[p].e & tailMask) == 0) pred.push_back(p);
    for (int j = 0; j < L.b; j++) pred.push_back(1000 + c + j);                 // entered later through new slot j
    for (int a = 0; a < K; a++) if ((L.A >> a & 1) && !apexUsed[a]) pred.push_back(2000 + a);
    vector<int> succ; if (!hasGE) succ.push_back(-2); succ.push_back(-3);
    for (int q = 0; q < m; q++) if (st.fr[q].s >= 0 && st.fr[q].s < 100 && (tailMask >> st.fr[q].s & 1)) succ.push_back(q);
    int newGroup = 0; for (int g : st.grp) newGroup = max(newGroup, g + 1);
    for (int pc : pred) for (int sc : succ) {
        if (pc >= 0 && pc < 1000 && sc == pc) continue;
        Frag nf; nf.s = (pc == -1) ? GS : (pc >= 2000 ? 100 + (pc - 2000) : (pc >= 1000 ? pc - 1000 : st.fr[pc].s));
        nf.e = (sc == -2) ? GEND : (sc == -3 ? (headMask | apexMask) : st.fr[sc].e);
        CutState ns; ns.c = c - __builtin_popcount(L.T) + L.b;
        for (int i = 0; i < m; i++) if (i != pc && i != sc) ns.fr.push_back(st.fr[i]); ns.fr.push_back(nf);
        bool ok = true; for (auto& f : ns.fr) if (f.s >= 0 && f.s < 100 && (tailMask >> f.s & 1)) { ok = false; break; }
        if (!ok) continue;
        // renumber slots in creation order: surviving old slots, then the new ones; groups renumbered by first appearance
        vector<int> newIdx(c + L.b, -1); int t = 0; for (int i = 0; i < c; i++) if (!(tailMask >> i & 1)) newIdx[i] = t++; for (int j = 0; j < L.b; j++) newIdx[c + j] = t++;
        vector<int> oldGrp(st.grp); for (int j = 0; j < L.b; j++) oldGrp.push_back(newGroup);
        ns.grp.assign(ns.c, -1); for (int i = 0; i < c + L.b; i++) if (newIdx[i] >= 0) ns.grp[newIdx[i]] = oldGrp[i];
        { map<int,int> rel; for (int& g : ns.grp) { auto it = rel.find(g); if (it == rel.end()) { int nid = (int)rel.size(); rel[g] = nid; g = nid; } else g = it->second; } }
        for (auto& f : ns.fr) { if (f.s >= 0 && f.s < 100) f.s = newIdx[f.s];
            if (f.e != GEND) { u64 mm = f.e & ~((1ULL << 8) - 1); for (int i = 0; i < c + L.b; i++) if ((f.e >> i & 1) && newIdx[i] >= 0) mm |= 1ULL << newIdx[i]; f.e = mm; } }
        sort(ns.fr.begin(), ns.fr.end());
        out[ns] += 1; }
    for (auto& kv : out) res.push_back({getId(kv.first), kv.second}); return res; }
static vector<pair<int, u64>> applyApex(const CutState& st, int a) {
    map<CutState, u64> out; int m = (int)st.fr.size(); bool hasGS = false, hasGE = false;
    for (auto& f : st.fr) { if (f.s == GS) hasGS = true; if (f.e == GEND) hasGE = true; }
    vector<int> pred; if (!hasGS) pred.push_back(-1); for (int p = 0; p < m; p++) if (st.fr[p].e != GEND && !(st.fr[p].e >> (8 + a) & 1)) pred.push_back(p);
    vector<int> succ; if (!hasGE) succ.push_back(-2); succ.push_back(-3); for (int q = 0; q < m; q++) if (st.fr[q].s == 100 + a) succ.push_back(q);
    for (int pc : pred) for (int sc : succ) { if (pc >= 0 && sc == pc) continue;
        Frag nf; nf.s = (pc == -1) ? GS : st.fr[pc].s; nf.e = (sc == -2) ? GEND : (sc == -3 ? 0 : st.fr[sc].e);
        CutState ns; ns.c = st.c; ns.grp = st.grp; for (int i = 0; i < m; i++) if (i != pc && i != sc) ns.fr.push_back(st.fr[i]); ns.fr.push_back(nf);
        bool ok = true; for (auto& f : ns.fr) if (f.s == 100 + a) { ok = false; break; } if (!ok) continue;
        for (auto& f : ns.fr) if (f.e != GEND) f.e &= ~(1ULL << (8 + a));
        sort(ns.fr.begin(), ns.fr.end()); out[ns] += 1; }
    vector<pair<int, u64>> res; for (auto& kv : out) res.push_back({getId(kv.first), kv.second}); return res; }
static int MODE = 0;
static bool letterAllowed(int c, int T, int b) {
    if (MODE == 0) return true;
    if (b > 1 || __builtin_popcount(T) > 1) return false;
    if (MODE == 2 && c == 2 && T == 2) return false;      // only the older slot (bit 0) may close
    if (MODE == 3 && c == 2 && T == 1) return false;      // only the newer slot (bit 1) may close
    return true; }
static vector<Letter> lettersFor(int c) { vector<Letter> L; for (int T = 0; T < (1 << c); T++) { int t = __builtin_popcount(T); for (int b = 0; c - t + b <= W; b++) for (int A = 0; A < (1 << K); A++) if (letterAllowed(c, T, b)) L.push_back({T, b, A}); } return L; }
static vector<vector<vector<pair<int,u64>>>> transCache; static vector<vector<vector<pair<int,u64>>>> apexCache; static vector<char> built;
static void ensureState(int id) {
    if ((int)transCache.size() <= id) { transCache.resize(id + 1); apexCache.resize(id + 1); built.resize(id + 1, 0); }
    if (built[id]) return; built[id] = 1;
    CutState st = states[id]; vector<Letter> L = lettersFor(st.c); transCache[id].resize(L.size());
    for (size_t i = 0; i < L.size(); i++) transCache[id][i] = applyLetter(st, L[i]);
    apexCache[id].resize(K); for (int a = 0; a < K; a++) apexCache[id][a] = applyApex(st, a); }
typedef vector<u64> Vec;
static u64 nodes = 0, recorded = 0, deduped = 0, runCut = 0, closerCut = 0; static vector<unsigned char> att; static vector<string> witness;
static vector<int> bestFa, bestM; static int curFa = 0;   // witness with the fewest reversed arcs (then shortest chain) per value
static unordered_map<u64, u64> seen; static vector<int> curWord;
static u64 hashVec(const Vec& v, u64 seed) { u64 h = seed; for (size_t i = 0; i < v.size(); i++) if (v[i]) { h ^= (v[i] + 0x9E3779B97F4A7C15ULL * (u64)(i + 1)); h *= 0xBF58476D1CE4E5B9ULL; h ^= h >> 31; } return h; }
static string describeWord() { string s; for (int x : curWord) s += to_string(x) + ","; return s; }
static Vec applyApexes(const Vec& v) {
    Vec cur = v; for (int a = 0; a < K; a++) { Vec nxt(states.size(), 0);
        for (size_t s = 0; s < cur.size(); s++) if (cur[s]) { ensureState((int)s); for (auto& pr : apexCache[s][a]) { if (pr.first >= (int)nxt.size()) nxt.resize(states.size(), 0); nxt[pr.first] = satAdd(nxt[pr.first], satMul(cur[s], pr.second)); } }
        cur.swap(nxt); }
    return cur; }
static u64 finishValue(const Vec& v) { Vec cur = applyApexes(v); u64 h = 0; for (size_t s = 0; s < cur.size(); s++) if (cur[s] && isComplete(states[s])) h = satAdd(h, cur[s]); return h; }
static void record(const Vec& v) { u64 h = finishValue(v); recorded++; if (h >= 1 && h <= X) { int m = (int)curWord.size();
    if (!att[h] || curFa < bestFa[h] || (curFa == bestFa[h] && m < bestM[h])) { att[h] = 1; bestFa[h] = curFa; bestM[h] = m; witness[h] = describeWord(); } } }
static Vec applyLetterVec(const Vec& v, int li) {
    Vec nxt(states.size(), 0);
    for (size_t s = 0; s < v.size(); s++) if (v[s]) { ensureState((int)s); for (auto& pr : transCache[s][li]) { if (pr.first >= (int)nxt.size()) nxt.resize(states.size(), 0); nxt[pr.first] = satAdd(nxt[pr.first], satMul(v[s], pr.second)); } }
    return nxt; }
static void dfs(int depth, const Vec& v, int sinkRun) {
    nodes++;
    int c = -1; for (size_t s = 0; s < v.size(); s++) if (v[s]) { c = states[s].c; break; } if (c < 0) return;
    if (c == 0 && depth > 0 && sinkRun == 0) record(v);
    if (depth >= MAXM) return;
    if (finishValue(v) > X) return;
    if (c > 0) {   // closer bound: the completion induces, on the prefix plus the first future tail, the prefix closed by one
                   // vertex (terminating a nonempty admissible slot set, any apex bits); the least such count bounds H below
        u64 lb = CAP; vector<Letter> Lc = lettersFor(c);
        for (size_t li = 0; li < Lc.size(); li++) { if (Lc[li].b != 0 || Lc[li].T == 0) continue;
            Vec nx = applyLetterVec(v, (int)li); bool any = false; for (u64 x : nx) if (x) { any = true; break; } if (!any) continue;
            u64 h = finishValue(nx); if (h < lb) lb = h; }
        if (lb > X) { closerCut++; return; } }
    if (DEDUPE) { u64 ha = hashVec(v, 1), hb = hashVec(v, 2); auto it = seen.find(ha); if (it != seen.end() && it->second == hb) { deduped++; return; } seen[ha] = hb; }
    vector<Letter> L = lettersFor(c); int allA = (1 << K) - 1;
    for (size_t li = 0; li < L.size(); li++) {
        const Letter& Lt = L[li]; int cNext = c - __builtin_popcount(Lt.T) + Lt.b;
        bool slotFree = (Lt.T == 0 && Lt.b == 0);
        bool isSink = slotFree && Lt.A == allA;
        bool isSource = slotFree && Lt.A == 0;
        if (depth == 0 && (K == 0 ? Lt.b == 0 : isSource)) continue;
        int nextRun = isSink ? sinkRun + 1 : 0;
        if (isSink && (nextRun >= 62 || (1ULL << nextRun) + 1 > X)) { runCut++; continue; }
        Vec nxt = applyLetterVec(v, (int)li);
        bool any = false; for (u64 x : nxt) if (x) { any = true; break; } if (!any) continue;
        curWord.push_back((int)li); curFa += __builtin_popcount(Lt.T);
        if (K == 0 && cNext == 0) record(nxt);
        else dfs(depth + 1, nxt, nextRun);
        curWord.pop_back(); curFa -= __builtin_popcount(Lt.T); } }
int main(int argc, char** argv) {
    if (argc < 4) { fprintf(stderr, "usage: grammar k w X [maxm] [nodedupe-flag] [witness-limit] [word]\n"); return 1; }
    if (getenv("SPINE_MODE")) MODE = atoi(getenv("SPINE_MODE"));
    K = atoi(argv[1]); W = atoi(argv[2]); X = strtoull(argv[3], 0, 10); MAXM = argc > 4 ? atoi(argv[4]) : 64; if (argc > 5 && atoi(argv[5])) DEDUPE = false;
    att.assign(X + 1, 0); witness.assign(X + 1, ""); att[1] = 1; witness[1] = "transitive"; bestFa.assign(X + 1, 1 << 30); bestM.assign(X + 1, 1 << 30); bestFa[1] = 0; bestM[1] = 0;
    CutState empty; empty.c = 0; int e0 = getId(empty); ensureState(e0);
    Vec v0(1, 0); v0[e0] = 1;
    if (argc > 7) {
        Vec cur = v0; string wd = argv[7]; size_t pos = 0; int depth = 0;
        while (pos < wd.size()) { size_t q = wd.find(',', pos); string tok = wd.substr(pos, q == string::npos ? string::npos : q - pos); pos = (q == string::npos) ? wd.size() : q + 1; if (tok.empty()) continue; int li = atoi(tok.c_str());
            int c = -1; for (size_t s = 0; s < cur.size(); s++) if (cur[s]) { c = states[s].c; break; }
            if (c < 0) { printf("dead word\n"); return 0; }
            vector<Letter> L = lettersFor(c); if (li < 0 || li >= (int)L.size()) { printf("letter %d out of range for c = %d\n", li, c); return 1; }
            cur = applyLetterVec(cur, li); depth++;
            printf("after letter %d = (T=%d,b=%d,A=%d) (depth %d): prefix H with apexes = %llu; states:", li, L[li].T, L[li].b, L[li].A, depth, finishValue(cur));
            for (size_t s = 0; s < cur.size(); s++) if (cur[s]) { printf(" [c=%d", states[s].c); for (auto& f : states[s].fr) { printf(" (%d,", f.s); if (f.e == GEND) printf("E)"); else printf("%llx)", f.e); } printf("]x%llu", cur[s]); } printf("\n"); }
        printf("H = %llu\n", finishValue(cur)); return 0; }
    if (MAXM < 0) {
        vector<int> queue; queue.push_back(e0); vector<char> seenS(1, 1);
        for (size_t qi = 0; qi < queue.size(); qi++) { int id = queue[qi]; ensureState(id);
            for (auto& lst : transCache[id]) for (auto& pr : lst) { if (pr.first >= (int)seenS.size()) seenS.resize(pr.first + 1, 0); if (!seenS[pr.first]) { seenS[pr.first] = 1; queue.push_back(pr.first); } }
            for (int a = 0; a < K; a++) for (auto& pr : apexCache[id][a]) { if (pr.first >= (int)seenS.size()) seenS.resize(pr.first + 1, 0); if (!seenS[pr.first]) { seenS[pr.first] = 1; queue.push_back(pr.first); } } }
        printf("STATES %zu\n", queue.size());
        for (int id : queue) { const CutState& st = states[id]; printf("S %d c=%d complete=%d fr=", id, st.c, (int)isComplete(st)); for (auto& f : st.fr) { printf("(%d,", f.s); if (f.e == GEND) printf("E)"); else printf("%llu)", f.e); } printf(" grp="); for (int g : st.grp) printf("%d", g); printf("\n");
            vector<Letter> L = lettersFor(st.c); for (size_t li = 0; li < L.size(); li++) { printf("T %d %zu %d %d %d :", id, li, L[li].T, L[li].b, L[li].A); for (auto& pr : transCache[id][li]) printf(" %d*%llu", pr.first, pr.second); printf("\n"); }
            for (int a = 0; a < K; a++) { printf("A %d %d :", id, a); for (auto& pr : apexCache[id][a]) printf(" %d*%llu", pr.first, pr.second); printf("\n"); } }
        return 0; }
    dfs(0, v0, 0);
    u64 cnt = 0; for (u64 x = 1; x <= X; x += 2) if (att[x]) cnt++;
    printf("G(k=%d, w=%d), X = %llu, maxm = %d, dedupe = %d: %zu automaton states, %llu DFS nodes, %llu dedup hits, %llu sink-run cuts, %llu closer cuts, %llu words completed; odd values attained directly: %llu of %llu\n", K, W, X, MAXM, (int)DEDUPE, states.size(), nodes, deduped, runCut, closerCut, recorded, cnt, (X + 1) / 2);
    vector<unsigned char> cl(att.begin(), att.end()); for (u64 a = 3; a <= X; a += 2) if (cl[a]) for (u64 b = 3; a * b <= X; b += 2) if (att[b]) cl[a * b] = 1;
    u64 cntcl = 0; vector<u64> miss; for (u64 x = 1; x <= X; x += 2) { if (cl[x]) cntcl++; else miss.push_back(x); }
    printf("with ordered-join products: %llu of %llu; %zu misses; first misses:", cntcl, (X + 1) / 2, miss.size()); for (size_t i = 0; i < miss.size() && i < 80; i++) printf(" %llu", miss[i]); printf("%s\n", miss.size() > 80 ? " ..." : "");
    printf("misses per dyadic scale (with products):"); for (int k = 2; (1ULL << k) <= X; k++) { u64 lo = 1ULL << k, hi = min(X, (1ULL << (k + 1)) - 1), cc = 0, tot = 0; for (u64 x = lo | 1; x <= hi; x += 2) { tot++; if (!cl[x]) cc++; } printf(" [2^%d:%llu/%llu]", k, cc, tot); } printf("\n");
    printf("ATTAINED_DIRECT:"); for (u64 x = 1; x <= X; x += 2) if (att[x]) printf(" %llu", x); printf("\n");
    { map<int, u64> byFa; for (u64 x = 1; x <= X; x += 2) if (att[x]) byFa[bestFa[x]]++; printf("values by minimal number of reversed arcs (over the words explored; exact minima only without dedupe):"); for (auto& kv : byFa) printf(" fa%d:%llu", kv.first, kv.second); printf("\n"); }
    if (argc > 6) { u64 lim = strtoull(argv[6], 0, 10); if (lim) { printf("direct witnesses (letter index per chain vertex; letters enumerated as (T,b,A) for the current slot count) for odd values <= %llu:\n", lim); for (u64 x = 1; x <= lim && x <= X; x += 2) if (att[x]) printf("  %llu = %s | fa=%d m=%d\n", x, witness[x].c_str(), bestFa[x], bestM[x]); } }
    return 0; }

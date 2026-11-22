// Count every one-arc neighbor of a backchain by its partial directed path cover.
// c++ -std=c++17 -O3 backchain_flip_sweep.cpp -lgmpxx -lgmp -o /tmp/backchain_flip_sweep
#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <optional>
#include <set>
#include <string>
#include <type_traits>
#include <vector>
#include <gmpxx.h>

using std::vector;
using Big = mpz_class;
static_assert(sizeof(unsigned long)==sizeof(uint64_t),"GMP low-word extraction requires 64-bit unsigned long");
constexpr int absent = 3, empty_state = 63, final_state = 60;

// A state has one two-bit field per possible path head: permanently forgotten,
// first active vertex, second active vertex. Its field is the tail or absent.
// Each forgotten head needs to start the final Hamilton path, so at most one
// component can have such a head. A zero tail beats every future vertex.
struct Transition {
    std::array<vector<int>, 64> destinations;
    Transition(int before, int back, int keep, bool add) {
        int rename[4] = {}, next_active = 0;
        for (int i = 1; i <= before; ++i)
            if (keep & (1 << (i-1))) rename[i] = ++next_active;
        if (add) rename[3] = ++next_active;
        assert(next_active <= 2);
        for (int state = 0; state < 64; ++state) {
            int heads[3], tails[3], count = 0;
            bool valid = true;
            for (int h = 0; h < 3; ++h) {
                int t = (state >> (2*h)) & 3;
                if (t == absent) continue;
                if (h > before || t > before) valid = false;
                heads[count] = h;
                tails[count++] = t;
            }
            if (!valid) continue;
            auto beats_new = [&](int t) { return !t || !(back & (1 << (t-1))); };
            for (int pred = -1; pred < count; ++pred) {
                if (pred >= 0 && !beats_new(tails[pred])) continue;
                for (int succ = -1; succ < count; ++succ) {
                    if (succ >= 0 && (succ == pred || beats_new(heads[succ]))) continue;
                    int dest = empty_state;
                    auto place = [&](int h, int t) {
                        const int shift = 2 * rename[h];
                        if (((dest >> shift) & 3) != absent) return false;
                        dest = (dest & ~(3 << shift)) | (rename[t] << shift);
                        return true;
                    };
                    bool good = place(pred < 0 ? 3 : heads[pred], succ < 0 ? 3 : tails[succ]);
                    for (int i = 0; good && i < count; ++i)
                        if (i != pred && i != succ) good = place(heads[i], tails[i]);
                    if (good) destinations[state].push_back(dest);
                }
            }
        }
    }
};

const Transition& transition(int before, int back, int keep, bool add) {
    static std::array<std::optional<Transition>, 128> cache;
    int key = (before << 5) | (back << 3) | (keep << 1) | add;
    if (!cache[key]) cache[key].emplace(before, back, keep, add);
    return *cache[key];
}

struct Host {
    vector<int> next;
    explicit Host(const vector<int>& widths) {
        int lo = 0;
        for (int w : widths) {
            next.resize(lo+w+2, -1);
            next[lo] = lo+w+1;
            lo += w+1;
        }
    }
    Big count(int lo = -1, int hi = -1) const {
        std::array<Big, 64> first{}, second{};
        auto *dp = &first, *out = &second;
        (*dp)[empty_state] = 1;
        int bag[2] = {}, size = 0;
        int last = lo < 0 ? -1 : (next[lo] == hi ? -1 : std::max(hi, next[lo]));
        auto end = [&](int u) { return u == lo ? last : next[u]; };
        for (int v = 0; v < int(next.size()); ++v) {
            int keep = 0, back = 0, future[2] = {}, new_size = 0;
            for (int i = 0; i < size; ++i) {
                const int u = bag[i];
                if ((next[u] == v) != (u == lo && v == hi)) back |= 1 << i;
                if (end(u) > v) {
                    future[new_size++] = u;
                    keep |= 1 << i;
                }
            }
            bool add = end(v) > v;
            if (add) {
                assert(new_size < 2);
                future[new_size++] = v;
            }
            const auto& op = transition(size, back, keep, add);
            for (int s = 0; s < 64; ++s) {
                if ((*dp)[s] == 0) continue;
                for (int t : op.destinations[s]) (*out)[t] += (*dp)[s];
                (*dp)[s] = 0;
            }
            std::swap(dp,out);
            std::copy(future, future+new_size, bag);
            size = new_size;
        }
        return (*dp)[final_state];
    }
    unsigned long direct(int lo, int hi) const {
        const int n = int(next.size());
        vector<vector<unsigned long>> dp(1<<n, vector<unsigned long>(n));
        for (int i = 0; i < n; ++i) dp[1<<i][i] = 1;
        for (int mask = 1; mask < (1<<n); ++mask) {
            for (int u = 0; u < n; ++u) {
                if (!dp[mask][u]) continue;
                for (int v = 0; v < n; ++v) if (!(mask & (1<<v))) {
                    int a = std::min(u,v), b = std::max(u,v);
                    bool forward = (next[a] != b) != (a == lo && b == hi);
                    if ((u < v) == forward) dp[mask | (1<<v)][v] += dp[mask][u];
                }
            }
        }
        unsigned long ans = 0;
        for (auto value : dp.back()) ans += value;
        return ans;
    }
    bool strong(int lo=-1,int hi=-1) const {
        const int n=int(next.size());
        for(bool reverse:{false,true}) {
            vector<bool> reached(n);
            vector<int> queue={0};
            reached[0]=true;
            for(std::size_t k=0;k<queue.size();++k) {
                int u=queue[k];
                for(int v=0;v<n;++v) if(!reached[v]) {
                    int a=std::min(u,v),b=std::max(u,v);
                    bool orientation=(next[a]!=b) != (a==lo && b==hi);
                    if(((u<v) == orientation) != reverse) {
                        reached[v]=true;
                        queue.push_back(v);
                    }
                }
            }
            if(int(queue.size())!=n) return false;
        }
        return true;
    }
};

void controls() {
    int tests = 0, unit_increases = 0;
    for (const vector<int>& word : vector<vector<int>>{{0},{1},{1,0,1},{0,0,0,0},{2,0,1,2},{1,0,0,0,1,1}}) {
        Host h(word);
        Big base = h.count();
        assert(base == h.direct(-1,-1));
        for (int i=0; i<int(h.next.size()); ++i)
            for (int j=i+1; j<int(h.next.size()); ++j) {
                Big val = h.count(i,j);
                assert(val == h.direct(i,j));
                ++tests;
                if (val == base+2) ++unit_increases;
            }
    }
    assert(unit_increases > 0);
    std::cout << "Subset-DP controls: " << tests << " edited graphs, "
              << unit_increases << " unit increases\n";
}

void sweep(int s, int t) {
    vector<int> word = {2,2,1,0,0,1,1,0,0,0,0,1,0,0,0,2,0,0,0,1,0,1,2,0,0,0,4,0,0,1,0,0,0,1,s,t};
    Host h(word);
    const Big target = (Big(1)<<60) * ((Big(1)<<61)-1) + 1;
    const Big lower = h.count();
    assert(h.next.size() == 133 && lower == target-2);
    std::array<int,3> hist{};
    int mod64 = 0, exact = 0, down = 0, equal = 0, up = 0;
    std::optional<Big> best;
    vector<std::pair<int,int>> nearest;
    for (int i=0; i<int(h.next.size()); ++i) for (int j=i+1; j<int(h.next.size()); ++j) {
        Big count = h.count(i,j);
        ++hist[mpz_fdiv_ui(count.get_mpz_t(),3)];
        if (mpz_congruent_2exp_p(count.get_mpz_t(),target.get_mpz_t(),6)) ++mod64;
        if (count == target) ++exact;
        if (count < lower) ++down;
        if (count == lower) ++equal;
        if (count <= lower) continue;
        ++up;
        Big jump = count-lower;
        if (!best || jump < *best) { best = jump; nearest.clear(); }
        if (jump == *best) nearest.emplace_back(i,j);
    }
    std::cout << "suffix " << s << ',' << t << " residues mod 3 " << hist[0] << ' '
              << hist[1] << ' ' << hist[2] << " pass upper mod 64 " << mod64
              << " exact " << exact << " down/equal/up " << down << '/' << equal << '/' << up;
    if (best) {
        std::cout << " nearest upward jump " << *best << " multiplicity " << nearest.size() << " arcs";
        for (auto [i,j] : nearest) std::cout << " (" << i << ',' << j << ')';
    }
    std::cout << std::endl;
}

struct Bag {
    int size=0;
    std::array<int,2> vertices{};
    bool contains(int u) const {
        return std::find(vertices.begin(), vertices.begin()+size, u) != vertices.begin()+size;
    }
    Bag with(int u) const {
        Bag out=*this;
        if (!out.contains(u)) {
            assert(out.size < 2);
            out.vertices[out.size++]=u;
            if(out.size==2 && out.vertices[0]>out.vertices[1]) std::swap(out.vertices[0],out.vertices[1]);
        }
        return out;
    }
};

using Counts = std::array<uint64_t,64>;
constexpr uint64_t odd_ogg=809482495054428195ULL;
constexpr std::array<int,15> ogg_primes={2,3,5,7,11,13,17,19,23,29,31,41,47,59,71};

template<bool odd> uint64_t add(uint64_t a,uint64_t b) {
    uint64_t c=a+b;
    if constexpr(odd) return c>=odd_ogg ? c-odd_ogg : c;
    return c;
}

const Transition& column(const Host& host, const Bag& left, const Bag& right, int v, int toggle=-1) {
    int back=0,keep=0;
    for(int i=0;i<left.size;++i) {
        int u=left.vertices[i];
        if((host.next[u]==v) != (u==toggle)) back|=1<<i;
        if(right.contains(u)) keep|=1<<i;
    }
    return transition(left.size,back,keep,right.contains(v));
}

template<bool odd> Counts forward(const Counts& counts, const Transition& op) {
    Counts answer{};
    for(int s=0;s<64;++s) if(counts[s])
        for(int t:op.destinations[s]) answer[t]=add<odd>(answer[t],counts[s]);
    return answer;
}

template<bool odd> Counts backward(const Counts& counts, const Transition& op) {
    Counts answer{};
    for(int s=0;s<64;++s)
        for(int t:op.destinations[s]) answer[s]=add<odd>(answer[s],counts[t]);
    return answer;
}

void print_mask(int mask) {
    for(int p:ogg_primes) {
        if(mask&1) std::cout<<' '<<p;
        mask>>=1;
    }
}

struct Example {
    int mask=0, i=0, j=0;
    vector<int> widths;
    Big value;
    bool strong_before=false,strong_after=false;
    void save(const Host& h,int left,int right,int new_mask) {
        mask=new_mask; i=left; j=right;
        value=h.count(i,j);
        const Big discrepancy=value-((Big(1)<<60)*((Big(1)<<61)-1)+1);
        for(int k=0;k<15;++k)
            assert(bool(mask&(1<<k))==bool(mpz_divisible_ui_p(discrepancy.get_mpz_t(),ogg_primes[k])));
        strong_before=h.strong();strong_after=h.strong(i,j);
        widths.clear();
        for(int v=0;v<int(h.next.size());++v)
            if(h.next[v]>=0) widths.push_back(h.next[v]-v-1);
    }
    void show(const char* label) const {
        std::cout<<label<<" primes";print_mask(mask);
        std::cout<<", arc "<<i<<','<<j<<", upper "<<value<<", strong "
                 <<strong_before<<'/'<<strong_after<<", widths";
        for(int k:widths) std::cout<<' '<<k;
        std::cout<<'\n';
    }
};

struct Phases {
    std::array<long long,32768> histogram{};
    std::array<uint64_t,32768> product{};
    Example by_count,by_product,by_count_three,by_product_three;
    Phases() {
        product[0]=1;
        for(unsigned i=1;i<product.size();++i)
            product[i]=product[i&(i-1)]*ogg_primes[__builtin_ctz(i)];
        assert(product.back()==2*odd_ogg);
    }
    void record(const Host& h,int i,int j,uint64_t remainder) {
        int mask=1;
        // Compilers can replace each of these fixed prime divisions with a multiply.
        for(int k=1;k<15;++k)
            if(remainder%ogg_primes[k]==0) mask|=1<<k;
        ++histogram[mask];
        auto better_count=[&](const Example& x){
            int diff=__builtin_popcount(unsigned(mask))-__builtin_popcount(unsigned(x.mask));
            return diff>0 || (diff==0 && product[mask]>product[x.mask]);
        };
        if(better_count(by_count)) by_count.save(h,i,j,mask);
        if(product[mask]>product[by_product.mask]) by_product.save(h,i,j,mask);
        if((mask&2) && better_count(by_count_three)) by_count_three.save(h,i,j,mask);
        if((mask&2) && product[mask]>product[by_product_three.mask]) by_product_three.save(h,i,j,mask);
    }
    void show() const {
        auto supersets=histogram;
        for(int k=0;k<15;++k) for(int m=0;m<32768;++m)
            if(!(m&(1<<k))) supersets[m]+=supersets[m|(1<<k)];
        int size=16, ways=0, cheapest=0;
        for(int m=1;m<32768;++m) if(!supersets[m]) {
            const int bits=__builtin_popcount(unsigned(m));
            if(bits<size) {size=bits;ways=0;cheapest=m;}
            if(bits==size) {
                ++ways;
                if(product[m]<product[cheapest]) cheapest=m;
            }
        }
        std::cout<<"all 15 pass "<<supersets.back()<<", total "<<supersets[0]
                 <<", min excluding cardinality "<<size<<", subsets "<<ways
                 <<", least product "<<product[cheapest]<<", primes";
        print_mask(cheapest);std::cout<<'\n';
        for(int k=1;k<15;++k) if(cheapest&(1<<k))
            std::cout<<"survivors deleting "<<ogg_primes[k]<<": "<<supersets[cheapest^(1<<k)]<<'\n';
        by_count.show("max-cardinality");
        by_product.show("max-product");
        by_count_three.show("max-cardinality with three");
        by_product_three.show("max-product with three");
    }
};

struct Screen {
    long long trials=0, mod64=0, mod2pow64=0, mod65521=0, exact=0;
};

template<bool odd> void screens(const Host& h, int end_left, int begin_right, const Big& target,
             Screen& stats, bool compare_all=false, Phases* phases=nullptr) {
    const int n=int(h.next.size());
    vector<Bag> ordinary(n+1);
    vector<Counts> prefix(n+1),suffix(n+1);
    prefix[0][empty_state]=1;
    for(int i=0;i<n;++i) {
        Bag next_bag;
        if(ordinary[i].size && h.next[ordinary[i].vertices[0]]>i) next_bag=ordinary[i];
        if(h.next[i]>i) next_bag=next_bag.with(i);
        ordinary[i+1]=next_bag;
        prefix[i+1]=forward<odd>(prefix[i],column(h,ordinary[i],ordinary[i+1],i));
    }
    suffix[n][final_state]=1;
    for(int i=n-1;i>=0;--i)
        suffix[i]=backward<odd>(suffix[i+1],column(h,ordinary[i],ordinary[i+1],i));
    assert(prefix[n][final_state]==suffix[0][empty_state]);
    auto residue=[](const Big& x){
        if constexpr(odd) return mpz_fdiv_ui(x.get_mpz_t(),odd_ogg);
        return mpz_get_ui(x.get_mpz_t());
    };
    const uint64_t wanted=residue(target);
    if(!compare_all) assert(prefix[n][final_state]==residue(target-2));
    for(int i=0;i<end_left;++i) {
        Bag live=ordinary[i+1].with(i);
        Counts dp=forward<odd>(prefix[i], column(h,ordinary[i],live,i));
        for(int j=i+1;j<n;++j) {
            if(j>=begin_right) {
                // Finish the edit now; all columns strictly after j are ordinary.
                const auto& last=column(h,live,ordinary[j+1],j,i);
                using Accumulator=std::conditional_t<odd,unsigned __int128,uint64_t>;
                Accumulator accum=0;
                for(int state=0;state<64;++state) if(dp[state]) {
                    uint64_t total=0;
                    for(int t:last.destinations[state]) total=add<odd>(total,suffix[j+1][t]);
                    // At most 64 products, each below odd_ogg squared, fit in 128 bits.
                    accum+=Accumulator(dp[state])*total;
                }
                uint64_t val;
                if constexpr(odd) val=accum%odd_ogg;
                else val=accum;
                if(compare_all) assert(val==residue(h.count(i,j)));
                ++stats.trials;
                if constexpr(odd) {
                    if(phases) phases->record(h,i,j,val>=wanted ? val-wanted : val+odd_ogg-wanted);
                } else if((val&63)==(wanted&63)) ++stats.mod64;
                if(!odd && val==wanted) {
                    ++stats.mod2pow64;
                    const Big value=h.count(i,j);
                    if(mpz_congruent_ui_p(value.get_mpz_t(),mpz_fdiv_ui(target.get_mpz_t(),65521),65521))
                        ++stats.mod65521;
                    if(value==target) {
                        ++stats.exact;
                        std::cout<<"EXACT arc "<<i<<','<<j<<" boundaries";
                        for(int v=0;v<n;++v) if(h.next[v]>=0) std::cout<<' '<<v;
                        std::cout<<std::endl;
                    }
                }
            }
            if(j+1<n) {
                Bag future=ordinary[j+1].with(i);
                dp=forward<odd>(dp,column(h,live,future,j));
                live=future;
            }
        }
    }
}

template<bool odd> void factor_bank(bool all_arcs) {
    const vector<vector<int>> factors={{2},{2},{1,0,0,1},
        {1,0,0,0,1,0,0,4,0,0,0,2},{1,0,0,0,0,1,0,0,0,2,0,0,0,1,0,1},{15},{61}};
    const Big target=(Big(1)<<60)*((Big(1)<<61)-1)+1;
    Screen checks;
    for(const vector<int>& w:vector<vector<int>>{{1,0,0,1},{2,0,1,2},{1,2,0,0,1,1}}) {
        Host h(w);
        screens<odd>(h,int(h.next.size()),0,Big(0),checks,true);
    }
    Screen result;
    Phases phases;
    std::set<vector<int>> distinct;
    long long representatives=0;
    for(int flag=0;flag<4;++flag) {
        auto words=factors;
        if(flag&1) std::reverse(words[3].begin(),words[3].end());
        if(flag&2) std::reverse(words[4].begin(),words[4].end());
        std::array<int,7> order={0,0,2,3,4,5,6};
        do {
            vector<int> w;
            for(int idx:order) w.insert(w.end(),words[idx].begin(),words[idx].end());
            assert(distinct.insert(w).second);
            if(std::lexicographical_compare(w.rbegin(),w.rend(),w.begin(),w.end())) continue;
            ++representatives;
            auto increase=[](const vector<int>& block) {
                int distance=0;
                for(int width:block) distance+=width+1;
                return distance;
            };
            int end_left=increase(words[order.front()])+1+words[order[1]].front();
            int end_from_right=increase(words[order.back()])+1+words[order[5]].back();
            Host h(w);
            screens<odd>(h,all_arcs ? int(h.next.size()) : end_left,
                    all_arcs ? 0 : int(h.next.size())-end_from_right,target,result,false,&phases);
        } while(std::next_permutation(order.begin(),order.end()));
    }
    assert(distinct.size()==10080 && representatives==5040);
    std::cout<<"factor words "<<distinct.size()<<", converse representatives "<<representatives
             <<", candidate arcs among representatives "<<result.trials;
    if constexpr(odd) {std::cout<<'\n';phases.show();}
    else std::cout<<", pass mod 64 "<<result.mod64
             <<", pass mod 2^64 "<<result.mod2pow64<<", also pass mod 65521 "<<result.mod65521
             <<", exact "<<result.exact<<std::endl;
}

int main(int argc,char** argv) {
    controls();
    if(argc>1 && std::string(argv[1])=="--factor-bank") factor_bank<false>(false);
    else if(argc>1 && std::string(argv[1])=="--factor-bank-all-arcs") factor_bank<false>(true);
    else if(argc>1 && std::string(argv[1])=="--ogg") factor_bank<true>(true);
    else { sweep(15,61); sweep(61,15); }
}

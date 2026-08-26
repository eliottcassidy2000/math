// Primary exact depth-nine closure for THM-4207.
//
// This translation unit reuses the literal joint-wall geometry and exact
// depth-eight deck from the maintained primary THM-4207 implementation.  It
// then enumerates every nine-subset of P.  The exceptional nine-subsets that
// hit every E_8(50,51) edge are closed by complement-local zeta transforms
// which produce explicit lawful depth-nine separator edges.

#define main thm4207_depth_eight_original_main
#include "lrc14_two_newcomer_depth_eight_obstruction_thm4207.cpp"
#undef main

namespace {

struct Separator {
    std::uint32_t repair = 0;
    i128 delta = 0;
    std::uint64_t score = 0;
};

Separator best_separator(
    std::uint32_t target,
    const std::vector<std::uint32_t>& bodies,
    const std::vector<unsigned char>& unresolved,
    const PairGeometry& geometry
) {
    std::array<int, 30> compressed{};
    compressed.fill(-1);
    std::array<int, 21> expanded_vertex{};
    int dimension = 0;
    for (int vertex = 0; vertex < 30; ++vertex) {
        if ((target & (std::uint32_t{1} << vertex)) == 0) {
            require(dimension < 21, "target is not a nine-set");
            compressed[vertex] = dimension;
            expanded_vertex[dimension++] = vertex;
        }
    }
    require(dimension == 21, "target complement dimension");

    std::vector<i64> mass(std::size_t{1} << dimension, 0);
    for (const PairAtom& atom : geometry.atoms) {
        if ((atom.failed_pool & target) != 0) continue;
        std::uint32_t key = 0;
        for (int vertex = 0; vertex < 30; ++vertex) {
            if ((atom.failed_pool & (std::uint32_t{1} << vertex)) != 0) {
                key |= std::uint32_t{1} << compressed[vertex];
            }
        }
        mass[key] += atom.length;
    }
    for (int bit = 0; bit < dimension; ++bit) {
        const std::uint32_t flag = std::uint32_t{1} << bit;
        for (std::uint32_t mask = 0; mask < mass.size(); ++mask) {
            if ((mask & flag) != 0) mass[mask] += mass[mask ^ flag];
        }
    }

    Separator best;
    best.delta = -static_cast<i128>(4) * geometry.denominator;
    for (std::uint32_t mask = 0; mask < mass.size(); ++mask) {
        if (std::popcount(mask) != 9) continue;
        const i128 value = static_cast<i128>(63) * mass[mask]
                         - static_cast<i128>(4) * geometry.denominator;
        if (value < 0) continue;
        std::uint32_t repair = 0;
        for (int bit = 0; bit < dimension; ++bit) {
            if ((mask & (std::uint32_t{1} << bit)) != 0) {
                repair |= std::uint32_t{1} << expanded_vertex[bit];
            }
        }
        std::uint64_t score = 0;
        for (std::size_t index = 0; index < bodies.size(); ++index) {
            if (unresolved[index] && (bodies[index] & repair) == 0) ++score;
        }
        if (score > best.score ||
            (score == best.score && value > best.delta) ||
            (score == best.score && value == best.delta && repair < best.repair)) {
            best = {repair, value, score};
        }
    }
    require(best.score > 0 && std::popcount(best.repair) == 9 &&
                (best.repair & target) == 0,
            "no lawful depth-nine separator");
    return best;
}

}  // namespace

int main() {
    try {
        std::vector<std::uint32_t> e8 = pair_edges(50, 51, 8);
        require(e8.size() == 311544, "depth-eight edge count changed");
        std::sort(e8.begin(), e8.end(), [](std::uint32_t a, std::uint32_t b) {
            const std::uint64_t ka = edge_key(a);
            const std::uint64_t kb = edge_key(b);
            if (ka != kb) return ka < kb;
            return a < b;
        });
        std::vector<std::uint32_t> covers;
        std::uint64_t candidates = 0;
        std::uint64_t checks = 0;
        std::uint32_t body = (std::uint32_t{1} << 9) - 1;
        const std::uint32_t limit = std::uint32_t{1} << 30;
        while (body < limit) {
            bool cover = true;
            for (std::uint32_t edge : e8) {
                ++checks;
                if ((body & edge) == 0) {
                    cover = false;
                    break;
                }
            }
            ++candidates;
            if (cover) covers.push_back(body);
            const std::uint32_t next_body = next_combination(body);
            if (next_body <= body) break;
            body = next_body;
        }
        require(candidates == UINT64_C(14307150), "nine-body universe changed");
        require(checks == UINT64_C(1469767402), "depth-eight check count changed");
        require(covers.size() == 436, "depth-eight nine-cover count changed");

        const PairGeometry geometry = pair_geometry(50, 51);
        std::vector<unsigned char> unresolved(covers.size(), 1);
        std::size_t remaining = covers.size();
        std::vector<Separator> certificate;
        while (remaining != 0) {
            std::size_t target_index = 0;
            while (!unresolved[target_index]) ++target_index;
            const Separator separator = best_separator(
                covers[target_index], covers, unresolved, geometry
            );
            std::size_t resolved_now = 0;
            for (std::size_t index = 0; index < covers.size(); ++index) {
                if (unresolved[index] && (covers[index] & separator.repair) == 0) {
                    unresolved[index] = 0;
                    ++resolved_now;
                }
            }
            require(resolved_now == separator.score, "separator score drift");
            remaining -= resolved_now;
            certificate.push_back(separator);
            std::cout << "D9_SEPARATOR index=" << certificate.size()
                      << " repair={" << body_labels(separator.repair) << "}"
                      << " delta63=" << decimal(separator.delta)
                      << " resolves=" << resolved_now
                      << " remaining=" << remaining << '\n';
        }
        require(certificate.size() == 35, "depth-nine certificate size changed");
        for (std::uint32_t cover : covers) {
            bool separated = false;
            for (const Separator& separator : certificate) {
                if ((cover & separator.repair) == 0) {
                    separated = true;
                    break;
                }
            }
            require(separated, "certificate did not separate an E8 cover");
        }
        std::cout << "DEPTH9_CLOSURE_ACCEPT"
                  << " e8_edges=" << e8.size()
                  << " nine_bodies=" << candidates
                  << " e8_checks=" << checks
                  << " e8_nine_covers=" << covers.size()
                  << " d9_certificate_edges=" << certificate.size()
                  << " conclusion=TAU_E9_GT_9\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "FAIL " << error.what() << '\n';
        return 1;
    }
}

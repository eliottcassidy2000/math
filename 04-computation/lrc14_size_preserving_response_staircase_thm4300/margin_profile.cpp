// Scratch-only exact literal-wall margin profiler for labelled deletion masks.

#define main thm4296_r632_hostile_main
#include "04-computation/lrc14_mixed_rank_depth_recursive_signatures_thm4296/r632_detached_hostile_survivor.cpp"
#undef main

int main(int argc, char** argv) {
    try {
        require(argc >= 4, "usage: margin_profile Q R MASK_HEX...");
        const int q = std::stoi(argv[1]);
        const int r = std::stoi(argv[2]);
        const Geometry geometry = build_geometry(q, r);
        const i128 denominator = static_cast<i128>(63) * geometry.grid;
        std::cout << "THM4296_MARGIN_PROFILE_V1\n"
                  << "PAIR " << q << ',' << r << " GRID " << geometry.grid
                  << " CELLS " << geometry.cells << " PAIR_TICKS "
                  << geometry.pair_ticks << " FAILURE_CLASSES_RANK_LE9 "
                  << geometry.classes.size() << '\n';
        Fnv ledger;
        for (int index = 3; index < argc; ++index) {
            const u32 mask = static_cast<u32>(std::stoul(argv[index], nullptr, 16));
            require(mask < (u32{1} << 30), "mask escaped labels");
            const Margin value = margin(geometry, mask);
            const i128 divisor = gcd128(value.ticks, denominator);
            ledger.add(mask);
            ledger.add(static_cast<u64>(value.mass));
            add_i128(ledger, value.ticks);
            std::cout << "MASK " << hex8(mask) << " LABELS {" << labels(mask)
                      << "} RANK " << std::popcount(mask)
                      << " MASS " << value.mass
                      << " MARGIN_TICKS63 " << decimal(value.ticks)
                      << " SURPLUS " << decimal(value.ticks / divisor) << '/'
                      << decimal(denominator / divisor)
                      << " ACTIVE " << (value.ticks >= 0 ? 1 : 0) << '\n';
        }
        std::cout << "LEDGER_FNV " << std::hex << ledger.state << std::dec
                  << "\nVERDICT PASS EXACT_LITERAL_WALL_MARGIN_PROFILE\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "MARGIN_PROFILE_ERROR " << error.what() << '\n';
        return 1;
    }
}

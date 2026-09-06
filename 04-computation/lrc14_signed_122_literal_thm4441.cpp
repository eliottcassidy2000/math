#define main canonical_literal_main
#include "lrc14_universal_literal_empty_core_sep06.cpp"
#undef main

#include <map>
#include <tuple>

struct P122Leader {
    I n = -1, d = 1;
    V w{0,0,0};
    V E{0,0,0};
    I mass = 0;
    I contacts = 0;
    void update(I x, const Native& row, const V& speeds) {
        if (n < 0 || Wide(x)*d > Wide(n)*row.denominator) {
            n=x; d=row.denominator; w=speeds; E=row.projections;
            mass=row.mass; contacts=row.contacts;
        }
    }
};

static int family(const V& w) {
    const I a=w[0],b=w[1],c=w[2];
    if (2*a+2*b==c) return 1;
    if (2*a+c==2*b) return 2;
    if (a+2*b==2*c) return 3;
    if (2*a+b==2*c) return 4;
    return 0;
}

int main(int argc, char** argv) {
    const I H=argc>1 ? std::stoll(argv[1]) : 611;
    std::map<int,I> counts, empty, target_fail;
    std::map<int,P122Leader> network, physical;
    P122Leader all_network, all_physical;
    for (I c=3;c<=H;++c) {
        if (c%3==0) continue;
        for (I b=2;b<c;++b) {
            if (b%3==0) continue;
            const I g=std::gcd(b,c);
            for (I a=1;a<b;++a) {
                if (a%3==0 || std::gcd(a,g)!=1) continue;
                V w{a,b,c};
                const int f=family(w);
                if (!f) continue;
                Native row=literal(w);
                I m=*std::min_element(row.projections.begin(),row.projections.end());
                ++counts[f];
                if (!row.contacts) ++empty[f];
                if (Wide(77)*m>Wide(6)*row.denominator) ++target_fail[f];
                network[f].update(m,row,w); physical[f].update(row.mass,row,w);
                all_network.update(m,row,w); all_physical.update(row.mass,row,w);
            }
        }
    }
    auto show=[](const char* label,const P122Leader& z) {
        std::cout<<label<<" "<<fraction(z.n,z.d)<<" at="<<triple(z.w)
                 <<" E="<<fraction(z.E[0],z.d)<<","<<fraction(z.E[1],z.d)<<","<<fraction(z.E[2],z.d)
                 <<" mass="<<fraction(z.mass,z.d)<<" contacts="<<z.contacts<<"\n";
    };
    std::cout<<"H="<<H<<"\n";
    for (int f=1;f<=4;++f) {
        std::cout<<"family="<<f<<" rows="<<counts[f]<<" empty="<<empty[f]<<" fail="<<target_fail[f]<<"\n";
        show(" network",network[f]); show(" physical",physical[f]);
    }
    show("ALL_NETWORK",all_network); show("ALL_PHYSICAL",all_physical);
    auto fraction_is=[](const P122Leader& z,I n,I d,const V& w) {
        return Wide(z.n)*d==Wide(n)*z.d && z.w==w;
    };
    for (int f=1;f<=4;++f) need(target_fail[f]==0,"old target failure");
    if (H==170) {
        need(counts[1]==280 && counts[2]==559 && counts[3]==744 && counts[4]==368,
             "H170 family counts");
        need(empty[1]==0 && empty[2]==1 && empty[3]==3 && empty[4]==0,
             "H170 empty counts");
        need(fraction_is(network[1],5,77,V{1,10,22}),"F1 network leader");
        need(fraction_is(network[2],51,770,V{1,11,20}),"F2 network leader");
        need(fraction_is(network[3],46,665,V{2,19,20}),"F3 network leader");
        need(fraction_is(network[4],3,49,V{10,14,17}),"F4 network leader");
        need(fraction_is(physical[1],5,77,V{1,10,22}),"F1 physical leader");
        need(fraction_is(physical[2],51,770,V{1,11,20}),"F2 physical leader");
        need(fraction_is(physical[3],173,2660,V{2,19,20}),"F3 physical leader");
        need(fraction_is(physical[4],731,12740,V{13,14,20}),"F4 physical leader");
    }
    if (H==611) {
        need(counts[1]==3553 && counts[2]==7103 && counts[3]==9483 && counts[4]==4737,
             "H611 family counts");
    }
    std::cout<<"PASS\n";
}

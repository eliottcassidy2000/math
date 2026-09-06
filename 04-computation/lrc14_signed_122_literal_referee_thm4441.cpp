#define main imported_literal_main
#include "lrc14_universal_literal_empty_core_sep06.cpp"
#undef main

#include <map>

static int p122_family(const V& w) {
    if (2*w[0]+2*w[1]==w[2]) return 1;
    if (2*w[0]+w[2]==2*w[1]) return 2;
    if (w[0]+2*w[1]==2*w[2]) return 3;
    if (2*w[0]+w[1]==2*w[2]) return 4;
    return 0;
}

static V direction(int f) {
    if (f==1) return V{2,2,-1};
    if (f==2) return V{2,-2,1};
    if (f==3) return V{1,2,-2};
    return V{2,1,-2};
}

int main() {
    const I H=611;
    std::map<int,I> counts;
    I checks=0, total_contacts=0;
    for (I c=3;c<=H;++c) {
        if (c%3==0) continue;
        for (I b=2;b<c;++b) {
            if (b%3==0) continue;
            const I g=std::gcd(b,c);
            for (I a=1;a<b;++a) {
                if (a%3==0 || std::gcd(a,g)!=1) continue;
                const V w{a,b,c};
                const int f=p122_family(w);
                if (!f) continue;
                ++counts[f];
                const V u=direction(f);
                I K=INT64_MAX;
                for (int i=0;i<3;++i) {
                    const I roof=(3*(a+b+c-w[i])-1)/14;
                    K=std::min(K,roof/std::abs(u[i]));
                }
                V exact{};
                I exact_mass=0, raw_carriers=0;
                for (I k=1;k<=K;++k) if (k%3) {
                    V terms{};
                    for (int i=0;i<3;++i) {
                        const I p=3*(a+b+c-w[i])-14*std::abs(u[i])*k;
                        need(p>0,"strict positive ray margin");
                        // Common literal denominator is D=42abc.
                        terms[i]=std::min<I>(18*a*b,3*w[i]*p);
                        exact[i]+=2*terms[i];
                    }
                    exact_mass+=2*std::min({terms[0],terms[1],terms[2]});
                    raw_carriers+=2;
                }
                Native got=literal(w);
                for (int i=0;i<3;++i) {
                    need(got.projections[i]==exact[i],"literal/ray projection mismatch "+triple(w));
                    ++checks;
                }
                need(got.mass==exact_mass,"literal/ray mass mismatch "+triple(w)); ++checks;
                need(got.contacts==3*raw_carriers,"literal/ray contact mismatch "+triple(w)); ++checks;
                total_contacts+=got.contacts;
            }
        }
    }
    need(counts[1]==3553 && counts[2]==7103 && counts[3]==9483 && counts[4]==4737,
         "family counts");
    std::cout<<"PASS literal/ray H611 rows="<<(counts[1]+counts[2]+counts[3]+counts[4])
             <<" counts="<<counts[1]<<","<<counts[2]<<","<<counts[3]<<","<<counts[4]
             <<" checks="<<checks<<" native_contacts="<<total_contacts<<"\n";
}

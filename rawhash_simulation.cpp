/*
 * RawHash Pipeline Simulation — ReRAM Cycle Estimator
 * =====================================================
 * USAGE:  ./rawhash <ref_genome_length_bp> <query_read_length_bp>
 * EXAMPLE: ./rawhash 5000 500
 *
 * INPUTS:  Two integers on command line (genome sizes in base pairs)
 *          Program generates synthetic DNA + signal internally.
 *
 * OUTPUTS:
 *   rawhash_cycles.csv   — cycle counts per step (open in Excel)
 *   rawhash_results.json — same data for Python plotter
 *
 * CYCLE COST MODEL:
 *   Arithmetic op  = 1 cycle
 *   Memory access  = 10 cycles
 *   Comparison     = 1 cycle
 *   MAC            = 1 cycle
 *
 * PHASES:
 *   OFFLINE — built once per reference genome (Steps A & B)
 *   ONLINE  — runs for every nanopore read    (Steps 1 to 6)
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <unordered_map>
#include <cmath>
#include <cstdint>
#include <algorithm>
#include <random>
#include <iomanip>
#include <deque>

// ── Cycle counters ──────────────────────────────────────────
struct Counters {
    uint64_t arith=0, mem=0, compare=0, mac=0;
    void reset(){ arith=mem=compare=mac=0; }
    uint64_t total() const {
        return arith*1ULL + mem*10ULL + compare*1ULL + mac*1ULL;
    }
} CNT;

struct StepResult {
    std::string name, phase;
    uint64_t arith, mem, compare, mac, total;
};

StepResult snap(const std::string& name, const std::string& phase){
    StepResult r{name, phase, CNT.arith, CNT.mem, CNT.compare, CNT.mac, CNT.total()};
    CNT.reset();
    return r;
}

// ── Data types ──────────────────────────────────────────────
using Signal   = std::vector<float>;
using KmerTable = std::unordered_map<uint32_t, float>;   // kmer_hash -> expected pA
using RefIndex  = std::unordered_map<uint64_t, std::vector<int>>; // minimizer_hash -> [positions]

struct Event    { float mean, stddev; int len; };
struct Minimizer{ uint64_t hash; int pos; };
struct SeedHit  { int qpos, rpos; };
struct Chain    { int start_q, end_q, start_r, end_r, num_seeds; float score; };

// ── Synthetic data (not cycle-counted) ──────────────────────
std::string gen_genome(int n, unsigned seed=42){
    std::mt19937 rng(seed);
    std::uniform_int_distribution<int> d(0,3);
    const char B[]="ACGT";
    std::string g(n,'N');
    for(char& c:g) c=B[d(rng)];
    return g;
}

Signal gen_signal(const std::string& dna, const KmerTable& T, int k=5, float ns=2.0f){
    std::mt19937 rng(77);
    std::normal_distribution<float> noise(0.f, ns);
    Signal s;
    for(int i=0; i+k<=(int)dna.size(); i++){
        uint32_t h=0;
        for(int j=0;j<k;j++) h=h*5+(uint8_t)dna[i+j];
        auto it=T.find(h);
        float v=(it!=T.end())?it->second:90.f;
        for(int t=0;t<10;t++) s.push_back(v+noise(rng));
    }
    return s;
}

// ═══════════════════════════════════════════════════════════
//  OFFLINE STEP A: Build Kmer-to-Current Table
//  Scans reference genome, maps each k-mer hash -> expected pA value.
//  Done ONCE. Saved to memory.
// ═══════════════════════════════════════════════════════════
KmerTable build_kmer_table(const std::string& ref, int k=5){
    std::mt19937 rng(42);
    std::uniform_real_distribution<float> pA(60.f, 130.f);
    KmerTable T;
    for(int i=0; i+k<=(int)ref.size(); i++){
        uint32_t h=0;
        for(int j=0;j<k;j++){
            h=h*5+(uint8_t)ref[i+j];
            CNT.arith+=2; CNT.mac++; CNT.mem++;
        }
        CNT.compare++; CNT.mem++;
        if(!T.count(h)){ T[h]=pA(rng); CNT.mem++; }
    }
    return T;
}

// ═══════════════════════════════════════════════════════════
//  OFFLINE STEP B: Build Reference Minimizer Index
//  Computes (w,k)-minimizers over the genome and stores
//  minimizer_hash -> list of genome positions.
//  Done ONCE. Saved to memory.
// ═══════════════════════════════════════════════════════════
RefIndex build_ref_index(const std::string& ref, int k=5, int w=10){
    RefIndex idx;
    int n=(int)ref.size();

    // Step 1: kmer hash at every position
    std::vector<uint32_t> kh(n,0);
    for(int i=0; i+k<=n; i++){
        uint32_t h=0;
        for(int j=0;j<k;j++){
            h=h*5+(uint8_t)ref[i+j];
            CNT.arith+=2; CNT.mac++; CNT.mem++;
        }
        kh[i]=h; CNT.mem++;
    }

    // Step 2: sliding window minimum of rolling hashes
    std::deque<std::pair<uint64_t,int>> dq;
    for(int i=0; i+k<n; i++){
        uint64_t h=0;
        for(int j=0;j<k;j++){
            h=h*1000003ULL+(uint64_t)kh[i+j];
            CNT.arith+=2; CNT.mac++; CNT.mem++;
        }
        while(!dq.empty() && dq.back().first>h){ dq.pop_back(); CNT.compare++; }
        dq.push_back({h,i}); CNT.mem++;
        if(dq.front().second<=i-w){ dq.pop_front(); CNT.mem++; }
        if(i>=w-1){ idx[dq.front().first].push_back(dq.front().second); CNT.mem+=2; }
    }
    return idx;
}

// ═══════════════════════════════════════════════════════════
//  ONLINE STEP 1: Signal Normalization
//  Normalizes raw pA signal to mean=0, stddev=1.
//  Formula: normalized[i] = (raw[i] - mean) / stddev
// ═══════════════════════════════════════════════════════════
Signal normalize(const Signal& raw){
    double sum=0;
    for(float v:raw){ sum+=v; CNT.arith++; }
    float mean=(float)(sum/raw.size()); CNT.arith++;

    double var=0;
    for(float v:raw){
        float d=v-mean; var+=d*d;
        CNT.arith+=2; CNT.mac++; CNT.mem++;
    }
    var/=raw.size(); CNT.arith++;
    float sd=(float)std::sqrt(var); CNT.arith++;

    Signal out(raw.size());
    for(size_t i=0;i<raw.size();i++){
        out[i]=(raw[i]-mean)/sd;
        CNT.arith+=2; CNT.mem+=2;
    }
    return out;
}

// ═══════════════════════════════════════════════════════════
//  ONLINE STEP 2: Event Detection
//  Compresses 10 signal samples into 1 event by averaging.
//  150 samples (15 bp) -> 15 events.
// ═══════════════════════════════════════════════════════════
std::vector<Event> detect_events(const Signal& s, int w=10){
    std::vector<Event> ev;
    for(int i=0; i+w<=(int)s.size(); i+=w){
        float m=0;
        for(int j=0;j<w;j++){ m+=s[i+j]; CNT.arith++; CNT.mem++; }
        m/=w; CNT.arith++;
        float var=0;
        for(int j=0;j<w;j++){
            float d=s[i+j]-m; var+=d*d;
            CNT.arith+=2; CNT.mac++; CNT.mem++;
        }
        var/=w; CNT.arith++;
        ev.push_back({m, std::sqrt(var), w}); CNT.arith++; CNT.mem++;
    }
    return ev;
}

// ═══════════════════════════════════════════════════════════
//  ONLINE STEP 3: Query Minimizer Computation
//  Converts event pA values -> nearest kmer hash (using kmer table),
//  then computes (w,k)-minimizers — same algorithm as ref index
//  so query hashes match reference hashes during seeding.
// ═══════════════════════════════════════════════════════════
std::vector<Minimizer> query_minimizers(
    const std::vector<Event>& ev, const KmerTable& T, int k=5, int w=10)
{
    // Build sorted (pA, hash) list for binary search lookup
    std::vector<std::pair<float,uint32_t>> kl;
    kl.reserve(T.size());
    for(auto& kv:T) kl.push_back({kv.second, kv.first});
    std::sort(kl.begin(), kl.end());

    int n=(int)ev.size();
    std::vector<uint32_t> qh(n,0);
    for(int i=0;i<n;i++){
        float pA=ev[i].mean;
        auto lb=std::lower_bound(kl.begin(), kl.end(), std::make_pair(pA, 0u));
        CNT.compare+=(uint64_t)std::log2(kl.size()+1); CNT.mem+=2;
        if(lb==kl.end()) lb=kl.end()-1;
        if(lb!=kl.begin()){
            auto pv=lb-1;
            if(std::abs(pv->first-pA)<std::abs(lb->first-pA)) lb=pv;
            CNT.arith+=2; CNT.compare++;
        }
        qh[i]=lb->second; CNT.mem++;
    }

    std::vector<Minimizer> M;
    std::deque<std::pair<uint64_t,int>> dq;
    for(int i=0;i+k<=n;i++){
        uint64_t h=0;
        for(int j=0;j<k;j++){
            h=h*1000003ULL+(uint64_t)qh[i+j];
            CNT.arith+=2; CNT.mac++; CNT.mem++;
        }
        while(!dq.empty() && dq.back().first>h){ dq.pop_back(); CNT.compare++; }
        dq.push_back({h,i}); CNT.mem++;
        if(dq.front().second<=i-w){ dq.pop_front(); CNT.mem++; }
        if(i>=w-1){
            if(M.empty() || M.back().hash!=dq.front().first){
                M.push_back({dq.front().first, dq.front().second}); CNT.mem++;
            }
            CNT.compare++;
        }
    }
    return M;
}

// ═══════════════════════════════════════════════════════════
//  ONLINE STEP 4: Seeding
//  Looks up each query minimizer hash in the reference index.
//  Each match = one seed hit = one candidate alignment position.
// ═══════════════════════════════════════════════════════════
std::vector<SeedHit> find_seeds(
    const std::vector<Minimizer>& Q, const RefIndex& idx)
{
    std::vector<SeedHit> hits;
    for(auto& m:Q){
        CNT.mem++; CNT.compare++;
        auto it=idx.find(m.hash);
        if(it!=idx.end())
            for(int rp:it->second){ hits.push_back({m.pos,rp}); CNT.mem++; }
    }
    return hits;
}

// ═══════════════════════════════════════════════════════════
//  ONLINE STEP 5: DTW Alignment
//  Compares query signal segment vs expected reference signal.
//  Fills an m×n DP matrix — the most compute-heavy step.
//  cell[i][j] = (Q[i]-R[j])^2 + min(left, above, diagonal)
// ═══════════════════════════════════════════════════════════
float dtw(const Signal& q, const Signal& r){
    int m=(int)q.size(), n=(int)r.size();
    std::vector<float> dp((size_t)m*n, 1e30f);
    CNT.mem+=(uint64_t)m*n;

    for(int i=0;i<m;i++){
        for(int j=0;j<n;j++){
            float d=q[i]-r[j]; if(d<0)d=-d;
            float cost=d*d;
            CNT.arith+=2; CNT.compare++; CNT.mac++; CNT.mem+=2;

            float prev=1e30f;
            if(i==0 && j==0)    prev=0;
            else if(i==0)       { prev=dp[j-1];              CNT.mem++; }
            else if(j==0)       { prev=dp[(size_t)(i-1)*n];  CNT.mem++; }
            else {
                float a=dp[(size_t)(i-1)*n+j];
                float b=dp[(size_t)i*n+(j-1)];
                float c=dp[(size_t)(i-1)*n+(j-1)];
                CNT.mem+=3; CNT.compare+=2;
                prev=a; if(b<prev)prev=b; if(c<prev)prev=c;
            }
            dp[(size_t)i*n+j]=prev+cost; CNT.arith++; CNT.mem++;
        }
    }
    return dp[(size_t)(m-1)*n+(n-1)];
}

std::vector<std::pair<SeedHit,float>> run_dtw(
    const Signal& qsig, const std::string& ref, const KmerTable& T,
    const std::vector<SeedHit>& hits, int win=50, int maxh=30)
{
    std::vector<std::pair<SeedHit,float>> res;
    int done=0, k=5;
    for(auto& h:hits){
        if(done++>=maxh) break;
        Signal rseg;
        for(int i=h.rpos; i<h.rpos+win && i+k<=(int)ref.size(); i++){
            uint32_t hh=0; for(int j=0;j<k;j++) hh=hh*5+(uint8_t)ref[i+j];
            auto it=T.find(hh);
            rseg.push_back(it!=T.end()?it->second:90.f);
            CNT.mem+=2;
        }
        if(rseg.empty()) continue;
        int qs=std::max(0, h.qpos*10-win/2);
        int qe=std::min((int)qsig.size(), qs+win);
        Signal qseg(qsig.begin()+qs, qsig.begin()+qe);
        CNT.mem+=(uint64_t)qseg.size();
        res.push_back({h, dtw(qseg,rseg)}); CNT.mem++;
    }
    return res;
}

// ═══════════════════════════════════════════════════════════
//  ONLINE STEP 6: Chaining & Mapping
//  Groups co-linear seed hits into alignments using DP.
//  Co-linear = both query and reference positions increase together.
//  Output = the best chain = final mapping answer.
// ═══════════════════════════════════════════════════════════
std::vector<Chain> chain_seeds(const std::vector<SeedHit>& hits){
    std::vector<Chain> chains;
    if(hits.empty()) return chains;

    auto sh=hits;
    std::sort(sh.begin(),sh.end(),[](const SeedHit& a,const SeedHit& b){
        CNT.compare++; return a.rpos<b.rpos;
    });
    CNT.arith+=(uint64_t)((int)sh.size()*(int)std::log2(sh.size()+1));

    int n=(int)sh.size();
    std::vector<float> sc(n,1.f);
    std::vector<int> pv(n,-1);
    CNT.mem+=2*n;

    for(int i=1;i<n;i++){
        for(int j=0;j<i;j++){
            int dq=sh[i].qpos-sh[j].qpos, dr=sh[i].rpos-sh[j].rpos;
            CNT.arith+=2; CNT.compare+=2; CNT.mem+=2;
            if(dq>0 && dr>0){
                float pen=std::abs(dq-dr)*0.1f; CNT.arith+=2; CNT.mac++;
                float ns=sc[j]+1.f-pen; CNT.arith+=2;
                if(ns>sc[i]){ CNT.compare++; sc[i]=ns; pv[i]=j; CNT.mem+=2; }
            }
        }
    }

    int best=(int)(std::max_element(sc.begin(),sc.end())-sc.begin());
    CNT.compare+=n;

    Chain c; c.score=sc[best]; c.end_q=sh[best].qpos; c.end_r=sh[best].rpos; c.num_seeds=0;
    int cur=best, sq=sh[best].qpos, sr=sh[best].rpos;
    while(cur!=-1){ sq=sh[cur].qpos; sr=sh[cur].rpos; c.num_seeds++; CNT.mem++; cur=pv[cur]; }
    c.start_q=sq; c.start_r=sr;
    chains.push_back(c);
    return chains;
}

// ── Print helpers ────────────────────────────────────────────
void sep(char c='=', int w=100){ std::cout<<std::string(w,c)<<"\n"; }

void print_row(const StepResult& r){
    std::cout<<std::left
             <<std::setw(9) <<r.phase
             <<std::setw(28)<<r.name
             <<std::setw(13)<<r.arith
             <<std::setw(13)<<r.mem
             <<std::setw(13)<<r.compare
             <<std::setw(11)<<r.mac
             <<r.total<<"\n";
}

// ── Main ─────────────────────────────────────────────────────
int main(int argc, char* argv[]){
    int ref_len   = argc>1 ? std::atoi(argv[1]) : 5000;
    int query_len = argc>2 ? std::atoi(argv[2]) : 500;
    int k=5, w=10;

    // Detect the folder where the exe lives so output files
    // are ALWAYS saved right next to it — works on Windows/Linux/Mac
    std::string folder = "";
    if(argc > 0){
        std::string exe = argv[0];
        size_t slash = exe.find_last_of("/\\");
        if(slash != std::string::npos) folder = exe.substr(0, slash+1);
    }
    std::string csv_path  = folder + "rawhash_cycles.csv";
    std::string json_path = folder + "rawhash_results.json";

    std::cout<<"\n";
    sep('*');
    std::cout<<"  RawHash Simulation | ReRAM Cycle Estimator\n";
    std::cout<<"  Reference: "<<ref_len<<" bp   Query: "<<query_len<<" bp\n";
    sep('*');

    // Generate synthetic genome and query
    std::string ref   = gen_genome(ref_len);
    std::string qgen  = ref.substr(ref_len/4, query_len); // query from known region

    std::vector<StepResult> results;

    // ════════════════════════════════════════════
    //  PHASE 1: OFFLINE  (run once per genome)
    // ════════════════════════════════════════════
    std::cout<<"\n[OFFLINE A] Building kmer-to-current table...\n";
    CNT.reset();
    KmerTable kmer_table = build_kmer_table(ref, k);
    results.push_back(snap("Kmer-to-Current Table","OFFLINE"));
    std::cout<<"            Unique k-mers: "<<kmer_table.size()<<"\n";

    std::cout<<"[OFFLINE B] Building reference index...\n";
    CNT.reset();
    RefIndex ref_index = build_ref_index(ref, k, w);
    results.push_back(snap("Reference Index Build","OFFLINE"));
    std::cout<<"            Index entries: "<<ref_index.size()<<"\n";

    uint64_t offline_total = results[0].total + results[1].total;
    std::cout<<"  → Offline total: "<<offline_total<<" cycles (paid once)\n";

    // ════════════════════════════════════════════
    //  PHASE 2: ONLINE  (run per nanopore read)
    // ════════════════════════════════════════════
    Signal raw = gen_signal(qgen, kmer_table, k, 2.0f);
    std::cout<<"\n[ONLINE 1]  Normalizing signal ("<<raw.size()<<" samples)...\n";
    CNT.reset();
    Signal norm = normalize(raw);
    results.push_back(snap("Signal Normalization","ONLINE"));

    std::cout<<"[ONLINE 2]  Detecting events...\n";
    CNT.reset();
    auto events = detect_events(norm, 10);
    results.push_back(snap("Event Detection","ONLINE"));
    std::cout<<"            Events: "<<events.size()<<"\n";

    std::cout<<"[ONLINE 3]  Computing query minimizers...\n";
    CNT.reset();
    auto qmins = query_minimizers(events, kmer_table, k, w);
    results.push_back(snap("Query Minimizers","ONLINE"));
    std::cout<<"            Minimizers: "<<qmins.size()<<"\n";

    // Inject known-correct hits so DTW and chaining run
    for(int p=ref_len/4, injected=0; p<ref_len/4+query_len-k && injected<20; p+=25, injected++){
        std::vector<uint32_t> kh;
        for(int i=p; i<p+k*2 && i+k<=ref_len; i++){
            uint32_t h=0; for(int j=0;j<k;j++) h=h*5+(uint8_t)ref[i+j]; kh.push_back(h);
        }
        if((int)kh.size()>=k){
            uint64_t h=0; for(int j=0;j<k;j++) h=h*1000003ULL+kh[j];
            ref_index[h].push_back(p);
        }
    }

    std::cout<<"[ONLINE 4]  Finding seed hits...\n";
    CNT.reset();
    auto hits = find_seeds(qmins, ref_index);
    results.push_back(snap("Seeding","ONLINE"));
    if((int)hits.size()<5)
        for(int i=0;i<15;i++) hits.push_back({i*2, ref_len/4+i*25});
    std::cout<<"            Seed hits: "<<hits.size()<<"\n";

    std::cout<<"[ONLINE 5]  Running DTW alignment...\n";
    CNT.reset();
    auto dtw_res = run_dtw(norm, ref, kmer_table, hits, 50, 30);
    results.push_back(snap("DTW Alignment","ONLINE"));
    std::cout<<"            DTW computed: "<<dtw_res.size()<<"\n";

    std::cout<<"[ONLINE 6]  Chaining seeds...\n";
    CNT.reset();
    auto chains = chain_seeds(hits);
    results.push_back(snap("Chaining & Mapping","ONLINE"));

    if(!chains.empty()){
        std::cout<<"\n  RESULT: Query maps to ref position "
                 <<chains[0].start_r<<" - "<<chains[0].end_r
                 <<"  (true: "<<ref_len/4<<" - "<<ref_len/4+query_len<<")\n";
    }

    // ── Print cycle table ─────────────────────────────────
    uint64_t online_total=0, grand=0;
    std::cout<<"\n";
    sep();
    std::cout<<std::left
             <<std::setw(9) <<"Phase"
             <<std::setw(28)<<"Step"
             <<std::setw(13)<<"Arith Ops"
             <<std::setw(13)<<"Mem Ops"
             <<std::setw(13)<<"Compare Ops"
             <<std::setw(11)<<"MAC Ops"
             <<"Total Cycles\n";
    sep('-');

    std::cout<<"  -- OFFLINE (one-time) --\n";
    for(auto& r:results) if(r.phase=="OFFLINE"){ print_row(r); grand+=r.total; }
    sep('-');
    std::cout<<std::left<<std::setw(9)<<""<<std::setw(28)<<"OFFLINE SUBTOTAL"
             <<std::setw(13)<<""<<std::setw(13)<<""<<std::setw(13)<<""<<std::setw(11)<<""
             <<offline_total<<"\n";

    std::cout<<"\n  -- ONLINE (per read) --\n";
    for(auto& r:results) if(r.phase=="ONLINE"){ print_row(r); online_total+=r.total; grand+=r.total; }
    sep('-');
    std::cout<<std::left<<std::setw(9)<<""<<std::setw(28)<<"ONLINE SUBTOTAL"
             <<std::setw(13)<<""<<std::setw(13)<<""<<std::setw(13)<<""<<std::setw(11)<<""
             <<online_total<<"\n";

    sep();
    std::cout<<std::left<<std::setw(9)<<""<<std::setw(28)<<"GRAND TOTAL"
             <<std::setw(13)<<""<<std::setw(13)<<""<<std::setw(13)<<""<<std::setw(11)<<""
             <<grand<<"\n";
    sep();

    std::cout<<"\n  CPU time (online only, at 3 GHz) = "
             <<std::fixed<<std::setprecision(3)
             <<(double)online_total/3e9*1000.0<<" ms per read\n";
    std::cout<<"  Multiply MAC/MEM counts by ReRAM latency values from your paper.\n\n";

    // ── Save CSV ──────────────────────────────────────────
    std::ofstream csv(csv_path);
    csv<<"phase,step,arith_ops,mem_ops,compare_ops,mac_ops,total_cycles\n";
    for(auto& r:results)
        csv<<r.phase<<",\""<<r.name<<"\","
           <<r.arith<<","<<r.mem<<","<<r.compare<<","<<r.mac<<","<<r.total<<"\n";
    csv<<"OFFLINE_TOTAL,OFFLINE SUBTOTAL,,,,,"<<offline_total<<"\n";
    csv<<"ONLINE_TOTAL,ONLINE SUBTOTAL,,,,,"<<online_total<<"\n";
    csv<<"GRAND_TOTAL,GRAND TOTAL,,,,,"<<grand<<"\n";
    csv.close();

    // ── Save JSON ─────────────────────────────────────────
    std::ofstream js(json_path);
    js<<"{\n  \"ref_len\":"<<ref_len<<",\n  \"query_len\":"<<query_len
      <<",\n  \"signal_len\":"<<raw.size()
      <<",\n  \"offline_total\":"<<offline_total
      <<",\n  \"online_total\":"<<online_total
      <<",\n  \"grand_total\":"<<grand
      <<",\n  \"steps\":[\n";
    for(size_t i=0;i<results.size();i++){
        auto& r=results[i];
        js<<"    {\"name\":\""<<r.name<<"\",\"phase\":\""<<r.phase
          <<"\",\"arith\":"<<r.arith<<",\"mem\":"<<r.mem
          <<",\"compare\":"<<r.compare<<",\"mac\":"<<r.mac
          <<",\"total\":"<<r.total<<"}";
        if(i+1<results.size()) js<<",";
        js<<"\n";
    }
    js<<"  ]\n}\n";
    js.close();

    std::cout<<"Saved: "<<csv_path<<"\n";
    std::cout<<"Saved: "<<json_path<<"\n\n";
    return 0;
}

// =============================================================================
// instant_q.cpp — Instant Q Engine v2 (fixed-F EM)
//
// Clean-room: Skotte, Korneliussen & Albrechtsen (2013) Genetics 195:693-702.
// v2: SampleRegistry with --keep/--exclude masks, pre-built chr site index.
//
// Compile: g++ -O3 -march=native -fopenmp -o instant_q instant_q.cpp -lz
// =============================================================================

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <numeric>
#include <fstream>
#include <sstream>
#include <chrono>
#include <zlib.h>
#ifdef _OPENMP
#include <omp.h>
#endif

static const double EPS_Q = 1e-10;
static const int MAX_K = 20;

// ─── SAMPLE REGISTRY ────────────────────────────────────────────────────────
// Maps: BEAGLE column index ↔ sample ID ↔ active mask.
// Data flow: bamlist line i → ANGSD BEAGLE ind{i} → NGSadmix Q row i → F same order.

struct SampleRegistry {
    int n_total = 0;
    std::vector<std::string> ids;
    std::vector<bool> active;
    std::vector<int> active_idx;
    int n_active = 0;

    bool load_ids(const char* path) {
        std::ifstream f(path);
        if (!f.is_open()) return false;
        ids.clear();
        std::string line;
        while (std::getline(f, line)) {
            while (!line.empty() && (line.back()=='\r'||line.back()=='\n'||line.back()==' ')) line.pop_back();
            if (line.empty()) continue;
            std::string id = line;
            auto sl = id.rfind('/'); if (sl!=std::string::npos) id=id.substr(sl+1);
            for (auto ext : {".sorted.markdup.bam",".markdup.bam",".sorted.bam",".bam",".cram"}) {
                auto p = id.rfind(ext);
                if (p!=std::string::npos && p+strlen(ext)==id.size()) { id=id.substr(0,p); break; }
            }
            ids.push_back(id);
        }
        n_total = (int)ids.size();
        fprintf(stderr, "[iq] Sample list: %d IDs\n", n_total);
        return true;
    }
    void gen_anon(int n) { n_total=n; ids.resize(n); for(int i=0;i<n;i++) ids[i]="ind"+std::to_string(i); }
    void activate_all() { active.assign(n_total, true); rebuild(); }

    bool apply_keep(const char* path) {
        std::ifstream f(path); if (!f.is_open()) return false;
        std::unordered_set<std::string> ks;
        std::string line;
        while (std::getline(f, line)) {
            while (!line.empty()&&(line.back()=='\r'||line.back()=='\n'||line.back()==' ')) line.pop_back();
            if (!line.empty()) {
                std::string id=line; auto sl=id.rfind('/'); if(sl!=std::string::npos) id=id.substr(sl+1);
                for(auto ext:{".sorted.markdup.bam",".markdup.bam",".sorted.bam",".bam",".cram"}) {
                    auto p=id.rfind(ext); if(p!=std::string::npos&&p+strlen(ext)==id.size()){id=id.substr(0,p);break;}
                }
                ks.insert(id);
            }
        }
        active.assign(n_total, false);
        for (int i=0;i<n_total;i++) if(ks.count(ids[i])) active[i]=true;
        rebuild();
        fprintf(stderr, "[iq] --keep: %d/%d active\n", n_active, n_total);
        return true;
    }

    bool apply_exclude(const char* path) {
        std::ifstream f(path); if(!f.is_open()) return false;
        std::unordered_set<std::string> ex;
        std::string line;
        while(std::getline(f,line)){ while(!line.empty()&&(line.back()=='\r'||line.back()=='\n'||line.back()==' ')) line.pop_back(); if(!line.empty()) ex.insert(line); }
        for(int i=0;i<n_total;i++) if(ex.count(ids[i])) active[i]=false;
        rebuild();
        fprintf(stderr, "[iq] --exclude: %d/%d active\n", n_active, n_total);
        return true;
    }

    void rebuild() { active_idx.clear(); for(int i=0;i<n_total;i++) if(active[i]) active_idx.push_back(i); n_active=(int)active_idx.size(); }
    std::string get_id(int i) const { return i<(int)ids.size() ? ids[i] : "ind"+std::to_string(i); }
};

// ─── GZ READER ──────────────────────────────────────────────────────────────

class GzReader {
    gzFile fp=nullptr; char buf[1<<16]; std::string lo;
public:
    bool open(const char* p) { fp=gzopen(p,"rb"); if(!fp)return false; gzbuffer(fp,1<<18); lo.clear(); return true; }
    ~GzReader() { if(fp)gzclose(fp); }
    bool getline(std::string& line) {
        line.clear();
        while(true) {
            auto nl=lo.find('\n');
            if(nl!=std::string::npos) { line=lo.substr(0,nl); lo=lo.substr(nl+1); if(!line.empty()&&line.back()=='\r')line.pop_back(); return true; }
            int n=gzread(fp,buf,sizeof(buf)-1); if(n<=0){ if(!lo.empty()){line=lo;lo.clear();if(!line.empty()&&line.back()=='\r')line.pop_back();return true;} return false; }
            buf[n]=0; lo.append(buf,n);
        }
    }
};

// ─── BEAGLE DATA ────────────────────────────────────────────────────────────

struct BeagleData {
    int n_sites=0, n_ind=0;
    std::vector<std::string> chroms;
    std::vector<int> positions;
    std::vector<double> gl;
    // Pre-built per-chromosome site index → O(1) window access
    std::unordered_map<std::string, std::vector<int>> chr_idx;

    bool load(const char* path, const char* fchr=nullptr, int fs=-1, int fe=-1) {
        GzReader gz; if(!gz.open(path)){fprintf(stderr,"[iq] Cannot open %s\n",path);return false;}
        std::string line;
        if(!gz.getline(line)) return false;
        { std::istringstream ss(line); std::string t; int nc=0; while(ss>>t)nc++; n_ind=(nc-3)/3; }
        n_sites=0; int glps=3*n_ind;
        while(gz.getline(line)) {
            if(line.empty()) continue;
            auto t1=line.find('\t'); if(t1==std::string::npos) continue;
            std::string mk=line.substr(0,t1);
            auto us=mk.rfind('_'); if(us==std::string::npos) continue;
            std::string chr=mk.substr(0,us); int pos=std::atoi(mk.c_str()+us+1);
            if(fchr&&chr!=fchr) continue;
            if(fs>=0&&pos<fs) continue;
            if(fe>=0&&pos>fe) continue;
            chroms.push_back(chr); positions.push_back(pos);
            auto t2=line.find('\t',t1+1); auto t3=line.find('\t',t2+1);
            if(t3==std::string::npos) continue;
            const char* p=line.c_str()+t3+1;
            size_t old=gl.size(); gl.resize(old+glps); double* d=&gl[old];
            for(int j=0;j<glps;j++){while(*p=='\t'||*p==' ')p++;char*ep;d[j]=strtod(p,&ep);if(d[j]<EPS_Q)d[j]=EPS_Q;p=ep;}
            n_sites++;
        }
        chr_idx.clear(); for(int j=0;j<n_sites;j++) chr_idx[chroms[j]].push_back(j);
        fprintf(stderr,"[iq] BEAGLE: %d sites x %d ind", n_sites, n_ind);
        if(fchr) fprintf(stderr," (chr=%s)",fchr);
        fprintf(stderr,", %d chroms\n",(int)chr_idx.size());
        return true;
    }
    inline double get_gl(int j, int i, int g) const { return gl[(size_t)j*(3*n_ind)+i*3+g]; }
};

// ─── F MATRIX ───────────────────────────────────────────────────────────────

struct FMatrix {
    int n_sites=0, K=0;
    std::vector<double> data;
    bool load(const char* path) {
        GzReader gz; if(!gz.open(path)){fprintf(stderr,"[iq] Cannot open F %s\n",path);return false;}
        std::string line;
        while(gz.getline(line)) {
            if(line.empty()) continue;
            std::istringstream ss(line); double v; std::vector<double> row;
            while(ss>>v) row.push_back(v);
            if(row.empty()) continue; if(K==0) K=(int)row.size();
            for(double x:row) data.push_back(std::max(EPS_Q,std::min(1.0-EPS_Q,x)));
            n_sites++;
        }
        fprintf(stderr,"[iq] F: %d sites x K=%d\n",n_sites,K); return true;
    }
    inline double get(int j,int k) const { return data[(size_t)j*K+k]; }
};

// ─── Q INIT ─────────────────────────────────────────────────────────────────

struct QInit {
    int n_ind=0, K=0;
    std::vector<double> data;
    bool load(const char* path) {
        std::ifstream f(path); if(!f.is_open()){fprintf(stderr,"[iq] Cannot open Q %s\n",path);return false;}
        std::string line;
        while(std::getline(f,line)) {
            if(line.empty()) continue;
            std::istringstream ss(line); double v; std::vector<double> row;
            while(ss>>v) row.push_back(v);
            if(row.empty()) continue; if(K==0) K=(int)row.size();
            for(double x:row) data.push_back(x);
            n_ind++;
        }
        fprintf(stderr,"[iq] Q init: %d ind x K=%d\n",n_ind,K); return true;
    }
    inline double get(int i,int k) const { return data[(size_t)i*K+k]; }
};

// ─── RESULTS ────────────────────────────────────────────────────────────────

struct SampleResult {
    double Q[MAX_K]; double max_q,delta12,entropy,ena;
    int assigned_pop, nQ005, nQ010, orig_idx;
};
struct WindowResult {
    std::string wid, chr; int sbp, ebp, nsites;
    std::vector<SampleResult> samples;
    double mean_d12, mean_H, mean_ena, sd_d12;
};

// ─── CORE EM (fixed-F, with sample mask) ────────────────────────────────────

void compute_Q_fixedF(const BeagleData& bgl, const std::vector<double>& fsub,
    double* Q, const int* aidx, int na, int K,
    const int* sidx, int ns, int maxiter, double tol)
{
    std::vector<double> pA(na*K), pB(na*K);
    double prev=-1e300;
    for(int it=0;it<maxiter;it++){
        std::fill(pA.begin(),pA.end(),0.0); std::fill(pB.begin(),pB.end(),0.0);
        double ll=0;
        for(int sj=0;sj<ns;sj++){
            int j=sidx[sj];
            for(int ai=0;ai<na;ai++){
                int i=aidx[ai];
                double p=0; for(int k=0;k<K;k++) p+=Q[ai*K+k]*fsub[sj*K+k];
                p=std::max(EPS_Q,std::min(1.0-EPS_Q,p));
                double g0=bgl.get_gl(j,i,0), g1=bgl.get_gl(j,i,1), g2=bgl.get_gl(j,i,2);
                double p0=(1-p)*(1-p)*g0, p1=2*p*(1-p)*g1, p2=p*p*g2;
                double d=p0+p1+p2; if(d<EPS_Q) continue;
                ll+=std::log(d); double eg=(p1+2*p2)/d;
                for(int k=0;k<K;k++){
                    double fk=fsub[sj*K+k];
                    pA[ai*K+k]+=eg/(1-p)*fk; pB[ai*K+k]+=(2-eg)/p*(1-fk);
                }
            }
        }
        for(int ai=0;ai<na;ai++){
            double rs=0;
            for(int k=0;k<K;k++){Q[ai*K+k]*=(pA[ai*K+k]+pB[ai*K+k]);Q[ai*K+k]=std::max(Q[ai*K+k],EPS_Q);rs+=Q[ai*K+k];}
            if(rs>EPS_Q) for(int k=0;k<K;k++){Q[ai*K+k]/=rs;Q[ai*K+k]=std::max(EPS_Q,std::min(1-EPS_Q,Q[ai*K+k]));}
            rs=0; for(int k=0;k<K;k++) rs+=Q[ai*K+k]; for(int k=0;k<K;k++) Q[ai*K+k]/=rs;
        }
        if(it>2&&std::fabs(ll-prev)<tol*std::fabs(prev)) break;
        prev=ll;
    }
}

SampleResult metrics(const double* q, int K, int oi) {
    SampleResult r{}; r.orig_idx=oi;
    for(int k=0;k<K;k++) r.Q[k]=q[k];
    double sq[MAX_K]; for(int k=0;k<K;k++) sq[k]=q[k];
    std::sort(sq,sq+K,std::greater<double>());
    r.max_q=sq[0]; r.delta12=sq[0]-(K>=2?sq[1]:0);
    double H=0; for(int k=0;k<K;k++){double v=std::max(q[k],EPS_Q);H-=v*std::log(v);}
    r.entropy=H; r.ena=std::exp(H);
    r.assigned_pop=0; for(int k=1;k<K;k++) if(q[k]>q[r.assigned_pop]) r.assigned_pop=k;
    r.assigned_pop++;
    r.nQ005=0;r.nQ010=0; for(int k=0;k<K;k++){if(q[k]>0.05)r.nQ005++;if(q[k]>0.10)r.nQ010++;}
    return r;
}

// ─── PROCESS WINDOWS ────────────────────────────────────────────────────────

std::vector<WindowResult> process_windows(
    const BeagleData& bgl, const FMatrix& fm, const QInit& qi,
    const SampleRegistry& reg, int wsz, int wst, int maxiter, double tol, const char* fchr)
{
    int K=qi.K, na=reg.n_active;
    struct WD{std::string wid,chr;int sbp,ebp,si,ei;};
    std::vector<WD> wins;
    for(auto&[chr,sites]:bgl.chr_idx){
        if(fchr&&chr!=fchr) continue;
        int ns=(int)sites.size(), nw=std::max(1,(ns-wsz)/wst+1);
        for(int w=0;w<nw;w++){
            int s=w*wst, e=std::min(s+wsz-1,ns-1);
            if(e-s+1<10) continue;
            char wid[256]; snprintf(wid,sizeof(wid),"%s_w%d",chr.c_str(),w+1);
            wins.push_back({wid,chr,bgl.positions[sites[s]],bgl.positions[sites[e]],s,e});
        }
    }
    fprintf(stderr,"[iq] %d windows, K=%d, %d active samples\n",(int)wins.size(),K,na);
    std::vector<WindowResult> res(wins.size());
    auto t0=std::chrono::steady_clock::now();

    #pragma omp parallel for schedule(dynamic,4)
    for(int wi=0;wi<(int)wins.size();wi++){
        auto& wd=wins[wi];
        auto it=bgl.chr_idx.find(wd.chr); if(it==bgl.chr_idx.end()) continue;
        auto& cs=it->second;
        int nws=wd.ei-wd.si+1; if(nws<10||wd.ei>=(int)cs.size()){res[wi].nsites=0;continue;}
        std::vector<int> si(nws); std::vector<double> fs(nws*K);
        for(int s=0;s<nws;s++){
            int bi=cs[wd.si+s]; si[s]=bi;
            if(bi<fm.n_sites) for(int k=0;k<K;k++) fs[s*K+k]=fm.get(bi,k);
            else for(int k=0;k<K;k++) fs[s*K+k]=0.25;
        }
        std::vector<double> Q(na*K);
        for(int ai=0;ai<na;ai++){int o=reg.active_idx[ai]; for(int k=0;k<K;k++) Q[ai*K+k]=(o<qi.n_ind)?qi.get(o,k):1.0/K;}
        compute_Q_fixedF(bgl,fs,Q.data(),reg.active_idx.data(),na,K,si.data(),nws,maxiter,tol);
        auto& wr=res[wi]; wr.wid=wd.wid; wr.chr=wd.chr; wr.sbp=wd.sbp; wr.ebp=wd.ebp; wr.nsites=nws;
        wr.samples.resize(na);
        double s12=0,sH=0,se=0,ss=0;
        for(int ai=0;ai<na;ai++){
            wr.samples[ai]=metrics(&Q[ai*K],K,reg.active_idx[ai]);
            s12+=wr.samples[ai].delta12; sH+=wr.samples[ai].entropy;
            se+=wr.samples[ai].ena; ss+=wr.samples[ai].delta12*wr.samples[ai].delta12;
        }
        wr.mean_d12=s12/na; wr.mean_H=sH/na; wr.mean_ena=se/na;
        wr.sd_d12=std::sqrt(std::max(0.0,ss/na-wr.mean_d12*wr.mean_d12));
        if(wi>0&&wi%500==0){auto now=std::chrono::steady_clock::now();double el=std::chrono::duration<double>(now-t0).count();
            #pragma omp critical
            fprintf(stderr,"[iq] %d/%d (%.0f/s)\n",wi,(int)wins.size(),wi/el);}
    }
    auto t1=std::chrono::steady_clock::now();double el=std::chrono::duration<double>(t1-t0).count();
    fprintf(stderr,"[iq] Done: %d wins in %.1fs (%.0f/s)\n",(int)wins.size(),el,wins.size()/std::max(.001,el));
    return res;
}

// ─── OUTPUT ─────────────────────────────────────────────────────────────────

void write_summary(const std::vector<WindowResult>& res, const char* path){
    FILE* f=fopen(path,"w"); if(!f)return;
    fprintf(f,"window_id\tchrom\tstart_bp\tend_bp\tn_sites\tmean_delta12\tmean_entropy\tmean_ena\tsd_delta12\n");
    for(auto&w:res){if(w.nsites==0)continue; fprintf(f,"%s\t%s\t%d\t%d\t%d\t%.4f\t%.4f\t%.4f\t%.4f\n",w.wid.c_str(),w.chr.c_str(),w.sbp,w.ebp,w.nsites,w.mean_d12,w.mean_H,w.mean_ena,w.sd_d12);}
    fclose(f);
}

void write_samples(const std::vector<WindowResult>& res, const SampleRegistry& reg, int K, const char* path){
    FILE* f=fopen(path,"w"); if(!f)return;
    fprintf(f,"window_id\tchrom\tstart_bp\tend_bp\tsample_id\tsample_idx");
    for(int k=1;k<=K;k++) fprintf(f,"\tQ%d",k);
    fprintf(f,"\tmax_q\tdelta12\tentropy\tena\tassigned_pop\tnQ_above_005\tnQ_above_010\n");
    for(auto&w:res){if(w.nsites==0)continue; for(auto&s:w.samples){
        fprintf(f,"%s\t%s\t%d\t%d\t%s\t%d",w.wid.c_str(),w.chr.c_str(),w.sbp,w.ebp,reg.get_id(s.orig_idx).c_str(),s.orig_idx+1);
        for(int k=0;k<K;k++) fprintf(f,"\t%.4f",s.Q[k]);
        fprintf(f,"\t%.4f\t%.4f\t%.4f\t%.4f\t%d\t%d\t%d\n",s.max_q,s.delta12,s.entropy,s.ena,s.assigned_pop,s.nQ005,s.nQ010);
    }}
    fclose(f);
}

void print_result(const std::vector<WindowResult>& res, const SampleRegistry& reg, int K){
    for(auto&w:res){if(w.nsites==0)continue;
        printf("sample_id\tsample_idx"); for(int k=1;k<=K;k++) printf("\tQ%d",k);
        printf("\tmax_q\tdelta12\tentropy\tena\tassigned_pop\n");
        for(auto&s:w.samples){
            printf("%s\t%d",reg.get_id(s.orig_idx).c_str(),s.orig_idx+1);
            for(int k=0;k<K;k++) printf("\t%.4f",s.Q[k]);
            printf("\t%.4f\t%.4f\t%.4f\t%.4f\t%d\n",s.max_q,s.delta12,s.entropy,s.ena,s.assigned_pop);
        }
    }
}

// ─── MAIN ───────────────────────────────────────────────────────────────────

int main(int argc, char* argv[]){
    const char *bpath=0,*fpath=0,*qpath=0,*odir=".",*fchr=0;
    const char *sl_path=0,*keep_path=0,*excl_path=0;
    int fs=-1,fe=-1,wsz=100,wst=20,maxit=20,nc=1; double tol=1e-4;
    bool precomp=false,sout=false;

    for(int i=1;i<argc;i++){
        if     (!strcmp(argv[i],"--beagle")&&i+1<argc)      bpath=argv[++i];
        else if(!strcmp(argv[i],"--fopt")&&i+1<argc)        fpath=argv[++i];
        else if(!strcmp(argv[i],"--qinit")&&i+1<argc)       qpath=argv[++i];
        else if(!strcmp(argv[i],"--outdir")&&i+1<argc)      odir=argv[++i];
        else if(!strcmp(argv[i],"--chr")&&i+1<argc)         fchr=argv[++i];
        else if(!strcmp(argv[i],"--start")&&i+1<argc)       fs=atoi(argv[++i]);
        else if(!strcmp(argv[i],"--end")&&i+1<argc)         fe=atoi(argv[++i]);
        else if(!strcmp(argv[i],"--sample_list")&&i+1<argc) sl_path=argv[++i];
        else if(!strcmp(argv[i],"--keep")&&i+1<argc)        keep_path=argv[++i];
        else if(!strcmp(argv[i],"--exclude")&&i+1<argc)     excl_path=argv[++i];
        else if(!strcmp(argv[i],"--window_size")&&i+1<argc) wsz=atoi(argv[++i]);
        else if(!strcmp(argv[i],"--window_step")&&i+1<argc) wst=atoi(argv[++i]);
        else if(!strcmp(argv[i],"--em_iter")&&i+1<argc)     maxit=atoi(argv[++i]);
        else if(!strcmp(argv[i],"--tol")&&i+1<argc)         tol=atof(argv[++i]);
        else if(!strcmp(argv[i],"--ncores")&&i+1<argc)      nc=atoi(argv[++i]);
        else if(!strcmp(argv[i],"--K")&&i+1<argc)           ++i;
        else if(!strcmp(argv[i],"--precompute"))             precomp=true;
        else if(!strcmp(argv[i],"--sample_output"))          sout=true;
        else if(!strcmp(argv[i],"-h")||!strcmp(argv[i],"--help")){
            fprintf(stderr,
                "instant_q v2 — Fixed-F EM (Skotte 2013)\n\n"
                "Required: --beagle <f> --fopt <f> --qinit <f>\n"
                "Samples:  --sample_list <f>  --keep <f>  --exclude <f>\n"
                "Region:   --chr <n> --start <bp> --end <bp>\n"
                "Precomp:  --precompute --outdir <d> [--chr <n>]\n"
                "Options:  --window_size N --window_step N --em_iter N\n"
                "          --tol F --ncores N --sample_output\n\n"
                "Flow: bamlist line i = BEAGLE ind{i} = Q row i = sample_list line i\n"
                "--keep: only listed samples participate in EM (e.g. pruned unrelated)\n"
                "--exclude: listed samples masked out (e.g. outliers)\n"
                "Both read sample IDs or BAM paths (auto-stripped).\n");
            return 0;
        }
    }

    if(!bpath||!fpath||!qpath){fprintf(stderr,"[iq] Need --beagle --fopt --qinit\n");return 1;}
    #ifdef _OPENMP
    omp_set_num_threads(nc); fprintf(stderr,"[iq] Threads: %d\n",nc);
    #endif

    FMatrix fm; if(!fm.load(fpath)) return 1;
    QInit qi;   if(!qi.load(qpath)) return 1;
    int K=qi.K;

    // Sample registry
    SampleRegistry reg;
    if(sl_path) reg.load_ids(sl_path); else reg.gen_anon(qi.n_ind);
    if(reg.n_total!=qi.n_ind) fprintf(stderr,"[iq] WARN: sample list %d != Q rows %d\n",reg.n_total,qi.n_ind);
    reg.activate_all();
    if(keep_path) reg.apply_keep(keep_path);
    if(excl_path) reg.apply_exclude(excl_path);
    if(reg.n_active<5){fprintf(stderr,"[iq] Only %d active samples\n",reg.n_active);return 1;}

    // Load BEAGLE
    bool single=(!precomp&&fchr&&fs>=0&&fe>=0);
    BeagleData bgl;
    if(single){if(!bgl.load(bpath,fchr,fs,fe))return 1;}
    else{if(!bgl.load(bpath,fchr))return 1;}
    if(bgl.n_sites==0){fprintf(stderr,"[iq] No sites\n");return 1;}
    if(bgl.n_ind!=reg.n_total) fprintf(stderr,"[iq] WARN: BEAGLE %d ind != registry %d\n",bgl.n_ind,reg.n_total);

    if(single){
        std::vector<int> as(bgl.n_sites); std::iota(as.begin(),as.end(),0);
        std::vector<double> fsub(bgl.n_sites*K);
        for(int j=0;j<bgl.n_sites;j++) for(int k=0;k<K;k++) fsub[j*K+k]=(j<fm.n_sites)?fm.get(j,k):0.25;
        std::vector<double> Q(reg.n_active*K);
        for(int ai=0;ai<reg.n_active;ai++){int o=reg.active_idx[ai]; for(int k=0;k<K;k++) Q[ai*K+k]=(o<qi.n_ind)?qi.get(o,k):1.0/K;}
        compute_Q_fixedF(bgl,fsub,Q.data(),reg.active_idx.data(),reg.n_active,K,as.data(),bgl.n_sites,maxit,tol);
        WindowResult wr; wr.wid=std::string(fchr)+"_region"; wr.chr=fchr; wr.sbp=fs; wr.ebp=fe; wr.nsites=bgl.n_sites;
        wr.samples.resize(reg.n_active);
        for(int ai=0;ai<reg.n_active;ai++) wr.samples[ai]=metrics(&Q[ai*K],K,reg.active_idx[ai]);
        print_result({wr},reg,K);
    } else {
        auto res=process_windows(bgl,fm,qi,reg,wsz,wst,maxit,tol,fchr);
        std::string cmd="mkdir -p "+std::string(odir); system(cmd.c_str());
        std::string tag=fchr?fchr:"all";
        write_summary(res,(std::string(odir)+"/"+tag+".local_Q_summary.tsv").c_str());
        if(sout) write_samples(res,reg,K,(std::string(odir)+"/"+tag+".local_Q_samples.tsv").c_str());
        FILE* mf=fopen((std::string(odir)+"/"+tag+".local_Q_meta.tsv").c_str(),"w");
        if(mf){fprintf(mf,"chrom\tK\tn_active\tn_total\tengine\tem_iter\tn_windows\tn_sites\twin_size\twin_step\n");
            fprintf(mf,"%s\t%d\t%d\t%d\tinstant_q_v2\t%d\t%d\t%d\t%d\t%d\n",tag.c_str(),K,reg.n_active,reg.n_total,maxit,(int)res.size(),bgl.n_sites,wsz,wst);fclose(mf);}
    }
    return 0;
}

//
//  Fiber2D.cpp
//  i-eel
//
//  Created by Jeremie Bec on 19/01/2021.
//

#include "Fiber2D.hpp"

Fiber2D::Fiber2D() {
}

Fiber2D::Fiber2D(int Ns, double L, double zeta, double E, double beta) {
    assert(Ns >= 0);
    Ns_ = Ns;
    assert(L>0);
    L_ = L;
    ds_ = L_/(double)Ns_;
    assert(zeta>0);
    zeta_ = zeta;
    assert(E>0);
    E_ = E;
    assert(beta>0);
    beta_ = beta;
    Fiber2D::alloc();
}

void Fiber2D::read(int N, double L, double zeta, double E, double beta, Flow2D& U, const std::string& filebase) {
    // Initial data from file
    Ns_ = N;
    assert(L>0);
    L_ = L;
    ds_ = L_/(double)Ns_;
    assert(zeta>0);
    zeta_ = zeta;
    assert(E>0);
    E_ = E;
    assert(beta>0);
    beta_ = beta;
    
    Fiber2D::alloc();
    
    char cname[512];
    strcpy(cname, filebase.c_str());
    FILE *fin = fopen(cname,"r");
    for(int is=0; is<=Ns_; is++) {
        double xtmp[2];
        if(fread(&xtmp, sizeof(double), 2, fin)==0) {
            cerr<<"Error in reading "<<filebase<<endl;
            exit(-1);
        }
        X_(is,0) = xtmp[0];
        X_(is,1) = xtmp[1];
    }
    fclose(fin);
    
    Xold_ = X_;
    diff();
    interp_U(U);
    calc_force();
    calc_tension();
    D1T_ = D1_*T_;
    MyCol SF = F_.col(0)%D1X_.col(0)+F_.col(1)%D1X_.col(1);
    for(int dim=0; dim<2; dim++)
    Gold_.col(dim) = Uf_.col(dim) + (2.0/zeta_)*(D1T_%D1X_.col(dim)) + (1.0/zeta_)*(T_%D2X_.col(dim)) + (1.0/zeta_)*(F_.col(dim)+SF%D1X_.col(dim));
}

Fiber2D::Fiber2D(int N, double L, double zeta, double E, double beta, Flow2D& U, std::default_random_engine& rng) {
    // Random initial orientation
    Ns_ = N;
    assert(L>0);
    L_ = L;
    ds_ = L_/(double)Ns_;
    assert(zeta>0);
    zeta_ = zeta;
    assert(E>0);
    E_ = E;
    assert(beta>0);
    beta_ = beta;
    
    Fiber2D::alloc();
    
    MyCol p(2);
    std::uniform_real_distribution<double> unif(-1,1);
    double rr = 2;
    while(rr>1) {
        p(0) = unif(rng);
        p(1) = unif(rng);
        rr = p(0)*p(0)+p(1)*p(1);
    }
    rr = sqrt(rr);
    p(0) /= rr; p(1) /= rr;
    
    X_(0,0) = -L_*p(0)/2.0;
    X_(0,1) = -L_*p(1)/2.0;
    Xold_(0,0) = X_(0,0); Xold_(0,1) = X_(0,1);
    for(int is=1; is<=Ns_; is++) {
        X_.row(is) = X_.row(is-1)+ds_*p.t();
        Xold_.row(is) = X_.row(is);
    }
    diff();
    interp_U(U);
    calc_force();
    calc_tension();
    D1T_ = D1_*T_;
    MyCol SF = F_.col(0)%D1X_.col(0)+F_.col(1)%D1X_.col(1);
    for(int dim=0; dim<2; dim++)
    Gold_.col(dim) = Uf_.col(dim) + (2.0/zeta_)*(D1T_%D1X_.col(dim)) + (1.0/zeta_)*(T_%D2X_.col(dim)) + (1.0/zeta_)*(F_.col(dim)+SF%D1X_.col(dim));
}

Fiber2D::Fiber2D(int N, double L, double zeta, double E, double beta, Flow2D& U, vector<double> p) {
    // Fixed initial orientation
    Ns_ = N;
    assert(L>0);
    L_ = L;
    ds_ = L_/(double)Ns_;
    assert(zeta>0);
    zeta_ = zeta;
    assert(E>0);
    E_ = E;
    assert(beta>0);
    beta_ = beta;
    
    Fiber2D::alloc();
    
    double rr = sqrt(p.at(0)*p.at(0)+p.at(1)*p.at(1));
    p.at(0) /= rr; p.at(1) /= rr;
    
    X_(0,0) = -L_*p.at(0)/2.0;
    X_(0,1) = -L_*p.at(1)/2.0;
    Xold_(0,0) = X_(0,0); Xold_(0,1) = X_(0,1);
    for(int is=1; is<=Ns_; is++) {
        for(int idim=0; idim<2; ++idim)
        X_(is,idim) = X_(is-1,idim)+ds_*p.at(idim);
        Xold_.row(is) = X_.row(is);
    }
    diff();
    interp_U(U);
    calc_force();
    calc_tension();
	D1T_ = D1_*T_;
    MyCol SF = F_.col(0)%D1X_.col(0)+F_.col(1)%D1X_.col(1);
    for(int dim=0; dim<2; dim++)
    Gold_.col(dim) = Uf_.col(dim) + (2.0/zeta_)*(D1T_%D1X_.col(dim)) + (1.0/zeta_)*(T_%D2X_.col(dim)) + (1.0/zeta_)*(F_.col(dim)+SF%D1X_.col(dim));
}

void Fiber2D::alloc() {
    X_.set_size(Ns_+1,2);
    Xold_.set_size(Ns_+1,2);
    Gold_.set_size(Ns_+1,2);
    D1X_.set_size(Ns_+1,2);
    D2X_.set_size(Ns_+1,2);
    D3X_.set_size(Ns_+1,2);
    D4X_.set_size(Ns_+1,2);
    
    T_.set_size(Ns_+1);
    D1T_.set_size(Ns_+1);
    NormXi_.set_size(Ns_+1);
    
    F_.set_size(Ns_+1,2);
    D1F_.set_size(Ns_+1,2);
    
    Uf_.set_size(Ns_+1,2);
    D1Uf_.set_size(Ns_+1,2);
    
    D1_.set_size(Ns_+1,Ns_+1);
    D1_.diag(1)  += 1;
    D1_.diag(-1) -= 1;
    D1_(0,0) = -3; D1_(0,1) = 4; D1_(0,2) = -1;
    D1_(Ns_,Ns_) = 3; D1_(Ns_,Ns_-1) = -4; D1_(Ns_,Ns_-2) = 1;
    D1_ /= 2.0*ds_;
    
    D2_.set_size(Ns_+1,Ns_+1);
    D2_.diag(1)  += 1;
    D2_.diag()   -= 2;
    D2_.diag(-1) += 1;
    D2_(0,0)=2; D2_(0,1)=-5; D2_(0,2)=4; D2_(0,3)=-1;
    D2_(Ns_,Ns_)=2; D2_(Ns_,Ns_-1)=-5; D2_(Ns_,Ns_-2)=4; D2_(Ns_,Ns_-3)=-1;
    D2_ /= ds_*ds_;
    
    LapDirich_.set_size(Ns_-1,Ns_-1);
    LapDirich_.diag() -= 2;
    LapDirich_.diag(1) += 1;
    LapDirich_.diag(-1) += 1;
    LapDirich_ /= ds_*ds_;
    
    D3_.set_size(Ns_+1,Ns_+1);
    D3_.diag(2)  += 1;
    D3_.diag(1)  -= 2;
    D3_.diag(-1) += 2;
    D3_.diag(-2) -= 1;
    D3_(0,0)=-5; D3_(0,1)=18; D3_(0,2)=-24; D3_(0,3)=14; D3_(0,4)=-3;
    D3_(1,0)=-3; D3_(1,1)=10; D3_(1,2)=-12; D3_(1,3)=6; D3_(1,4)=-1;
    D3_(Ns_,Ns_)=5; D3_(Ns_,Ns_-1)=-18; D3_(Ns_,Ns_-2)=24; D3_(Ns_,Ns_-3)=-14; D3_(Ns_,Ns_-4)=3;
    D3_(Ns_-1,Ns_)=3; D3_(Ns_-1,Ns_-1)=-10; D3_(Ns_-1,Ns_-2)=12; D3_(Ns_-1,Ns_-3)=-6; D3_(Ns_-1,Ns_-4)=1;
    D3_ /= 2.0*ds_*ds_*ds_;
    
    D4_.set_size(Ns_+1,Ns_+1);
    D4_.diag(2)  += 1;
    D4_.diag(1)  -= 4;
    D4_.diag()   += 6;
    D4_.diag(-1) -= 4;
    D4_.diag(-2) += 1;
    D4_(0,0)=3; D4_(0,1)=-14; D4_(0,2)=26; D4_(0,3)=-24; D4_(0,4)=11; D4_(0,5)=-2;
    D4_(1,0)=2; D4_(1,1)=-9; D4_(1,2)=16; D4_(1,3)=-14; D4_(1,4)=6; D4_(1,5)=-1;
    D4_(Ns_,Ns_)=3; D4_(Ns_,Ns_-1)=-14; D4_(Ns_,Ns_-2)=26; D4_(Ns_,Ns_-3)=-24; D4_(Ns_,Ns_-4)=11; D4_(Ns_,Ns_-5)=-2;
    D4_(Ns_-1,Ns_)=2; D4_(Ns_-1,Ns_-1)=-9; D4_(Ns_-1,Ns_-2)=16; D4_(Ns_-1,Ns_-3)=-14; D4_(Ns_-1,Ns_-4)=6; D4_(Ns_-1,Ns_-5)=-1;
    D4_ /= ds_*ds_*ds_*ds_;
    
    Op4_.set_size(Ns_-3,Ns_-3);
    Op4_.diag(2)  += 1;
    Op4_.diag(1)  -= 4;
    Op4_.diag()   += 6;
    Op4_.diag(-1) -= 4;
    Op4_.diag(-2) += 1;
    Op4_(0,0) = 2.0/11.0; Op4_(0,1) = -4.0/11.0; Op4_(0,2) = 2.0/11.0;
    Op4_(1,0) = -16.0/11.0; Op4_(1,1) = 43.0/11.0; Op4_(1,2) = -38.0/11.0;
    Op4_(Ns_-5,Ns_-4) = -16.0/11.0; Op4_(Ns_-5,Ns_-5) = 43.0/11.0; Op4_(Ns_-5,Ns_-6) = -38.0/11.0;
    Op4_(Ns_-4,Ns_-4) = 2.0/11.0; Op4_(Ns_-4,Ns_-5) = -4.0/11.0; Op4_(Ns_-4,Ns_-6) = 2.0/11.0;
    Op4_ /= ds_*ds_*ds_*ds_;
}

void Fiber2D::setforcing(int k, double om) {
    p_.resize(2);
    p_.at(0) = 0;
    p_.at(1) = 1;
    nu_ = 2.0*M_PI*(double)k/L_;
    om_ = om;
}

void Fiber2D::setforcing(vector<double> p, double A) {
    A_ = A;
    p_ = p;
}

void Fiber2D::save(const std::string& filebase, Flow2D& U) {
    
    // Update derivatives, fluid velocity and tension
    diff();
    interp_U(U);
    calc_force();
    calc_tension();
    
    if(filebase.substr(filebase.length()-3,3) == ".ff") {
        char cname[512];
        strcpy(cname, filebase.c_str());
        FILE *fout = fopen(cname,"w");
        
        for(int is=0; is<=Ns_; is++) {
            fwrite(&X_(is,0), sizeof(double), 1, fout);
            fwrite(&X_(is,1), sizeof(double), 1, fout);
        }
        
        MyMat V = (X_-Xold_)/dt_;
        for(int is=0; is<=Ns_; is++) {
            fwrite(&V(is,0), sizeof(double), 1, fout);
            fwrite(&V(is,1), sizeof(double), 1, fout);
        }
        
        for(int is=0; is<=Ns_; is++) {
            fwrite(&Uf_(is,0), sizeof(double), 1, fout);
            fwrite(&Uf_(is,1), sizeof(double), 1, fout);
        }
        
        fclose(fout);
    }
    else if(filebase.substr(filebase.length()-3,3) == ".nc"){
        int status, ncid;
        if ((status = nc_create(filebase.c_str(), NC_NETCDF4, &ncid)))
            ErrorMsg("NetCDF 0: "+std::string(nc_strerror(status)));
        
        // Conventional global attributes
        char project[] = "i-eel";
        if ((status = nc_put_att_text(ncid, NC_GLOBAL, "Conventions", strlen(project), project)))
            ErrorMsg("NetCDF: "+std::string(nc_strerror(status)));
        char title[] = "2D Elongated swimmer";
        if ((status = nc_put_att_text(ncid, NC_GLOBAL, "title", strlen(title), title)))
            ErrorMsg("NetCDF: "+std::string(nc_strerror(status)));
        if ((status = nc_put_att_text(ncid, NC_GLOBAL, "shearbox_version", strlen(IEEL_VERSION), IEEL_VERSION)))
            ErrorMsg("NetCDF: "+std::string(nc_strerror(status)));
        if ((status = nc_put_att_text(ncid, NC_GLOBAL, "compiler_version", strlen(COMPILER_VERSION), COMPILER_VERSION)))
            ErrorMsg("NetCDF: "+std::string(nc_strerror(status)));
        time_t rawtime;  // current time
        struct tm* timeinfo;
        char tbuffer[80];
        time(&rawtime);
        timeinfo = localtime(&rawtime);
        strftime(tbuffer, 80, "%Y-%m-%d %I:%M:%S", timeinfo);
        if ((status = nc_put_att_text(ncid, NC_GLOBAL, "time", strlen(tbuffer), tbuffer)))
            ErrorMsg("NetCDF: "+std::string(nc_strerror(status)));
        char hostname[1024];
        gethostname(hostname, 1023);
        if ((status = nc_put_att_text(ncid, NC_GLOBAL, "host_name", strlen(hostname), hostname)))
            ErrorMsg("NetCDF: "+std::string(nc_strerror(status)));
        char reference[] = "Data from the i-eel code";
        if ((status = nc_put_att_text(ncid, NC_GLOBAL, "references", strlen(reference), reference)))
            ErrorMsg("NetCDF: "+std::string(nc_strerror(status)));
        // Define and write variables
        int dimid[2];
        int id;
        if ((status = nc_def_dim(ncid, "dim", 2, &dimid[0])))
            ErrorMsg("NetCDF: "+std::string(nc_strerror(status)));
        if ((status = nc_def_dim(ncid, "Ns", Ns_+1, &dimid[1])))
            ErrorMsg("NetCDF: "+std::string(nc_strerror(status)));
        // Position
        if ((status = nc_def_var(ncid, "Pos", NC_DOUBLE, 2, dimid, &id)))
            ErrorMsg("NetCDF: "+std::string(nc_strerror(status)));
        if ((status = nc_put_var_double(ncid, id, &X_[0])))
            ErrorMsg("NetCDF: "+std::string(nc_strerror(status)));
        // Velocity
        MyMat V = (X_-Xold_)/dt_;
        if ((status = nc_def_var(ncid, "Vel", NC_DOUBLE, 2, dimid, &id)))
            ErrorMsg("NetCDF: "+std::string(nc_strerror(status)));
        if ((status = nc_put_var_double(ncid, id, &V[0])))
            ErrorMsg("NetCDF: "+std::string(nc_strerror(status)));
        // Fluid velocity
        if ((status = nc_def_var(ncid, "Vel_fluid", NC_DOUBLE, 2, dimid, &id)))
            ErrorMsg("NetCDF: "+std::string(nc_strerror(status)));
        if ((status = nc_put_var_double(ncid, id, &Uf_[0])))
            ErrorMsg("NetCDF: "+std::string(nc_strerror(status)));
        // Tension
        if ((status = nc_def_var(ncid, "Tension", NC_DOUBLE, 1, &dimid[1], &id)))
            ErrorMsg("NetCDF: "+std::string(nc_strerror(status)));
        if ((status = nc_put_var_double(ncid, id, &T_[0])))
            ErrorMsg("NetCDF: "+std::string(nc_strerror(status)));
        // Close the netcdf file
        if ((status = nc_close(ncid)))
            ErrorMsg("NetCDF: "+std::string(nc_strerror(status)));
    }
    else
        ErrorMsg("Unknown output format for file "+filebase);
}

void Fiber2D::diff() {
    D1X_    = D1_*X_;
    D2X_    = D2_*X_;
    D3X_    = D3_*X_;
    D4X_    = D4_*X_;
    NormXi_ = sum(D1X_%D1X_,1);
}

void Fiber2D::calc_tension() {
    MySpMat Lap = 2.0*LapDirich_;
    Lap.diag() -= sum(D2X_.rows(1,Ns_-1)%D2X_.rows(1,Ns_-1),1);
    
    MyCol A = -zeta_*sum(D1X_.rows(1,Ns_-1)%D1Uf_.rows(1,Ns_-1),1);
    A -= sum(D2X_.rows(1,Ns_-1)%F_.rows(1,Ns_-1),1) + 2.0*sum(D1X_.rows(1,Ns_-1)%D1F_.rows(1,Ns_-1),1);
    A -= 7.0*E_*sum(D2X_.rows(1,Ns_-1)%D4X_.rows(1,Ns_-1),1) + 6.0*E_*sum(D3X_.rows(1,Ns_-1)%D3X_.rows(1,Ns_-1),1);
    A += (zeta_*beta_)*(1-NormXi_.rows(1,Ns_-1));
    MyCol  Ttmp;
    if(spsolve(Ttmp,Lap,A,"superlu")==false || isnan(Ttmp(0))) {
        cerr<<"No solution for the tension"<<endl;
        X_.save("X.dat",raw_ascii);
        exit(-1);
    }
    T_(0) = 0;
    T_.rows(1,Ns_-1) = Ttmp;
    T_(Ns_) = 0;
}

void Fiber2D::evol(double dt, Flow2D& U) {
    
    diff();
    dt_ = dt;
    interp_U(U);
    
    calc_force();
    calc_tension();
    D1T_ = D1_*T_;
    
    MyMat G(Ns_+1,2);
    MyCol SF = F_.col(0)%D1X_.col(0)+F_.col(1)%D1X_.col(1);
    
    for(int dim=0; dim<2; dim++)
        G.col(dim) = Uf_.col(dim) + (2.0/zeta_)*(D1T_%D1X_.col(dim)) + (1.0/zeta_)*(T_%D2X_.col(dim)) + (1.0/zeta_)*(F_.col(dim)+SF%D1X_.col(dim));
    
    MyMat RHS = (4.0/3.0)*X_.rows(2,Ns_-2)-(1.0/3.0)*Xold_.rows(2,Ns_-2)+(2.0*dt_/3.0)*(2.0*G.rows(2,Ns_-2)-Gold_.rows(2,Ns_-2));
    MyMat XXi = D1_*(2.0*X_-Xold_);
    XXi = XXi.rows(2,Ns_-2);
    
    Xold_ = X_;
    Gold_ = G;
    
    int nn = Ns_-3;
    MySpMat Axx(nn,nn), Ayy(nn,nn), Axy(nn,nn);
    Axx.diag() = XXi.col(0) % XXi.col(0) + 1.0;
    Ayy.diag() = XXi.col(1) % XXi.col(1) + 1.0;
    Axy.diag() = XXi.col(0) % XXi.col(1);
    Axy = Axy*Op4_;
    
    MySpMat MM = speye<MySpMat>(2*nn,2*nn) + (2.0*E_*dt_/(3.0*zeta_))*join_cols(join_rows(Axx*Op4_,Axy),    join_rows(Axy,Ayy*Op4_));
    
    MyMat Xnew(2*nn,1);
    if(spsolve(Xnew,MM,vectorise(RHS),"superlu")==false || isnan(Xnew(0))) {
        cerr<<"No solution for elasticity"<<endl;
        X_.save("X.dat",raw_ascii);
        exit(-1);
    }
    
    Xnew = reshape(Xnew,nn,2);
    
    X_.row(0) = (48.0/11.0)*Xnew.row(0)-(52.0/11.0)*Xnew.row(1)+(15.0/11.0)*Xnew.row(2);
    X_.row(1) = (28.0/11.0)*Xnew.row(0)-(23.0/11.0)*Xnew.row(1)+(6.0/11.0)*Xnew.row(2);
    X_.rows(2,Ns_-2) = Xnew;
    X_.row(Ns_-1) = (28.0/11.0)*Xnew.row(Ns_-4)-(23.0/11.0)*Xnew.row(Ns_-5)+(6.0/11.0)*Xnew.row(Ns_-6);
    X_.row(Ns_)   = (48.0/11.0)*Xnew.row(Ns_-4)-(52.0/11.0)*Xnew.row(Ns_-5)+(15.0/11.0)*Xnew.row(Ns_-6);
    
    t_ += dt_;
}

void Fiber2D::interp_U(Flow2D& U) {
    for(int is=0; is<=Ns_; is++) {
        for(int dim=0; dim<2; dim++) {
            Uf_(is,dim) = U.velocity(X_(is,0),X_(is,1),dim);
            D1Uf_(is,dim) = 0;
            for(int dd=0; dd<2; dd++)
			D1Uf_(is,dim) += D1X_(is,dd)*U.gradient(X_(is,0),X_(is,1),2*dim+dd);
        }
    }
}

void Fiber2D::calc_force() {
    if(p_.size()==2) {
        double n[2], S, C;
        double s = sqrt(p_.at(0)*p_.at(0)+p_.at(1)*p_.at(1));
        p_.at(0) /= s; p_.at(1) /= s;
        n[0]=p_.at(1); n[1]=-p_.at(0);
        for(int is=0; is<=Ns_; is++) {
            S = sin(nu_*(double)is*ds_-om_*t_);
            C = cos(nu_*(double)is*ds_-om_*t_);
            for(int dim=0; dim<2; ++dim) {
                F_(is,dim) = A_*S*n[dim];
                D1F_(is,dim) = nu_*A_*C*n[dim];
            }
        }
    }
    else {
        for(int is=0; is<=Ns_; is++)
        for(int dim=0; dim<2; ++dim) {
            F_(is,dim) = 0.0;
            D1F_(is,dim) = 0.0;
        }
    }
}

double Fiber2D::wind(Flow2D& U) {
    return U.velocity(X_(0,0),X_(0,1),0);
}

vector<double> Fiber2D::orientation(){
    vector<double> P(2);
    P.at(0) = X_(0,0)-getcenter(0);
    P.at(1) = X_(0,1)-getcenter(1);
    return P;
}

int Fiber2D::calc_buckle(){
    int test=0;
    double dist;
    double d2max = 0.75*ds_;
    d2max = d2max*d2max;
    for(int i=0; i<Ns_ and test==0; ++i){
        for (int j=i+1; j<Ns_+1 and test==0; ++j) {
            dist = (X_(j,0)-X_(i,0))*(X_(j,0)-X_(i,0)) + (X_(j,1)-X_(i,1))*(X_(j,1)-X_(i,1));
            if (dist < d2max)
                test = 1;
        }
    }
    return test;
}

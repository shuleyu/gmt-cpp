#include<iostream>
#include<fstream>
#include<iterator>
#include<vector>
#include<cstdlib>
#include<cstring>
#include<chrono>
#include<random>
#include<algorithm>
#include<cmath>

extern "C" {
// Current version: GMT-5.4.4
#include<gmt/gmt.h>
}


/***********************************************************
 * This C++ template convert input longitude (in deg) to
 * [0,360] deg.
 *
 * input(s):
 * const double &lon  ----  Longitude.
 *
 * return(s):
 * double ans  ----  Longitude in 0 ~ 360.
 *
 * Shule Yu
 * Dec 29 2017
 *
 * Key words: coordinates, longitude.
***********************************************************/

inline double Lon2360(const double &lon){

    double ans=lon;
    if (ans>=0) ans-=360.0*((int)(ans/360));
    else ans+=360.0*(1+(int)(-ans/360));

    if (ans>=360) ans=0;
    return ans;
}


/***********************************************************
 * This C++ template returns the 1D grid meshing results.
 *
 * input(s):
 * const double &lowerbound   ----  Grid lower bound.
 * const double &upperbound   ----  Grid upper bound.
 * const double &Para         ----  meaning-changing parameter.
 *                                  it's meaning depends on which mode:
 *
 * const int    &mode         ----  select mode.
 *                                  Define the meaning of Para and return vector.
 *              
 *                                  0: Para means number of gird points.
 *                                     Return a vector (ans) of size Para:
 *              
 *                                     ans.size()=Para;
 *                                     ans[0]=lowerbound;
 *                                     ans[Para-1]=upperbound;
 *              
 *                                  1: Para means grid increment.
 *                                     Last grid point is less equal to higherbound.
 *                                     Return a vector (ans) of calculated size:
 *              
 *                                     ans[0]=lowerbound;
 *                                     ans[1]-ans[0]=Para;
 *                                     upperbound-Para<ans.back();
 *                                     ans.back()<=upperbound;
 *              
 *                                 -1: Same as 1. Will only return the grid property:
 *              
 *                                     ans.size()=2;
 *                                     ans[0] = calculated grid size.
 *                                     ans[1] = adjusted upperbound.
 *              
 *              
 *                                  2: Para means an estimation of grid increment.
 *                                     The calculated grid increment (Para*) is (possibly)
 *                                     sightly decreased such that the higherbound is meet.
 *                                     Return a vector (ans) of calculated size:
 *              
 *                                     ans[0]=lowerbound;
 *                                     ans[1]-ans[0]=Para* (<=Para);
 *                                     ans.back()=upperbound;
 *              
 *                                 -2: Same as 2. Will only return the grid property:
 *              
 *                                     ans.size()=2;
 *                                     ans[0] = calculated grid size.
 *                                     ans[1] = adjusted Para (Para*).
 *
 * return(s):
 * vector<double> ans  ----  Created grid or grid properties, depending on mode.
 *
 * Shule Yu
 * Jan 23 2018
 *
 * Key words: creating grid.
***********************************************************/

std::vector<double> CreateGrid(const double &lowerbound, const double &upperbound,
                               const double &Para, const int mode){

    // check lower and upper bound.
    if (upperbound<lowerbound) {
        std::cerr <<  "Error in " << __func__ << ": lower bound > upper bound ..." << std::endl;
        return {};
    }

    if (mode==0){

        int N=Para;
        double Inc=1.0*(upperbound-lowerbound)/(N-1);

        std::vector<double> ans(N,lowerbound);
        for (int i=1;i<N;++i) ans[i]=ans[i-1]+Inc;
        return ans;
    }
    if (mode==1 || mode==-1){

        double Inc=Para;

        if (mode==-1) {
            double N=1+floor(1.0*(upperbound-lowerbound)/Inc);
            if (fabs((upperbound-lowerbound-N*Inc))<Inc*1e-6) // float number rounding error?
                N+=1;
            return {N,lowerbound+(N-1)*Inc};
        }

        std::vector<double> ans;

        double Cur=lowerbound;

        while (Cur<=upperbound || Cur-upperbound<Inc*1e-6) {
            ans.push_back(Cur);
            Cur+=Inc;
        }
        return ans;
    }
    if (mode==2 || mode==-2){

        double Inc=Para;
        int N=1+(int)ceil((upperbound-lowerbound)/Inc);
        Inc=1.0*(upperbound-lowerbound)/(N-1);

        if (mode==-2) {
            return {1.0*N,Inc};
        }

        std::vector<double> ans;
        double Cur=lowerbound;

        for (int i=0;i<N;++i) {
            ans.push_back(Cur);
            Cur+=Inc;
        }
        return ans;
    }

    std::cerr <<  "Error in " << __func__ << ": mode error ..." << std::endl;
    return {};
}


/*************************************************
 * This is a wrapper for GMT API.
 *
 *
 * Shule Yu
 * Aug 13 2018
 *
 * Key words: gmt plotting c++
*************************************************/

namespace GMT { // the order of the function definition matters: dependencies should appear first.

    // basic operations (no dependencies needed)

    // gmt set.
    void set(const std::string &cmd){
        void *API=GMT_Create_Session(__func__,2,0,NULL);
        char *command=strdup(cmd.c_str());
        GMT_Call_Module(API,"set",GMT_MODULE_CMD,command);

        delete [] command;
        GMT_Destroy_Session(API);
    }


    // gmt pscoast.
    void pscoast(const std::string &outfile, const std::string &cmd){
        void *API=GMT_Create_Session(__func__,2,0,NULL);
        char *command=strdup((cmd+" ->>"+outfile).c_str());
        GMT_Call_Module(API,"pscoast",GMT_MODULE_CMD,command);

        delete [] command;
        GMT_Destroy_Session(API);
    }


    // gmt makecpt.
    void makecpt(const std::string &cmd){
        void *API=GMT_Create_Session(__func__,2,0,NULL);
        char *command=strdup(cmd.c_str());
        GMT_Call_Module(API,"makecpt",GMT_MODULE_CMD,command);

        delete [] command;
        GMT_Destroy_Session(API);
    }


    // gmt psscale.
    void psscale(const std::string &outfile, const std::string &cmd){
        void *API=GMT_Create_Session(__func__,2,0,NULL);
        char *command=strdup((cmd+" ->>"+outfile).c_str());
        GMT_Call_Module(API,"psscale",GMT_MODULE_CMD,command);

        delete [] command;
        GMT_Destroy_Session(API);
    }


    // gmt psbasemap.
    void psbasemap(const std::string &outfile, const std::string &cmd){
        void *API=GMT_Create_Session(__func__,2,0,NULL);
        char *command=strdup((cmd+" ->>"+outfile).c_str());
        GMT_Call_Module(API,"psbasemap",GMT_MODULE_CMD,command);

        delete [] command;
        GMT_Destroy_Session(API);
    }


    // gmt pstext.
    class Text {
        public:
        std::string Justify,Font,Content;
        double X,Y;
        int Size;
        Text(const double &x, const double &y, const std::string &c) : Text(x,y,c,12) {}
        Text(const double &x, const double &y, const std::string &c, const int &s) :
             Text(x,y,c,s,"CB") {}
        Text(const double &x, const double &y, const std::string &c, const int &s,
             const std::string &j) : Text(x,y,c,s,j,"Helvetica") {}
        Text(const double &x, const double &y, const std::string &c, const int &s,
             const std::string &j, const std::string &f) {
                 X=x;Y=y;Content=c;Size=s;Justify=j;Font=f;
             }
    };
    void pstext(const std::string &outfile, const std::vector<GMT::Text> &texts,
                const std::string &cmd){

        std::size_t n=texts.size();
        if (n==0) return;

        void *API=GMT_Create_Session(__func__,2,0,NULL);

        uint64_t par[]={1,1,1};
        GMT_TEXTSET *txt=(GMT_TEXTSET *)GMT_Create_Data(API,GMT_IS_TEXTSET,GMT_IS_NONE,
                                                        0,par,NULL,NULL,0,-1,NULL);

        char filename[20];
        for (const auto &item:texts) {
            std::string s=std::to_string(item.X)+" "+std::to_string(item.Y)+" "+item.Content;

            txt->table[0]->segment[0]->n_rows = 1;
            txt->table[0]->segment[0]->data[0] = new char [s.size()+1];
            strcpy(txt->table[0]->segment[0]->data[0],s.c_str());

            GMT_Open_VirtualFile(API,GMT_IS_TEXTSET,GMT_IS_NONE,GMT_IN,txt,filename);
            char *command=strdup(("-<"+std::string(filename)+" "+cmd+" -F+j"+item.Justify+"+f"+
                                  std::to_string(item.Size)+"p,"+item.Font+" ->>"+outfile).c_str());
            GMT_Call_Module(API,"pstext",GMT_MODULE_CMD,command);

            delete [] command;
            delete [] txt->table[0]->segment[0]->data[0];
            txt->table[0]->segment[0]->data[0]=nullptr;
            GMT_Close_VirtualFile(API,filename);
        }

        GMT_Destroy_Data(API,txt);
        GMT_Destroy_Session(API);
        return;
    }


    // gmt psxy.
    template <typename T>
    void psxy(const std::string &outfile, const T XBegin, const T XEnd,
              const T YBegin, const T YEnd, const std::string &cmd){

        // Check array size.
        std::size_t n=std::distance(XBegin,XEnd),m=std::distance(YBegin,YEnd);
        if (n==0) return;

        if (m!=n) throw std::runtime_error("In "+std::string(__func__)+", input x,y size don't match.");

        void *API=GMT_Create_Session(__func__,2,0,NULL);


        // Set vector dimensions..
        uint64_t par[2];
        par[0]=2;  // x,y value (2 columns).
        par[1]=n;  // npts (n rows).

        // Create plot data.
        GMT_VECTOR *vec=(GMT_VECTOR *)GMT_Create_Data(API,GMT_IS_VECTOR,GMT_IS_POINT,
                                                      0,par,NULL,NULL,0,-1,NULL);

        // Inject data.
        double *X = new double [n], *Y = new double [n];
        std::size_t i=0;
        for (auto it=XBegin;it!=XEnd;++it) X[i++]=*it;
        i=0;
        for (auto it=YBegin;it!=YEnd;++it) Y[i++]=*it;

        GMT_Put_Vector(API,vec,0,GMT_DOUBLE,X);
        GMT_Put_Vector(API,vec,1,GMT_DOUBLE,Y);

        // Get the virtual file.
        char filename[20];
        GMT_Open_VirtualFile(API,GMT_IS_VECTOR,GMT_IS_POINT,GMT_IN,vec,filename);

        // Plot.
        char *command=strdup(("-<"+std::string(filename)+" "+cmd+" ->>"+outfile).c_str());
        GMT_Call_Module(API,"psxy",GMT_MODULE_CMD,command);

        // free spaces (X,Y seem to be deleted by GMT_Destroy_Session? very confusing).
        delete [] command;
        GMT_Close_VirtualFile(API,filename);
        GMT_Destroy_Data(API,vec);
        GMT_Destroy_Session(API);

        return;
    }

    // psxy with multiple coloms.
    template <typename T>
    void psxy(const std::string &outfile, const std::vector<std::vector<T>> &Data,
              const std::string &cmd){


        // Check array size.
        if (Data.empty()) return;
        std::size_t n=Data.size(),m=Data[0].size(); // m rows, n columns
        for (std::size_t i=0;i<n;++i)
            if (Data[i].size()!=m)
                throw std::runtime_error("In "+std::string(__func__)
                                         +", input data column size don't match.");

        void *API=GMT_Create_Session(__func__,2,0,NULL);

        // Set vector dimensions..
        uint64_t par[2];
        par[0]=n;  // x,y,z,... value (n columns).
        par[1]=m;  // (m of rows).

        // Create plot data.
        GMT_VECTOR *vec=(GMT_VECTOR *)GMT_Create_Data(API,GMT_IS_VECTOR,GMT_IS_POINT,
                                                      0,par,NULL,NULL,0,-1,NULL);

        // Inject data.
        for (std::size_t i=0;i<n;++i) {
            double *X = (double *)malloc(m*sizeof(double));
            for (std::size_t j=0;j<m;++j)
                X[j]=Data[i][j];
            GMT_Put_Vector(API,vec,i,GMT_DOUBLE,X);
        }

        // Get the virtual file.
        char filename[20];
        GMT_Open_VirtualFile(API,GMT_IS_VECTOR,GMT_IS_POINT,GMT_IN,vec,filename);

        // Plot.
        char *command=strdup(("-<"+std::string(filename)+" "+cmd+" ->>"+outfile).c_str());
        GMT_Call_Module(API,"psxy",GMT_MODULE_CMD,command);

        // free spaces (X's seem to be deleted by GMT_Destroy_Session? very confusing).
        delete [] command;
        GMT_Close_VirtualFile(API,filename);
        GMT_Destroy_Data(API,vec);
        GMT_Destroy_Session(API);

        return;
    }

    // gmt pshistogram.
    template <typename T>
    void pshistogram(const std::string &outfile, const T XBegin, const T XEnd,
                     const std::string &cmd){

        // Check array size.
        std::size_t n=std::distance(XBegin,XEnd);
        if (n==0) return;

        void *API=GMT_Create_Session(__func__,2,0,NULL);


        // Set vector dimensions..
        uint64_t par[2];
        par[0]=1;  // value (1 columns).
        par[1]=n;  // npts (n rows).

        // Create plot data.
        GMT_VECTOR *vec=(GMT_VECTOR *)GMT_Create_Data(API,GMT_IS_VECTOR,GMT_IS_POINT,
                                                      0,par,NULL,NULL,0,-1,NULL);

        // Inject data.
        double *X = new double [n];
        std::size_t i=0;
        for (auto it=XBegin;it!=XEnd;++it) X[i++]=*it;

        GMT_Put_Vector(API,vec,0,GMT_DOUBLE,X);

        // Get the virtual file.
        char filename[20];
        GMT_Open_VirtualFile(API,GMT_IS_VECTOR,GMT_IS_POINT,GMT_IN,vec,filename);

        // Plot.
        char *command=strdup(("-<"+std::string(filename)+" "+cmd+" ->>"+outfile).c_str());
        GMT_Call_Module(API,"pshistogram",GMT_MODULE_CMD,command);

        // free spaces.
        delete [] command;
        GMT_Close_VirtualFile(API,filename);
        GMT_Destroy_Data(API,vec);
        GMT_Destroy_Session(API);

        return;
    }

    // gmt grdimage.

    template<typename T>
    void grdimage(const std::string &outfile, const std::vector<std::vector<T>> &G,
                  const double &xinc, const double &yinc, const std::string &cmd){

        if (G.empty()) return;

        // Check array size.
        for (const auto &item:G)
            if (item.size()!=3)
                throw std::runtime_error("In "+std::string(__func__)+", input is not x,y,z data");

        void *API=GMT_Create_Session(__func__,2,0,NULL);

        // Set grid limits.
        double MinVal=std::numeric_limits<double>::max(),MaxVal=-MinVal;
        double wesn[]={MinVal,MaxVal,MinVal,MaxVal};
        for (const auto &item:G) {
            wesn[0]=std::min(wesn[0],item[0]);
            wesn[1]=std::max(wesn[1],item[0]);
            wesn[2]=std::min(wesn[2],item[1]);
            wesn[3]=std::max(wesn[3],item[1]);
            MinVal=std::min(MinVal,item[2]);
            MaxVal=std::max(MaxVal,item[2]);
        }


        // Set grid increments.
        double inc[]={xinc,yinc};

        // Set up grid size.
        auto res=CreateGrid(wesn[0],wesn[1],inc[0],-1);
        std::size_t m=(std::size_t)res[0]+4;
        res=CreateGrid(wesn[2],wesn[3],inc[1],-1);
        std::size_t n=(std::size_t)res[0]+4;


        // Create plot data.
        GMT_GRID *grid=(GMT_GRID *)GMT_Create_Data(API,GMT_IS_GRID,GMT_IS_SURFACE,
                                                   GMT_CONTAINER_ONLY,NULL,wesn,inc,
                                                   GMT_GRID_NODE_REG,-1,NULL);

        // Inject data.
        float *aux_data = new float [m*n];
        for (std::size_t i=0;i<m*n;++i) aux_data[i]=0.0/0.0;
        for (const auto &item:G) {
            std::size_t X=(std::size_t)round((item[0]-wesn[0])/xinc);
            std::size_t Y=(std::size_t)round((item[1]-wesn[2])/yinc);
            // swap X,Y position and flip along y-axis. 
            aux_data[(n-3-Y)*m+X+2]=item[2];
        }
        grid->data=aux_data;


        // Adjust something in the header. (magic)
        grid->header->z_min=MinVal;
        grid->header->z_max=MaxVal;

        grid->header->grdtype=3;
        grid->header->gn=1;
        grid->header->gs=1;
        grid->header->BC[2]=3;
        grid->header->BC[3]=3;

        // periodic along longitude (x) direction.
        grid->header->BC[0]=2;
        grid->header->BC[1]=2;
        for (std::size_t Y=0;Y<n;++Y) {
            aux_data[(n-1-Y)*m+0]=aux_data[(n-1-Y)*m+m-4];
            aux_data[(n-1-Y)*m+1]=aux_data[(n-1-Y)*m+m-3];
            aux_data[(n-1-Y)*m+m-1]=aux_data[(n-1-Y)*m+3];
            aux_data[(n-1-Y)*m+m-2]=aux_data[(n-1-Y)*m+2];
        }

        // Get the virtual file.
        char filename[20];
        GMT_Open_VirtualFile(API,GMT_IS_GRID,GMT_IS_SURFACE,GMT_IN,grid,filename);

        // Plot.
        char *command=strdup(("-<"+std::string(filename)+" "+cmd+" ->>"+outfile).c_str());
        GMT_Call_Module(API,"grdimage",GMT_MODULE_CMD,command);

        // Free spaces.
        delete [] command;
        GMT_Close_VirtualFile(API,filename);
        GMT_Destroy_Data(API,grid);
        GMT_Destroy_Session(API);

        return;
    }

    // Functions(templates) dependent on others.

    // Move reference point.
    void MoveReferencePoint(const std::string &outfile, const std::string &cmd){
        GMT::psbasemap(outfile,"-J -R -Bwens -O -K "+cmd);
        return;
    }

    // Plot a time stamp at left bottom corner.
    void timestamp(const std::string &outfile){
        auto it=std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
        std::string date(ctime(&it));
        std::vector<GMT::Text> texts;
        texts.push_back(GMT::Text(0,0,date,4,"LB"));
        GMT::MoveReferencePoint(outfile,"-Xf0i -Yf0i");
        GMT::pstext(outfile,texts,"-JX1i/1i -R0/1/0/1 -N -O -K");
        return;
    }


    template<typename T>
    void psxy(const std::string &outfile, const std::vector<T> &X,
              const std::vector<T> &Y, const std::string &cmd){
        GMT::psxy(outfile,X.begin(),X.end(),Y.begin(),Y.end(),cmd);
        return;
    }

    template<typename T>
    void pshistogram(const std::string &outfile, const std::vector<T> &X, const std::string &cmd){
        GMT::pshistogram(outfile,X.begin(),X.end(),cmd);
        return;
    }

    void BeginPlot(const std::string &outfile){
        remove(outfile.c_str());
        GMT::psbasemap(outfile,"-JX1i/1i -R-1/1/-1/1 -Bwens -P -K");
        GMT::timestamp(outfile);
        return;
    }

    void SealPlot(const std::string &outfile){
        GMT::psbasemap(outfile,"-JX1i/1i -R-1/1/-1/1 -Bwens -O");
        remove("gmt.conf");
        remove("gmt.history");
        return;
    }
}

using namespace std;

int main(){

    bool PlotGrid=false;

    string outfile="GMT.ps";
    size_t NRow=3,NCol=3;
    size_t row,col;
    double Len=5,SpaceRatio=0.1,XSIZE=(NCol+2*SpaceRatio)*Len,YSIZE=(NRow+SpaceRatio)*Len+1,xp,yp;

    // config, define media ("gmt set PS_MEDIA 8.5ix8.5i").
    GMT::set("PS_MEDIA "+to_string(XSIZE)+"ix"+to_string(YSIZE)+"i");

    GMT::set("MAP_FRAME_PEN 0.4p,black");
//     GMT::set("MAP_TICK_PEN_PRIMARY 0.4p,black");
//     GMT::set("MAP_TICK_LENGTH_PRIMARY 0.007i");
//     GMT::set("MAP_TICK_LENGTH_SECONDARY 0.004i");
//     GMT::set("MAP_ANNOT_OFFSET_PRIMARY 0.02i");
//     GMT::set("FONT_LABEL 2.5p");
//     GMT::set("MAP_GRID_PEN_PRIMARY 0.25p,gray,.");
//     GMT::set("MAP_TICK_PEN_PRIMARY 0.25p,black");
    GMT::set("FONT_ANNOT 8p");

    // begin plot.
    GMT::BeginPlot(outfile);

    // pstitle.
    GMT::MoveReferencePoint(outfile,"-Xc -Yf"+to_string(YSIZE-1)+"i");
    if (PlotGrid) GMT::psbasemap(outfile,"-JX"+to_string(Len)+"i -R-10/10/-10/10 -Bxa10f5 -Bya10f5 -BWSne -O -K");

    vector<GMT::Text> texts;
    texts.push_back(GMT::Text(0,0,"GMT API Examples",24,"CB"));
    GMT::pstext(outfile,texts,"-JX1i/1i -R-10/10/-10/10 -N -O -K");

    // pscoast.
    row=1,col=1;
    xp=(col-1+SpaceRatio)*Len,yp=YSIZE-1-row*Len;
    GMT::MoveReferencePoint(outfile,"-Xf"+to_string(xp)+"i -Yf"+to_string(yp)+"i");
    if (PlotGrid) GMT::psbasemap(outfile,"-JX"+to_string(Len)+"i -R-10/10/-10/10 -Bxa10f5 -Bya10f5 -BWSne -O -K");

    GMT::pscoast(outfile,"-JR30/"+to_string(Len*(1-SpaceRatio))+"i -Rg -Bxa60g60 -Bya60g60 -BWSne -Glightgray -A10000 -O -K");

    // psxy.
    row=1,col=2;
    xp=(col-1+SpaceRatio)*Len,yp=YSIZE-1-row*Len;
    GMT::MoveReferencePoint(outfile,"-Xf"+to_string(xp)+"i -Yf"+to_string(yp)+"i");
    if (PlotGrid) GMT::psbasemap(outfile,"-JX"+to_string(Len)+"i -R-10/10/-10/10 -Bxa10f5 -Bya10f5 -BWSne -O -K");
    GMT::pscoast(outfile,"-JG-10/20/"+to_string(Len*(1-SpaceRatio))+"i -Rg -Bxa60g60 -Bya60g60 -BWSne -Glightgray -A10000 -O -K");

    vector<double> lon{-60,0,40},lat{-20,10,30};
    GMT::psxy(outfile,lon,lat,"-J -R -W0.5p,black,- -L -O -K");
    GMT::psxy(outfile,lon,lat,"-J -R -W0.5p,red -O -K");
    GMT::psxy(outfile,lon,lat,"-J -R -Gblue -Sc0.05i -W1p -O -K");

    texts.clear();
    texts.push_back(GMT::Text(lon[0],lat[0],"P1",10,"RM"));
    texts.push_back(GMT::Text(lon[1],lat[1],"P2",10,"LM"));
    GMT::pstext(outfile,texts,"-J -R -N -O -K");

    // pstext.
    row=1,col=3;
    xp=(col-1+SpaceRatio)*Len,yp=YSIZE-1-row*Len;
    GMT::MoveReferencePoint(outfile,"-Xf"+to_string(xp)+"i -Yf"+to_string(yp)+"i");
    if (PlotGrid) GMT::psbasemap(outfile,"-JX"+to_string(Len)+"i -R-5/5/-5/5 -Bxa5g1 -Bya10g1 -BWSne -O -K");

    texts.clear();
    texts.push_back(GMT::Text(0,0,"BigTitle",24,"CB"));
    texts.push_back(GMT::Text(0,-1,"@;red;SmallTitle@;;",12,"CB","Times-Bold"));
    GMT::pstext(outfile,texts,"-JX"+to_string(Len)+"i -R-5/5/-5/5 -N -O -K");

    // grdimage 1.
    row=2,col=1;
    xp=(col-1+SpaceRatio)*Len,yp=YSIZE-1-row*Len;
    GMT::MoveReferencePoint(outfile,"-Xf"+to_string(xp)+"i -Yf"+to_string(yp)+"i");
    if (PlotGrid) GMT::psbasemap(outfile,"-JX"+to_string(Len)+"i -R-10/10/-10/10 -Bxa5g1 -Bya10g1 -BWSne -O -K");

    vector<vector<double>> grid; 

    ifstream fpin("/home/shule/Research/Fun.Bash.c001/ritsema.2880");
    double la,lo,val,xinc=1,yinc=1,minval=numeric_limits<double>::max(),maxval=-minval;
    while (fpin >> la >> lo >> val) {
        minval=min(minval,val);
        maxval=max(maxval,val);
        grid.push_back({Lon2360(lo),la,val});
    }
    fpin.close();
    minval=-2.5;maxval=2.5;
    GMT::makecpt("-Cpolar -T"+to_string(minval)+"/"+to_string(maxval)+"/0.5 -I -Z > tmp.cpt");
    yp+=1.1;
    GMT::MoveReferencePoint(outfile,"-Xf"+to_string(xp)+"i -Yf"+to_string(yp)+"i");
    GMT::grdimage(outfile,grid,xinc,yinc,"-JR180/"+to_string(Len*(1-SpaceRatio))+"i -Rg -Ctmp.cpt -O -K");
    // -D scale center position relative to last reference piont.
    GMT::psscale(outfile,"-Ctmp.cpt -D"+to_string(Len*(1-SpaceRatio)/2)+"i/-0.5i/2i/0.1ih -O -K -B0.5:dVs(%):");
    remove("tmp.cpt");
    GMT::pscoast(outfile,"-JR180/"+to_string(Len*(1-SpaceRatio))+"i -Rg -Bxa60g60 -Bya60g60 -BWSne -W0.5p,black -A10000 -O -K");


    // grdimage 2.
    row=2,col=2;
    xp=(col-1+SpaceRatio)*Len,yp=YSIZE-1-row*Len;
    GMT::MoveReferencePoint(outfile,"-Xf"+to_string(xp)+"i -Yf"+to_string(yp)+"i");
    if (PlotGrid) GMT::psbasemap(outfile,"-JX"+to_string(Len)+"i -R-10/10/-10/10 -Bxa5g1 -Bya10g1 -BWSne -O -K");

    grid.clear();
    xinc=1,yinc=2,val=0;
    for (double i=0;i<3;i+=xinc)
        for (double j=0;j<3;j+=yinc)
            grid.push_back({i,j,val+=1});

    GMT::grdimage(outfile,grid,xinc,yinc,"-JX"+to_string(Len)+"i -R-10/10/-10/10 -Chot -O -K");


    // pshistogram.
    row=2,col=3;
    xp=(col-1+SpaceRatio)*Len,yp=YSIZE-1-row*Len;
    GMT::MoveReferencePoint(outfile,"-Xf"+to_string(xp)+"i -Yf"+to_string(yp)+"i");
    if (PlotGrid) GMT::psbasemap(outfile,"-JX"+to_string(Len)+"i -R-10/10/-10/10 -Bxa5g1 -Bya10g1 -BWSne -O -K");

    std::default_random_engine generator;
    std::normal_distribution<double> distribution(0,2.0);
    vector<double> gaussian_val;
    for (size_t i=0;i<3000;++i) gaussian_val.push_back(distribution(generator));

    GMT::pshistogram(outfile,gaussian_val,"-JX"+to_string(Len*(1-SpaceRatio))+"i -R-10/10/0/40 "
                                          "-Z1 -W1 -L0.5p -G50/50/250 -Bxa5f1 -Bya10f2+u\"%\" -BWS -O -K");

    // psxy with different colors.
    row=3,col=1;
    xp=(col-1+SpaceRatio)*Len,yp=YSIZE-1-row*Len;
    GMT::MoveReferencePoint(outfile,"-Xf"+to_string(xp)+"i -Yf"+to_string(yp)+"i");
    if (PlotGrid) GMT::psbasemap(outfile,"-JX"+to_string(Len)+"i -R-10/10/-10/10 -Bxa5g1 -Bya10g1 -BWSne -O -K");

    vector<vector<double>> data={{1.0,1.0,2.0},{3,-6,9.0},{0.3,0.6,0.9}};
    GMT::makecpt("-Cgray -T0/1 -I > tmp.cpt");
    GMT::psxy(outfile,data,"-J -R-10/10/-10/10 -Sc0.1i -Ctmp.cpt -W1p,black -Bxa5g1 -Bya10g1  -N -O -K");
    remove("tmp.cpt");


    // Seal the page.
    GMT::SealPlot(outfile);

    return 0;
}

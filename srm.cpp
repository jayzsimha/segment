#include <Magick++.h>
#include <math.h>
#include <stdlib.h>
#include <string>
#include <iostream> 
using namespace std; 
using namespace Magick;

/*
* Class to keep track of pixel pairs
*/
class Rmpair {
public:
        int r1,r2;
        int diff;
public:
        Rmpair() {
	        r1=0; 
	        r2=0; 
	        diff=0;
        }
        
        ostream& operator<<(ostream &out) {
                cout << this->r1 <<" "<< this->r2 <<"\n";
                return out;
        }
};

/*
* Union Find Data Structure of Tarjan for Disjoint Sets
*/
class UnionFind {
public:
        int *rank; 
        int *parent; 
        
        UnionFind(){cout << "Creating a UF DS for 0 elements";}
        UnionFind(int n);
        int Find(int k);
        int UnionRoot(int x, int y);
        ~UnionFind();
};

//
// Create a UF for n elements
//
UnionFind::UnionFind(int n) {
        int k;
        
        try {
                parent=new int[n]; 
                rank=new int[n];
        } catch(bad_alloc&) {
                cout << "Out of memory !!!" << endl;
        }
        cout << "Creating a UF DS for "<< n <<" elements"<< endl;

        for (k = 0; k < n; k++) {
                parent[k] = k;
                rank[k] = 0;     
        }
}

UnionFind::~UnionFind() {
        delete[] rank;
        delete[] parent;
}

//
// Find procedures
//
int UnionFind::Find(int k) {
        while (parent[k]!=k ) 
                k=parent[k];
        return k;
}

//
// Assume x and y being roots
//
int UnionFind::UnionRoot(int x, int y) {  
        if ( x == y ) 
                return -1;

        if (rank[x] > rank[y]) {
                parent[y]=x; 
                return x;
        }
        else {
                parent[x]=y;
                if (rank[x]==rank[y]) 
                        rank[y]++;
                return y;
        }
}

//
// Class SRM which is the main class consisting of all
// methods and variables required for statistical region
// merging. 
//
class SRM {
public:
//
//      Constructor
//      ifile : Input filename
//      ofile : Output filename
//      threshq : Threshold value for Q.
//
        SRM(string ifile,double threshq);
        
//
//      Public methods for performing segmentation.
//
        void InitializeSegmentation();
        Image* FullSegmentation(bool drawBorder);
        void GetRegion(int x,int y,PixelPacket* region);
        void copyDress(PixelPacket* dress, PixelPacket* target,int tw,int th);
        
//      Destructor
        ~SRM();
        
private:
//
//      Second level constructor to avoid exception during
//      object construction.
//
        void init();

//
//      Private methods for performing the actual segmentation
//        
        void Segmentation();
        bool MergePredicate(int reg1, int reg2);
        void MergeRegions(int C1, int C2);
        void MergeSmallRegion();
        void OutputSegmentation();
        void DrawBorder();        
        
//
//      Helper functions to get the min, max elements
//      and sorting the pixel pairs.
//
        double min(double a, double b);
        double max(double a, double b);
        double max3(double a,double b,double c);
        Rmpair* BucketSort(Rmpair a[], int n);
        
private:
        Image *img; // Input image
        Image *imgseg;  // Target image
        string filenameimg;  // input filename
        unsigned int depth;  // Color depth
        
        //for srm
        PixelPacket *raster;  // Pixel cache for input image
        PixelPacket *rastero; // Pixel cache for output image.
        double aspectratio;     
        
        double Q;             // Threshold Q.  
        UnionFind *UF;
        double g;  
        double logdelta;
        
        // Auxilliary buffers for union-find operations
        int *N;
	double *Ravg;
	double *Gavg;
	double *Bavg;
	int *C;
	
	// number of pixels that define a "small region" to be collapsed
	int smallregion;
	int borderthickness;
	
public:
        int w,h,n;
};

//
//      Constructor
//      ifile : Input filename
//      ofile : Output filename
//      threshq : Threshold value for Q.
//
SRM::SRM(string ifile, double threshq) {
        filenameimg = ifile;  // Input file
        Q = threshq;      // Threshold Q
        
        init(); // Second level construction.
}

//      Destructor
SRM::~SRM() {
        delete img;
        delete UF;
        delete[] N;
        delete[] Ravg;
        delete[] Gavg;
        delete[] Bavg;
        delete[] C;
}

void SRM::init() {
        
	borderthickness=1; // border thickness of regions
	
	// Color depth of the image
	depth = pow(2,(sizeof(Quantum)*8));
	g=(double)depth;
        
	img = new Image();
	img->read(filenameimg);
	//img->quantizeColorSpace(YCbCrColorspace);
	//img->display();
	
	int columns = img->columns();
	int rows = img->rows();
	
	// create pixel cache for performing pixel level
	// operations.
	raster = img->getPixels(0,0,columns,rows);
	
        w = columns;
        h = rows;
        aspectratio = (double)h/(double)w;
        n = w*h;
        logdelta = 2.0*log(6.0*n);
        
        // small regions are less than 1% of image pixels
        smallregion=(int)(0.005*n);
        
        // Create the required auxillary buffers.
        UF=new UnionFind(n);
        Ravg=new double[n];
        Gavg=new double[n];
        Bavg=new double[n];
        N=new int[n];
        C=new int[n];
}

double SRM::min(double a, double b) {
        if (a<b) return a; else return b;	
}

double SRM::max(double a, double b) {
        if (a>b) return a; else return b;	
}

double SRM::max3(double a,double b,double c) {
        return max(a,max(b,c));
}

Rmpair* SRM::BucketSort(Rmpair *a, int n)
{
        int i;
        int *nbe;
        int  *cnbe;
        Rmpair *b;
        
        nbe=new int[depth];
        cnbe=new int[depth];

        b=new Rmpair[n];
        
        for(i=0;i<depth;i++) nbe[i]=0;
        // class all elements according to their family
        for(i=0;i<n;i++) nbe[a[i].diff]++;
        // cumulative histogram
        cnbe[0]=0;
        for(i=1;i<depth;i++)
                cnbe[i]=cnbe[i-1]+nbe[i-1]; // index of first element of category i

        // allocation
        for(i=0;i<n;i++) {
                b[cnbe[a[i].diff]++]=a[i];
        }
        delete[] nbe;
        delete[] cnbe;
        delete[] a;
        return b;
}

//
// Initialization
//
void SRM::InitializeSegmentation() {
        int x,y,index;
        Quantum red,green,blue;
        
        for (y = 0; y < h; y++) {
                for (x = 0; x < w ; x++) {
                        index=y*w+x;

                        red = raster[index].red ;   
                        green = raster[index].green;
                        blue = raster[index].blue;

                        Ravg[index]=red;
                        Gavg[index]=green;
                        Bavg[index]=blue;
                        N[index]=1;
                        C[index]=index;
                }
        }
}

//
// Function to perform SRM segmentation
// drawBorder : Boolean value to suggest whether to draw region borders.
// value : Image*
// Note: Transfers ownership of the image.
// 
Image* SRM::FullSegmentation(bool drawBorder) {
        Segmentation();
        MergeSmallRegion();
        OutputSegmentation();
        if(drawBorder)
                DrawBorder();
        
        imgseg->syncPixels();
        return imgseg;
}
	        
bool SRM::MergePredicate(int reg1, int reg2) {
        double dR, dG, dB;
        double logreg1, logreg2;
        double dev1, dev2, dev;

        dR=(Ravg[reg1]-Ravg[reg2]); 
        dR*=dR;

        dG=(Gavg[reg1]-Gavg[reg2]); 
        dG*=dG;

        dB=(Bavg[reg1]-Bavg[reg2]); 
        dB*=dB;

        logreg1 = min(g,N[reg1])*log(1.0+N[reg1]);
        logreg2 = min(g,N[reg2])*log(1.0+N[reg2]);

        dev1=((g*g)/(2.0*Q*N[reg1]))*(logreg1 + logdelta);
        dev2=((g*g)/(2.0*Q*N[reg2]))*(logreg2 + logdelta);

        dev=dev1+dev2;

        return ( (dR<dev) && (dG<dev) && (dB<dev) );
}

//
// Merge two regions
//
void SRM::MergeRegions(int C1, int C2) {
        int reg,nreg;
        double ravg,gavg,bavg;

	reg=UF->UnionRoot(C1,C2);
	
	nreg=N[C1]+N[C2]; 
	ravg=(N[C1]*Ravg[C1]+N[C2]*Ravg[C2])/nreg;
	gavg=(N[C1]*Gavg[C1]+N[C2]*Gavg[C2])/nreg;
	bavg=(N[C1]*Bavg[C1]+N[C2]*Bavg[C2])/nreg;
	
	N[reg]=nreg;
	
	Ravg[reg]=ravg;
	Gavg[reg]=gavg;
	Bavg[reg]=bavg;
}

//
// Main segmentation procedure here
//        
void SRM::Segmentation() {
        int i,j,index;
        int reg1,reg2;
        Rmpair *order;
        int npair;
        int cpair=0;
        int C1,C2;
        Quantum r1,g1,b1; 
        Quantum r2,g2,b2;

        // Consider C4-connectivity here
        npair=2*(w-1)*(h-1)+(h-1)+(w-1);
        order = new Rmpair[npair];
	
        cout << "Building the initial image RAG (" << npair << " edges)"<< endl;
	
        for(i=0;i<h-1;i++) {
	        for(j=0;j<w-1;j++) {
	                index=i*w+j;

	                // C4  left
	                order[cpair].r1=index;
	                order[cpair].r2=index+1;
	
	
	                r1=raster[index].red;
	                g1=raster[index].green;
	                b1=raster[index].blue;
	
	                r2=raster[index+1].red;
	                g2=raster[index+1].green;
	                b2=raster[index+1].blue;
	
	                order[cpair].diff=(int)max3(abs(r2-r1),abs(g2-g1),abs(b2-b1));
	                cpair++;
	

	                // C4 below
                        order[cpair].r1=index;
	                order[cpair].r2=index+w;
	
	                r2=raster[index+w].red;
	                g2=raster[index+w].green;
	                b2=raster[index+w].blue;
		
                        order[cpair].diff=(int)max3(abs(r2-r1),abs(g2-g1),abs(b2-b1));
	                cpair++;
                }
        }
        //
        // The two border lines
        //
        for(i=0;i<h-1;i++) {
                index=i*w+w-1;
                order[cpair].r1=index;
	        order[cpair].r2=index+w;
	
	        r1=raster[index].red;
	        g1=raster[index].green;
	        b1=raster[index].blue;
	        
	        r2=raster[index+w].red;
	        g2=raster[index+w].green;
	        b2=raster[index+w].blue;
	        
                order[cpair].diff=(int)max3(abs(r2-r1),abs(g2-g1),abs(b2-b1));
	        cpair++;
        }

        for(j=0;j<w-1;j++) {
	        index=(h-1)*w+j;

                order[cpair].r1=index;
	        order[cpair].r2=index+1;
	
	        r1=raster[index].red;
	        g1=raster[index].green;
	        b1=raster[index].blue;
	        
	        r2=raster[index+1].red;
	        g2=raster[index+1].green;
	        b2=raster[index+1].blue;
	        	
                order[cpair].diff=(int)max3(abs(r2-r1),abs(g2-g1),abs(b2-b1));
	
	        cpair++;
        }

        cout<<"Sorting all "<< cpair <<" edges using BucketSort"<< endl;

        //
        // Sort the edges according to the maximum color channel difference
        //
        order = BucketSort(order,npair);

        // Main algorithm is here!!!
        cout << "Testing the merging predicate in a single loop"<< endl;

        for(i=0;i<npair;i++) {
	        reg1=order[i].r1;
                C1=UF->Find(reg1);
	        reg2=order[i].r2;
	        C2=UF->Find(reg2);
	        if ((C1!=C2)&&(MergePredicate(C1,C2))) 
	                MergeRegions(C1,C2);			
        }
        delete[] order;        
}

//
// Output the segmentation: Average color for each region
//
void SRM::OutputSegmentation() {
        int i,j, index, indexb;
        int r,g,b;
        
        imgseg = new Image(Geometry(w,h), Color(MaxRGB, MaxRGB, MaxRGB, 0));
        rastero = imgseg->getPixels(0,0,w,h);
        
        index=0;

        for(i=0;i<h;i++) { // for each row
	        for(j=0;j<w;j++) { // for each column
                        index=i*w+j;
                        indexb=UF->Find(index); // Get the root index 

                        //
                        // average color choice in this demo
                        //
                        r = (Quantum)Ravg[indexb];
                        g = (Quantum)Gavg[indexb];
                        b = (Quantum)Bavg[indexb];

                        rastero[index].red=r;
                        rastero[index].green=g;
                        rastero[index].blue=b;
	        }
        }
}

//
// Merge small regions
//
void SRM::MergeSmallRegion() {
        int i,j, C1,C2, index;

        index=0;

        for(i=0;i<h;i++) // for each row
	        for(j=1;j<w;j++) { // for each column
	                index=i*w+j;
	                C1=UF->Find(index);
	                C2=UF->Find(index-1);
	                if (C2!=C1) {
	                        if ((N[C2]<smallregion)||(N[C1]<smallregion)) 
	                                MergeRegions(C1,C2);
                        }	
                }
}

// Draw white borders delimiting perceptual regions
void SRM::DrawBorder() {
        int i, j, k, l, C1,C2, reg,index;

        for(i=1;i<h;i++) // for each row
                for(j=1;j<w;j++) { // for each column
                        index=i*w+j;

                        C1=UF->Find(index);
                        C2=UF->Find(index-1-w);

                        if (C2!=C1) {
                                for(k=-borderthickness;k<=borderthickness;k++)
                                        for(l=-borderthickness;l<=borderthickness;l++) {
                                                index=(i+k)*w+(j+l);
                                                if ((index>=0)&&(index<w*h)) {
                                                        rastero[index]=	Color("white");
                                                }
                                        }
                        }
                }
}

void SRM::GetRegion(int x,int y,PixelPacket* region) {
        
        int i,j,C1,C2,index;
        C2=UF->Find(x*w+y);
        
        for(i=1;i<h;i++) // for each row
                for(j=1;j<w;j++) { // for each column
                        index=i*w+j;

                        C1=UF->Find(index);
                        
                        if (C2==C1) {
                                region[index]=	raster[index];
                        }
                }
}

void SRM::copyDress(PixelPacket* dress, PixelPacket* target,int tw,int th) {
        int i,j,index;
        
        for(i=1;i<th;i++) // for each row
                for(j=1;j<tw;j++) { // for each column
                        index=i*tw+j;
                        if(dress[index].red > 0 ||
                           dress[index].blue > 0 ||
                           dress[index].green > 0) {
                                target[index] = dress[index];
                        }
                }
}

int main(int argc,const char** argv) {
        cout<<"starting now..."<<endl;
        SRM *srm = new SRM(argv[1],atof(argv[3]));
        bool showBorders = (bool)atoi(argv[4]);
        // Initialize various parameters of image magick
	InitializeMagick(*argv);
	
        srm->InitializeSegmentation();
        Image *outimg = srm->FullSegmentation(showBorders);
        
        cout << "writing output to file" << argv[2] << endl;
        //outimg->write(argv[2]);
        
        Image *dress = new Image(Geometry(srm->w,srm->h), Color(0, 0, 0, 0));
        PixelPacket *region = dress->getPixels(0,0,srm->w,srm->h);
        int x = atoi(argv[5]);
        int y = atoi(argv[6]);
        srm->GetRegion(x,y,region);
        dress->syncPixels();
        dress->write(argv[2]);
        
        /*Image *target_img = new Image();
        string target_image = argv[7];
        target_img->read(target_image);
        int width = target_img->columns();
        int height = target_img->rows();
        PixelPacket *target = target_img->getPixels(0,0,width,height);
        
        //dress->scale(Geometry(width,height));
        //dress->display();
        //dress->gaussianBlur(0,0.5);
        PixelPacket *region_mod = dress->getPixels(0,0,width,height);
        srm->copyDress(region_mod,target,width,height);
        target_img->syncPixels();
        //target_img->display();
        target_img->write(argv[2]);*/
        
        delete srm;
        delete outimg;
}

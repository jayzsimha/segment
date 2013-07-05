#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>
#include <math.h>
#include <stdarg.h>
#include <mm_malloc.h>
#include <emmintrin.h>
#include <highgui.h>
#include <cv.h>
#include <cvaux.h>
#include <iostream>

//#define MAX(x,y) x > y ? return x : return y
//#define MAX(x,y) x > y ? return y : return x

//#define __DEBUG_
#define _INLINE inline

using namespace std;
using namespace cv;

#ifdef __DEBUG_
const string winName = "image";
#endif


namespace RTTL {
  // It is an empty implementation for stand-alone app.
  class AtomicCounter  {
  private:
    int m_counter;

  public:
    AtomicCounter() {
      m_counter = -1;
    }

    _INLINE int inc() {
      return ++m_counter;
    }
  };
}
typedef unsigned char colortype;
struct rgba {
  colortype r, g, b, a;
};

#include "MLAA.hxx"
int pixelmode = 0x40;

/** 
* getBinMask
* Function to get a binary mask from the mask returned by grabCut.
* 0 indicates background and 1 indicates a foreground pixel. 
* 
* @param comMask Input mask.
* @param binMask Output binary mask.
*
**/
void getBinMask( const Mat& comMask, Mat& binMask )
{
    if( comMask.empty() || comMask.type()!=CV_8UC1 )
        CV_Error( CV_StsBadArg, "comMask is empty or has incorrect type (not CV_8UC1)" );
    if( binMask.empty() || binMask.rows!=comMask.rows || binMask.cols!=comMask.cols )
        binMask.create( comMask.size(), CV_8UC1 );
    binMask = comMask & 1;
}

/**
* saveMask
* Function to save the mask into a file outfilename.
* 
* @param outfilename, the file to write the mask into.
* @param comMask Input Mask to be written.
*
**/
void saveMask(const string outfilename, const Mat& comMask)
{
    Mat binMask;
    if( comMask.empty() || comMask.type()!=CV_8UC1 )
        CV_Error( CV_StsBadArg, "comMask is empty or has incorrect type (not CV_8UC1)" );
    if( binMask.empty() || binMask.rows!=comMask.rows || binMask.cols!=comMask.cols )
        binMask.create( comMask.size(), CV_8UC1 );
    binMask = comMask & 1;
    imwrite( outfilename, binMask);
}

/** 
* changeMask
* Method to change the mask based on the foreground and 
* background pixels.
* 
* @param mask output mask created using params bgdPixels & fgdPixels
* @param bgdPixels list of background pixels.
* @param fgdPixels list of foreground pixels.
* @param trueMask boolean value indicating whether bgdPixels & fgdPixels are true or probabilistic.
*
**/
void changeMask( Mat& mask, vector<Point> bgdPixels, vector<Point> fgdPixels , bool trueMask=true)
{
    // GC_BGD - Background pixel
    // GC_FGD - Foreground pixel
    // GC_PR_BGD - Possibly a background pixel.
    // GC_PR_FGD - Possibly a foreground pixel.
    vector<Point>::const_iterator it, itEnd;
    if(!bgdPixels.empty()) {
        it = bgdPixels.begin(), itEnd = bgdPixels.end();
        for( ; it != itEnd; ++it )
            if(trueMask)
                mask.at<uchar>(*it) = GC_BGD;
            else
                mask.at<uchar>(*it) = GC_PR_BGD;
    }
    
    if(!fgdPixels.empty()) {
        it = fgdPixels.begin(), itEnd = fgdPixels.end();
        for( ; it != itEnd; ++it )
            if(trueMask)
                mask.at<uchar>(*it) = GC_FGD;
            else
                mask.at<uchar>(*it) = GC_PR_FGD;
    }
}

#ifdef __DEBUG_
/** 
* showImage
* Diagnostic method to display images. A wrapper for imshow.
* 
* @param _img image to display
* @param _mask mask value as returned by grabCut.
* @param bgdPixels list of background pixels.
* @param _fgdPixels list of foreground pixels.
*
**/
void showImage( Mat& _img, Mat& _mask, vector<Point>& _bgdPxls, vector<Point>& _fgdPxls )
{
    const Scalar RED = Scalar(0,0,255);
    const Scalar BLUE = Scalar(255,0,0);
    
    Mat res = Mat::zeros(_img.size(),_img.type());
    Point pts[4];
    pts[0] = Point(0,0);
    pts[1] = Point(_img.cols,0);
    pts[2] = Point(_img.cols,_img.rows);
    pts[3] = Point(0,_img.rows);
    fillConvexPoly(res,pts,4,CV_RGB(255,0,255));
    
    Mat binMask;
    if( _mask.empty() )
        _img.copyTo( res );
    else
    {
        getBinMask( _mask, binMask );
        _img.copyTo( res, binMask );
    }

    vector<Point>::const_iterator it;
    for( it = _bgdPxls.begin(); it != _bgdPxls.end(); ++it )
        circle( res, *it, 1, BLUE );
    for( it = _fgdPxls.begin(); it != _fgdPxls.end(); ++it )
        circle( res, *it, 1, RED );
    imshow( winName, res );
}
#endif

void imwrite_withalpha(string outfilename,Mat& res) {
    Mat target = Mat::zeros(res.size(),CV_8UC4);
    for (int ix=0; ix<res.cols; ix++) {
        for (int iy=0; iy<res.rows; iy++) {
            uchar b = res.at<Vec3b>( iy, ix )[0];
            uchar g = res.at<Vec3b>( iy, ix )[1];
            uchar r = res.at<Vec3b>( iy, ix )[2];

            // skip zeros
            if( b == 255 && g == 0 && r == 255 ){}
            else {

                    cv::Vec4b& rgba = target.at<cv::Vec4b>(iy, ix);

                    rgba[0] = b;
                    rgba[1] = g;
                    rgba[2] = r;
                    rgba[3] = UCHAR_MAX;
            }
        }
    }
    std::vector<int> compression_params;
    compression_params.push_back(CV_IMWRITE_PNG_COMPRESSION);
    compression_params.push_back(9);
    imwrite(outfilename.c_str(), target, compression_params);
}

void antiAlias(Mat& target, Mat& _matrix) {
    int width = _matrix.rows;
    int height = _matrix.cols;
    int width0 = width;
    rgba* data = (rgba*)_mm_malloc(sizeof(rgba) * width*height, 16);
    for (int y = height - 1; y >= 0; y--) {
      for (int x = 0; x < width0; x++ ) {
        if (x >= width) continue;
        int index = x + y * width;
        Vec3b pixel = _matrix.at<Vec3b>(x,y);
        data[index].r = pixel[2];
        data[index].g = pixel[1];
        data[index].b = pixel[0];
        data[index].a = 0;
      }
    }
    RTTL::AtomicCounter job;
    MLAA((unsigned int*)data, NULL, width, height, job);
    for (int y = height - 1; y >= 0; y-- ) {
        for (int x = 0; x < width; x++) {
          int index = x + y* width;
          target.at<Vec3b>(x,y)[2] = data[index].r;
          target.at<Vec3b>(x,y)[1] = data[index].g;
          target.at<Vec3b>(x,y)[0] = data[index].b;
        }
    }
}

void morph(Mat& _mask) {
    Mat opened = Mat::zeros(_mask.size(), _mask.type());
    Mat closed = Mat::zeros(_mask.size(), _mask.type());
    int erosion_type = MORPH_ELLIPSE;
    Mat element = getStructuringElement( erosion_type,
                                       Size( 11,11 ),
                                       Point(0, 0));

    /// Apply the open operation
    morphologyEx(_mask,opened,MORPH_OPEN,element);
    morphologyEx(opened, closed, MORPH_CLOSE, element);
    closed.copyTo(_mask);    
}

/** 
* saveSegmentedImage
* method to save image to a file. A wrapper for imwrite.
* 
* @param outfilename output file name
* @param _img image to display
* @param _mask mask value as returned by grabCut.
* @param bgdPixels list of background pixels.
* @param _fgdPixels list of foreground pixels.
*
**/
void saveSegmentedImage( string outfilename, Mat& _img, Mat& _mask, vector<Point>& _bgdPxls, vector<Point>& _fgdPxls )
{
    Mat res = Mat::zeros(_img.size(),_img.type());
    Point pts[4];
    pts[0] = Point(0,0);
    pts[1] = Point(_img.cols,0);
    pts[2] = Point(_img.cols,_img.rows);
    pts[3] = Point(0,_img.rows);
    fillConvexPoly(res,pts,4,CV_RGB(255,0,255));

    Mat target = Mat::zeros(_img.size(),_img.type());
    Mat binMask;

    if( _mask.empty() )
        _img.copyTo( res );
    else {
        getBinMask( _mask, binMask );
        morph(binMask);
        _img.copyTo( res, binMask );
    }
    antiAlias(target,res);
    imwrite_withalpha( outfilename, res);
    //imwrite( outfilename, res);
}

Rect bounding_box(Mat _img,int y0,int y1) {   
    Mat img_gray,edge_image;
    cvtColor(_img, img_gray, CV_BGR2GRAY );                                                                                     
    Canny(img_gray,edge_image,180,10);
    imwrite("edgeImage.png",edge_image);
    int left = edge_image.cols;
    int right = 0;
    int x = 0;
    int y = 0;
    
    while(x < edge_image.cols) {
            y = y0;
            while(y < y1) {
                    int pix = (int)edge_image.at<uchar>(y,x);
                    if(pix > 0) {
                            if(left > x) { 
                                left = x;
                                #ifdef __DEBUG_
                                    cout << x <<" "<< y << " "<<pix << endl;
                                #endif
                            }
                            if(right < x) { right = x; }
                    }
                    y++;
            }
            x++;
    }
    return Rect(Point(left,y0),Point(right,y1));
}

/*Rect bounding_box(Mat _img)
{
    Mat img_gray,edge_image;
    vector<vector<Point> > contours;
    vector<Point> flattened_contours;
    vector<Point> approx_curve;
    vector<Vec4i> hierarchy;

    cvtColor(_img, img_gray, CV_BGR2GRAY );
    Canny(img_gray,edge_image,200,255,5,true);
    imwrite("edge_image.png",edge_image);
    findContours(edge_image,contours,hierarchy,CV_RETR_EXTERNAL,CV_CHAIN_APPROX_SIMPLE);
    
    vector<Point>::const_iterator it, itEnd;
    for(int i=0;i<contours.size();i++) {
        for(int j=0;j<contours[i].size();j++)
            flattened_contours.push_back(contours[i][j]);
    }

    /*vector<vector<Point> >hull( contours.size() );
    for( int i = 0; i < contours.size(); i++ ) {
        convexHull( Mat(contours[i]), hull[i], false ); 
    }

    Mat dst = Mat::zeros(_img.size(),_img.type());
    int idx = 0;
    for( ; idx >= 0; idx = hierarchy[idx][0] )
    {    
        Scalar color( rand()&255, rand()&255, rand()&255 );        
        drawContours( dst, contours, idx, color, CV_FILLED, 8, hierarchy );
        //drawContours( dst, hull, idx, Scalar(255,255,255), 1, 8, vector<Vec4i>());
    }

    RotatedRect rect = minAreaRect(flattened_contours);
    Point2f vertices[4];
    rect.points(vertices);
    for (int i = 0; i < 4; i++)
        line(dst, vertices[i], vertices[(i+1)%4], Scalar(0,255,0));
    
    //approxPolyDP(flattened_contours,approx_curve,1.2,true);
    Rect rect = boundingRect(flattened_contours);
    imwrite("contours.png",dst);
    //return Rect(Point(0,0),Point(_img.rows,_img.cols));
    return rect;
}*/

bool detectFace(Mat& _img, Rect& r) {
    String haar_face_cascade_name = "/usr/local/share/OpenCV/haarcascades/haarcascade_frontalface_alt.xml";
    String haar_profile_face_cascade_name = "/usr/local/share/OpenCV/haarcascades/haarcascade_profileface.xml";
    CascadeClassifier face_cascade, face_profile_cascade;
    CascadeClassifier hand_cascade;
    std::vector<Rect> faces;
    Mat img_gray;
    
    if( !face_cascade.load( haar_face_cascade_name ) ) { return false; } 
    if( !face_profile_cascade.load( haar_profile_face_cascade_name ) ){ return false; }
    
    cvtColor( _img, img_gray, CV_BGR2GRAY );
    equalizeHist( img_gray, img_gray );
    face_cascade.detectMultiScale( img_gray, faces, 1.1, 3, 0|CV_HAAR_SCALE_IMAGE, Size(5,5));
    if(faces.size() < 0) {
        face_profile_cascade.detectMultiScale( img_gray, faces, 1.1, 3, 0|CV_HAAR_SCALE_IMAGE,Size(5,5));
    }
    
    int i=faces.size()-1;
    if(faces.size() > 0) {
        #ifdef __DEBUG_
        cout << "Yippee !! face found " << faces.size() <<endl;
        #endif
        r = faces[i];
        return true;
    }
    else {return false;}
}

/** 
* segmentClothing
* segment clothing using graph mincut based segmentation.
* Uses grabCut method in OpenCV to do the segmentation.
* 
* @param _img input image.
* @param _rect rectangle which encompasses the object.
*
**/
void segmentImage(Mat& _img, Rect& _rect,string outfilename)
{
    int max_iter = 10;
    Mat mask = Mat::zeros(_img.size(),CV_8UC1);
    Mat bgdModel, fgdModel;
    vector<Point> fgdPxls, bgdPxls;
    vector<Point> edgePxls;
    int iterCount = 0;
    
#ifdef __DEBUG_
    long start = clock();
#endif
    Rect face;
    int y0 = 10;
    int y1 = _img.rows;

    if(detectFace(_img,face))
        _rect.y = face.y-(_img.rows/5);
    
    _rect = bounding_box(_img,y0,y1);

    Mat img_gray,edge_image;
    cvtColor(_img, img_gray, CV_BGR2GRAY );                                                                                     
    Canny(img_gray,edge_image,180,10);

    vector<Vec4i> lines;
    HoughLinesP(edge_image, lines, 1, CV_PI/180, 50, 50, 10 );
    Mat str_edges = Mat::zeros(_img.size(),CV_8UC1);
    for( size_t i = 0; i < lines.size(); i++ ) {
        Vec4i l = lines[i];
        line( str_edges, Point(l[0], l[1]), Point(l[2], l[3]), 255, 3);
    }
    imwrite("edgeImage.png",str_edges);

    for(int x=0;x<str_edges.cols;x++) {
        for(int y=0;y<str_edges.rows;y++) {
            int pix = (int)str_edges.at<uchar>(y,x);
            Point pt = Point(x,y);
            if(pix > 0) {//|| (faceDetected && pt.inside(face))) {
                bgdPxls.push_back(pt);
            }
        }
    }
    changeMask(mask,bgdPxls,fgdPxls,true);
    bgdPxls.clear(); fgdPxls.clear();
    grabCut( _img, mask, _rect, bgdModel, fgdModel, 0, GC_INIT_WITH_RECT|GC_INIT_WITH_RECT);
#ifdef __DEBUG_    
    showImage( _img, mask, bgdPxls, fgdPxls );
#endif
    assert( bgdPxls.empty() && fgdPxls.empty());
    
    while(iterCount < max_iter) {
        changeMask(mask,bgdPxls,fgdPxls);
        bgdPxls.clear(); fgdPxls.clear();
        grabCut( _img, mask, _rect, bgdModel, fgdModel, 1 );
        iterCount++;
    }
    
    saveSegmentedImage(outfilename, _img, mask, bgdPxls, fgdPxls );

#ifdef __DEBUG_
    long end = (clock() - start)/CLOCKS_PER_SEC;
    cout << "Entire thing took " << end << "clock ticks" << endl;
    showImage( _img, mask, bgdPxls, fgdPxls );
    cvWaitKey(0);
#endif
}

/** 
* temporary main method.
* 
* @usage grabcut <input image> x1 y1 x2 y2 (rectangle co-ordinates.)
* @param input image.
* @param top left and bottom right co-ordinates for rectangle which encompasses the object.
*
**/
int main( int argc, char** argv)
{
    Mat image;
    Rect rect;
    Point p1, p2;
    
    // argv[1] is the filename
    string filename = argv[1];
    if( filename.empty() )
        return 1;
        
    // read image file.
    image = imread(filename,1);
    if( image.empty())
        return 1;

    string outfilename = argv[2];
    
#ifdef __DEBUG_    
    cvNamedWindow( winName.c_str(), CV_WINDOW_AUTOSIZE );
#endif
    
    p1 = Point(0,0);
    p2 = Point((image.cols),(image.rows));
    rect = Rect(p1,p2);
    
    segmentImage(image,rect,outfilename);
#ifdef __DEBUG_      
    cvDestroyWindow( winName.c_str() );
#endif

    return 0;
}

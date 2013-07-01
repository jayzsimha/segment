Graph cut based segmention
==========================

Instructions for making it work:

Dependency:
OpenCV-2.4.x or higher

Build:
cmake .
make

Run:
./segment <inputfilename> <outputfilename>

Algorithmic parameters:
No manual parameters.

Advantages:
    A fairly generic solution. Works reliably even if the person is wearing quite colourful clothing.
    Binary classification, hence "ideally" no manual input would be necessary
Disadvantages:
    Initialisation of the algorithm. Either a smallest rectangle, or a mask is needed. This makes life a little difficult.

Statistical region merging (SRM):
================================

Instructions for making it work:

Dependency:
ImageMagick-6.x.x or higher

Build:
g++ srm.cpp `Magick++-config --cppflags --cxxflags --ldflags --libs` -o srm

Run:
./srm <inputfilename> <outputfilename> <sigma> <draw_borders> <x> <y>

Algorithmic parameters:
sigma : The complexity parameter. Lesser the complexity, we have simpler segments. The main problem here is that, since the images are quite varied statistically, we cannot determine this parameter empirically. Hence needs to be selected interactively by the user. Something like a sliding scale might be helpful. 
<x> & <y> are the seed pixels for region selection. Since this is not binary segmentation, there are going to be multiple regions in the same image. The co-ordinates allows us to return a particular region. Again needs to be interactively obtained from the user.

Advantages:
    Simple and fast. Very little processing power is used. 
Disadvantages:
    Not so reliable in the case of complex images.
    parameters need to be selected for individual images. Hence, user interaction is vital.

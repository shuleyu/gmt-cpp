# gmt-cpp
c++ gmt plotting simple wrapper (see https://github.com/shuleyu/CPP-Library/blob/master/GMT.hpp)

Example programs see https://github.com/shuleyu/CPP-Library/blob/master/Examples/GMT.cpp

The files above will be updated in the future. The example in this repository (CMT.cpp) will not be updated.

Usage:

$ c++ GMT.cpp -lgmt # -I${PathToGMTIncludeFiles} -L${PathToGMTLibFiles} # if GMT is not installed in the default directory.

$ ./a.out ; gs GMT.ps # or ps2pdf GMT.ps GMT.pdf

![alt text](https://github.com/shuleyu/gmt-cpp/blob/master/Examples.png)

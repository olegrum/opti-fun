#include <iostream>
#include <cmath>
/*
#include <optitms.h>
#include <optidfp.h>
#include <optipowell.h>
*/
using namespace std;
inline double f(double x,double y)
{
    return sin(x+y)+(x-y)*(x-y)-1.5*x+2.5*y+1;
}
int main()
{
    //optitms(x1,x2,y1,y2);
    //optidfp(x1,x2,y1,y2);
    //optipowell(x1,x2,y1,y2);
    cout<<f(-0.54719,-1.54719);
    return 0;
}

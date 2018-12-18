#ifndef _POINT2D_H_
#define _POINT2D_H_

#include <iostream>
using namespace std;
/**
 * Primitive 2D point class used as input for the LineParamEstimator.
 *
 * Author: Ziv Yaniv (zivy@cs.huji.ac.il)
 */
class Point2D {
public:
	Point2D(double px, double py) : x(px), y(py) {}
	double x;
	double y;
};

inline ostream &operator<<(ostream &output,const Point2D &pnt)
{
	output<<pnt.x<<' '<<pnt.y;
	return(output);
}

#endif //_POINT2D_H_
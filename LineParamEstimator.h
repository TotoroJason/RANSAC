#ifndef _LINE_PARAM_ESTIMATOR_H_
#define _LINE_PARAM_ESTIMATOR_H_

#include "ParameterEsitmator.h"
#include "Point2D.h"

/**
 * This class estimates the parameters of 2D lines.
 * A 2D line is represented as: (*) dot(n,p-a)=0 
 *                              where n is the line normal (|n| = 1) and 'a' is a 
 *                              point on the line. 
 * All points 'p' which satisfy equation (*) are on the line.
 *
 * The reason for choosing this line parametrization is that it can represent
 * all lines, including vertical and horizontal, unlike the slope intercept (y=ax+b)
 * parametrization.
 *
 * Author: Ziv Yaniv (zivy@cs.huji.ac.il)
 */

class LineParamEstimator : public ParameterEsitmator<Point2D,double> {
public:
	LineParamEstimator(double delta);

	/**
	 * Compute the line defined by the given data points.
	 * @param data A vector containing two 2D points.
	 * @param This vector is cleared and then filled with the computed parameters.
	 *        The parameters of the line passing through these points [n_x,n_y,a_x,a_y]
	 *        where ||(n_x,ny)|| = 1.
	 *        If the vector contains less than two points then the resulting parameters
	 *        vector is empty (size = 0).
	 */
	virtual void estimate(std::vector<Point2D *> &data, std::vector<double> &parameters);

	/**
	 * Compute a least squares estimate of the line defined by the given points.
	 * This implementation is of an orthogonal least squares error.
	 *
	 * @param data The line should minimize the least squares error to these points.
	 * @param parameters This vector is cleared and then filled with the computed parameters.
	 *                   Fill this vector with the computed line parameters [n_x,n_y,a_x,a_y]
	 *                   where ||(n_x,ny)|| = 1.
	 *                   If the vector contains less than two points then the resulting parameters
	 *                   vector is empty (size = 0).
	 */
	virtual void leastSquaresEstimate(std::vector<Point2D *> &data, std::vector<double> &parameters);

	/**
	 * Return true if the distance between the line defined by the parameters and the
	 * given point is smaller than 'delta' (see constructor).
	 * @param parameters The line parameters [n_x,n_y,a_x,a_y].
	 * @param data Check that the distance between this point and the line is smaller than 'delta'.
	 */
	virtual bool agree(std::vector<double> &parameters, Point2D &data);

	/**
	 * Test the class methods.
	 */
	static void debugTest(ostream &out);

private:
	double m_deltaSquared; //given line L and point P, if dist(L,P)^2 < m_delta^2 then the point is on the line
};

#endif //_LINE_PARAM_ESTIMATOR_H_
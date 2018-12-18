#include <math.h>
#include "LineParamEstimator.h"

LineParamEstimator::LineParamEstimator(double delta) : m_deltaSquared(delta*delta) {}
/*****************************************************************************/
/*
 * Compute the line parameters  [n_x,n_y,a_x,a_y]
 */
void LineParamEstimator::estimate(std::vector<Point2D *> &data, 
																	std::vector<double> &parameters)
{
	parameters.clear();
	if(data.size()<2)
		return;
	double nx = data[1]->y - data[0]->y;
	double ny = data[0]->x - data[1]->x;
	double norm = sqrt(nx*nx + ny*ny);
	
	parameters.push_back(nx/norm);
	parameters.push_back(ny/norm);
	parameters.push_back(data[0]->x);
	parameters.push_back(data[0]->y);		
}
/*****************************************************************************/
/*
 * Compute the line parameters  [n_x,n_y,a_x,a_y]
 */
void LineParamEstimator::leastSquaresEstimate(std::vector<Point2D *> &data, 
																							std::vector<double> &parameters)
{
	double meanX, meanY, nx, ny, norm;
	double covMat11, covMat12, covMat21, covMat22; // The entries of the symmetric covarinace matrix
	int i, dataSize = data.size();

	parameters.clear();
	if(data.size()<2)
		return;

	meanX = meanY = 0.0;
	covMat11 = covMat12 = covMat21 = covMat22 = 0;
	for(i=0; i<dataSize; i++) {
		meanX +=data[i]->x;
		meanY +=data[i]->y;

		covMat11	+=data[i]->x * data[i]->x;
		covMat12	+=data[i]->x * data[i]->y;
		covMat22	+=data[i]->y * data[i]->y;
	}

	meanX/=dataSize;
	meanY/=dataSize;

	covMat11 -= dataSize*meanX*meanX;
  covMat12 -= dataSize*meanX*meanY;
	covMat22 -= dataSize*meanY*meanY;
	covMat21 = covMat12;

	if(covMat11<1e-12) {
		nx = 1.0;
	  ny = 0.0;
	}
	else {	    //lamda1 is the largest eigen-value of the covariance matrix 
	           //and is used to compute the eigne-vector corresponding to the smallest
	           //eigenvalue, which isn't computed explicitly.
		double lamda1 = (covMat11 + covMat22 + sqrt((covMat11-covMat22)*(covMat11-covMat22) + 4*covMat12*covMat12)) / 2.0;
		nx = -covMat12;
		ny = lamda1 - covMat22;
		norm = sqrt(nx*nx + ny*ny);
		nx/=norm;
		ny/=norm;
	}
	parameters.push_back(nx);
	parameters.push_back(ny);
	parameters.push_back(meanX);
	parameters.push_back(meanY);
}
/*****************************************************************************/
/*
 * Given the line parameters  [n_x,n_y,a_x,a_y] check if
 * [n_x, n_y] dot [data.x-a_x, data.y-a_y] < m_delta
 */
bool LineParamEstimator::agree(std::vector<double> &parameters, Point2D &data)
{
	double signedDistance = parameters[0]*(data.x-parameters[2]) + parameters[1]*(data.y-parameters[3]); 
	return ((signedDistance*signedDistance) < m_deltaSquared);
}
/*****************************************************************************/
void LineParamEstimator::debugTest(ostream &out)
{
	std::vector<double> lineParameters;
	LineParamEstimator lpEstimator(0.5);
	std::vector<Point2D *> pointData;

	pointData.push_back(new Point2D(7,7));
	pointData.push_back(new Point2D(-1,-1));
	lpEstimator.estimate(pointData,lineParameters);
	out<<"[n_x,n_y,a_x,a_y] [ "<<lineParameters[0]<<", "<<lineParameters[1]<<", ";
	out<<lineParameters[2]<<", "<<lineParameters[3]<<" ]"<<endl;
	lpEstimator.leastSquaresEstimate(pointData,lineParameters);
	out<<"[n_x,n_y,a_x,a_y] [ "<<lineParameters[0]<<", "<<lineParameters[1]<<", ";
	out<<lineParameters[2]<<", "<<lineParameters[3]<<" ]"<<endl;
	delete pointData[0];
	delete pointData[1];
	pointData.clear();

	
	pointData.push_back(new Point2D(6,12));
	pointData.push_back(new Point2D(6,6));
	lpEstimator.estimate(pointData,lineParameters);
	out<<"[n_x,n_y,a_x,a_y] [ "<<lineParameters[0]<<", "<<lineParameters[1]<<", ";
	out<<lineParameters[2]<<", "<<lineParameters[3]<<" ]"<<endl;
	lpEstimator.leastSquaresEstimate(pointData,lineParameters);
	out<<"[n_x,n_y,a_x,a_y] [ "<<lineParameters[0]<<", "<<lineParameters[1]<<", ";
	out<<lineParameters[2]<<", "<<lineParameters[3]<<" ]"<<endl;
	delete pointData[0];
	delete pointData[1];
	pointData.clear();


	pointData.push_back(new Point2D(7,9));
	pointData.push_back(new Point2D(10,9));
	lpEstimator.estimate(pointData,lineParameters);
	out<<"[n_x,n_y,a_x,a_y] [ "<<lineParameters[0]<<", "<<lineParameters[1]<<", ";
	out<<lineParameters[2]<<", "<<lineParameters[3]<<" ]"<<endl;
	lpEstimator.leastSquaresEstimate(pointData,lineParameters);
	out<<"[n_x,n_y,a_x,a_y] [ "<<lineParameters[0]<<", "<<lineParameters[1]<<", ";
	out<<lineParameters[2]<<", "<<lineParameters[3]<<" ]"<<endl;
	delete pointData[0];
	delete pointData[1];
	pointData.clear();
}

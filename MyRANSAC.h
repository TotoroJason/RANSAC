
#ifndef MY_RANSAC_H_
#define MY_RANSAC_H_

#include <set>
#include <vector>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "ParameterEsitmator.h"
#include "Point2D.h"

/*
		这个实现，只处理二维点，提取直线


		 Compute the k parameter (k=log(1-p)/log(1-w^n))
		 k，迭代次数
		 p，probability_表示一些迭代过程中从数据集内随机选取出的点均为局内点的概率
		 w = 局内点的数目 / 数据集的数目
		 n，数据集的点数


		 model_parameters: 模型参数，这里是直线方程的参数。直线方程使用 dot(n, p-a) = 0来表示
							n为法线，|n|=1，
							a为直线上一个点
							所以这个参数存放顺序是， [n_x,n_y,a_x,a_y]

		param_estimator：	模型估计对象，有两种方法，一种是直接取两点估计直线方程，另一种是最小二乘法

		data：	数据集

		num_for_estimate：估计模型所需的点数，直线的估计，随机取两个点就行
								
		
		最小二乘法，网上讲得比较好的资料
		https://www.zhihu.com/question/37031188

		随机采样一致性，网上资料
		http://www.cnblogs.com/xrwang/archive/2011/03/09/ransac-1.html
									
*/




class MyRANSAC 
{
public:
	MyRANSAC();


	double compute(std::vector<double> &model_parameters, 
		ParameterEsitmator<Point2D,double> *param_estimator , 
		std::vector<Point2D> &data, 
		int num_for_estimate		//估计模型时，需要的点数。此处直处理直线，设为2, 初始采样2个点得到直线方程
		);


private:

	class SubSetIndexComparator 
	{
	private:
		int m_length;
	public:
		SubSetIndexComparator(int arrayLength) : m_length(arrayLength){}
		bool operator()(const int *arr1, const int *arr2) const 
		{
			//大小都要比较，不然就会有二义性
			for(int i=0; i<m_length; i++)
			{	if(arr1[i] < arr2[i])
				{
					return true;
				}
				if (arr1[i] > arr2[i])
				{
					return false;
				}
			}
			return false;			
		}
	};

	double probability_;	//p表示一些迭代过程中从数据集内随机选取出的点均为局内点的概率

	int max_iterations_;	//最大迭代次数
};


#endif
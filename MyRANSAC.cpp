#include "MyRANSAC.h"
#include <math.h>
#include <algorithm>

MyRANSAC::MyRANSAC()
	:probability_(0.99)
	, max_iterations_ (1000)
{

}

double MyRANSAC::compute(std::vector<double> &model_parameters, 
	ParameterEsitmator<Point2D,double> *param_estimator , 
	std::vector<Point2D> &data, 
	int num_for_estimate	//适用于模型的最少数据个数
	)
{

	int num_data = data.size();

	int iterations = 0;	
	double k = max_iterations_;	//算法的迭代次数


	double log_probability  = log (1.0 - probability_);
	double one_over_indices = 1.0 / static_cast<double> (num_data);
	
	

	short *best_votes = new short[num_data]; //最好模型的点集索引标志, data[i] 符合最好模型置1，否则0，用于记录符合模型的点的索引，以及统计点数
	short *cur_votes = new short[num_data];  //当前模型的点集索引标志，data[i] 符合当前模型置1，否则0

	SubSetIndexComparator sub_set_index_comparator(num_for_estimate);		//比较器，用于判定初始样本是否一样
	std::set<int *, SubSetIndexComparator > chosen_sub_sets(sub_set_index_comparator);	//用于记录已经用过的子集，避免重复


	int* cur_inti_sub_set_indexs = NULL;	//初始选择的样本
	

	int best_model_num = -1;	//最好模型中的点数
	int maybe_inliers_num = 0;	//可能模型中的点数
	std::vector<double> maybe_model;	//当前找到的模型

	std::vector<int> shuffled_indices(num_data);//用于取初始样本

	while (iterations < k)
	{

		maybe_inliers_num = 0;

		//从数据集中随机选择n个点
		cur_inti_sub_set_indexs = new int[num_for_estimate];
					

		//当前找到的模型参数
		maybe_model.clear();

		//重置		
		for (int i=0;i < (int)shuffled_indices.size(); i++)
		{
			shuffled_indices[i]=i;
		}

		//随机选择两个点
		int max_index = num_data-1;
		for (int i=0; i<num_for_estimate; i++)
		{
			std::swap(shuffled_indices[i], shuffled_indices[i + rand() % (data.size() - i)]);
		}
		
		memset(cur_votes, 0, num_data*sizeof(short));

		for (int i=0; i<num_for_estimate; i++)
		{
			cur_inti_sub_set_indexs[i] = shuffled_indices[i];
			cur_votes[shuffled_indices[i]] = 1;
		}
		maybe_inliers_num = num_for_estimate;


		//查看是否已经用过这个子集
		std::pair< std::set<int *, SubSetIndexComparator >::iterator, bool > res = chosen_sub_sets.insert(cur_inti_sub_set_indexs);

		if (res.second)//true,表示插入成功，第一次用到这个子集
		{
			vector<Point2D*> exactEstimateData;
			for (int i=0; i<num_for_estimate; i++)
			{
				exactEstimateData.push_back(&(data[cur_inti_sub_set_indexs[i]]));
			}
			//根据两点得到直线方程
			param_estimator->estimate(exactEstimateData,maybe_model);

			//判定剩余的点是否符合模型
			for(int i=0; i<num_data; i++)
			{
				if(0 == cur_votes[i] && 
					param_estimator->agree(maybe_model, data[i]))
				{
					cur_votes[i] = 1;
					maybe_inliers_num++;					
				}
			}
			//比之前更好？
			if (maybe_inliers_num > best_model_num)
			{
				best_model_num = maybe_inliers_num;
				memcpy(best_votes, cur_votes, num_data*sizeof(short));
			}

			//重新计算k, k=log(1-p)/log(1-pow(w,n))
			double w = static_cast<double> (best_model_num) * one_over_indices;
			double p_no_outliers = 1.0 - std::pow(w, static_cast<double> (maybe_inliers_num));

			p_no_outliers = (std::max) (std::numeric_limits<double>::epsilon (), p_no_outliers);       // Avoid division by -Inf
			p_no_outliers = (std::min) (1.0 - std::numeric_limits<double>::epsilon (), p_no_outliers);   // Avoid division by 0.
			k = log_probability / log(p_no_outliers);
			
		}
		else
		{
			delete [] cur_inti_sub_set_indexs;
			--iterations;	//这次迭代不算数
		}

		++iterations;
	}

	//清理
	std::set<int *, SubSetIndexComparator >::iterator it = chosen_sub_sets.begin();
	std::set<int *, SubSetIndexComparator >::iterator chosenSubSetsEnd = chosen_sub_sets.end();
	while(it!=chosenSubSetsEnd) {
		delete [] (*it);
		it++;
	}
	chosen_sub_sets.clear();


	//对找到的点集，用最小二乘法重新估计模型参数
	std::vector<Point2D*> leastSquaresEstimateData;
	for(int j=0; j<num_data; j++) {
		if(best_votes[j])
			leastSquaresEstimateData.push_back(&(data[j]));
	}
	param_estimator->leastSquaresEstimate(leastSquaresEstimateData,model_parameters);

	delete []best_votes;
	delete []cur_votes;


	return (double)best_model_num/(double)num_data;
}
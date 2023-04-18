#pragma once
#include "OpticalLinearSensorModel.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

#include <Eigen/Core>
#include <Eigen/StdVector>
#include <Eigen/Geometry>
#include <Eigen/Sparse>
typedef Eigen::Triplet<double> T;

#include "g2o/types/slam3d/vertex_se3.h"
#include "g2o/types/slam3d/edge_se3.h"
#include "g2o/stuff/sampler.h"
#include "g2o/core/factory.h"

#include "g2o/core/sparse_optimizer.h"
#include "g2o/core/block_solver.h"
#include "g2o/solvers/dense/linear_solver_dense.h"
#include "g2o/core/optimization_algorithm_levenberg.h"

#include "g2o/types/sba/types_six_dof_expmap.h"



//三线阵立体模型 用于做匹配时，确定定位区域
class StereoSensorModel
{
public:
	StereoSensorModel(string sAttPathNd, string sEtPathNd, string sEpPathNd,string sInnerParaNd,string sOffsetParaNd ,int iCameraIndexNd,bool existCalibND, \
		string sAttPathS1, string sEtPathS1, string sEpPathS1, string sInnerParaS1, string sOffsetParaS1, int iCameraIndexS1, bool existCalibS1, \
		string sAttPathS2, string sEtPathS2, string sEpPathS2, string sInnerParaS2, string sOffsetParaS2, int iCameraIndexS2, bool existCalibS2, \
		string sCamerainstallPath);
	~StereoSensorModel();

	//初始化影像
	bool initImg(string sImgNd, string sImgS1, string sImgS2);

	//初始化DEM，由于三幅影像共用一套DEM，存在四次DEM读取
	bool initDem(string demfile);	

	//严密模型计算的地理区间采用vector存储四个角点的形式；在立体交会时采用起点+向量的形式
	GeoRange transtoGeoRange(vector<rpcGCP>& vGeoRange); 
	
	//计算任意两个区域的重叠
	GeoRange calOverlapArea2Img(GeoRange& model1Range, GeoRange& model2Range);
	
	//三景影像重叠区域
	GeoRange calOverlapArea3Img(); 

	//划分格网+sift匹配（两景影像） 
	//参数：range指定区域 resolution指定格网分辨率，单位：°，例如1/6°表示格网约10*10km2，等于300*300pix
	//格网会存储在m_GridsModel1和m_GridsModel2中
	bool segmentGrid(OpticalLinearSensorModel &model1, OpticalLinearSensorModel &model2, GeoRange &range, double resolution); //明确使用那两个几何定位模型
	bool gridMatch(OpticalLinearSensorModel& model1, OpticalLinearSensorModel& model2, string sTiePointsFile);//明确使用那两个几何定位模型
	//三景影像
	//格网会存储在 m_GridsModel1、m_GridsModel2和m_GridsModel3中
	bool segmentGrid3Img(GeoRange& range, double resolution);
	bool gridMatch(string sTiePointsFile);

	//读取连接点文件，注意无论是立体还是三视，都是读取三组像点坐标
	bool readtiepoints(string sTiepointsFilepath);
	
	//立体模型交会，指明使用哪两个立体定位模型
	bool fromtiepoint2XYZ(OpticalLinearSensorModel& model1, OpticalLinearSensorModel& model2, double x1, double y1, double x2, double y2, double& X, double& Y, double& Z);

	//计算所有连接点 XYZ 经纬度高程+残差,消除残差大的无匹配点
	bool caltieGeopose(bool eliminateTiepoints);
	//更新地理坐标后重新计算重投影
	bool updataResidual();
	//三个立体定位模型交会，单个连接点
	//两两交会+优化地理坐标
	tiePoint3Geo fromtiepoint2XYZ(double x1, double y1, double x2, double y2, double x3, double y3);
	//迭代优化三组地理坐标，得到唯一解
	bool caltiepointBestXYZ(Eigen::Matrix3d R1, Eigen::Vector3d XYZs1, Eigen::Vector3d XYZc1, Eigen::Matrix3d R2, Eigen::Vector3d XYZs2, Eigen::Vector3d XYZc2, Eigen::Matrix3d R3, Eigen::Vector3d XYZs3, Eigen::Vector3d XYZc3, double &X_ori, double &Y_ori, double &Z_ori);
	//输出连接点像素坐标
	bool savetiePointsPixel(string sPixelfilePath);

	//输出连接点坐标
	bool savetiePointsGeopose(string sfilePath);

	//计算单个连接点前方交会的像方残差
	bool caltieGeoresidual(OpticalLinearSensorModel& model, double x, double y, double lat, double lon, double alt, double& x_res, double& y_res);
	//输出连接点残差
	bool savetiePointsResival(string sResivalfilePath);

	//根据连接点残差,剔除误匹配点
	bool eliminat_tiePoints_residual(double x_ResidualThreshold = 100.0, double y_ResidualThreshold = 100.0);

	//外定标（三景立体连接点） 
	//定标参数为phi,omega,kappa
	bool extCalibrate();
	
	//读取DEM高程
	double getHvaluefromDEM(double lat, double lon);

public:
	OpticalLinearSensorModel m_cModelNd;	//考虑可以用数组元素 分别表示ND S1 S2相机定位模型
	OpticalLinearSensorModel m_cModelS1;
	OpticalLinearSensorModel m_cModelS2;

	vector<DRect> m_GridsModel1;	//存储格网的角点像素坐标和宽高,配合匹配使用 
	vector<DRect> m_GridsModel2;    //注意并不与S1 S2 ND对应
	vector<DRect> m_GridsModel3;    //当同时求三景重叠时，与S1 S2 ND对应

	GeoImageIO m_dem;
	TransCoordAPI m_coord;

	vector<tie3point> m_tiePointsPix; // s1 s2 nd
	vector<tiePoint3Geo> m_tiePointsGeo; // s1s2 s1nd nds2 best四类交会结果
	vector<tieResidual> m_tiePointsRes;
};


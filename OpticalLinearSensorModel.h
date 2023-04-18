#pragma once
#include <iostream>
#include <cmath>
#include <thread>
#include <Eigen/Core>
#include <Eigen/Dense>
#include "CommonFunc.h"
#include "CommonData.h"
#include "GeoImageIO.h"
#include "TransCoordAPI.h"
#include "siftgpu.h"

using namespace std;

class OpticalLinearSensorModel
{
public:

	//模型初始化：卫星辅助数据、相机参数、是否存在定标参数
	OpticalLinearSensorModel(string attfile, string etfile, string epfile, string camerafile,string innercalibrationfile,string offsetcalibrationfile, int cameraindex,bool calibratPara);

	~OpticalLinearSensorModel();

	bool initauxModel();//拟合姿轨 多项式系数
	bool initCamera(bool exist_calibratParameter);	//初始化相机安装、指向角

	//初始化影像 计算 影像的行列数、分辨率、仿射系数
	bool initImg(string imgfile);
	bool initDem(string demfile);//初始化DEM

	//计算三个相机指向角模型系数，三阶多项式来拟合每个探元位置（a0 a1 a2 a3 b0 b1 b2 b3）,并输出指向文件
	CoefPointAngle calCoePointAngle(int index = 0);

	//按照默认值（相机参数、单位阵）设置指向角参数和偏置矩阵
	bool setPointAngle();
	bool setOffsetAngle();

	//添加已存在的定标参数
	bool addOffsetAngle(string filePath);
	bool addInnerMdl(string filePath);

	//返回Ru矩阵
	Eigen::Matrix3d get_RuMatrix();

	//返回相机安装矩阵
	Eigen::Matrix3d get_InstallMatrix();

	//基于DEM迭代计算高程
	bool calHeightUsDEM(double x,double y,double &h,double dh);

	//读取DEM高程
	double getHvaluefromDEM(double lat, double lon);
	
	//像素坐标计算成像时间
	double fromxyCalET(double x, double y);
	
	//成像时间计算卫星位置
	bool frometCalpose(double et, double& pose_X, double& pose_Y, double& pose_Z);

	//成像时间计算卫星姿态
	Eigen::Matrix3d frometCalatt(double et);
	
	//坐标正算，已知高程值
	bool fromxyh2XYZ(double x, double y, double h, double& X, double& Y, double& Z);
	bool fromxyh2latlon(double x, double y, double h, double& lat, double& lon);

	//基于DEM迭代计算地理坐标，高程未知
	bool fromxy2XYZUsDEM(double x, double y, double &X, double &Y, double &Z);
	bool fromxy2latlonaltUsDEM(double x, double y, double& lat, double& lon, double& alt);

	//计算给定范围 影像四个角点的地理范围 使用DEM高程耗时，不使用高程设置为0
	vector<rpcGCP> getGeoRange(double x, double y, double width, double height, bool useDEM);

	//坐标反算
	bool fromlatlonh2xy(double lat, double lon, double h, double& x, double& y);

	//计算影像的像素范围
	vector<rpcGCP> getPixRange(GeoRange Area);

	//读取高程值
	double readHeight(double lat, double lon);

	//归一化相机模型，默认采用下视ND模型 编号0；立体相机S1模型 编号1,立体相机S2模型 编号2
	bool fromxy2XYZc(double x, double y, double& X, double& Y, double& Z, int index = 0);
	
	//指向角模型计算相机坐标
	Eigen::Vector3d fromxy2Camerapose(double x,double y);
	
	//计算影像分辨率
	double computeImgLambda();

	//计算影像仿射系数
	affineCoefficient computeImgAffineCoe();
	
	//正射纠正 整景或局部区域 有地理编码 
	//开始想的是：影像第一个波段是灰度，第二个波段是高程。但是高程值必须从DEM中读取，所以没意义了
	//基于严密几何模型，反复坐标反算，坐标正算 速度太慢了
	bool OrthoCorrection(vector<rpcGCP> GeoRange, string sOrthImgfilePath);

	//返回X' Y' Z'，用于定标
	Eigen::Vector3d fromxyXYZ_XYZpie(double x, double y, double X, double Y, double Z);

public:

	//影像参数
	GeoImageIO m_Img;
	int m_Img_Width;
	int m_Img_Height;
	double m_Img_Lambda;//分辨率 单位：度/像素
	affineCoefficient m_Img_Affine;//仿射系数

	//DEM参数
	GeoImageIO m_Dem;
	TransCoordAPI m_Dem_Coord;//DEM坐标系转换关系

	//辅助数据路径
	string m_sAttfilePath;
	string m_sEpfilePath;
	string m_sEtfilePath;
	string m_sCamerafilePath;
	string m_sCameraInnerParameterPath;
	string m_sCameraOffserParameterPath;

	//辅助数据参数
	vector<Attitude> m_att;
	vector<Orbit> m_ep;
	vector<CameraTime> m_et;
	int m_num_line;//数量
	//存储姿轨数据多项式系数（时间）
	CoeAttitude m_coeAtt;
	CoeOrbit m_coeOrbit;
	
	double m_et_ori; //起始时刻
	
	//相机参数
	Eigen::Matrix3d m_cameraInstall;
	double m_camera_f = 0.175;
	int m_camera_index;

	//ND相机参数
	double m_camera_ND_a = 0.000007;//探元尺寸
	double m_camera_ND_x0 = 2588.0;//像主点
	int m_camera_ND_numccd = 5176;//探元个数
	double m_camera_ND_y0 = 0.;//归一化相机坐标系下Y坐标
	

	//S1相机参数
	double m_camera_S1_a = 0.000014;//探元尺寸
	double m_camera_S1_x0 = 1292.0;//像主点
	int m_camera_S1_numccd = 2584;//探元个数
	double m_camera_S1_y0 = tan(18.9 * M_PI / 180.0);

	//S2相机参数
	double m_camera_S2_a = 0.000014;//探元尺寸
	double m_camera_S2_x0 = 1292.0;//像主点
	int m_camera_S2_numccd = 2584;//探元个数
	double m_camera_S2_y0 = tan(-18.9 * M_PI / 180.0);

	//内定标参数
	bool m_exist_CamPointAngle = false; //默认没有指向角信息
	CoefPointAngle m_camera_coePointAngle;//相机指向角系数
	//外定标参数
	bool m_exist_CamOffsetAngleMatrix = false;//默认没有偏置矩阵
	OffsetAngleMatrix m_camera_OffsetAngleMatrix;//相机的外定标参数

};


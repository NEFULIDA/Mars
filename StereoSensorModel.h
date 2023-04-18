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



//����������ģ�� ������ƥ��ʱ��ȷ����λ����
class StereoSensorModel
{
public:
	StereoSensorModel(string sAttPathNd, string sEtPathNd, string sEpPathNd,string sInnerParaNd,string sOffsetParaNd ,int iCameraIndexNd,bool existCalibND, \
		string sAttPathS1, string sEtPathS1, string sEpPathS1, string sInnerParaS1, string sOffsetParaS1, int iCameraIndexS1, bool existCalibS1, \
		string sAttPathS2, string sEtPathS2, string sEpPathS2, string sInnerParaS2, string sOffsetParaS2, int iCameraIndexS2, bool existCalibS2, \
		string sCamerainstallPath);
	~StereoSensorModel();

	//��ʼ��Ӱ��
	bool initImg(string sImgNd, string sImgS1, string sImgS2);

	//��ʼ��DEM����������Ӱ����һ��DEM�������Ĵ�DEM��ȡ
	bool initDem(string demfile);	

	//����ģ�ͼ���ĵ����������vector�洢�ĸ��ǵ����ʽ�������彻��ʱ�������+��������ʽ
	GeoRange transtoGeoRange(vector<rpcGCP>& vGeoRange); 
	
	//������������������ص�
	GeoRange calOverlapArea2Img(GeoRange& model1Range, GeoRange& model2Range);
	
	//����Ӱ���ص�����
	GeoRange calOverlapArea3Img(); 

	//���ָ���+siftƥ�䣨����Ӱ�� 
	//������rangeָ������ resolutionָ�������ֱ��ʣ���λ���㣬����1/6���ʾ����Լ10*10km2������300*300pix
	//������洢��m_GridsModel1��m_GridsModel2��
	bool segmentGrid(OpticalLinearSensorModel &model1, OpticalLinearSensorModel &model2, GeoRange &range, double resolution); //��ȷʹ�����������ζ�λģ��
	bool gridMatch(OpticalLinearSensorModel& model1, OpticalLinearSensorModel& model2, string sTiePointsFile);//��ȷʹ�����������ζ�λģ��
	//����Ӱ��
	//������洢�� m_GridsModel1��m_GridsModel2��m_GridsModel3��
	bool segmentGrid3Img(GeoRange& range, double resolution);
	bool gridMatch(string sTiePointsFile);

	//��ȡ���ӵ��ļ���ע�����������廹�����ӣ����Ƕ�ȡ�����������
	bool readtiepoints(string sTiepointsFilepath);
	
	//����ģ�ͽ��ᣬָ��ʹ�����������嶨λģ��
	bool fromtiepoint2XYZ(OpticalLinearSensorModel& model1, OpticalLinearSensorModel& model2, double x1, double y1, double x2, double y2, double& X, double& Y, double& Z);

	//�����������ӵ� XYZ ��γ�ȸ߳�+�в�,�����в�����ƥ���
	bool caltieGeopose(bool eliminateTiepoints);
	//���µ�����������¼�����ͶӰ
	bool updataResidual();
	//�������嶨λģ�ͽ��ᣬ�������ӵ�
	//��������+�Ż���������
	tiePoint3Geo fromtiepoint2XYZ(double x1, double y1, double x2, double y2, double x3, double y3);
	//�����Ż�����������꣬�õ�Ψһ��
	bool caltiepointBestXYZ(Eigen::Matrix3d R1, Eigen::Vector3d XYZs1, Eigen::Vector3d XYZc1, Eigen::Matrix3d R2, Eigen::Vector3d XYZs2, Eigen::Vector3d XYZc2, Eigen::Matrix3d R3, Eigen::Vector3d XYZs3, Eigen::Vector3d XYZc3, double &X_ori, double &Y_ori, double &Z_ori);
	//������ӵ���������
	bool savetiePointsPixel(string sPixelfilePath);

	//������ӵ�����
	bool savetiePointsGeopose(string sfilePath);

	//���㵥�����ӵ�ǰ��������񷽲в�
	bool caltieGeoresidual(OpticalLinearSensorModel& model, double x, double y, double lat, double lon, double alt, double& x_res, double& y_res);
	//������ӵ�в�
	bool savetiePointsResival(string sResivalfilePath);

	//�������ӵ�в�,�޳���ƥ���
	bool eliminat_tiePoints_residual(double x_ResidualThreshold = 100.0, double y_ResidualThreshold = 100.0);

	//�ⶨ�꣨�����������ӵ㣩 
	//�������Ϊphi,omega,kappa
	bool extCalibrate();
	
	//��ȡDEM�߳�
	double getHvaluefromDEM(double lat, double lon);

public:
	OpticalLinearSensorModel m_cModelNd;	//���ǿ���������Ԫ�� �ֱ��ʾND S1 S2�����λģ��
	OpticalLinearSensorModel m_cModelS1;
	OpticalLinearSensorModel m_cModelS2;

	vector<DRect> m_GridsModel1;	//�洢�����Ľǵ���������Ϳ��,���ƥ��ʹ�� 
	vector<DRect> m_GridsModel2;    //ע�Ⲣ����S1 S2 ND��Ӧ
	vector<DRect> m_GridsModel3;    //��ͬʱ�������ص�ʱ����S1 S2 ND��Ӧ

	GeoImageIO m_dem;
	TransCoordAPI m_coord;

	vector<tie3point> m_tiePointsPix; // s1 s2 nd
	vector<tiePoint3Geo> m_tiePointsGeo; // s1s2 s1nd nds2 best���ཻ����
	vector<tieResidual> m_tiePointsRes;
};


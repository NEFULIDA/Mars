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

	//ģ�ͳ�ʼ�������Ǹ������ݡ�����������Ƿ���ڶ������
	OpticalLinearSensorModel(string attfile, string etfile, string epfile, string camerafile,string innercalibrationfile,string offsetcalibrationfile, int cameraindex,bool calibratPara);

	~OpticalLinearSensorModel();

	bool initauxModel();//����˹� ����ʽϵ��
	bool initCamera(bool exist_calibratParameter);	//��ʼ�������װ��ָ���

	//��ʼ��Ӱ�� ���� Ӱ������������ֱ��ʡ�����ϵ��
	bool initImg(string imgfile);
	bool initDem(string demfile);//��ʼ��DEM

	//�����������ָ���ģ��ϵ�������׶���ʽ�����ÿ��̽Ԫλ�ã�a0 a1 a2 a3 b0 b1 b2 b3��,�����ָ���ļ�
	CoefPointAngle calCoePointAngle(int index = 0);

	//����Ĭ��ֵ�������������λ������ָ��ǲ�����ƫ�þ���
	bool setPointAngle();
	bool setOffsetAngle();

	//����Ѵ��ڵĶ������
	bool addOffsetAngle(string filePath);
	bool addInnerMdl(string filePath);

	//����Ru����
	Eigen::Matrix3d get_RuMatrix();

	//���������װ����
	Eigen::Matrix3d get_InstallMatrix();

	//����DEM��������߳�
	bool calHeightUsDEM(double x,double y,double &h,double dh);

	//��ȡDEM�߳�
	double getHvaluefromDEM(double lat, double lon);
	
	//��������������ʱ��
	double fromxyCalET(double x, double y);
	
	//����ʱ���������λ��
	bool frometCalpose(double et, double& pose_X, double& pose_Y, double& pose_Z);

	//����ʱ�����������̬
	Eigen::Matrix3d frometCalatt(double et);
	
	//�������㣬��֪�߳�ֵ
	bool fromxyh2XYZ(double x, double y, double h, double& X, double& Y, double& Z);
	bool fromxyh2latlon(double x, double y, double h, double& lat, double& lon);

	//����DEM��������������꣬�߳�δ֪
	bool fromxy2XYZUsDEM(double x, double y, double &X, double &Y, double &Z);
	bool fromxy2latlonaltUsDEM(double x, double y, double& lat, double& lon, double& alt);

	//���������Χ Ӱ���ĸ��ǵ�ĵ���Χ ʹ��DEM�̺߳�ʱ����ʹ�ø߳�����Ϊ0
	vector<rpcGCP> getGeoRange(double x, double y, double width, double height, bool useDEM);

	//���귴��
	bool fromlatlonh2xy(double lat, double lon, double h, double& x, double& y);

	//����Ӱ������ط�Χ
	vector<rpcGCP> getPixRange(GeoRange Area);

	//��ȡ�߳�ֵ
	double readHeight(double lat, double lon);

	//��һ�����ģ�ͣ�Ĭ�ϲ�������NDģ�� ���0���������S1ģ�� ���1,�������S2ģ�� ���2
	bool fromxy2XYZc(double x, double y, double& X, double& Y, double& Z, int index = 0);
	
	//ָ���ģ�ͼ����������
	Eigen::Vector3d fromxy2Camerapose(double x,double y);
	
	//����Ӱ��ֱ���
	double computeImgLambda();

	//����Ӱ�����ϵ��
	affineCoefficient computeImgAffineCoe();
	
	//������� ������ֲ����� �е������ 
	//��ʼ����ǣ�Ӱ���һ�������ǻҶȣ��ڶ��������Ǹ̡߳����Ǹ߳�ֵ�����DEM�ж�ȡ������û������
	//�������ܼ���ģ�ͣ��������귴�㣬�������� �ٶ�̫����
	bool OrthoCorrection(vector<rpcGCP> GeoRange, string sOrthImgfilePath);

	//����X' Y' Z'�����ڶ���
	Eigen::Vector3d fromxyXYZ_XYZpie(double x, double y, double X, double Y, double Z);

public:

	//Ӱ�����
	GeoImageIO m_Img;
	int m_Img_Width;
	int m_Img_Height;
	double m_Img_Lambda;//�ֱ��� ��λ����/����
	affineCoefficient m_Img_Affine;//����ϵ��

	//DEM����
	GeoImageIO m_Dem;
	TransCoordAPI m_Dem_Coord;//DEM����ϵת����ϵ

	//��������·��
	string m_sAttfilePath;
	string m_sEpfilePath;
	string m_sEtfilePath;
	string m_sCamerafilePath;
	string m_sCameraInnerParameterPath;
	string m_sCameraOffserParameterPath;

	//�������ݲ���
	vector<Attitude> m_att;
	vector<Orbit> m_ep;
	vector<CameraTime> m_et;
	int m_num_line;//����
	//�洢�˹����ݶ���ʽϵ����ʱ�䣩
	CoeAttitude m_coeAtt;
	CoeOrbit m_coeOrbit;
	
	double m_et_ori; //��ʼʱ��
	
	//�������
	Eigen::Matrix3d m_cameraInstall;
	double m_camera_f = 0.175;
	int m_camera_index;

	//ND�������
	double m_camera_ND_a = 0.000007;//̽Ԫ�ߴ�
	double m_camera_ND_x0 = 2588.0;//������
	int m_camera_ND_numccd = 5176;//̽Ԫ����
	double m_camera_ND_y0 = 0.;//��һ���������ϵ��Y����
	

	//S1�������
	double m_camera_S1_a = 0.000014;//̽Ԫ�ߴ�
	double m_camera_S1_x0 = 1292.0;//������
	int m_camera_S1_numccd = 2584;//̽Ԫ����
	double m_camera_S1_y0 = tan(18.9 * M_PI / 180.0);

	//S2�������
	double m_camera_S2_a = 0.000014;//̽Ԫ�ߴ�
	double m_camera_S2_x0 = 1292.0;//������
	int m_camera_S2_numccd = 2584;//̽Ԫ����
	double m_camera_S2_y0 = tan(-18.9 * M_PI / 180.0);

	//�ڶ������
	bool m_exist_CamPointAngle = false; //Ĭ��û��ָ�����Ϣ
	CoefPointAngle m_camera_coePointAngle;//���ָ���ϵ��
	//�ⶨ�����
	bool m_exist_CamOffsetAngleMatrix = false;//Ĭ��û��ƫ�þ���
	OffsetAngleMatrix m_camera_OffsetAngleMatrix;//������ⶨ�����

};


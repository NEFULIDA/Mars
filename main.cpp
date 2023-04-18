#pragma once
#include <iostream>
#include <omp.h>
#include <sciplot/sciplot.hpp>
#include <filesystem>

#include "gdal.h"
#include "spice.h"
#include "OpticalLinearSensorModel.h"
#include "StereoSensorModel.h"
#include "siftgpu.h"
#include "rsgeodll.h"
#include "RpcModel.h"
#include "eigenspice.h"
#include "MOLA.h"

using namespace std;
using namespace sciplot;

// ����һ HRSC�����װ����
void test_CameraInstall() 
{	
	string filecameraInstall = "E:\\date_keyanrenwu\\Mars\\MEX\\HRSC\\cameraInstall_HEAD.txt";
	MatrixcameraInstall(filecameraInstall);
}

// ���Զ� ���ݽ���
void test_Predata_MEX_HRSC()
{
	//��ʱ·�� �趨����·��
	string sEtPath = "../data/i864/it/hi864_0000_s12.txt";
	//��̬ ��� �������� �趨����·�����ҵĴ���ͨ�ø�ʽ��
	string spiceatt = "../data/i864/auxiliary/my/hi864_0000_s12.att";
	string spiceep = "../data/i864/auxiliary/my/hi864_0000_s12.eph";

	Predata_MEX_HRSC(spiceatt, spiceep, sEtPath);

	//ת�� ��׼��ʽ���·�� ���������������ͨ�ø�ʽ��
	string stdEPfile = "E:\\leeDa\\ˮ�ֺ�Ͽ�����򡪹���߶�500km����\\��������\\std\\hg733_0000_nd.eph";
	string stdETfile = "E:\\leeDa\\ˮ�ֺ�Ͽ�����򡪹���߶�500km����\\��������\\std\\hg733_0000_nd.it";
	string stdATTfile = "E:\\leeDa\\ˮ�ֺ�Ͽ�����򡪹���߶�500km����\\��������\\std\\hg733_0000_nd.att";
	string stdlineNumfile = "E:\\leeDa\\ˮ�ֺ�Ͽ�����򡪹���߶�500km����\\��������\\std\\hg733_0000_nd.txt";

	int startlineNum = 0; int endlineNum = 0;
	double startET, exposure;

	FILE* fp = fopen(sEtPath.c_str(), "r");
	char c_read[1024];
	fgets(c_read, 1024, fp);
	sscanf(c_read, "%d%lf%lf", &startlineNum, &startET, &exposure);

	while (!feof(fp)) 
	{
		fgets(c_read, 1024, fp);
		endlineNum++;
	}
	fclose(fp);

	fp = fopen(stdlineNumfile.c_str(), "w");
	fprintf(fp, "%d \t%d", startlineNum, endlineNum - 1);
	fclose(fp);

	char utcstr[1024];
	et2utc_c(startET, "C", 3, 35, utcstr);

	calcu_stdEPfile(utcstr, spiceep, stdEPfile);

	calcu_stdATTfile(utcstr, spiceatt, stdATTfile);

	calcu_stdITfile(utcstr, sEtPath, stdETfile);

}

// ������ �������������λģ�� ��DEM���������������ַ������������� 1Ӱ��·�� 2��ʱ·�� 3��������� 4��̬���ݱ���·�� 5������ݱ���·��
void test_GeoImg_DEM() 
{
	// //s1Ӱ��
	//string sImgPath = "D:\\Code\\Mars\\data\\i864\\img\\hi864_0000_s12.img";
	// //��������
	//string sEtPath = "D:\\Code\\Mars\\data\\i864\\auxiliary\\hi864_0000_s12.it";
	//string spiceatt = "D:\\Code\\Mars\\data\\i864\\auxiliary\\hi864_0000_s12.att";
	//string spiceep = "D:\\Code\\Mars\\data\\i864\\auxiliary\\hi864_0000_s12.eph";
	// //���������
	//int CameraIndex = 1;

	 //s2Ӱ��
	string sImgPath = "../data/i864/img/hi864_0000_s22.img";
	 //��������
	string sEtPath = "../data/i864/auxiliary/hi864_0000_s22.it";
	string spiceatt = "../data/i864/auxiliary/hi864_0000_s22.att";
	string spiceep = "../data/i864/auxiliary/hi864_0000_s22.eph";
	string sInnerParameterS2 = "../data/i864/calibrationParameters/hi864_0000_s22.int";
	string sOffsetParameterS2 = "../data/i864/calibrationParameters/hi864_0000_s22.ext";
	bool bExist_calibrationS2 = false;
	//���������
	int CameraIndex = 2;

	//�����װ
	string sCamerainstallPath = "../data/camera/cameraInstall_HEAD.txt";

	//�������������ģ��ʵ����
	OpticalLinearSensorModel cModel(spiceatt, sEtPath, spiceep, sCamerainstallPath, sInnerParameterS2, sOffsetParameterS2, CameraIndex, bExist_calibrationS2);

	//��ʼ��
	string sDEMPath = "../data/mola/megr00n270hb.lbl"; //DEMͶӰ���ݣ�����ͶӰ����
	cModel.initDem(sDEMPath);
	cModel.initImg(sImgPath);
	
	//������������
	/*
	vector<rpcGCP> points;
	FILE* fp = fopen("../data/i864/match/Gcp_Img_s2.txt", "r");
	FILE* fp_write = fopen("../data/i864/match/Gcp_Img_s2_toGeo.txt", "w");
	fprintf(fp_write, "\n");

	int time_index = 1;
	char c_read[1024];
	while (!feof(fp)) {
		if (fgets(c_read, 1024, fp) == NULL)
			continue;
		rpcGCP point;
		if (fscanf(fp, "%lf %lf %lf %lf %lf", &point.x, &point.y, &point.lat, &point.lon, &point.h) != 5)
			continue;
		cModel.fromxy2latlonaltUsDEM(point.x, point.y, point.lat, point.lon, point.h);
		fprintf(fp_write, "%lf\t%lf\t%lf\t%lf\t%lf\n", point.x, point.y, point.lat, point.lon, point.h);
		points.push_back(point);
		cout << time_index << endl;
		time_index++;
	}
	fclose(fp);
	fclose(fp_write);
	*/

	//�������귴��
	/*
	FILE* fp = fopen("../data/i864/match/Gcp_Img_s2_toGeo.txt", "r");
	FILE* fp_write = fopen("../data/i864/match/Gcp_Img_s2_topixel.txt", "w");
	fprintf(fp_write, "\n");
	int time_index = 1;

	char c_read[1024];
	while (!feof(fp)) {
		if (fgets(c_read, 1024, fp) == NULL)
			continue;
		rpcGCP point;
		if (fscanf(fp, "%lf %lf %lf %lf %lf", &point.x, &point.y, &point.lat, &point.lon, &point.h) != 5)	//sample line γ�� ���� �߳�  
			continue;
		cModel.fromlatlonh2xy(point.lat, point.lon, point.h, point.x, point.y);
		fprintf(fp_write, "%lf\t%lf\t%lf\t%lf\t%lf\n", point.x, point.y, point.lat, point.lon, point.h);
		cout << time_index << endl;
		time_index++;
	}
	fclose(fp);
	fclose(fp_write);
	*/
	
	//������� ������ָ���ֲ���Χ��
	/*
	vector<rpcGCP> GeoRange = cModel.getGeoRange(500, 500, 400, 400, true);//ȫ����Χ
	string sOrthImgfilePath = "../data/i864/orth/hi864_0000_s22.tif";
	cModel.OrthoCorrection(GeoRange, sOrthImgfilePath);
	*/
}

// ������ �������彻��ģ��
void test_GeoStereo() {
	
	string sAttPathNd = "../data/i864/auxiliary/hi864_0000_nd2.att";
	string sEtPathNd = "../data/i864/auxiliary/hi864_0000_nd2.it";
	string sEpPathNd = "../data/i864/auxiliary/hi864_0000_nd2.eph";
	string sInnerParameterNd = "../data/i864/calibrationParameters/hi864_0000_nd2.int";
	string sOffsetParameterNd = "../data/i864/calibrationParameters/hi864_0000_nd2.ext";
	bool bExist_calibrationNd = true;
	int iCameraIndexNd = 0;

	string sAttPathS1 = "../data/i864/auxiliary/hi864_0000_s12.att";
	string sEtPathS1 = "../data/i864/auxiliary/hi864_0000_s12.it"; 
	string sEpPathS1 = "../data/i864/auxiliary/hi864_0000_s12.eph";
	string sInnerParameterS1 = "../data/i864/calibrationParameters/hi864_0000_s12.int";
	string sOffsetParameterS1 = "../data/i864/calibrationParameters/hi864_0000_s12.ext";
	bool bExist_calibrationS1 = true;
	int iCameraIndexS1 = 1;

	string sAttPathS2 = "../data/i864/auxiliary/hi864_0000_s22.att"; 
	string sEtPathS2 = "../data/i864/auxiliary/hi864_0000_s22.it"; 
	string sEpPathS2 = "../data/i864/auxiliary/hi864_0000_s22.eph";
	string sInnerParameterS2 = "../data/i864/calibrationParameters/hi864_0000_s22.int";
	string sOffsetParameterS2 = "../data/i864/calibrationParameters/hi864_0000_s22.ext";
	bool bExist_calibrationS2 = true;
	int iCameraIndexS2 = 2;

	string sCamerainstallPath = "../data/camera/cameraInstall_HEAD.txt";

	StereoSensorModel cModel(sAttPathNd, sEtPathNd, sEpPathNd, sInnerParameterNd, sOffsetParameterNd, iCameraIndexNd, bExist_calibrationNd, \
		sAttPathS1, sEtPathS1, sEpPathS1, sInnerParameterS1, sOffsetParameterS1,iCameraIndexS1, bExist_calibrationS1, \
		sAttPathS2, sEtPathS2, sEpPathS2, sInnerParameterS2, sOffsetParameterS2, iCameraIndexS2, bExist_calibrationS2, \
		sCamerainstallPath);

	string sDEMPath = "../data/mola/megr00n270hb.lbl";
	cModel.initDem(sDEMPath);

	string sImgNd = "../data/i864/img/hi864_0000_nd2.img"; 
	string sImgS1 = "../data/i864/img/hi864_0000_s12.img"; 
	string sImgS2 = "../data/i864/img/hi864_0000_s22.img";

	cModel.initImg(sImgNd, sImgS1, sImgS2);

	//ƥ������������ص����� �������� �ֿ�ƥ��
	/*
	GeoRange overRange = cModel.calOverlapArea3Img();//�ص�����
	double resolution = 1.0 / 6.0;//�ֱ���1/6��=10*10km2=300*300pix 
	cModel.segmentGrid3Img(overRange, resolution); //��������
	string tiepointfilepath = "../data/i864/match/s1s2nd.pts";
	cModel.gridMatch(tiepointfilepath); //ƥ�����ӵ�
	*/

	//���彻��+�в��޳�
	string tiepointfilepath = "../data/i864/match/s1s2nd.pts";
	cModel.readtiepoints(tiepointfilepath);
	cModel.caltieGeopose(true);
	cModel.savetiePointsPixel("../data/i864/match/");

	//����ǰ��λ���+�в�
	/*string reportSavefile = "../data/i864/calibrationParameters/before/";
	cModel.savetiePointsGeopose(reportSavefile);
	cModel.savetiePointsResival(reportSavefile);*/
	
	//�ⶨ��
	cModel.extCalibrate();
	string reportSavefile_after = "../data/i864/calibrationParameters/after/";
	cModel.updataResidual();
	cModel.savetiePointsGeopose(reportSavefile_after);
	cModel.savetiePointsResival(reportSavefile_after);

	//�ڶ���
}

// ������ ������������ ����ƽ����
void test_MOLA() {

	//Ԥ����Щ�������DEM
	//predictOrbitNum();
	//���ƹ����������
	//drawSatalliteorbit(molafile, "..\\..\\..\\date_keyanrenwu\\Mars\\MGS\\MOLA\\PEDR\\i864_mola\\output\\10819.pdf");
	
	//string demfile = "E:\\IMG_i863_controlPoints\\i863_mola\\run-DEM.tif";//���� ASP������DEM����
	//string molafile = "E:\\IMG_i863_controlPoints\\i863_mola\\AP10652L.TAB";//���� MOLA����

	//string InitGcpfile = "E:\\IMG_i863_controlPoints\\i863_mola\\��ͼ���\\AP10652L.txt";//��� 
	//readMolaGpoints(molafile, InitGcpfile, demfile);//������ȡMOLA�ľ�γ�ȸ̣߳��Լ���Ӧλ�õ�DEM�̺߳���������

	//string demfile = "E:\\IMG_i863_controlPoints\\i863_mola\\run-DEM.tif";//���� ASP������DEM����
	//string InitGcpfile = "E:\\IMG_i863_controlPoints\\i863_mola\\��ͼ���\\AP1011.txt";//����

	//string slopePoints = "E:\\IMG_i863_controlPoints\\i863_mola\\��ͼ���\\slopePoints.txt";
	//string Nccfile = "E:\\IMG_i863_controlPoints\\i863_mola\\��ͼ���\\offsetH.dat";
	//string BestGcpfile = "E:\\IMG_i863_controlPoints\\i863_mola\\��ͼ���\\BestH.dat";

	//calNccOffset(InitGcpfile, demfile, slopePoints, Nccfile);//����һ�����ϵ����ѵ�ƫ��λ��
	//calBestOffset(InitGcpfile, Nccfile, demfile, BestGcpfile);//���յ�MOLA���ݺ�DEM����

	//i864���Ƶ����
	//string demfile = "E:\\IMG_i864_All\\TiePoints_calibaration\\calibration\\sensorCorrect_result\\run-DEM.tif";//���� ASP������DEM����
	////string demfile = "E:\\IMG_i864_All\\dem_calibration\\test\\sensorCorrect_result\\run-DEM.tif";//���� ASP������DEM����(DEMԼ������)
	//
	//string molafile = "E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\PEDR\\i864_mola\\AP10228L.TAB";//���� MOLA����

	//string InitGcpfile = "E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\PEDR\\i864_mola\\output\\AP10228L.txt";//��� 
	//readMolaGpoints(molafile, InitGcpfile, demfile);//������ȡMOLA�ľ�γ�ȸ̣߳��Լ���Ӧλ�õ�DEM�̺߳���������

	string InitGcpfile = "E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\PEDR\\i864_mola\\output\\AP10228L.txt";
	string fileCritePoints = "E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\PEDR\\i864_mola\\output\\critePoints\\AP10228L.txt";
	saveCritePoints(InitGcpfile, fileCritePoints);

	//string demfile = "E:\\IMG_i864_All\\TiePoints_calibaration\\calibration\\sensorCorrect_result\\run-DEM.tif";//���� ASP������DEM����
	//string demfile = "E:\\IMG_i864_All\\dem_calibration\\test\\sensorCorrect_result\\run-DEM.tif";//���� ASP������DEM����(DEMԼ������)
	//string InitGcpfile = "E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\PEDR\\i864_mola\\output\\AP10228L.txt";//����

	//string slopePoints = "E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\PEDR\\i864_mola\\output\\slopePoints.txt";
	//string BestGcpfile = "E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\PEDR\\i864_mola\\output\\bestH\\AP10228_best";

	//string Nccfile = "E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\PEDR\\i864_mola\\output\\offset\\AP10228_offset";
	//calNccOffset(InitGcpfile, demfile, slopePoints, Nccfile, BestGcpfile);//����һ�����ϵ����ѵ�ƫ��λ��

	//string Nccfile = "E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\PEDR\\i864_mola\\output\\offset\\offset1.txt";
	//calBestOffset(InitGcpfile, Nccfile, demfile, BestGcpfile);//���յ�MOLA���ݺ�DEM����
	
}

// ������ dem�������� ����ƽ����
void test_MOLA_grid() {

	string NCCfile = "E:\\IMG_i864_All\\molaGridMatch\\ncc\\";
	string BestHfile = "E:\\IMG_i864_All\\molaGridMatch\\bestH\\";
	
	//����DEM��Χ�޶�MOLA��������
	string MOLA = "E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\MEG128\\radius\\00n270value.tif";//���Ƿ�Χ��,ͶӰ���룬��Чֵ-32768.0
	string DEM = "E:\\IMG_i864_All\\TiePoints_calibaration\\calibration\\sensorCorrect_result\\run-DEM.tif";
	GeoImageIO dem;
	dem.open(DEM.c_str());
	double demTrans[4];
	dem.getGeoTransform(demTrans);
	int demWid = dem.getImageWid();
	int demHei = dem.getImageHei();
	double ori_Dem_lat = demTrans[2], ori_Dem_lon = demTrans[0];
	double end_Dem_lat = demTrans[2] + demTrans[3] * demHei, end_Dem_lon = demTrans[0] + demTrans[1] * demWid;
	dem.readToBuffer_proj(0, ori_Dem_lon, end_Dem_lon, ori_Dem_lat, end_Dem_lat);

	GeoImageIO mola;
	mola.open(MOLA.c_str());
	double molaTrans[4];
	mola.getGeoTransform(molaTrans);

	//����DEM��MOLA�ص�������MOLA����ϵ�µ�����
	TransCoordAPI coord;
	//coord.init(dem.getProjectionRef().c_str(), GcsMarsDTMWTK.c_str());//��֤�˶��߱任һ��
	// �����һ��ʹ�ý��ȴ���������
	//coord.init(dem.getProjectionRef().c_str(), mola.getProjectionRef().c_str());
	// ���ֻ������һ��
	coord.init(GcsMarsDTMWTK.c_str(), mola.getProjectionRef().c_str());
	cout <<"DEMͶӰ��Ϣ��\n" << dem.getProjectionRef() << endl;
	cout <<"���������Ϣ��\n" << GcsMarsDTMWTK << endl;
	
	double ori_Mola_west, ori_Mola_north;
	double end_Mola_east, end_Mola_south;
	
	//ע����ΪҪ�ٽ�����ƽ����������mola����ӦС��dem���� dem�ֱ���*���=35
	//����dem������Χ 160*160 ����
	const int gridnumX = 160;//lon
	const int gridnumY = 40;//lat
	coord.transCoord(ori_Dem_lat + gridnumY * demTrans[3], ori_Dem_lon + gridnumX * demTrans[1], ori_Mola_north, ori_Mola_west);
	coord.transCoord(end_Dem_lat - gridnumY * demTrans[3], end_Dem_lon - gridnumX * demTrans[1], end_Mola_south, end_Mola_east);
	mola.readToBuffer_proj(0, ori_Mola_west, end_Mola_east, ori_Mola_north, end_Mola_south);

	//�ص�����mola������������
	int molaSample, molaLine;
	molaSample = (end_Mola_east - ori_Mola_west) / molaTrans[1];
	molaLine = (end_Mola_south - ori_Mola_north) / molaTrans[3];

	//����MOLA���굽��������
	TransCoordAPI coord_geo;
	coord_geo.init(mola.getProjectionRef().c_str(), GcsMarsDTMWTK.c_str());

	//��MOLA�����л��ָ���10*10��������Χ�Ը������ĵ���Χ200*200
	int molaGrid = 10;
	int molaGridNums_wid, molaGridNums_hei;
	molaGridNums_wid = molaSample / molaGrid;
	molaGridNums_hei = molaLine / molaGrid;
	cout << "��mola���ݻ���Ϊ�� " << molaGridNums_hei << " �У� " << molaGridNums_wid << " �еĸ���" << endl;

	//�������¶ȵĵ㼯
	FILE* fp_mola = fopen("E:\\IMG_i864_All\\molaGridMatch\\slope\\slope.txt", "w");
	fprintf(fp_mola, "# MOLA���� MOLAγ�� MOLA�߳� DEM�߳� DEM_sample DEM_line\n");

	string bestfile_all = BestHfile + "all.txt";
	FILE* fpBestH_all = fopen(bestfile_all.c_str(), "w");
	fprintf(fpBestH_all, "# MOLA���� MOLAγ�� MOLA�߳� DEM���� DEMγ�� DEM�߳� DEM_sample DEM_line DEM��ʼ�߳� \n");

	//ѭ��ÿһ��MOLA����
	for (int i = 0; i < molaGridNums_hei; i++)
	{
		for (int j = 0; j < molaGridNums_wid; j++)
		{
			cout << "ִ�еڣ� " << i << " �У� ��" << j << " �е�MOLA����" << endl;
			//��������mola��
			vector<GeoPoint> molaGeoPoints;
			for (int m = 0; m < molaGrid; m++)
			{
				for (int n = 0; n < molaGrid; n++)
				{
					GeoPoint molapoint;
					double x = ori_Mola_west + (j * molaGrid + n) * molaTrans[1];
					double y = ori_Mola_north + (i * molaGrid + m) * molaTrans[3];
					mola.getBufferValue(x,y,2,&molapoint.h);
					if (molapoint.h)
					{
						coord_geo.transCoord(y, x, molapoint.lat, molapoint.lon);
						molaGeoPoints.push_back(molapoint);
					}
				}
			}

			//����̷߳���
			double variance = 0.;
			double means = 0.;
			for (int t = 0; t < molaGeoPoints.size(); t++)
			{
				means += molaGeoPoints[t].h / molaGeoPoints.size();
			}
			for (int t = 0; t < molaGeoPoints.size(); t++)
			{
				variance += pow(molaGeoPoints[t].h - means, 2) / molaGeoPoints.size();
			}
			cout << "����Ϊ�� " << variance << endl;

			//����������
			if (variance>1000)
			{
				for (int t = 0; t < molaGeoPoints.size(); t++)
				{
					fprintf(fp_mola, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", molaGeoPoints[t].lon, molaGeoPoints[t].lat, molaGeoPoints[t].h, 0., 0., 0.);
				}

				//����ƫ��
				double offsetH[gridnumY][gridnumX];//�洢ÿ��ƫ�Ƶ�ĸ߳�ƫ���Ƥ��ѷϵ��
				double coefficientNcc = 0.0;
				double bestCoefficientNcc = 0.0;
				double sumCoefficientNcc = 0.0;
				int xOffset = 0;
				int yOffset = 0;
				for (int a = -1 * gridnumY / 2; a < gridnumY / 2; a++) //Y�� γ��ƫ�� -80��79
				{
					for (int b = -1 * gridnumX / 2; b < gridnumX / 2; b++) //X�� ����ƫ�� -80��79
					{
						vector<double> DEMvalue, MOLAvalue;//��������߳�
						for (int k = 0; k < molaGeoPoints.size(); k++)
						{
							double h_value = 0.;
							dem.getBufferValue(molaGeoPoints[k].lon + b * demTrans[1], molaGeoPoints[k].lat + a * demTrans[3], 2, &h_value);
							if (h_value) //��0���������㣻
							{
								DEMvalue.push_back(h_value);
								MOLAvalue.push_back(molaGeoPoints[k].h);
							}
						}

						//���ƥ��ϵ��
						coefficientNcc = computerNCC(DEMvalue, MOLAvalue); //�������ϵ��(-1��1֮��)
						sumCoefficientNcc += (coefficientNcc > 0.0 ? coefficientNcc : 0.0); // �ۼ�ϵ����
						if (coefficientNcc > bestCoefficientNcc)
						{
							bestCoefficientNcc = coefficientNcc;
							xOffset = b;
							yOffset = a;
						}

						//�洢ϵ��ֵ
						offsetH[a + gridnumY / 2][b + gridnumX / 2] = coefficientNcc;

					}
				}
				if (sumCoefficientNcc)
				{
					double threshold = bestCoefficientNcc * gridnumX * gridnumY / sumCoefficientNcc;
					cout << "���ƥ��ϵ����" << bestCoefficientNcc << endl;
					cout << "��ֵ��" << threshold << endl;
					int index = j + i * molaGridNums_wid;

					if (bestCoefficientNcc > 0.0 && threshold > 2.0)
					{

						//�洢 ƫ���������ڵĸ߳�ƫ��ͼ dat�ļ�
						string nccfile = NCCfile + to_string(index) + ".txt";
						FILE* fpoffsetH = fopen(nccfile.c_str(), "w");
						fprintf(fpoffsetH, "# γ��ƫ�� ����ƫ�� NCCһ����\n");
						for (int i = 0; i < gridnumY; i++)
						{
							for (int j = 0; j < gridnumX; j++)
							{
								fprintf(fpoffsetH, "%d\t%d\t%lf\n", i - gridnumY / 2, j - gridnumX / 2, offsetH[i][j]);
							}
							fprintf(fpoffsetH, "\n");
						}
						fclose(fpoffsetH);

						//�洢ƫ�ƺ������
						string bestfile = BestHfile + to_string(index) + ".txt";
						FILE* fpBestH = fopen(bestfile.c_str(), "w");
						fprintf(fpBestH, "# MOLA���� MOLAγ�� MOLA�߳� DEM���� DEMγ�� DEM�߳� DEM_sample DEM_line DEM��ʼ�߳� \n");
						for (int k = 0; k < molaGeoPoints.size(); k++)
						{
							double h_value = 0.;
							double lon = molaGeoPoints[k].lon + xOffset * demTrans[1];
							double lat = molaGeoPoints[k].lat + yOffset * demTrans[3];
							dem.getBufferValue(lon, lat, 2, &h_value);

							if (h_value)
							{
								double sample = 0.0;
								double line = 0.0;
								dem.getBufferPixel(lon, lat, sample, line);
								fprintf(fpBestH, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", molaGeoPoints[k].lon, molaGeoPoints[k].lat, molaGeoPoints[k].h, lon, lat, h_value, sample, line, 0);
								fprintf(fpBestH_all, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", molaGeoPoints[k].lon, molaGeoPoints[k].lat, molaGeoPoints[k].h, lon, lat, h_value, sample, line, 0);
							}
						}
						fclose(fpBestH);

					}
				}
				else
				{
					cout << "���ƥ��ϵ����" << bestCoefficientNcc << endl;
				}

			}

		}
	}
	fclose(fp_mola);
	fclose(fpBestH_all);

	//��ȡMOLA�����ڸ߳�ֵ������DEM��Ӧ�ĸ߳�ֵ

	//�̵߳㼯ƥ�䣬�������ϵ����ƽ����

	//�޳������Ե͵ĵ㼯�Լ������Բ����Եĵ㼯
}

// ������ ƴ��ȫ��DEM��ת��lbl��ʽΪtif��ִ��һ�鼴�ɣ�������ִ�У�
void test_DEM() {
	//ƴ��ȫ��DEM��ת��lbl��ʽΪtif
	//read_Img_Information("E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\MEG128\\radius\\megr00n090hb.lbl");
	GeoImageIO mars_Dem;
	string mars_Dem_file = "E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\MEG128\\radius\\mars_dem_global.tif";
	int mars_Dem_wid, mars_Dem_hei;
	string mars_Dem_wkt;
	GeoImageIO mars_dem_megr88n000hb;
	mars_dem_megr88n000hb.open("E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\MEG128\\radius\\megr88n000hb.lbl");
	mars_Dem_wid = mars_dem_megr88n000hb.getImageWid() * 4;
	mars_Dem_hei = mars_dem_megr88n000hb.getImageHei() * 4;
	mars_Dem_wkt = mars_dem_megr88n000hb.getProjectionRef();
	double geoTransform_megr88n000hb[4];
	mars_dem_megr88n000hb.getGeoTransform(geoTransform_megr88n000hb);
	double desTransform[] = { geoTransform_megr88n000hb[0],geoTransform_megr88n000hb[1],0,geoTransform_megr88n000hb[2],0.0,geoTransform_megr88n000hb[3] };
	mars_Dem.create(mars_Dem_file.c_str(), "GTiff", GDT_Int16, mars_Dem_wid, mars_Dem_hei, 1, desTransform, mars_Dem_wkt);

	int wid_single = mars_dem_megr88n000hb.getImageWid();
	int hei_single = mars_dem_megr88n000hb.getImageHei();

	//��һ��
	void* pData = new int16_t[__int64(wid_single) * hei_single];
	mars_dem_megr88n000hb.readBlock(0, 0, wid_single, hei_single, 0, pData, GDT_Int16);
	mars_Dem.writeBlock(0, 0, wid_single, hei_single, 0, pData);

	//�ڶ���
	GeoImageIO mars_dem_megr88n090hb;
	mars_dem_megr88n090hb.open("E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\MEG128\\radius\\megr88n090hb.lbl");
	mars_dem_megr88n090hb.readBlock(0, 0, wid_single, hei_single, 0, pData, GDT_Int16);
	mars_Dem.writeBlock(wid_single - 1, 0, wid_single, hei_single, 0, pData);

	//������
	GeoImageIO mars_dem_megr88n180hb;
	mars_dem_megr88n180hb.open("E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\MEG128\\radius\\megr88n180hb.lbl");
	mars_dem_megr88n180hb.readBlock(0, 0, wid_single, hei_single, 0, pData, GDT_Int16);
	mars_Dem.writeBlock(wid_single * 2 - 1, 0, wid_single, hei_single, 0, pData);

	//���Ŀ�
	GeoImageIO mars_dem_megr88n270hb;
	mars_dem_megr88n270hb.open("E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\MEG128\\radius\\megr88n270hb.lbl");
	mars_dem_megr88n270hb.readBlock(0, 0, wid_single, hei_single, 0, pData, GDT_Int16);
	mars_Dem.writeBlock(wid_single * 3 - 1, 0, wid_single, hei_single, 0, pData);

	//�����
	GeoImageIO mars_dem_megr44n000hb;
	mars_dem_megr44n000hb.open("E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\MEG128\\radius\\megr44n000hb.lbl");
	mars_dem_megr44n000hb.readBlock(0, 0, wid_single, hei_single, 0, pData, GDT_Int16);
	mars_Dem.writeBlock(0, hei_single - 1, wid_single, hei_single, 0, pData);

	//������
	GeoImageIO mars_dem_megr44n090hb;
	mars_dem_megr44n090hb.open("E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\MEG128\\radius\\megr44n090hb.lbl");
	mars_dem_megr44n090hb.readBlock(0, 0, wid_single, hei_single, 0, pData, GDT_Int16);
	mars_Dem.writeBlock(wid_single - 1, hei_single - 1, wid_single, hei_single, 0, pData);

	//���߿�
	GeoImageIO mars_dem_megr44n180hb;
	mars_dem_megr44n180hb.open("E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\MEG128\\radius\\megr44n180hb.lbl");
	mars_dem_megr44n180hb.readBlock(0, 0, wid_single, hei_single, 0, pData, GDT_Int16);
	mars_Dem.writeBlock(wid_single * 2 - 1, hei_single - 1, wid_single, hei_single, 0, pData);

	//�ڰ˿�
	GeoImageIO mars_dem_megr44n270hb;
	mars_dem_megr44n270hb.open("E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\MEG128\\radius\\megr44n270hb.lbl");
	mars_dem_megr44n270hb.readBlock(0, 0, wid_single, hei_single, 0, pData, GDT_Int16);
	mars_Dem.writeBlock(wid_single * 3 - 1, hei_single - 1, wid_single, hei_single, 0, pData);

	//�ھſ�
	GeoImageIO mars_dem_megr00n000hb;
	mars_dem_megr00n000hb.open("E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\MEG128\\radius\\megr00n000hb.lbl");
	mars_dem_megr00n000hb.readBlock(0, 0, wid_single, hei_single, 0, pData, GDT_Int16);
	mars_Dem.writeBlock(0, hei_single * 2 - 1, wid_single, hei_single, 0, pData);

	//��ʮ��
	GeoImageIO mars_dem_megr00n090hb;
	mars_dem_megr00n090hb.open("E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\MEG128\\radius\\megr00n090hb.lbl");
	mars_dem_megr00n090hb.readBlock(0, 0, wid_single, hei_single, 0, pData, GDT_Int16);
	mars_Dem.writeBlock(wid_single - 1, hei_single * 2 - 1, wid_single, hei_single, 0, pData);

	//��ʮһ��
	GeoImageIO mars_dem_megr00n180hb;
	mars_dem_megr00n180hb.open("E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\MEG128\\radius\\megr00n180hb.lbl");
	mars_dem_megr00n180hb.readBlock(0, 0, wid_single, hei_single, 0, pData, GDT_Int16);
	mars_Dem.writeBlock(wid_single * 2 - 1, hei_single * 2 - 1, wid_single, hei_single, 0, pData);

	//��ʮ����
	GeoImageIO mars_dem_megr00n270hb;
	mars_dem_megr00n270hb.open("E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\MEG128\\radius\\megr00n270hb.lbl");
	mars_dem_megr00n270hb.readBlock(0, 0, wid_single, hei_single, 0, pData, GDT_Int16);
	mars_Dem.writeBlock(wid_single * 3 - 1, hei_single * 2 - 1, wid_single, hei_single, 0, pData);

	//��ʮ����
	GeoImageIO mars_dem_megr44s000hb;
	mars_dem_megr44s000hb.open("E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\MEG128\\radius\\megr44s000hb.lbl");
	mars_dem_megr44s000hb.readBlock(0, 0, wid_single, hei_single, 0, pData, GDT_Int16);
	mars_Dem.writeBlock(0, hei_single * 3 - 1, wid_single, hei_single, 0, pData);

	//��ʮ�Ŀ�
	GeoImageIO mars_dem_megr44s090hb;
	mars_dem_megr44s090hb.open("E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\MEG128\\radius\\megr44s090hb.lbl");
	mars_dem_megr44s090hb.readBlock(0, 0, wid_single, hei_single, 0, pData, GDT_Int16);
	mars_Dem.writeBlock(wid_single - 1, hei_single * 3 - 1, wid_single, hei_single, 0, pData);

	//��ʮ���
	GeoImageIO mars_dem_megr44s180hb;
	mars_dem_megr44s180hb.open("E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\MEG128\\radius\\megr44s180hb.lbl");
	mars_dem_megr44s180hb.readBlock(0, 0, wid_single, hei_single, 0, pData, GDT_Int16);
	mars_Dem.writeBlock(wid_single * 2 - 1, hei_single * 3 - 1, wid_single, hei_single, 0, pData);

	//��ʮ����
	GeoImageIO mars_dem_megr44s270hb;
	mars_dem_megr44s270hb.open("E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\MEG128\\radius\\megr44s270hb.lbl");
	mars_dem_megr44s270hb.readBlock(0, 0, wid_single, hei_single, 0, pData, GDT_Int16);
	mars_Dem.writeBlock(wid_single * 3 - 1, hei_single * 3 - 1, wid_single, hei_single, 0, pData);

	delete[]pData, pData = nullptr;
	mars_Dem.destroy();
}

// ���԰� ��ȡָ��λ�ø߳�
void test_readGeoH() {
	string demFile = "../data/mola/megr00n270hb.lbl";
	string in_lonlatFile = "../data/mola/s1s2nd_Geo.pts"; 
	string out_lonlataltFile = "../data/mola/i864_alt.txt";

	// ��ȡ�ļ��ľ��ȣ�γ��
	FILE* fp = fopen(in_lonlatFile.c_str(), "r");
	vector<pair<double, double>> orbit;
	char c_read[1024];
	while (!feof(fp)) {
		if (fgets(c_read, 1024, fp) == NULL)
			continue;
		double a1, a2;
		if (fscanf(fp, "%lf%lf", &a1, &a2) != 2)	//1γ��2����
			continue;
		orbit.push_back(make_pair(a1, a2));
	}
	fclose(fp);

	// ��ȡDEM�߳�ֵ
	GeoImageIO dem;
	vector<double> h_dem;
	dem.open(demFile);
	//��ȡӰ�񸲸Ƿ�Χ
	double leftlon, rightlon, uplat, downlat;
	double trans[4];
	dem.getGeoTransform(trans);
	leftlon = trans[0];
	uplat = trans[2];
	rightlon = trans[0] + trans[1] * dem.getImageWid();
	downlat = trans[2] + trans[3] * dem.getImageHei();
	dem.readToBuffer_proj(0, leftlon, rightlon, uplat, downlat);
	
	TransCoordAPI coord;
	coord.init(GcsMarsDTMWTK.c_str(), dem.getProjectionRef().c_str());

	for (size_t i = 0; i < orbit.size(); i++)
	{
		double north = 0.; double east = 0.;
		coord.transCoord(orbit[i].first, orbit[i].second, north, east);
		double h_value = 0.;
		dem.getBufferValue(east, north, 0, &h_value);

		h_dem.push_back(h_value);
	}

	//�洢�߳�ֵ
	FILE* fpwrite = fopen(out_lonlataltFile.c_str(), "w");
	fprintf(fpwrite, "#γ�� ���� �߳�\n");
	for (size_t i = 0; i < orbit.size(); i++)
	{
		fprintf(fpwrite, "%lf\t%lf\t%lf\n", orbit[i].first, orbit[i].second, h_dem[i]);
	}
	fclose(fpwrite);
}

void test_IMG_Mola() {
	string demFile = "E:\\IMG_I977\\test2\\sensorCorrect_result\\run-DEM.tif";//DEMӰ��
	string MolaFile = "E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\PEDR\\i977_mola\\input\\MolaLists.txt";//MOLA�����б�
	string filePoints_out = "E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\PEDR\\i977_mola\\output\\MolaPoints.txt";//MOLA Points ����γ�ȸ߳���ʽ������·��
	string filePoints_out_xyz = "E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\PEDR\\i818_mola\\output\\MolaPoints_xyz.txt";//MOLA Points ��XYZ��ʽ������·��
	read_IMG_Mola(demFile, MolaFile, filePoints_out, filePoints_out_xyz);
}

// ���Ծ� �Զ���ȡ���彻��������DEM�з���ϴ������
void test_dem_slopePoints() {
	// �����С�޶�100*100����
    // ��������10������
	string DEM = "E:\\IMG_i864_All\\TiePoints_calibaration\\calibration\\sensorCorrect_result\\run-DEM.tif";
	GeoImageIO dem;
	dem.open(DEM.c_str());
	double demTrans[4];
	dem.getGeoTransform(demTrans);
	int demWid = dem.getImageWid();
	int demHei = dem.getImageHei();
	double ori_Dem_lat = demTrans[2], ori_Dem_lon = demTrans[0];
	double end_Dem_lat = demTrans[2] + demTrans[3] * demHei, end_Dem_lon = demTrans[0] + demTrans[1] * demWid;
	dem.readToBuffer_proj(0, ori_Dem_lon, end_Dem_lon, ori_Dem_lat, end_Dem_lat);

	int windows_size = 100;
	int slip_size = 10;

	vector<rpcGCP> GeoWindow;//������������
	vector<double> GeoWindow_variance;//�������򷽲�
	//ȫͼѭ��
	for (int line_index = 0; line_index + windows_size + slip_size < demHei; line_index += slip_size)
	{
		for (int sample_index = 0; sample_index + windows_size + slip_size < demWid; sample_index += slip_size)
		{
			//�������ͳ�Ƹ߳�ֵ
			rpcGCP GeoPoints;//��������
			GeoPoints.lat = ori_Dem_lat + (line_index + windows_size/2) * demTrans[3];
			GeoPoints.lon = ori_Dem_lon + (sample_index + windows_size/2) * demTrans[1];
			dem.getBufferValue(GeoPoints.lon, GeoPoints.lat, 0, &GeoPoints.h);
			GeoPoints.y = line_index + windows_size / 2;
			GeoPoints.x = sample_index + windows_size / 2;

			vector<double> H_window;
			for (int i = 0; i < windows_size; i++)
			{
				for (int j = 0; j < windows_size; j++)
				{
					double lon = ori_Dem_lon + (sample_index + j) * demTrans[1];
					double lat = ori_Dem_lat + (line_index + i) * demTrans[3];
					double h = 0.;
					dem.getBufferValue(lon, lat, 0, &h);
					if (h > -100000.0)//�޳���Чֵ
					{
						H_window.push_back(h);
					}
				}
			}

			if (H_window.size())//�߳�ֵ��Ч��ִ��
			{
				//����̷߳���
				double variance = 0.;
				double means = 0.;
				for (int t = 0; t < H_window.size(); t++)
				{
					means += H_window[t] / H_window.size();
				}
				for (int t = 0; t < H_window.size(); t++)
				{
					variance += (pow(H_window[t] - means, 2) / H_window.size());
				}
				//cout << "����Ϊ�� " << variance << endl;
				if (variance > 3000.0)
				{
					GeoWindow.push_back(GeoPoints);
					GeoWindow_variance.push_back(variance);
				}
			}
			
		}
	}

	//��������
	FILE* fp = fopen("E:\\IMG_i864_All\\dem_grid\\slope\\slope3000.txt", "w");
	fprintf(fp, "sample line lon lat variance > 3000.0\n");
	for (int t = 0; t < GeoWindow.size(); t++)
	{
		fprintf(fp, "%lf\t%lf\t%lf\t%lf\t%lf\n", GeoWindow[t].x, GeoWindow[t].y, GeoWindow[t].lon, GeoWindow[t].lat,GeoWindow_variance[t]);
	}
	fclose(fp);
}

// ����ʮ ��ȡtab��ʽ��MOLA����
void read_Tab_Mola() 
{
	// �޶�һ������ γ�� ��Χ
	// �����ڸ÷�Χ�����м����ľ��� γ�� �߳�
	
	//MC-17��Χ
	/*string filePoints_out = "E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\PEDR\\DATA3\\MC-17\\MOLA17.txt";
	double leftlon = -135.0;
	double rightlon = -90.0;
	double uplat = 0.0;
	double downlat = -30.0;*/

	//MC-18��Χ
	//�㼯����·��
	/*string filePoints_out = "E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\PEDR\\DATA3\\MC-18\\MOLA18.txt";
	double leftlon = -90.0;
	double rightlon = -45.0;
	double uplat = 0.0;
	double downlat = -30.0;*/

	//MC-19��Χ
	/*string filePoints_out = "E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\PEDR\\DATA3\\MC-19\\MOLA19.txt";
	double leftlon = -45.0;
	double rightlon = 0.0;
	double uplat = 0.0;
	double downlat = -30.0;*/

	// MC-20��Χ
	/*string filePoints_out = "E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\PEDR\\DATA3\\MC-20\\MOLA20.txt";
	double leftlon = 0.0;
	double rightlon = 45.0;
	double uplat = 0.0;
	double downlat = -30.0;*/

	// MC-23��Χ
	string filePoints_out = "E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\PEDR\\DATA3\\MC-23\\MOLA23.txt";
	double leftlon = 135.0;
	double rightlon = 180.0;
	double uplat = 0.0;
	double downlat = -30.0;

	//MC-24��Χ
	/*string filePoints_out = "E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\PEDR\\DATA3\\MC-24\\MOLA24.txt";
	double leftlon = -180.0;
	double rightlon = -120.0;
	double uplat = -30.0;
	double downlat = -65.0;*/

	//MC-25��Χ
	/*string filePoints_out = "E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\PEDR\\DATA3\\MC-25\\MOLA25.txt";
	double leftlon = -120.0;
	double rightlon = -60.0;
	double uplat = -30.0;
	double downlat = -65.0;*/

	//MC-26��Χ
	/*string filePoints_out = "E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\PEDR\\DATA3\\MC-26\\MOLA26.txt";
	double leftlon = -60.0;
	double rightlon = 0.0;
	double uplat = -30.0;
	double downlat = -65.0;*/

	//MC-27��Χ
	/*string filePoints_out = "E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\PEDR\\DATA3\\MC-27\\MOLA27.txt";
	double leftlon = 0.0;
	double rightlon = 60.0;
	double uplat = -30.0;
	double downlat = -65.0;*/

	//MC-28��Χ
	/*string filePoints_out = "E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\PEDR\\DATA3\\MC-28\\MOLA28.txt";
	double leftlon = 60.0;
	double rightlon = 120.0;
	double uplat = -30.0;
	double downlat = -65.0;*/

	//mc-29��Χ
	/*string filePoints_out = "E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\PEDR\\DATA3\\MC-29\\MOLA29.txt";
	double leftlon = 120.0;
	double rightlon = 180.0;
	double uplat = -30.0;
	double downlat = -65.0;*/

	string tabListFile = "E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\PEDR\\DATA\\list.txt";//771��MOLA�ļ�·������
	vector<string> tabLists;
	FILE* fp_tablist = fopen(tabListFile.c_str(), "r");
	string tab_file;
	char c_read[1024];
	while (!feof(fp_tablist)) {
		if (fgets(c_read, 1024, fp_tablist) == NULL)
			continue;
		tabLists.push_back(c_read);
	}
	fclose(fp_tablist);

	FILE* fpMolaPoints = fopen(filePoints_out.c_str(), "w");
	fprintf(fpMolaPoints, "# MOLA���� MOLAγ�� MOLA�߳� \n");

	for (int tab_index = 0; tab_index < tabLists.size(); tab_index++)
	{
		// ��ȡ�ļ��ľ��� γ�� �߳�
		string tabfile = tabLists[tab_index];
		cout << tabfile << endl;
		int n = tabfile.find_last_not_of("\n");
		tabfile.erase(n + 1, tabfile.size() - n);
		FILE* fp = fopen(tabfile.c_str(), "r");
		/*vector<pair<double, double>> orbit;
		vector<double> h;*/
		char c_read2[1024];
		while (!feof(fp)) {
			if (fgets(c_read2, 1024, fp) == NULL)
				continue;
			double a1, a2, a3, a4, a5, a6, a7;
			if (fscanf(fp, "%lf%lf%lf%lf%lf%lf%lf", &a1, &a2, &a3, &a4, &a5, &a6, &a7) != 7)	//1���� 2γ�� 3topgra 4��� 5dem+����뾶 
				continue;
			//ת�����ȵı�ʾ����
			if (a1 >= 180.0) {
				a1 -= 360.0;
			}
			if (a1>leftlon&&a1<rightlon&&a2>downlat&&a2<uplat)
			{
				//orbit.push_back(make_pair(a1, a2));
				double h_a5 = a5 - 3396000.0;
				//h.push_back(h_a5);*/
				fprintf(fpMolaPoints, "%lf\t%lf\t%lf\n", a1, a2, h_a5);
			}
			
		}
		fclose(fp);
	}

	fclose(fpMolaPoints);
}

// ����ʮ�� ��ȡMOLA���ݣ����չ����,һ���洢
vector<vector<GeoPoint>> read_orbit_Mola(string fileIMG_MOLA) {

	vector<vector<GeoPoint>> points_IMG_MOLA;//ȫ����MOLA�㼯
	vector<GeoPoint> points_orbit_MOLA;//�����MOLA�㼯
	GeoPoint point_single;//MOLA��
	double lat_old = 0.;//��һʱ�̵�MOLA���γ��ֵ

	FILE* fp = fopen(fileIMG_MOLA.c_str(), "r");
	char c_read[1024];
	while (!feof(fp))
	{
		if (fgets(c_read, 1024, fp) == NULL)
			continue;
		if (fscanf(fp, "%lf %lf %lf", &point_single.lon, &point_single.lat, &point_single.h) != 3)
			continue;
		//��������
		else
		{
			// �����������֮��γ�Ȳ�ֵ����0.1�ȣ��洢һ���µ�����
			if (abs(point_single.lat - lat_old) >= 0.1 && points_orbit_MOLA.size() != 0)
			{
				points_IMG_MOLA.push_back(points_orbit_MOLA);//�洢һ��
				vector<GeoPoint>().swap(points_orbit_MOLA);//���
			}

			points_orbit_MOLA.push_back(point_single);
			
			lat_old = point_single.lat;
		}
	}
	fclose(fp);

	cout << "������" << points_IMG_MOLA.size() << "������" << endl;
	return points_IMG_MOLA;
}

// ����ʮ�� �ֶ�ƥ�䣬����CCC�����治һ�µĵ㼯
void test_ccc_points() {
	string fileIMG_MOLA = "E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\PEDR\\i977_mola\\output\\MolaPoints.txt";
	string demFile = "E:\\IMG_I977\\test2\\sensorCorrect_result\\run-DEM.tif";
	string fileMOLA_filter = "E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\PEDR\\i977_mola\\output\\molapointsfilter.txt";
	saveNone_CCC_points(read_orbit_Mola(fileIMG_MOLA), demFile, fileMOLA_filter);
}

// ����ʮ�� �������ƽ����
void test_match_Mola2Dem() {

	// offset save picture
	string script = "C:\\Users\\leeda\\Documents\\offsetsavepicture.plt";
	FILE* fp_script = fopen(script.c_str(), "w");
	fprintf(fp_script, "set term pdfcairo\ncd 'H:\leeda\date_keyanrenwu\Mars\MGS\MOLA\PEDR\i864_mola\output\offset'\nset xlabel 'lat (pixel)'\nset ylabel 'lon (pixel)'\nset zlabel 'offset(m)'\nset size ratio 4\nunset key\nset pm3d map\n");

	string molafile = "H:\\leeda\\date_keyanrenwu\\Mars\\MGS\\MOLA\\PEDR\\i864_mola\\output\\mola_kmeans\\";
	string DEMfile = "D:\\TEST\\IMG_i864_All\\TiePoints_calibaration\\calibration\\sensorCorrect_result\\run-DEM.tif";  
	string offsetfile = "H:\\leeda\\date_keyanrenwu\\Mars\\MGS\\MOLA\\PEDR\\i864_mola\\output\\offset\\";
	string DEMGeoCoordinatefile = "H:\\leeda\\date_keyanrenwu\\Mars\\MGS\\MOLA\\PEDR\\i864_mola\\output\\ƫ�ƺ��DEM��������\\";
	string picturefile = "../../../TEST2/IMG_I864/picture/";

	string DEMGeoCoordinatefile_all = "H:\\leeda\\date_keyanrenwu\\Mars\\MGS\\MOLA\\PEDR\\i864_mola\\output\\ƫ�ƺ��DEM��������\\all.txt";
	FILE* fpDEMGeoCoordinatefile_all = fopen(DEMGeoCoordinatefile_all.c_str(), "w");
	fprintf(fpDEMGeoCoordinatefile_all, "#molalon molalat molah demlon demlat demh no no no\n");
	string ccc = "H:\\leeda\\date_keyanrenwu\\Mars\\MGS\\MOLA\\PEDR\\i864_mola\\output\\offset\\ccc.txt"; // i864

	/*string molafile = "E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\PEDR\\i863_mola\\output\\mola_kmeans\\";
	string DEMfile = "E:\\IMG_i863\\test_tiecalibration\\sensorCorrect_result\\run-DEM.tif";
	string offsetfile = "E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\PEDR\\i863_mola\\output\\offset\\";
	string DEMGeoCoordinatefile = "E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\PEDR\\i863_mola\\output\\ƫ�ƺ��DEM��������\\";
	string picturefile = "../../../date_keyanrenwu/Mars/MGS/MOLA/PEDR/i863_mola/output/picture/";

	string DEMGeoCoordinatefile_all = "E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\PEDR\\i863_mola\\output\\ƫ�ƺ��DEM��������\\all.txt";
	FILE* fpDEMGeoCoordinatefile_all = fopen(DEMGeoCoordinatefile_all.c_str(), "w");
	fprintf(fpDEMGeoCoordinatefile_all, "#molalon molalat molah demlon demlat demh no no no\n");
	string ccc = "E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\PEDR\\i863_mola\\output\\offset\\ccc.txt";*/	//i863

	/*string molafile = "E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\PEDR\\i818_mola\\output\\mola_kmeans\\";
	string DEMfile = "E:\\IMG_I818\\test_2\\sensorCorrect_result\\run-DEM.tif";
	string offsetfile = "E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\PEDR\\i818_mola\\output\\offset\\";
	string DEMGeoCoordinatefile = "E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\PEDR\\i818_mola\\output\\ƫ�ƺ��DEM��������\\";
	string picturefile = "../../../date_keyanrenwu/Mars/MGS/MOLA/PEDR/i818_mola/output/picture/";

	string DEMGeoCoordinatefile_all = "E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\PEDR\\i818_mola\\output\\ƫ�ƺ��DEM��������\\all.txt";
	FILE* fpDEMGeoCoordinatefile_all = fopen(DEMGeoCoordinatefile_all.c_str(), "w");
	fprintf(fpDEMGeoCoordinatefile_all, "#molalon molalat molah demlon demlat demh no no no\n");
	string ccc = "E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\PEDR\\i818_mola\\output\\offset\\ccc.txt";*/

	/*string molafile = "E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\PEDR\\i977_mola\\output\\mola_kmeans\\";
	string DEMfile = "E:\\IMG_I977\\test2\\sensorCorrect_result\\run-DEM.tif";
	string offsetfile = "E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\PEDR\\i977_mola\\output\\offset\\";
	string DEMGeoCoordinatefile = "E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\PEDR\\i977_mola\\output\\ƫ�ƺ��DEM��������\\";
	string picturefile = "../../../date_keyanrenwu/Mars/MGS/MOLA/PEDR/i977_mola/output/picture/";

	string DEMGeoCoordinatefile_all = "E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\PEDR\\i977_mola\\output\\ƫ�ƺ��DEM��������\\all.txt";
	FILE* fpDEMGeoCoordinatefile_all = fopen(DEMGeoCoordinatefile_all.c_str(), "w");
	fprintf(fpDEMGeoCoordinatefile_all, "#molalon molalat molah demlon demlat demh no no no\n");
	string ccc = "E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\PEDR\\i977_mola\\output\\offset\\ccc.txt";*/ //i977

	FILE* fpccc = fopen(ccc.c_str(), "w");

	for (int i = 139; i < 140; i++)
	{
		cout << "ƥ���" << i << "��" << endl;
		string indexMolafile = molafile + to_string(i) + ".txt";
		string indexOffsetfile = offsetfile + to_string(i) + ".txt"; 
		string indexDEMGeoCoordinatefile = DEMGeoCoordinatefile + to_string(i) + ".txt"; 
		string indexPicturefile = picturefile + to_string(i) + ".pdf";
		double initCoefficientNcc = 0.0;
		double bestCoefficientNcc = 0.0;
		match_Mola2Dem(indexMolafile, DEMfile, indexOffsetfile, indexDEMGeoCoordinatefile, indexPicturefile, initCoefficientNcc, bestCoefficientNcc);
		fprintf(fpccc, "%d\t%lf\t%lf\t%lf\n", i, initCoefficientNcc, bestCoefficientNcc, bestCoefficientNcc / initCoefficientNcc);
		if (bestCoefficientNcc>0.9) // ��ȡƫ�ƺ��DEM��������
		{
			FILE* fp_mola = fopen(indexDEMGeoCoordinatefile.c_str(), "r");
			char c_read1[1024];
			while (!feof(fp_mola))
			{
				if (fgets(c_read1, 1024, fp_mola) == NULL)
					continue;
				double a1, a2, a3, a4, a5, a6, a7, a8, a9;
				if (fscanf(fp_mola, "%lf%lf%lf%lf%lf%lf%lf%lf%lf", &a1, &a2, &a3, &a4, &a5, &a6, &a7, &a8, &a9) != 9)
					continue;
				fprintf(fpDEMGeoCoordinatefile_all, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", a1, a2, a3, a4, a5, a6, a7, a8, a9);
			}
			fclose(fp_mola);

			fprintf(fp_script, "set output '%d.pdf'\n", i);
			fprintf(fp_script, "splot '%d.txt'\n", i);
			fprintf(fp_script, "set output\n");
		}
	}
	fclose(fpDEMGeoCoordinatefile_all);
	fclose(fpccc);
	fclose(fp_script);
}

// ����ʮ�� 
void make_gnuplot_script()
{
	// 1 ���Ƽ�ͷ
	/*//��ȡ�㼯
	string offset = "E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\PEDR\\i864_mola\\output\\offset.txt";
	vector<rpcGCP> Points;
	FILE* fp = fopen(offset.c_str(), "r");
	char c_read[1024];
	while (!feof(fp)) {
		if (fgets(c_read, 1024, fp) == NULL)
			continue;
		rpcGCP Point;
		if (fscanf(fp, "%lf%lf%lf%lf%lf", &Point.lon, &Point.lat, &Point.x, &Point.y, &Point.h) != 5)	 
			continue;
		Points.push_back(Point);
	}
	fclose(fp);

	//�ű�·��
	string script = "C:\\Users\\lee\\Documents\\drawarrow.plt";
	FILE* fp_script = fopen(script.c_str(), "w");
	fprintf(fp_script, "cd 'E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\PEDR\\i864_mola\\output'\n");
	fprintf(fp_script, "set xlabel 'X-LON (PIXEL)'\n");
	fprintf(fp_script, "set ylabel 'Y-LAT (PIXEL)'\n");

	double ori_lat = -22.999944910792745;
	double ori_lon = -63.687646394017960;
	double resolution_lon = 0.00057747454488099592;
	double resolution_lat = -0.00057747454488099592;

	for (int i = 0; i < Points.size(); i++)
	{
		fprintf(fp_script, "set arrow %d from %lf,%lf rto %lf,%lf nofilled\n", i + 1, (Points[i].lon - ori_lon) / resolution_lon, (Points[i].lat - ori_lat) / resolution_lat, Points[i].x, Points[i].y);
	}
	fprintf(fp_script, "plot 'offset.txt' u (($1+63.687646394017960)/0.00057747454488099592):(($2+22.999944910792745)/-0.00057747454488099592) pt 7 ps 1title 'MOLA Points'\n");
	fclose(fp_script);
	*/
	// 2 ��������offsetͼ
	string script = "C:\\Users\\lee\\Documents\\offseti863.plt";
	FILE* fp_script = fopen(script.c_str(), "w");
	fprintf(fp_script, "set term pdfcairo\ncd 'E:\date_keyanrenwu\Mars\MGS\MOLA\PEDR\i863_mola\output\offset'\nset xlabel 'lat (pixel)'\nset ylabel 'lon (pixel)'\nset zlabel 'offset(m)'\nset size ratio 4\nunset key\nset pm3d map\n");
	for (int i = 0; i < 280; i++)
	{
		fprintf(fp_script, "set output '%d.pdf'\n", i);
		fprintf(fp_script, "splot '%d.txt'\n", i);
		fprintf(fp_script, "set output\n");
	}
	fclose(fp_script);
}

// ����ʮ�� ��ȡSTEREO DEM����,תΪXYZ����ʽtxt
void test_stereo_dem2xyz() {
	//DEMӰ��
	string demFile = "E:\\IMG_i864_All\\TiePoints_calibaration\\calibration\\sensorCorrect_result\\run-DEM.tif";
	GeoImageIO dem;
	dem.open(demFile);
	//��ȡӰ�񸲸Ƿ�Χ
	double leftlon, rightlon, uplat, downlat;
	double trans[4];
	dem.getGeoTransform(trans);
	int demWid = dem.getImageWid();
	int demHei = dem.getImageHei();
	double ori_Dem_lat = trans[2], ori_Dem_lon = trans[0];
	double end_Dem_lat = trans[2] + trans[3] * demHei, end_Dem_lon = trans[0] + trans[1] * demWid;
	dem.readToBuffer_proj(0, ori_Dem_lon, end_Dem_lon, ori_Dem_lat, end_Dem_lat);

	//�ֿ������ά��
	//�ļ���·��
	string filePoints_out = "E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\PEDR\\i864_mola\\output\\stereo_dem\\";

	int index_wid_nums = 2;
	int index_hei_nums = ceil(dem.getImageHei() / dem.getImageWid() * index_wid_nums);
	int index_nums = index_wid_nums * index_hei_nums;

	for (int index = 0; index < index_nums; index++)
	{
		string filePath = filePoints_out + to_string(index) + "I864_stereo_dem_PointsXYZ.txt";
		FILE* fpMolaPoints = fopen(filePath.c_str(), "w");

		//��ʼ�к� �к�
		int ori_line = (index / index_wid_nums) * (dem.getImageWid() / index_wid_nums);
		int ori_sample = (index % index_wid_nums) * (dem.getImageWid() / index_wid_nums);

		for (int i=ori_line; i < ori_line + dem.getImageWid() / index_wid_nums ||i > dem.getImageHei(); i++)
		{
			for (int j = ori_sample; j < ori_sample + dem.getImageWid() / index_wid_nums || j > dem.getImageWid(); j++)
			{
				double h_value = 0.;
				dem.getBufferValue(ori_Dem_lon + j * trans[1], ori_Dem_lat + i * trans[3], 0, &h_value);
				if (h_value > -20000.0)
				{
					double X = 0.0;
					double Y = 0.0;
					double Z = 0.0;
					fromlatlonalt2XYZ(ori_Dem_lat + i * trans[3], ori_Dem_lon + j * trans[1], h_value, X, Y, Z);
					fprintf(fpMolaPoints, "%lf\t%lf\t%lf\n", X, Y, Z);
				}

			}
		}
		fclose(fpMolaPoints);
	}
	
}

void test_MOLA_dem2xyz() {
	//MOLA DEMӰ��
	string demFile = "E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\MEG128\\radius\\resample\\overrange_i864_136.tif";
	GeoImageIO dem;
	dem.open(demFile);
	//��ȡӰ�񸲸Ƿ�Χ
	double leftlon, rightlon, uplat, downlat;
	double trans[4];
	dem.getGeoTransform(trans);
	int demWid = dem.getImageWid();
	int demHei = dem.getImageHei();
	double ori_Dem_lat = trans[2], ori_Dem_lon = trans[0];
	double end_Dem_lat = trans[2] + trans[3] * demHei, end_Dem_lon = trans[0] + trans[1] * demWid;
	dem.readToBuffer_proj(0, ori_Dem_lon, end_Dem_lon, ori_Dem_lat, end_Dem_lat);

	//��ҪתΪ��������ϵ
	TransCoordAPI coord;
	coord.init(dem.getProjectionRef().c_str(), GcsMarsDTMWTK.c_str());
	//�������Χ
	double ori_Dem_lat_Geo, ori_Dem_lon_Geo;
	coord.transCoord(ori_Dem_lat, ori_Dem_lon, ori_Dem_lat_Geo, ori_Dem_lon_Geo);
	double end_Dem_lat_Geo, end_Dem_lon_Geo;
	coord.transCoord(end_Dem_lat, end_Dem_lon, end_Dem_lat_Geo, end_Dem_lon_Geo);
	double resolution_lat = (end_Dem_lat_Geo - ori_Dem_lat_Geo) / demHei;
	double resolution_lon = (end_Dem_lon_Geo - ori_Dem_lon_Geo) / demWid;

	//�ļ���·��
	string filePoints_out = "E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\PEDR\\i864_mola\\output\\mola_dem\\";

	int index_wid_nums = 2;
	int index_hei_nums = ceil(dem.getImageHei() / dem.getImageWid() * index_wid_nums);
	int index_nums = index_wid_nums * index_hei_nums;

	for (int index = 0; index < index_nums; index++)
	{
		string filePath = filePoints_out + to_string(index) + "I864_mola_dem_PointsXYZ.txt";
		FILE* fpMolaPoints = fopen(filePath.c_str(), "w");

		//��ʼ�к� �к�
		int ori_line = (index / index_wid_nums) * (dem.getImageWid() / index_wid_nums);
		int ori_sample = (index % index_wid_nums) * (dem.getImageWid() / index_wid_nums);
		
		for (int i = ori_line; i < ori_line + dem.getImageWid() / index_wid_nums || i > dem.getImageHei(); i++)
		{
			for (int j = ori_sample; j < ori_sample + dem.getImageWid() / index_wid_nums || j > dem.getImageWid(); j++)
			{
				double h_value = 0.;
				dem.getBufferValue(ori_Dem_lon + j * trans[1], ori_Dem_lat + i * trans[3], 0, &h_value);
				if (h_value > -20000.0)
				{
					double X = 0.0, Y = 0.0, Z = 0.0;
					double lat = ori_Dem_lat_Geo + i * resolution_lat; 
					double lon = ori_Dem_lon_Geo + j * resolution_lon;
					fromlatlonalt2XYZ(lat, lon, h_value, X, Y, Z);
					fprintf(fpMolaPoints, "%lf\t%lf\t%lf\n", X, Y, Z);
				}

			}
		}
		fclose(fpMolaPoints);
	}
}

// ����ʮ�� ����STEREO DEM��MOLA DEM���ص��������ڲü�
void test_overrange() {
	//STEREO DEMӰ��
	string demFile = "E:\\IMG_i864_All\\TiePoints_calibaration\\calibration\\sensorCorrect_result\\run-DEM.tif";
	GeoImageIO dem;
	dem.open(demFile);
	//��ȡӰ�񸲸Ƿ�Χ
	double leftlon, rightlon, uplat, downlat;
	double trans[4];
	dem.getGeoTransform(trans);
	int demWid = dem.getImageWid();
	int demHei = dem.getImageHei();
	double ori_Dem_lat = trans[2], ori_Dem_lon = trans[0];
	double end_Dem_lat = trans[2] + trans[3] * demHei, end_Dem_lon = trans[0] + trans[1] * demWid;

	//MOLA DEMӰ��
	string MOLAdemFile = "E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\MEG128\\radius\\mars_dem_global.tif";
	GeoImageIO MOLAdem;
	MOLAdem.open(MOLAdemFile);
	double MOLAtrans[4];
	MOLAdem.getGeoTransform(MOLAtrans);
	
	TransCoordAPI coord;
	// �����һ��ʹ�ý��ȴ���������
	//coord.init(dem.getProjectionRef().c_str(), MOLAdem.getProjectionRef().c_str());
	// ���ֻ������һ��
	coord.init(GcsMarsDTMWTK.c_str(), MOLAdem.getProjectionRef().c_str());
	double ori_MOLA_lat, ori_MOLA_lon, end_MOLA_lat, end_MOLA_lon;
	coord.transCoord(ori_Dem_lat, ori_Dem_lon, ori_MOLA_lat, ori_MOLA_lon);
	coord.transCoord(end_Dem_lat, end_Dem_lon, end_MOLA_lat, end_MOLA_lon);

	int startx = (ori_MOLA_lon - MOLAtrans[0]) / MOLAtrans[1];
	int starty = (ori_MOLA_lat - MOLAtrans[2]) / MOLAtrans[3];
	int sizex = (end_MOLA_lon - ori_MOLA_lon) / MOLAtrans[1];
	int sizey = (end_MOLA_lat - ori_MOLA_lat) / MOLAtrans[3];

	string cutfile = "E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\MEG128\\radius\\cutrange\\overrange_i864.tif";
	
	segment_image(MOLAdemFile, cutfile, sizex, sizey, startx, starty);
	
}

// ����ʮ�� ��ȡVICAR����
void test15() {
	string fileimg = "../data/i864/img/hi864_0000_nd2.img";
	read_VICAR_img(fileimg);
}


// ����ʮ�� ��ȡDTM��Ӧλ�ã�һ�켤��㣩�߳�ֵ
void test_readDTM_singleMola() {
	string file_mola = "H:\\leeda\\date_keyanrenwu\\Mars\\MGS\\MOLA\\PEDR\\DATA2\\AP10228L.TAB";
	string file_DTM_beforeCalibration = "D:\\TEST\\IMG_i864_All\\TiePoints_calibaration\\calibration\\sensorCorrect_result\\run-DEM.tif";
	string file_altitude_mola_dtm_beforeCalibration = "D:\\TEST2\\mola\\i864\\10228_before.txt";
	readMolaGpoints(file_mola, file_altitude_mola_dtm_beforeCalibration, file_DTM_beforeCalibration);

	string file_DTM_afterCalibration = "H:\\leeda\\shiyan\\sensorCorrect_result_I864_success\\run-DEM.tif";
	string file_altitude_mola_dtm_afterCalibration = "D:\\TEST2\\mola\\i864\\10228_after.txt";
	readMolaGpoints(file_mola, file_altitude_mola_dtm_afterCalibration, file_DTM_afterCalibration);
}

//804�������װ+�ڷ�λ��������
void makeInnerPara() {
	FILE* fp = fopen("C:\\Users\\leeda\\Desktop\\inner.txt", "w");
	double ccdlamda = 14e-6;
	double f = 1.6;

	int NumofCCDs = 33;
	fprintf(fp, "NumofCCDs = %d\n", NumofCCDs);
	int pixelsPerCCD = 1024;
	fprintf(fp, "pixelsPerCCD = %d\n", pixelsPerCCD);
	double ori = 219.9528e-3;

	int CCDID = 0;
	fprintf(fp, "CCDID = =%d\n", CCDID);
	double ori_X = -219.9528e-3;
	double end_X = -205.6168e-3;
	double Y = -207.5220e-3;
	vector<cv::Point2f> points_tanx;
	vector<cv::Point2f> points_tany;
	int n_x = 3;	//����ʽϵ������
	int n_y = 3;	//����ʽϵ������
	for (int i = 0; i < pixelsPerCCD; i++)
	{
		double tanx = -Y / f;
		double tany = -(ori_X - ori + i * ccdlamda) / f;
		points_tanx.push_back(cv::Point2f(i, tanx));
		points_tany.push_back(cv::Point2f(i, tany));
	}
	cv::Mat xMat = polyfit(points_tanx, n_x);
	cv::Mat yMat = polyfit(points_tany, n_y);
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", xMat.at<double>(0, 0), xMat.at<double>(1, 0), xMat.at<double>(2, 0), xMat.at<double>(3, 0));
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", yMat.at<double>(0, 0), yMat.at<double>(1, 0), yMat.at<double>(2, 0), yMat.at<double>(3, 0));
	points_tanx.clear();
	points_tany.clear();

	CCDID = 1;
	fprintf(fp, "CCDID = =%d\n", CCDID);
	ori_X = -208.0366e-3;
	end_X = -193.7006e-3;
	Y = -193.7754e-3;
	for (int i = 0; i < pixelsPerCCD; i++)
	{
		double tanx = -Y / f;// X���췽��
		double tany = -(ori_X - ori + i * ccdlamda) / f;
		points_tanx.push_back(cv::Point2f(i, tanx));
		points_tany.push_back(cv::Point2f(i, tany));
	}
	xMat = polyfit(points_tanx, n_x);
	yMat = polyfit(points_tany, n_y);
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", xMat.at<double>(0, 0), xMat.at<double>(1, 0), xMat.at<double>(2, 0), xMat.at<double>(3, 0));
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", yMat.at<double>(0, 0), yMat.at<double>(1, 0), yMat.at<double>(2, 0), yMat.at<double>(3, 0));
	points_tanx.clear();
	points_tany.clear();

	CCDID = 2;
	fprintf(fp, "CCDID = =%d\n", CCDID);
	ori_X = -195.5563e-3;
	end_X = -181.2203e-3;
	Y = -207.7522e-3;
	for (int i = 0; i < pixelsPerCCD; i++)
	{
		double tanx = -Y / f;// X���췽��
		double tany = -(ori_X - ori + i * ccdlamda) / f;
		points_tanx.push_back(cv::Point2f(i, tanx));
		points_tany.push_back(cv::Point2f(i, tany));
	}
	xMat = polyfit(points_tanx, n_x);
	yMat = polyfit(points_tany, n_y);
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", xMat.at<double>(0, 0), xMat.at<double>(1, 0), xMat.at<double>(2, 0), xMat.at<double>(3, 0));
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", yMat.at<double>(0, 0), yMat.at<double>(1, 0), yMat.at<double>(2, 0), yMat.at<double>(3, 0));
	points_tanx.clear();
	points_tany.clear();

	CCDID = 3;
	fprintf(fp, "CCDID = =%d\n", CCDID);
	ori_X = -183.0670e-3;
	end_X = -168.7310e-3;
	Y = -193.9835e-3;
	for (int i = 0; i < pixelsPerCCD; i++)
	{
		double tanx = -Y / f;// X���췽��
		double tany = -(ori_X - ori + i * ccdlamda) / f;
		points_tanx.push_back(cv::Point2f(i, tanx));
		points_tany.push_back(cv::Point2f(i, tany));
	}
	xMat = polyfit(points_tanx, n_x);
	yMat = polyfit(points_tany, n_y);
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", xMat.at<double>(0, 0), xMat.at<double>(1, 0), xMat.at<double>(2, 0), xMat.at<double>(3, 0));
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", yMat.at<double>(0, 0), yMat.at<double>(1, 0), yMat.at<double>(2, 0), yMat.at<double>(3, 0));
	points_tanx.clear();
	points_tany.clear();

	CCDID = 4;
	fprintf(fp, "CCDID = =%d\n", CCDID);
	ori_X = -170.1485e-3;
	end_X = -155.8125e-3;
	Y = -207.9633e-3;
	for (int i = 0; i < pixelsPerCCD; i++)
	{
		double tanx = -Y / f;// X���췽��
		double tany = -(ori_X - ori + i * ccdlamda) / f;
		points_tanx.push_back(cv::Point2f(i, tanx));
		points_tany.push_back(cv::Point2f(i, tany));
	}
	xMat = polyfit(points_tanx, n_x);
	yMat = polyfit(points_tany, n_y);
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", xMat.at<double>(0, 0), xMat.at<double>(1, 0), xMat.at<double>(2, 0), xMat.at<double>(3, 0));
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", yMat.at<double>(0, 0), yMat.at<double>(1, 0), yMat.at<double>(2, 0), yMat.at<double>(3, 0));
	points_tanx.clear();
	points_tany.clear();

	CCDID = 5;
	fprintf(fp, "CCDID = =%d\n", CCDID);
	ori_X = -157.2622e-3;
	end_X = -142.9262e-3;
	Y = -194.1729e-3;
	for (int i = 0; i < pixelsPerCCD; i++)
	{
		double tanx = -Y / f;// X���췽��
		double tany = -(ori_X - ori + i * ccdlamda) / f;
		points_tanx.push_back(cv::Point2f(i, tanx));
		points_tany.push_back(cv::Point2f(i, tany));
	}
	xMat = polyfit(points_tanx, n_x);
	yMat = polyfit(points_tany, n_y);
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", xMat.at<double>(0, 0), xMat.at<double>(1, 0), xMat.at<double>(2, 0), xMat.at<double>(3, 0));
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", yMat.at<double>(0, 0), yMat.at<double>(1, 0), yMat.at<double>(2, 0), yMat.at<double>(3, 0));
	points_tanx.clear();
	points_tany.clear();

	CCDID = 6;
	fprintf(fp, "CCDID = =%d\n", CCDID);
	ori_X = -144.2951e-3;
	end_X = -129.9591e-3;
	Y = -208.1500e-3;
	for (int i = 0; i < pixelsPerCCD; i++)
	{
		double tanx = -Y / f;// X���췽��
		double tany = -(ori_X - ori + i * ccdlamda) / f;
		points_tanx.push_back(cv::Point2f(i, tanx));
		points_tany.push_back(cv::Point2f(i, tany));
	}
	xMat = polyfit(points_tanx, n_x);
	yMat = polyfit(points_tany, n_y);
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", xMat.at<double>(0, 0), xMat.at<double>(1, 0), xMat.at<double>(2, 0), xMat.at<double>(3, 0));
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", yMat.at<double>(0, 0), yMat.at<double>(1, 0), yMat.at<double>(2, 0), yMat.at<double>(3, 0));
	points_tanx.clear();
	points_tany.clear();

	CCDID = 7;
	fprintf(fp, "CCDID = =%d\n", CCDID);
	ori_X = -131.0922e-3;
	end_X = -116.7562e-3;
	Y = -194.3370e-3;
	for (int i = 0; i < pixelsPerCCD; i++)
	{
		double tanx = -Y / f;// X���췽��
		double tany = -(ori_X - ori + i * ccdlamda) / f;
		points_tanx.push_back(cv::Point2f(i, tanx));
		points_tany.push_back(cv::Point2f(i, tany));
	}
	xMat = polyfit(points_tanx, n_x);
	yMat = polyfit(points_tany, n_y);
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", xMat.at<double>(0, 0), xMat.at<double>(1, 0), xMat.at<double>(2, 0), xMat.at<double>(3, 0));
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", yMat.at<double>(0, 0), yMat.at<double>(1, 0), yMat.at<double>(2, 0), yMat.at<double>(3, 0));
	points_tanx.clear();
	points_tany.clear();

	CCDID = 8;
	fprintf(fp, "CCDID = =%d\n", CCDID);
	ori_X = -117.6013e-3;
	end_X = -103.2653e-3;
	Y = -208.3114e-3;
	for (int i = 0; i < pixelsPerCCD; i++)
	{
		double tanx = -Y / f;// X���췽��
		double tany = -(ori_X - ori + i * ccdlamda) / f;
		points_tanx.push_back(cv::Point2f(i, tanx));
		points_tany.push_back(cv::Point2f(i, tany));
	}
	xMat = polyfit(points_tanx, n_x);
	yMat = polyfit(points_tany, n_y);
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", xMat.at<double>(0, 0), xMat.at<double>(1, 0), xMat.at<double>(2, 0), xMat.at<double>(3, 0));
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", yMat.at<double>(0, 0), yMat.at<double>(1, 0), yMat.at<double>(2, 0), yMat.at<double>(3, 0));
	points_tanx.clear();
	points_tany.clear();

	CCDID = 9;
	fprintf(fp, "CCDID = =%d\n", CCDID);
	ori_X = -104.1480e-3;
	end_X = -89.8120e-3;
	Y = -194.4750e-3;
	for (int i = 0; i < pixelsPerCCD; i++)
	{
		double tanx = -Y / f;// X���췽��
		double tany = -(ori_X - ori + i * ccdlamda) / f;
		points_tanx.push_back(cv::Point2f(i, tanx));
		points_tany.push_back(cv::Point2f(i, tany));
	}
	xMat = polyfit(points_tanx, n_x);
	yMat = polyfit(points_tany, n_y);
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", xMat.at<double>(0, 0), xMat.at<double>(1, 0), xMat.at<double>(2, 0), xMat.at<double>(3, 0));
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", yMat.at<double>(0, 0), yMat.at<double>(1, 0), yMat.at<double>(2, 0), yMat.at<double>(3, 0));
	points_tanx.clear();
	points_tany.clear();

	CCDID = 10;
	fprintf(fp, "CCDID = =%d\n", CCDID);
	ori_X = -90.4539e-3;
	end_X = -76.1179e-3;
	Y = -208.4415e-3;
	for (int i = 0; i < pixelsPerCCD; i++)
	{
		double tanx = -Y / f;// X���췽��
		double tany = -(ori_X - ori + i * ccdlamda) / f;
		points_tanx.push_back(cv::Point2f(i, tanx));
		points_tany.push_back(cv::Point2f(i, tany));
	}
	xMat = polyfit(points_tanx, n_x);
	yMat = polyfit(points_tany, n_y);
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", xMat.at<double>(0, 0), xMat.at<double>(1, 0), xMat.at<double>(2, 0), xMat.at<double>(3, 0));
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", yMat.at<double>(0, 0), yMat.at<double>(1, 0), yMat.at<double>(2, 0), yMat.at<double>(3, 0));
	points_tanx.clear();
	points_tany.clear();

	CCDID = 11;
	fprintf(fp, "CCDID = =%d\n", CCDID);
	ori_X = -76.7763e-3;
	end_X = -62.4403e-3;
	Y = -194.5820e-3;
	for (int i = 0; i < pixelsPerCCD; i++)
	{
		double tanx = -Y / f;// X���췽��
		double tany = -(ori_X - ori + i * ccdlamda) / f;
		points_tanx.push_back(cv::Point2f(i, tanx));
		points_tany.push_back(cv::Point2f(i, tany));
	}
	xMat = polyfit(points_tanx, n_x);
	yMat = polyfit(points_tany, n_y);
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", xMat.at<double>(0, 0), xMat.at<double>(1, 0), xMat.at<double>(2, 0), xMat.at<double>(3, 0));
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", yMat.at<double>(0, 0), yMat.at<double>(1, 0), yMat.at<double>(2, 0), yMat.at<double>(3, 0));
	points_tanx.clear();
	points_tany.clear();

	CCDID = 12;
	fprintf(fp, "CCDID = =%d\n", CCDID);
	ori_X = -62.9214e-3;
	end_X = -48.5854e-3;
	Y = -208.5371e-3;
	for (int i = 0; i < pixelsPerCCD; i++)
	{
		double tanx = -Y / f;// X���췽��
		double tany = -(ori_X - ori + i * ccdlamda) / f;
		points_tanx.push_back(cv::Point2f(i, tanx));
		points_tany.push_back(cv::Point2f(i, tany));
	}
	xMat = polyfit(points_tanx, n_x);
	yMat = polyfit(points_tany, n_y);
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", xMat.at<double>(0, 0), xMat.at<double>(1, 0), xMat.at<double>(2, 0), xMat.at<double>(3, 0));
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", yMat.at<double>(0, 0), yMat.at<double>(1, 0), yMat.at<double>(2, 0), yMat.at<double>(3, 0));
	points_tanx.clear();
	points_tany.clear();

	CCDID = 13;
	fprintf(fp, "CCDID = =%d\n", CCDID);
	ori_X = -49.3248e-3;
	end_X = -34.9888e-3;
	Y = -194.6547e-3;
	for (int i = 0; i < pixelsPerCCD; i++)
	{
		double tanx = -Y / f;// X���췽��
		double tany = -(ori_X - ori + i * ccdlamda) / f;
		points_tanx.push_back(cv::Point2f(i, tanx));
		points_tany.push_back(cv::Point2f(i, tany));
	}
	xMat = polyfit(points_tanx, n_x);
	yMat = polyfit(points_tany, n_y);
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", xMat.at<double>(0, 0), xMat.at<double>(1, 0), xMat.at<double>(2, 0), xMat.at<double>(3, 0));
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", yMat.at<double>(0, 0), yMat.at<double>(1, 0), yMat.at<double>(2, 0), yMat.at<double>(3, 0));
	points_tanx.clear();
	points_tany.clear();

	CCDID = 14;
	fprintf(fp, "CCDID = =%d\n", CCDID);
	ori_X = -35.3225e-3;
	end_X = -20.9865e-3;
	Y = -208.5957e-3;
	for (int i = 0; i < pixelsPerCCD; i++)
	{
		double tanx = -Y / f;// X���췽��
		double tany = -(ori_X - ori + i * ccdlamda) / f;
		points_tanx.push_back(cv::Point2f(i, tanx));
		points_tany.push_back(cv::Point2f(i, tany));
	}
	xMat = polyfit(points_tanx, n_x);
	yMat = polyfit(points_tany, n_y);
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", xMat.at<double>(0, 0), xMat.at<double>(1, 0), xMat.at<double>(2, 0), xMat.at<double>(3, 0));
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", yMat.at<double>(0, 0), yMat.at<double>(1, 0), yMat.at<double>(2, 0), yMat.at<double>(3, 0));
	points_tanx.clear();
	points_tany.clear();

	CCDID = 15;
	fprintf(fp, "CCDID = =%d\n", CCDID);
	ori_X = -21.2895e-3;
	end_X = -6.9535e-3;
	Y = -194.6926e-3;
	for (int i = 0; i < pixelsPerCCD; i++)
	{
		double tanx = -Y / f;// X���췽��
		double tany = -(ori_X - ori + i * ccdlamda) / f;
		points_tanx.push_back(cv::Point2f(i, tanx));
		points_tany.push_back(cv::Point2f(i, tany));
	}
	xMat = polyfit(points_tanx, n_x);
	yMat = polyfit(points_tany, n_y);
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", xMat.at<double>(0, 0), xMat.at<double>(1, 0), xMat.at<double>(2, 0), xMat.at<double>(3, 0));
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", yMat.at<double>(0, 0), yMat.at<double>(1, 0), yMat.at<double>(2, 0), yMat.at<double>(3, 0));
	points_tanx.clear();
	points_tany.clear();

	CCDID = 16;
	fprintf(fp, "CCDID = =%d\n", CCDID);
	ori_X = -7.1680e-3;
	end_X = 7.1680e-3;
	Y = -208.6169e-3;
	for (int i = 0; i < pixelsPerCCD; i++)
	{
		double tanx = -Y / f;// X���췽��
		double tany = -(ori_X - ori + i * ccdlamda) / f;
		points_tanx.push_back(cv::Point2f(i, tanx));
		points_tany.push_back(cv::Point2f(i, tany));
	}
	xMat = polyfit(points_tanx, n_x);
	yMat = polyfit(points_tany, n_y);
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", xMat.at<double>(0, 0), xMat.at<double>(1, 0), xMat.at<double>(2, 0), xMat.at<double>(3, 0));
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", yMat.at<double>(0, 0), yMat.at<double>(1, 0), yMat.at<double>(2, 0), yMat.at<double>(3, 0));
	points_tanx.clear();
	points_tany.clear();

	CCDID = 17;
	fprintf(fp, "CCDID = =%d\n", CCDID);
	ori_X = 6.9535e-3;
	end_X = 21.2895e-3;
	Y = -194.6926e-3;
	for (int i = 0; i < pixelsPerCCD; i++)
	{
		double tanx = -Y / f;// X���췽��
		double tany = -(ori_X - ori + i * ccdlamda) / f;
		points_tanx.push_back(cv::Point2f(i, tanx));
		points_tany.push_back(cv::Point2f(i, tany));
	}
	xMat = polyfit(points_tanx, n_x);
	yMat = polyfit(points_tany, n_y);
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", xMat.at<double>(0, 0), xMat.at<double>(1, 0), xMat.at<double>(2, 0), xMat.at<double>(3, 0));
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", yMat.at<double>(0, 0), yMat.at<double>(1, 0), yMat.at<double>(2, 0), yMat.at<double>(3, 0));
	points_tanx.clear();
	points_tany.clear();

	CCDID = 18;
	fprintf(fp, "CCDID = =%d\n", CCDID);
	ori_X = 20.9865e-3;
	end_X = 35.3225e-3;
	Y = -208.5957e-3;
	for (int i = 0; i < pixelsPerCCD; i++)
	{
		double tanx = -Y / f;// X���췽��
		double tany = -(ori_X - ori + i * ccdlamda) / f;
		points_tanx.push_back(cv::Point2f(i, tanx));
		points_tany.push_back(cv::Point2f(i, tany));
	}
	xMat = polyfit(points_tanx, n_x);
	yMat = polyfit(points_tany, n_y);
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", xMat.at<double>(0, 0), xMat.at<double>(1, 0), xMat.at<double>(2, 0), xMat.at<double>(3, 0));
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", yMat.at<double>(0, 0), yMat.at<double>(1, 0), yMat.at<double>(2, 0), yMat.at<double>(3, 0));
	points_tanx.clear();
	points_tany.clear();

	CCDID = 19;
	fprintf(fp, "CCDID = =%d\n", CCDID);
	ori_X = 34.9888e-3;
	end_X = 49.3248e-3;
	Y = -194.6547e-3;
	for (int i = 0; i < pixelsPerCCD; i++)
	{
		double tanx = -Y / f;// X���췽��
		double tany = -(ori_X - ori + i * ccdlamda) / f;
		points_tanx.push_back(cv::Point2f(i, tanx));
		points_tany.push_back(cv::Point2f(i, tany));
	}
	xMat = polyfit(points_tanx, n_x);
	yMat = polyfit(points_tany, n_y);
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", xMat.at<double>(0, 0), xMat.at<double>(1, 0), xMat.at<double>(2, 0), xMat.at<double>(3, 0));
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", yMat.at<double>(0, 0), yMat.at<double>(1, 0), yMat.at<double>(2, 0), yMat.at<double>(3, 0));
	points_tanx.clear();
	points_tany.clear();

	CCDID = 20;
	fprintf(fp, "CCDID = =%d\n", CCDID);
	ori_X = 48.5854e-3;
	end_X = 62.9214e-3;
	Y = -208.5371e-3;
	for (int i = 0; i < pixelsPerCCD; i++)
	{
		double tanx = -Y / f;// X���췽��
		double tany = -(ori_X - ori + i * ccdlamda) / f;
		points_tanx.push_back(cv::Point2f(i, tanx));
		points_tany.push_back(cv::Point2f(i, tany));
	}
	xMat = polyfit(points_tanx, n_x);
	yMat = polyfit(points_tany, n_y);
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", xMat.at<double>(0, 0), xMat.at<double>(1, 0), xMat.at<double>(2, 0), xMat.at<double>(3, 0));
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", yMat.at<double>(0, 0), yMat.at<double>(1, 0), yMat.at<double>(2, 0), yMat.at<double>(3, 0));
	points_tanx.clear();
	points_tany.clear();

	CCDID = 21;
	fprintf(fp, "CCDID = =%d\n", CCDID);
	ori_X = 62.4403e-3;
	end_X = 76.7763e-3;
	Y = -194.5820e-3;
	for (int i = 0; i < pixelsPerCCD; i++)
	{
		double tanx = -Y / f;// X���췽��
		double tany = -(ori_X - ori + i * ccdlamda) / f;
		points_tanx.push_back(cv::Point2f(i, tanx));
		points_tany.push_back(cv::Point2f(i, tany));
	}
	xMat = polyfit(points_tanx, n_x);
	yMat = polyfit(points_tany, n_y);
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", xMat.at<double>(0, 0), xMat.at<double>(1, 0), xMat.at<double>(2, 0), xMat.at<double>(3, 0));
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", yMat.at<double>(0, 0), yMat.at<double>(1, 0), yMat.at<double>(2, 0), yMat.at<double>(3, 0));
	points_tanx.clear();
	points_tany.clear();

	CCDID = 22;
	fprintf(fp, "CCDID = =%d\n", CCDID);
	ori_X = 76.1179e-3;
	end_X = 90.4539e-3;
	Y = -208.4415e-3;
	for (int i = 0; i < pixelsPerCCD; i++)
	{
		double tanx = -Y / f;// X���췽��
		double tany = -(ori_X - ori + i * ccdlamda) / f;
		points_tanx.push_back(cv::Point2f(i, tanx));
		points_tany.push_back(cv::Point2f(i, tany));
	}
	xMat = polyfit(points_tanx, n_x);
	yMat = polyfit(points_tany, n_y);
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", xMat.at<double>(0, 0), xMat.at<double>(1, 0), xMat.at<double>(2, 0), xMat.at<double>(3, 0));
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", yMat.at<double>(0, 0), yMat.at<double>(1, 0), yMat.at<double>(2, 0), yMat.at<double>(3, 0));
	points_tanx.clear();
	points_tany.clear();

	CCDID = 23;
	fprintf(fp, "CCDID = =%d\n", CCDID);
	ori_X = 89.8120e-3;
	end_X = 104.1480e-3;
	Y = -194.4750e-3;
	for (int i = 0; i < pixelsPerCCD; i++)
	{
		double tanx = -Y / f;// X���췽��
		double tany = -(ori_X - ori + i * ccdlamda) / f;
		points_tanx.push_back(cv::Point2f(i, tanx));
		points_tany.push_back(cv::Point2f(i, tany));
	}
	xMat = polyfit(points_tanx, n_x);
	yMat = polyfit(points_tany, n_y);
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", xMat.at<double>(0, 0), xMat.at<double>(1, 0), xMat.at<double>(2, 0), xMat.at<double>(3, 0));
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", yMat.at<double>(0, 0), yMat.at<double>(1, 0), yMat.at<double>(2, 0), yMat.at<double>(3, 0));
	points_tanx.clear();
	points_tany.clear();

	CCDID = 24;
	fprintf(fp, "CCDID = =%d\n", CCDID);
	ori_X = 103.2653e-3;
	end_X = 117.6013e-3;
	Y = -208.3114e-3;
	for (int i = 0; i < pixelsPerCCD; i++)
	{
		double tanx = -Y / f;// X���췽��
		double tany = -(ori_X - ori + i * ccdlamda) / f;
		points_tanx.push_back(cv::Point2f(i, tanx));
		points_tany.push_back(cv::Point2f(i, tany));
	}
	xMat = polyfit(points_tanx, n_x);
	yMat = polyfit(points_tany, n_y);
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", xMat.at<double>(0, 0), xMat.at<double>(1, 0), xMat.at<double>(2, 0), xMat.at<double>(3, 0));
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", yMat.at<double>(0, 0), yMat.at<double>(1, 0), yMat.at<double>(2, 0), yMat.at<double>(3, 0));
	points_tanx.clear();
	points_tany.clear();

	CCDID = 25;
	fprintf(fp, "CCDID = =%d\n", CCDID);
	ori_X = 116.7562e-3;
	end_X = 131.0922e-3;
	Y = -194.3370e-3;
	for (int i = 0; i < pixelsPerCCD; i++)
	{
		double tanx = -Y / f;// X���췽��
		double tany = -(ori_X - ori + i * ccdlamda) / f;
		points_tanx.push_back(cv::Point2f(i, tanx));
		points_tany.push_back(cv::Point2f(i, tany));
	}
	xMat = polyfit(points_tanx, n_x);
	yMat = polyfit(points_tany, n_y);
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", xMat.at<double>(0, 0), xMat.at<double>(1, 0), xMat.at<double>(2, 0), xMat.at<double>(3, 0));
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", yMat.at<double>(0, 0), yMat.at<double>(1, 0), yMat.at<double>(2, 0), yMat.at<double>(3, 0));
	points_tanx.clear();
	points_tany.clear();

	CCDID = 26;
	fprintf(fp, "CCDID = =%d\n", CCDID);
	ori_X = 129.9591e-3;
	end_X = 144.2951e-3;
	Y = -208.1500e-3;
	for (int i = 0; i < pixelsPerCCD; i++)
	{
		double tanx = -Y / f;// X���췽��
		double tany = -(ori_X - ori + i * ccdlamda) / f;
		points_tanx.push_back(cv::Point2f(i, tanx));
		points_tany.push_back(cv::Point2f(i, tany));
	}
	xMat = polyfit(points_tanx, n_x);
	yMat = polyfit(points_tany, n_y);
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", xMat.at<double>(0, 0), xMat.at<double>(1, 0), xMat.at<double>(2, 0), xMat.at<double>(3, 0));
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", yMat.at<double>(0, 0), yMat.at<double>(1, 0), yMat.at<double>(2, 0), yMat.at<double>(3, 0));
	points_tanx.clear();
	points_tany.clear();

	CCDID = 27;
	fprintf(fp, "CCDID = =%d\n", CCDID);
	ori_X = 142.9262e-3;
	end_X = 157.2622e-3;
	Y = -194.1729e-3;
	for (int i = 0; i < pixelsPerCCD; i++)
	{
		double tanx = -Y / f;// X���췽��
		double tany = -(ori_X - ori + i * ccdlamda) / f;
		points_tanx.push_back(cv::Point2f(i, tanx));
		points_tany.push_back(cv::Point2f(i, tany));
	}
	xMat = polyfit(points_tanx, n_x);
	yMat = polyfit(points_tany, n_y);
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", xMat.at<double>(0, 0), xMat.at<double>(1, 0), xMat.at<double>(2, 0), xMat.at<double>(3, 0));
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", yMat.at<double>(0, 0), yMat.at<double>(1, 0), yMat.at<double>(2, 0), yMat.at<double>(3, 0));
	points_tanx.clear();
	points_tany.clear();

	CCDID = 28;
	fprintf(fp, "CCDID = =%d\n", CCDID);
	ori_X = 155.8125e-3;
	end_X = 170.1485e-3;
	Y = -207.9633e-3;
	for (int i = 0; i < pixelsPerCCD; i++)
	{
		double tanx = -Y / f;// X���췽��
		double tany = -(ori_X - ori + i * ccdlamda) / f;
		points_tanx.push_back(cv::Point2f(i, tanx));
		points_tany.push_back(cv::Point2f(i, tany));
	}
	xMat = polyfit(points_tanx, n_x);
	yMat = polyfit(points_tany, n_y);
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", xMat.at<double>(0, 0), xMat.at<double>(1, 0), xMat.at<double>(2, 0), xMat.at<double>(3, 0));
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", yMat.at<double>(0, 0), yMat.at<double>(1, 0), yMat.at<double>(2, 0), yMat.at<double>(3, 0));
	points_tanx.clear();
	points_tany.clear();

	CCDID = 29;
	fprintf(fp, "CCDID = =%d\n", CCDID);
	ori_X = 168.7310e-3;
	end_X = 183.0670e-3;
	Y = -193.9835e-3;
	for (int i = 0; i < pixelsPerCCD; i++)
	{
		double tanx = -Y / f;// X���췽��
		double tany = -(ori_X - ori + i * ccdlamda) / f;
		points_tanx.push_back(cv::Point2f(i, tanx));
		points_tany.push_back(cv::Point2f(i, tany));
	}
	xMat = polyfit(points_tanx, n_x);
	yMat = polyfit(points_tany, n_y);
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", xMat.at<double>(0, 0), xMat.at<double>(1, 0), xMat.at<double>(2, 0), xMat.at<double>(3, 0));
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", yMat.at<double>(0, 0), yMat.at<double>(1, 0), yMat.at<double>(2, 0), yMat.at<double>(3, 0));
	points_tanx.clear();
	points_tany.clear();

	CCDID = 30;
	fprintf(fp, "CCDID = =%d\n", CCDID);
	ori_X = 181.2203e-3;
	end_X = 195.5563e-3;
	Y = -207.7522e-3;
	for (int i = 0; i < pixelsPerCCD; i++)
	{
		double tanx = -Y / f;// X���췽��
		double tany = -(ori_X - ori + i * ccdlamda) / f;
		points_tanx.push_back(cv::Point2f(i, tanx));
		points_tany.push_back(cv::Point2f(i, tany));
	}
	xMat = polyfit(points_tanx, n_x);
	yMat = polyfit(points_tany, n_y);
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", xMat.at<double>(0, 0), xMat.at<double>(1, 0), xMat.at<double>(2, 0), xMat.at<double>(3, 0));
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", yMat.at<double>(0, 0), yMat.at<double>(1, 0), yMat.at<double>(2, 0), yMat.at<double>(3, 0));
	points_tanx.clear();
	points_tany.clear();

	CCDID = 31;
	fprintf(fp, "CCDID = =%d\n", CCDID);
	ori_X = 193.7006e-3;
	end_X = 208.0366e-3;
	Y = -193.7754e-3;
	for (int i = 0; i < pixelsPerCCD; i++)
	{
		double tanx = -Y / f;// X���췽��
		double tany = -(ori_X - ori + i * ccdlamda) / f;
		points_tanx.push_back(cv::Point2f(i, tanx));
		points_tany.push_back(cv::Point2f(i, tany));
	}
	xMat = polyfit(points_tanx, n_x);
	yMat = polyfit(points_tany, n_y);
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", xMat.at<double>(0, 0), xMat.at<double>(1, 0), xMat.at<double>(2, 0), xMat.at<double>(3, 0));
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", yMat.at<double>(0, 0), yMat.at<double>(1, 0), yMat.at<double>(2, 0), yMat.at<double>(3, 0));
	points_tanx.clear();
	points_tany.clear();

	CCDID = 32;
	fprintf(fp, "CCDID = =%d\n", CCDID);
	ori_X = 205.6168e-3;
	end_X = 219.9528e-3;
	Y = -207.5220e-3;
	for (int i = 0; i < pixelsPerCCD; i++)
	{
		double tanx = -Y / f;// X���췽��
		double tany = -(ori_X - ori + i * ccdlamda) / f;
		points_tanx.push_back(cv::Point2f(i, tanx));
		points_tany.push_back(cv::Point2f(i, tany));
	}
	xMat = polyfit(points_tanx, n_x);
	yMat = polyfit(points_tany, n_y);
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", xMat.at<double>(0, 0), xMat.at<double>(1, 0), xMat.at<double>(2, 0), xMat.at<double>(3, 0));
	fprintf(fp, "%.9lf\t%.9lf\t%.9lf\t%.9lf\n", yMat.at<double>(0, 0), yMat.at<double>(1, 0), yMat.at<double>(2, 0), yMat.at<double>(3, 0));
	points_tanx.clear();
	points_tany.clear();

	fclose(fp);
}
void makeCameraINstall() {
	double seta11 = (3.0 / 60.0 + 53.0 / 60.0 / 60.0) * M_PI / 180.0;
	double seta12 = (90.0 + 2.0 / 60.0 + 36.0 / 60.0 / 60.0) * M_PI / 180.0;
	double seta13 = (89.0 + 57.0 / 60.0 + 6.0 / 60.0 / 60.0) * M_PI / 180.0;
	double seta21 = (89.0 + 57.0 / 60.0 + 24.0 / 60.0 / 60.0) * M_PI / 180.0;
	double seta22 = (2.0 / 60.0 + 37.0 / 60.0 / 60.0) * M_PI / 180.0;
	double seta23 = (90.0 + 18.0 / 60.0 / 60.0) * M_PI / 180.0;
	double seta31 = (90.0 + 2.0 / 60.0 + 54.0 / 60.0 / 60.0) * M_PI / 180.0;
	double seta32 = (89.0 + 59.0 / 60.0 + 42.0 / 60.0 / 60.0) * M_PI / 180.0;
	double seta33 = (2.0 / 60.0 + 55.0 / 60.0 / 60.0) * M_PI / 180.0;
	
	FILE* fp = fopen("C:\\Users\\leeda\\Desktop\\cameraInstall.txt", "w");//�����body��
	fprintf(fp, "%lf\t%lf\t%lf\n", cos(seta11), cos(seta21), cos(seta31));
	fprintf(fp, "%lf\t%lf\t%lf\n", cos(seta12), cos(seta22), cos(seta32));
	fprintf(fp, "%lf\t%lf\t%lf\n", cos(seta13), cos(seta23), cos(seta33));

	fclose(fp);
}

//804����ŷ���Ǳպ�
void testeul() {
	//1�����������תJ2000ϵ�µ�λ���ٶ�ʸ��
	/*Orbit ep_J2000;
	ep_J2000.et = 413593224.356578;
	ep_J2000.X = -6501497.423;
	ep_J2000.Y = -198485.176;
	ep_J2000.Z = 3663475.00;
	ep_J2000.Xv = 3191.535;
	ep_J2000.Yv = -3655.109;
	ep_J2000.Zv = 5467.207;*/
	EpSixElements six;
	six.a = 7471114; six.e = 0.00079; six.i = 63.444; six.AscendingNode = 165.399; six.omega = 27.0; six.M = 6.262329;
	//Orbit ep_J2000;
	//fromElments2PosVel_M(six, &ep_J2000);
	//2�����ϵ��J2000����ת����
	/*SpiceDouble state[6];
	state[0] = ep_J2000.X;
	state[1] = ep_J2000.Y;
	state[2] = ep_J2000.Z;
	state[3] = ep_J2000.Xv;
	state[4] = ep_J2000.Yv;
	state[5] = ep_J2000.Zv;*/
	SpiceDouble rotationOrbit2J2000[3][3];
	//ortationFromOrbit2J2000(state, rotationOrbit2J2000);
	eul2m_c(-1.0 * rpd_c() * six.AscendingNode, -1.0 * rpd_c() * six.i, -1.0 * rpd_c() * six.omega, 3, 1, 3, rotationOrbit2J2000);//˳ʱ����ת
	//3��Body��J2000����Ԫ��
	SpiceDouble		Qbody2J2000[4];
	Qbody2J2000[1] = 239408848 * pow(2, -31);
	Qbody2J2000[2] = -1737613440 * pow(2, -31);
	Qbody2J2000[3] = 430516512 * pow(2, -31);
	Qbody2J2000[0] = 1161776384 * pow(2, -31);

	SpiceDouble rotationBody2J2000[3][3];
	q2m_c(Qbody2J2000, rotationBody2J2000);
	//4��Body�����ϵ Rbody2orbit=RJ20002orbit*Rbody2J2000
	SpiceDouble rotationBody2Orbit[3][3];
	mtxm_c(rotationOrbit2J2000, rotationBody2J2000, rotationBody2Orbit);
	//5��תΪŷ����
	SpiceDouble roll, pitch, yaw;
	m2eul_c(rotationBody2Orbit, 3, 2, 1, &yaw, &pitch, &roll);//��ת��˳�� x y z
	cout << "ŷ���ǣ� " <<"��ת��" << dpr_c() * roll << "  ������" << dpr_c() * pitch << "  ƫ����" << dpr_c() * yaw << endl;
}

//804������Ԫ���պ�
void testque() {

	//1�� ��Body�����ϵ��ŷ���ǣ��õ���ת����Rbody2orbit��
	SpiceDouble rotationBody2Orbit[3][3];
	SpiceDouble ang3_yaw = -35286e-4, ang2_pitch = -809e-4, ang1_roll = 936e-4;//��λ/��
	eul2m_c(rpd_c() * ang3_yaw, rpd_c() * ang2_pitch, rpd_c() * ang1_roll, 3, 2, 1, rotationBody2Orbit);
	//2�����ϵ��J2000����ת����Rorbit2J2000��
	EpSixElements six;
	six.a = 7471114; six.e = 0.00079; six.i = 63.444; six.AscendingNode = 165.399; six.omega = 27.0; six.M = 6.262329;
	SpiceDouble rotationOrbit2J2000[3][3];
	eul2m_c(-1.0 * rpd_c() * six.AscendingNode, -1.0 * rpd_c() * six.i, -1.0 * rpd_c() * six.omega, 3, 1, 3, rotationOrbit2J2000);//˳ʱ����ת
	//3��Body��J2000����ת����Rbody2J2000 = Rorbit2J2000 * Rbody2orbit
	SpiceDouble rotationBody2J2000[3][3];
	mxm_c(rotationOrbit2J2000, rotationBody2Orbit, rotationBody2J2000);
	//4��Body��J2000����Ԫ��
	SpiceDouble		Qbody2J2000[4];
	m2q_c(rotationBody2J2000, Qbody2J2000);
	//5��������Ԫ��
	SpiceDouble q[4];
	q[1] = 239408848 * pow(2, -31);
	q[2] = -1737613440 * pow(2, -31);
	q[3] = 430516512 * pow(2, -31);
	q[0] = -1161776384 * pow(2, -31);
	SpiceDouble Q2M[3][3];
	q2m_c(q, Q2M);
	cout << "������Ԫ����" << " q1:" << 239408848 * pow(2, -31) << " q2:" << -1737613440 * pow(2, -31) << " q3:" << 430516512 * pow(2, -31) << " q4:" << -1161776384 * pow(2, -31) << endl;

	cout << "ת����Ԫ����" << " q0:" << Qbody2J2000[0] << " q2:" << Qbody2J2000[1] << " q3:" << Qbody2J2000[2] << " q4:" << Qbody2J2000[3] << endl;
	// -0.1203(pitch) -3.4885(roll) -13.6599(yaw)
}

int main() {
	
	cout << "hello Lee, welcome to Mars!\n";

	//makeCameraINstall();
	//makeInnerPara();
	//testeul();
	//testque();
	// ����һ HRSC�����װ����
	//test_CameraInstall();

	// ���Զ� ���ݽ���
	//test_Predata_MEX_HRSC();

	// ������ �������������λģ�� ��DEM������
	//test_GeoImg_DEM();

	// ������ ���彻��ģ��
	//test_GeoStereo();

	// ������ �����ѵļ������� 
	//test_MOLA();
	//test_MOLA_grid();
	//test_dem_slopePoints();

	// ������ ƴ��ȫ��DEM��ת��lbl��ʽΪtif��ִ��һ�鼴�ɣ�������ִ�У�
	//test_DEM();

	// ������ ��ȡָ��λ�ø߳�
	//test_readGeoH();
	//read_Img_Information("E:\\IMG_i864_All\\TiePoints_calibaration\\calibration\\sensorCorrect_result\\run-DEM.tif");
	//read_Img_Information("E:\\date_keyanrenwu\\Mars\\MGS\\MOSAIC\\mc18-256.img");

	// ���԰� ��ȡ�޶�����ľ�γ�ȸ߳�����
	//read_Tab_Mola();
	// ���԰� ����dem�����mola��
	//test_IMG_Mola();

	// ���԰� ��ȡ�����չ���洢
	//read_orbit_Mola("E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\PEDR\\i818_mola\\output\\MolaPoints.txt");

	// ���԰� �ֶμ���һ����
	//test_ccc_points();

	// ���Ծ� ƥ��MOLA��DEM�ĸ߳�ֵ
	//test_match_Mola2Dem();

	// ����ʮ ����gnuplot�ű�
	//make_gnuplot_script();

	// ����ʮһ DEM����LON LAT H תΪX Y Z 
	//test_stereo_dem2xyz();

	// ����ʮ�� GDALӰ��ü�,�ȼ����ص�����
	//test_overrange();

	// ����ʮ�� GDALӰ���ز���
	//ResampleGDAL("E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\MEG128\\radius\\cutrange\\overrange_i864.tif","E:\\date_keyanrenwu\\Mars\\MGS\\MOLA\\MEG128\\radius\\resample\\overrange_i864_136.tif",13.6,13.6);
	
	//test_MOLA_dem2xyz();

	//down_sample("E:\\IMG_i871\\img\\hi871_0000_s12_resample.img", "E:\\IMG_i871\\img\\hi871_0000_s12_resample2.img", 0.1);

	// ����ʮ�� ������������
	/*predictOrbitNum();
	string filereadpath = "H:\\leeda\\date_keyanrenwu\\Mars\\MGS\\MOLA\\PEDR\\DATA2\\AP10228L.TAB";
	string filewritepath = "D:\\TEST2\\mola\\i864\\10228.txt";
	string Filepath_dem = "D:\\TEST\\IMG_i864_All\\TiePoints_calibaration\\calibration\\sensorCorrect_result\\run-DEM.tif";
	readMolaGpoints(filereadpath, filewritepath, Filepath_dem);*/

	//����ʮ�� �������Զ�λ������֤����ȡdtm�߳�
	//test_readDTM_singleMola();

	//����ʮ��
	//test15();

	return 0; 
}

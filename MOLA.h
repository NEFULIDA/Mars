#pragma once
#include <iostream>
#include <fstream>
#include <sciplot/sciplot.hpp>
#include "GeoImageIO.h"
#include "FrechetDistance.h"

using namespace sciplot;
using namespace std;

//����MOLA����

//����Ŀ���������ĵ�γ�� ���ȣ� ����Ԥ��Ĺ����
bool predictOrbitNum();

//���� �������ݵ�ɨ�踲�Ƿ�Χ
bool drawSatalliteorbit(string filepath,string outFilepath);

//��ȡDEM��Χ��MOLA���Ƶ�
//������MOLA�ļ�����·�������Ƶ����·����DEM·����
//����MOLA���ݣ�ԭʼMOLA���ݣ�
void readMolaGpoints(string filereadpath, string filewritepath, string Filepath_dem);

// ���ĳһ��DEM�����Ҷ�Ӧ��MOLA���
// ������1��DEMӰ��2��MOLA�����б�;3��MOLA Points ����γ�ȸ߳���ʽ������·��;4��MOLA Points ��XYZ��ʽ������·��
void read_IMG_Mola(string demFile, string MolaFile, string filePoints_out, string filePoints_out_xyz);

//���������Բ����ײ���ӵ㼯
void saveCritePoints(string filePoints_in, string filePoints_out);

//����NCC׼��ͳ�Ƶ�������
void calNccOffset(string file,string demfile, string slopefile, string NCCfile, string bestfile);

//��ȡNCC������������ƫ������
void calBestOffset(string MolaFile, string Nccfile, string demfile, string bestfile);

// �ֶ�ƥ�䣬����CCC�����治һ�µĵ㼯
void saveNone_CCC_points(vector<vector<GeoPoint>> points_IMG_MOLA, string demFile, string fileMOLA_filter);

/* MOLA��ɢ��������ƥ��DEMͼ��
* string molafile                   MOLA��ɢ��������
* string DEMfile                    DEMͼ��
* string offsetfile                 ��ƫ������һ�������ϵ��
* string DEMGeoCoordinatefile       ƽ�ƺ��DEM��������
* string picturefile                �̶߳Ա�ʾ��ͼ�����·����
* double& initCoefficientNcc        ��ֵһ�������ϵ��
* double& bestCoefficientNcc        ����һ�������ϵ��
*/ 
void match_Mola2Dem(string molafile, string DEMfile, string offsetfile, string DEMGeoCoordinatefile, string picturefile, double& initCoefficientNcc, double& bestCoefficientNcc);
#pragma once
#include <iostream>
#include <fstream>
#include <sciplot/sciplot.hpp>
#include "GeoImageIO.h"
#include "FrechetDistance.h"

using namespace sciplot;
using namespace std;

//解析MOLA数据

//输入目标区域中心的纬度 经度， 返回预测的轨道号
bool predictOrbitNum();

//绘制 激光数据的扫描覆盖范围
bool drawSatalliteorbit(string filepath,string outFilepath);

//读取DEM范围内MOLA控制点
//参数（MOLA文件读入路径，控制点输出路径，DEM路径）
//单条MOLA数据（原始MOLA数据）
void readMolaGpoints(string filereadpath, string filewritepath, string Filepath_dem);

// 针对某一景DEM，查找对应的MOLA点对
// 参数：1、DEM影像；2、MOLA数据列表;3、MOLA Points （经纬度高程形式）保存路径;4、MOLA Points （XYZ形式）保存路径
void read_IMG_Mola(string demFile, string MolaFile, string filePoints_out, string filePoints_out_xyz);

//保存有明显差异的撞击坑点集
void saveCritePoints(string filePoints_in, string filePoints_out);

//按照NCC准则统计点相似性
void calNccOffset(string file,string demfile, string slopefile, string NCCfile, string bestfile);

//读取NCC结果，查找最佳偏移数据
void calBestOffset(string MolaFile, string Nccfile, string demfile, string bestfile);

// 分段匹配，计算CCC，保存不一致的点集
void saveNone_CCC_points(vector<vector<GeoPoint>> points_IMG_MOLA, string demFile, string fileMOLA_filter);

/* MOLA离散点云数据匹配DEM图像
* string molafile                   MOLA离散点云数据
* string DEMfile                    DEM图像
* string offsetfile                 各偏移量的一致相关性系数
* string DEMGeoCoordinatefile       平移后的DEM地理坐标
* string picturefile                高程对比示意图（相对路径）
* double& initCoefficientNcc        均值一致相关性系数
* double& bestCoefficientNcc        最终一致相关性系数
*/ 
void match_Mola2Dem(string molafile, string DEMfile, string offsetfile, string DEMGeoCoordinatefile, string picturefile, double& initCoefficientNcc, double& bestCoefficientNcc);
#include "route.h"
#include "lib_record.h"
#include <stdio.h>
#include "time.h"
#include<map>
#include "GA.h"
#include <sstream>
#include<iostream>
#include<string>
#include<vector>
#include<stdlib.h>
#define MAX 2100
#define INF (~(0x1<<31))
#define DMAX 101 
using namespace std;
//你要完成的功能总入口
Matrix matrix[MAX][MAX];//存储顶点连接矩阵
void search_route(char *topo[MAX_EDGE_NUM], int edge_num, char *demand[MAX_DEMAND_NUM], int demand_num)
{
	//定义开始时间
//-------------------------------------------------------------------------------------------
	clock_t begin;
	clock_t end;
	begin=clock();
//-------------------------------------------------------------------------------------------
	//读topo
//------------------------------------------------------------------------------------
	//map<int,vector<int> > designated_point;//存储必经点
	int designated_point[demand_num][DMAX];
	int source[demand_num],dest[demand_num];
	double largest;//存储最大编号
	stringstream ss;	//定义字符串流
	string str="";
	char trash;			//读取到的非数字字符
	int from, to, e_index, e_weight;	//读取到的数字字符

	for (int i = 0; i < edge_num; i++)
	{
		ss.str(topo[i]);
		ss >> e_index >> trash >> from >> trash >> to >> trash >> e_weight;
		matrix[from][to].weight = e_weight;				//初始化连接矩阵
		//cout<<from<<" "<<to<<" "<<e_weight<<endl;
		matrix[from][to].num = e_index;
		if(from>largest)
			largest=from;
		if(to>largest)
			largest=to;
		ss.clear();
	}
	//largest++;					//从0到最大的顶点编号之间的个数
	//cout<<largest<<endl;
	int  designate_len[demand_num];
	for(int i=0;i<demand_num;i++)
		 designate_len[i]=0;
	int demand_index = -1;
	for ( int i = 0; i < demand_num; i++)
	{
		ss.str(demand[i]);
		ss>>demand_index>>trash>>source[i]>>trash>>dest[i]>>trash>>str;
		ss.clear();ss.str(str);
		while (getline(ss,str,'|'))
		{
			designated_point[i][designate_len[i]++]=atoi(str.c_str());
		}
		ss.clear();
	}/**/
	GA* test;
	test = new GA();
	for(int i = 0; i<largest+1; i++)
			for(int j = 0; j<largest+1; j++)
			{
				if(matrix[i][j].weight==0)
					test->matrix[i][j].weight=INF;
				else
					test->matrix[i][j]=matrix[i][j];
			}
		//初始化必经点，写入遗传类中
	for(int i=0;i<demand_num;i++)
	{
		for(int j = 0; j < designate_len[i]; j++)
			test->designated_point_tran[i].push_back(designated_point[i][j]);
		test->source_tran[i]=source[i];
		test->dest_tran[i]=dest[i];
		test->designate_len[i]=designate_len[i];
	}
	test->largest=largest;
	test->whole_begin=begin;
	test->whole_end=end;
	//开始运行遗传算法
	test->runGA();	
	
    unsigned short result1[] = {0, 1, 2};//P'路径
    unsigned short result2[] = {5, 6, 2};//P''路径

    for (int i = 0; i < 3; i++)
    {
        record_result(WORK_PATH, result1[i]);
        record_result(BACK_PATH, result2[i]);
    }
}

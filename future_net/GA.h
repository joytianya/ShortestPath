#ifndef GA_H
#define GA_H
#include <iostream>
#include <fstream>
#include <limits>
#include <ctime>
#include <cmath>
#include <map>
#include <vector>
#include <string>
#include <algorithm>
#include <stdlib.h>
#include<math.h>
#include<queue>
#include<stack>
#include "time.h"
using namespace std;
#define MAX 2100
#define INF    (~(0x1<<31))        // 无穷大(即0X7FFFFFFF)
#define DMAX 101 
struct Matrix
{
	int num=0;
	int weight=INF;
public:
	Matrix &operator=(Matrix &other)
	{
		num=other.num;
		weight=other.weight;
	}
};
struct module_point
{
	int receive_vex=-1;
	int send_vex=-1;
	module_point(int i, int j ):receive_vex(i),send_vex(j){}
};
//在dfs去环时链表初始路径用到
struct point
{
	int next;
	point(int i):next(i){}
};
struct designated_node
{
	int value;
	int pre;
	int next;
	designated_node(int i,int j, int k):value(i),pre(j),next(k){}
};
//染色体结构
struct individual
{
public:
	bool operator<( const individual &other)
const
	{
		//return gene_designed_num>other.gene_designed_num;
		return length < other.length;
	}
	friend ostream& operator<<( ostream& out,const individual &other)
	{
		out<<"length:"<<other.length<<"\tpath:";
		for(int i=0;i<other.path.size();i++)
			out<<other.path[i]<<" ";
		return out;
	}
	double length;
	double gene_designed_num;
	vector<int> path;
};

//遗传算法
class GA
{
public:
	GA(double pm=0.15,double ps=0.75,int individual_num=50,int gen_max=1000);

	~GA();
	void runGA();//开始运算法 初始化种群 并进化 遗传算法过程
	void DFS();	//深度搜索遍历图
	int dijkstra(int vs,int vd,int prev[],int dist[]);	//寻找起点和终点之间的最短路径
	int dijkstra_search(vector<int> designated_point_para); //寻找下一个最近必经点的路径
	int dijkstra_search_combine(vector<int> designated_point_para,vector<int> limit_vec,int largest_temp);//为了出现多个集团时进行切割用
	int path_break_to_combine(map<int,vector<int> > all_path,int num_path);//出现两个集团时进行切割
	int dijkstra_repair(vector<int> d,int vs, int vd,int prev[],int dist[]);//补全首位
	void Shake_var(int k);//抖动变异
private:
	void Init_Graph();//初始化
	void Init_Groups(int k);//初始华一个种群
	void Evolution();//进化过程 调用下面三个进化操作
	//void Select();
	//void Cross();
	//void Varation();
	//void evaluate_individual();//评估函数
private://功能辅助函数
	void CalCulate_length(individual &p,int num_path);//个体的代价计算 即路径长度计算
	int repass_edge(vector<int> path_1,vector<int> path_2);//计算重边个数
	void Save_best_individual();//存储最好的个体 暂时最优解
	//int Find(const individual&, int);//交叉操作的辅助函数 查找子代是否已经存在某节点	
	void DFS(int i,int *visited);
	int firstVertex(int v);
	int nextVertex(int v,int w);
public:
	short source, dest;			//源节点和目的节点
	//vector<int> designated_point_1,designated_point_2; //存储必经点
	vector<int> designated_point;
	map<int, vector<int> >designated_point_tran;
	int source_tran[2],dest_tran[2];
	int  designate_len[2];
	double largest; 	//节点最大编号
	Matrix matrix[MAX][MAX];
	int path_all[2][MAX];
	vector<int> designate_point_order;//出现一个集团时，存储
	vector<int>	designate_point_order_break;//出现多个集团时，存储
	int pre_refresh[MAX]={0},dist_refresh[MAX]={0};
	map<int,vector<individual> > groups; //群体 即部分解空间
	vector<int> last_node_piece;
	clock_t whole_begin;//程序运行总时间记录
	clock_t whole_end;
	clock_t search_begin;
	clock_t search_end;
	clock_t var_begin;
	clock_t var_end;
	clock_t begin_stop;//可行解间隔时间记录
	clock_t end_stop;
	double duration;
	//time_t begin_first_stop;
	//time_t end_first_stop;
private: //相关属性
	double pm;
	double ps;
	int individual_num;//群体大小
	int gen_max;	//繁殖代数
	int designated_point_num;//必经点个数
	//剪枝时使用
	double path_point_num;//路径的节点数
	double ratio;//当前的必经点与所有点之间的比值
	int continue_no_designed=0;//连续没有必经点的个数
	int pre[MAX];//保存上一节点
	int num;//染色体的数目
	int visited[MAX];//标志节点是否访问过，1代表访问过，0代表没有
	//选择概率
	double change_individual[MAX];
	map<int,individual> best_individual;//保留最优解
	
	int pre_designed_num=0;
	int Min;
	int Min_temp[MAX];
};

#endif	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
		

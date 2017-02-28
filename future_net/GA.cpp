#include"GA.h"	//本cpp的头文件
//为遗传算法的初始华参数列表
GA::GA(double pm,double ps,int individual_num,int gen_max)://构造函数 初始化列表
pm(pm),
ps(ps),
individual_num(individual_num),
gen_max(gen_max)
{
	for(int i=0;i<2;i++)
		best_individual[i].length=numeric_limits<double>::max();//最大值
}
//---------------------------------------------------------------------------------
/*//DFS递归过程
int GA::firstVertex(int v)
{
	int i;
	if(v<0||v>largest-1)
		return -1;
	for(i=0;i<largest;i++)
		if(matrix[v][i].weight!=INF)
			return i;
	return -1;
}

int GA::nextVertex(int v,int w)
{
	int i;
	if(v<0||v>largest-1||w<0||w>largest-1)
		return -1;
	for(i=w+1; i<largest;i++)
		if(matrix[v][i].weight!=INF)
			return i;
	return -1;
}
void GA::DFS(int i,int *visited)
{
	
	int w;
	vector<int> temp,temp_pos;
	visited[i]=1;
	if(matrix[i][dest].weight!=INF)
	{
		pre[dest]=i;
	   	int pre_temp=dest,num_designated_path=0;
		while(pre_temp!=-1)
		{
			for(int i=0;i<designated_point_num;i++)
			{
				if(pre_temp==designated_point[i])
					num_designated_path++;
			}
			pre_temp=pre[pre_temp];
		}
		
		if(num_designated_path==designated_point_num)
		{
			pre_temp=dest;
			int min=0;
			while(pre[pre_temp]!=-1)
			{
				min+=matrix[pre[pre_temp]][pre_temp].weight;
				pre_temp=pre[pre_temp];
			}
			if(min<Min)
			{
				Min=min;
				pre_temp=dest;
				while(pre_temp!=-1)
				{
					temp.push_back(pre_temp);
					pre_temp=pre[pre_temp];
				}
				while(!temp.empty())
				{
					temp_pos.push_back(temp.back());
					temp.pop_back();
				}
				individual p;
				p.path=temp_pos;
				CalCulate_length(p);//计算代价
				groups.push_back(p);//存储在群体向量中
			}
		}
	}
	w=firstVertex(i);
	while(w>=0)
	{
		if(!visited[w]&&w!=dest)
		{	
			Min_temp[w]=Min_temp[i]+matrix[i][w].weight;
			if(Min_temp[w]<Min)
			{
				
				pre[w]=i;
				DFS(w,visited);
			}
		}
		w=nextVertex(i,w);
		if(w==-1)
		{
			visited[i]=0;
		}
	}
}

void GA::DFS()
{
	int i;
	for(i=0;i<MAX;i++)
	{
		visited[i]=0;
		pre[i]=-1;
		Min_temp[i]=0;
	}
	Min=INF;
	DFS(source,visited);
}*/

//----------------------------------------------------------------------------------------

/**
	运行算法 初始化种群 进化 输出最优个体 最优解
**/
void GA::runGA() //运行算法
{
	//designated_point_num=designated_point.size();
	//DFS();		//DFS暴搜加剪枝
	search_begin=clock();
	designated_point=designated_point_tran[1];
	source=source_tran[1];
	dest=dest_tran[1];
	Init_Groups(1);//初始化群体
	cout<<"下一条路径"<<endl;
	vector<int> weight_temp;
	for(int k=0;k<groups[1].size();k++)
	{
		for(int i=0;i<groups[1][k].path.size()-1;i++)
		{
			//weight_temp.push_back(matrix[groups[1][k].path[i]][groups[1][k].path[i+1]].weight);
			matrix[groups[1][k].path[i]][groups[1][k].path[i+1]].weight=INF;
		}
		search_begin=clock();
		designated_point=designated_point_tran[0];
		source=source_tran[0];
		dest=dest_tran[0];
		Init_Groups(0);//初始化群体
		//for(int i=0;i<groups[1][k].path.size()-1;i++)
		//{
			//matrix[groups[1][k].path[i]][groups[1][k].path[i+1]].weight=weight_temp[i];
		//}
		//time(&var_begin);
		//Shake_var(0);
		//while(!weight_temp.empty())
			//weight_temp.pop_back();
		//cout<<endl;
	}
	
	int min_repass_edge=INF;
	for(int i=0;i<groups[0].size();i++)
		for(int j=0;j<groups[1].size();j++)
			if(repass_edge(groups[0][i].path,groups[1][j].path)<min_repass_edge)
			{
				min_repass_edge=repass_edge(groups[0][0].path,groups[1][0].path);
				best_individual[0]=groups[0][0];
				best_individual[1]=groups[1][0];
			}
	cout<<"重边最少"<<min_repass_edge<<endl;
	cout<<"第一条路径权值"<<best_individual[0].length<<endl;
	cout<<"第二条路径权值"<<best_individual[1].length<<endl;
	cout<<"总权值："<<best_individual[0].length+best_individual[1].length<<endl;
	//cout<<"最优个体"<<endl;
	//for(int i=0;i<2;i++)
	//{
		//cout<<best_individual[i]<<endl;
	//}/**/
	
	//Evolution();//进化过程
	//cout<<best_individual;//输出最有个体
}
//------------------------------------------------------------------------------
int GA::dijkstra_search(vector<int> designated_point_para)
{
	int i=0,j=0,k=0,n=0;
	//新建动态二维数组
	int** prev = new int*[designated_point_para.size()];
	int** dist = new int*[designated_point_para.size()];
	int** flag = new int*[designated_point_para.size()];
	for(i = 0; i < designated_point_para.size(); i++)
	{
		prev[i]=new int[MAX];
		dist[i]=new int[MAX];
		flag[i]=new int[MAX];
	}
	//新建一维动态数组
	int* min = new int[designated_point_para.size()];
	int* tmp = new int[designated_point_para.size()];
	int* point_found = new int[designated_point_para.size()];
	//初始化访问标志，必经点作为其他所有点的前驱顶点，最短长度为必经点到其他所有点的长度
	//所有必经点并行进行
	for(i = 0; i < designated_point_para.size(); i++)
		for(j=0;j<largest+1;j++)
		{
			flag[i][j]=0;
			prev[i][j]=designated_point_para[i];		//顶点i的前驱顶点为0；
			dist[i][j]=matrix[designated_point_para[i]][j].weight;
		}
	//对“顶点vs”自身进行初始化，起点和终点不能访问
	for(i = 0; i < designated_point_para.size(); i++)
	{
		flag[i][designated_point_para[i]]=1;
		dist[i][designated_point_para[i]]=0;
		flag[i][source]=1;
		flag[i][dest]=1;
	}
	int path_all[MAX];
	//路径中已经访问过的点
	vector<int> pre_point;
	//初始化整体路径记录前驱数组
	for(i=0;i<MAX;i++)
		path_all[i]=-1;
	//整体路径中已经找到伙伴的必经点集合
	vector<int> visited_origin;
	//循环中的标记
	int flag1,flag2,flag4=0,flag_j=0,flag3=0;
	//存储上一个必经点和下一个必经点位置，为了寻找整体路径的头和尾
	int pre_path_point[MAX],next_path_point[MAX];
	for(i=0;i< designated_point_para.size(); i++)
	{
		pre_path_point[i]=-1;
		next_path_point[i]=-1;
	}
	//最多遍历largest次；每次找到出一个顶点的最短路径
	for(i=0;i<largest;i++)
	{
		//每个必经点寻找自己的小伙伴
		for(int m = 0; m < designated_point_para.size(); m++)
		{
		//奇数和偶数位置的必经点步数不同，增加多样性
			int n_bushu=0;
			if(m%2==0)
				n_bushu=1;
			else
				n_bushu=2;
		//cout<<m<<" ";
		//进入寻找，走一步或者两步
			while(n_bushu--)
			{
				//如果必经点已经找到自己的伙伴，将其位置存入visited_origin容器中，标记不再访问
				flag2=0;
				for(int x=0; x<visited_origin.size(); x++)
					if(m==visited_origin[x])
						flag2=1;
				//如果该必经点还没有找到自己的伙伴
				if(flag2!=1)
				{
					flag1=0;
					min[m]=INF;//初始化dj最短路径
					flag_j=0;//标志着能找到dj的下一个点
					//遍历所有点，找到此次最小距离点
					for(j=0;j<largest+1;j++)
					{	
						//if(designated_point_para[m]==1923)
								//cout<<"1923的下一个值"<<dist[m][1244]<<endl;
						if(flag[m][j]==0&&dist[m][j]<min[m])
						{
							
							min[m]=dist[m][j];
							point_found[m]=j;
							flag_j=1;
						}
					}
					
					//如果能找到下一点
					if(flag_j==1)
					{
						flag3=0;
						flag[m][point_found[m]]=1;//标记顶点point_found[m]为已经获取到的最短路径
						//判断这次的最近点是否是必经点中的点
						for(int p = 0; p < designated_point_para.size(); p++)
						{
							//如果是必经点
							if(point_found[m]==designated_point_para[p])
							{
								pre_path_point[p]=m;//该必经点位置p的前驱必经点位置是m
								next_path_point[m]=p;//必经点位置m的后继必经点位置是p
								//查看此段是否存在之前路径中的重复点
								//cout<<"发散点"<<designated_point_para[m];
								//cout<<"伙伴点"<<designated_point_para[p]<<endl;
								int next=designated_point_para[p];
								flag4=0;//标志有无重复点
								while(next!=designated_point_para[m])
								{
									for(int q=0; q< largest+1; q++)
										if(prev[m][next]==path_all[q])
										{   
											flag4=1;
											break;
										}
									if(flag4==1)
										break;
									next=prev[m][next];
								}
								//如果有重复点则进行dj修复
								int dj_flag=0;//标记这如果有重复点，是否可以dj修复
								if(flag4==1)
								{	
								//-------------------------------------------------------------------
									for(int i=0;i<designated_point_para.size();i++)//bug更改点
										if(i!=p)
											pre_point.push_back(designated_point_para[i]);
									if(dijkstra_repair(pre_point,designated_point_para[m], designated_point_para[p], pre_refresh, dist_refresh)==1)	
									{	
										dj_flag=1;
										//存储dj修复后的路径段
										int first=designated_point_para[p];
										//-------------------------------------------------------
										while(first!=designated_point_para[m])
										{
											//cout<<first<<" ";
											path_all[first]=pre_refresh[first];
											pre_point.push_back(first);
											first=pre_refresh[first];
										}
										//cout<<first<<" "<<endl;
									}
									else
										flag3=1;
								}
								//如果没有重复点，则按照正常处理
								else
								{
									dj_flag=1;
									//储存路径段
									next=designated_point_para[p];
									while(next!=designated_point_para[m])
									{
										//cout<<next<<" ";
										path_all[next]=prev[m][next];
										pre_point.push_back(next);
										next=prev[m][next];
									}
									//cout<<next<<" "<<endl;
								 }
								//如果路径段确定存储后
								if(dj_flag==1)
								{
									visited_origin.push_back(m);//将该寻找伙伴的必经点关闭
									//同理，将关闭其他必经点寻找这个伙伴（必经点）
									for(int q = 0; q < designated_point_para.size(); q++)
										flag[q][point_found[m]]=1;
									//开始寻找路径的头和尾
									//寻找路径的尾
									int last_point=p;
									while(pre_path_point[last_point]!=-1)
									{
										last_point=pre_path_point[last_point];
									}
									//寻找路径的头
									int first_point=m;
									while(next_path_point[first_point]!=-1)
									{
										first_point=next_path_point[first_point];
									}
									//找到之后，来使得头不能访问尾
									flag[first_point][designated_point_para[last_point]]=1;
									flag1=1;
									n++;//记录有多少段
								}
								//如果找不到这段路径则重置前驱和后继必经点
								else
								{
									pre_path_point[p]=-1;
									next_path_point[m]=-1;
								}
							}
						}
						
						if(flag1==1)
							continue;
						if(flag3==0)//如果最短距离点是必经点，有充分点，但是没有dj修复成功//或者最短距离点不是必经点
						{
							for(j=0;j<largest+1;j++)
							{	
								tmp[m]=(matrix[point_found[m]][j].weight==INF?INF:(min[m]+matrix[point_found[m]][j].weight));
								if(flag[m][j]==0&&(tmp[m]<dist[m][j]))
								{
									dist[m][j]=tmp[m];
									prev[m][j]=point_found[m];
								}
							}
						}
					}
				}
			}
		}
	}
	//释放动态二维数组
	//先释放第二维数组
	for(int i=0; i < designated_point_para.size(); i++)
	{
		delete[] prev[i];
		delete[] dist[i];
		delete[] flag[i];
	}
	//再释放第一维数组
		delete[] prev;
		delete[] dist;
		delete[] flag;
		delete[] min;
		delete[] tmp;
		delete[] point_found;
	//寻找终点必经点
	vector<int> path_last_node;
	int last_node=-1;
	for(i=0;i<designated_point_para.size();i++)
	{	
		int flag_last=0;
		for(j=0;j<largest+1;j++)
		{
			if(path_all[j]==designated_point_para[i])
			{	
				flag_last=1;
				break;
			}
		}
		if(flag_last==0)
		{
			last_node=designated_point_para[i];
			path_last_node.push_back(last_node);
			//break;
		}
	}
	
	cout<<"几个集团"<<path_last_node.size()<<" "<<endl;
	//将完整必经点路径段存入designate_point_order中
	if(path_last_node.size()==1)
	{
		last_node=path_last_node.back();
		while(last_node!=-1)
		{	
			designate_point_order.push_back(last_node);
			//cout<<last_node<<" ";
			last_node=path_all[last_node];
		}
	}
	else
	{
		last_node=-1;
		map<int,vector<int> > all_path_temp,all_path;
		//vector<int> aaa;
		for(i=0;i<path_last_node.size();i++)
		{
			last_node=path_last_node[i];
			while(last_node!=-1)
			{	
				all_path_temp[i].push_back(last_node);
				//cout<<last_node<<" ";
				last_node=path_all[last_node];
			}
			while(!all_path_temp[i].empty())
			{
				all_path[i].push_back(all_path_temp[i].back());
				//aaa.push_back(all_path_temp[i].back());
				//cout<<all_path_temp[i].back()<<" ";
				all_path_temp[i].pop_back();
			}
			//cout<<"结束"<<endl;
		}
		/*int prev[MAX]={0},dist[MAX]={0};
		vector<int> bbb;
		for(auto &cc:aaa)
			if(cc!=all_path[1].front())
				bbb.push_back(cc);
		cout<<dijkstra_repair(bbb, all_path[0].back(), all_path[1].front(),prev,dist)<<endl;
		int prev_temp= all_path[1].front();
		while(prev_temp!=all_path[0].back())
		{
			cout<<prev_temp<<" ";
			prev_temp=prev[prev_temp];
		}*/
		//cout<<all_path[0].back()<<endl;
		int flag_combine=path_break_to_combine(all_path,path_last_node.size());
	}/**/
		
}

int GA::path_break_to_combine(map<int,vector<int> > all_path,int num_path)
{
	//vector<int> path_1(all_path[0]),path_2(all_path[1]);//初始化需要切割的两条路径
	
	map<int, vector<int> > path_designated_location;
	for(int p=0;p<num_path;p++)
	{
		for(int i=0;i<all_path[p].size();i++)
		{
			for(int j=0;j<designated_point.size();j++)
			{
				if(all_path[p][i]==designated_point[j])
					path_designated_location[p].push_back(i);
			}
		}
	}
	
	double* length_all=new double[num_path];
	//double length_all[2]={0,0};
	//cout<<"集团数"<<num_path;
	for(int p=0;p<num_path;p++)
	{
		length_all[p]=0;
		for(int i=0;i<all_path[p].size()-1;i++)
			length_all[p]+=matrix[all_path[p][i]][all_path[p][i+1]].weight;
			//cout<<matrix[all_path[p][i]][all_path[p][i+1]].weight<<" ";
		//cout<<"集团的长度"<<length_all[p]<<" "<<endl;
	}
	//cout<<"KK"<<length_all[0];
	double** change_temp=new double* [num_path];
	//double** change=new double* [num_path];
	for(int i=0;i<num_path;i++)
	{
			change_temp[i]=new double [path_designated_location[i].size()-1];
		//change[i]=new double [path_designated_location[i].size()-1];
	}
	for(int i=0;i<num_path;i++)
		for(int j=0;j<path_designated_location[i].size()-1;j++)
			change_temp[i][j]=0;
	for(int i=0;i<num_path;i++)
	{	if(length_all[i]!=0)
		{	
			for(int j=0;j<path_designated_location[i].size()-1;j++)
			{
				double length_temp=0;
				for(int k=path_designated_location[i][j];k<path_designated_location[i][j+1];k++)
				{
					length_temp+=matrix[all_path[i][k]][all_path[i][k+1]].weight;
					//cout<<matrix[all_path[i][k]][all_path[i][k+1]].weight<<" ";
				}
				change_temp[i][j]=length_temp/length_all[i];
				//cout<<change_temp[i][j]<<" ";
			}
			//cout<<"结束"<<endl;
		}
	}
	//double temp=0;
	//for(int i=0;i<path_designated_location[0].size()-1;i++)
		//temp+=change_temp[0][i];
		//cout<<"总的概率"<<temp<<" "<<endl;
	int* break_point=new int[num_path];
	double* r=new double[num_path];  
	for(int i=0;i<num_path;i++)
	{
		break_point[i]=-1;
		r[i]=0;
	}
	//r=rand()/(RAND_MAX+1.0);
	//for(int i=0;i<num_path;i++)
		
	for(int i=0;i<num_path;i++)
	{	
		if(length_all[i]!=0)
		{	
			while(break_point[i]==-1)
			{
				r[i]=rand()/(RAND_MAX+1.0);
				double m=0;
				for(int j=0;j<path_designated_location[i].size()-1;j++)
				{
					m=m+change_temp[i][j];
					if(r[i]<=m)
					{
						break_point[i]=j;
						//cout<<j<<" ";
						break;
					}
				}
			}
			//cout<<endl;
		}
	}
	
	map <int ,vector<int> > module_path;
	vector<int> limit_vec;
	int num_module=0;
	for(int i=0;i<num_path;i++)
	{
		if(length_all[i]!=0)
		{	
			for(int j=0;j<=path_designated_location[i][break_point[i]];j++)
			{
				module_path[num_module].push_back(all_path[i][j]);
				limit_vec.push_back(all_path[i][j]);
			}
			num_module++;
			for(int j=path_designated_location[i][break_point[i]+1];j<all_path[i].size();j++)
			{
				module_path[num_module].push_back(all_path[i][j]);
				limit_vec.push_back(all_path[i][j]);
			}
			num_module++;
		}
		else
		{
			//cout<<"该集团只有一个点"<<all_path[i][0];
			module_path[num_module].push_back(all_path[i][0]);
			limit_vec.push_back(all_path[i][0]);
			num_module++;
		}
	}
	//cout<<"模块数"<<num_module<<" ";
	int* add_point=new int [num_module];
	vector<int> search_point;
	for(int i=0;i<num_module;i++)
	{
		add_point[i]=largest+1+i;
		search_point.push_back(add_point[i]);
	}
	for(int i=0;i<num_module;i++)
	{
		for(int j=0;j<largest+num_module+1;j++)
		{
			matrix[j][add_point[i]]=matrix[j][module_path[i].front()];
			//if(matrix[j][add_point[i]].weight!=INF)
				//cout<<j<<" ";
		}
		//cout<<endl;
	}
	//cout<<endl;
	for(int i=0;i<num_module;i++)
	{
		for(int j=0;j<largest+num_module+1;j++)
		{
			matrix[add_point[i]][j]=matrix[module_path[i].back()][j];
			//if(matrix[add_point[i]][j].weight!=INF)
				//cout<<j<<" ";
		}
		//cout<<endl;
	}

	random_shuffle(search_point.begin(),search_point.end());
	//cout<<largest+num_module<<endl;
	//vector<int> limit_vec1;
	int flag_combine=dijkstra_search_combine(search_point,limit_vec,largest+num_module);
	cout<<"修复成功"<<flag_combine<<" "<<endl;
		
	
	int* add_path_location=new int [num_module];
	vector<int> add_order;
	vector<int> designate_point_order_break_pos;
	vector<int> path_all;
	map<int, vector<int> > combine_path;
	int n_order_add=0;
	if(flag_combine)
	{	
		while(!designate_point_order_break.empty())
		{
			designate_point_order_break_pos.push_back(designate_point_order_break.back());
			//cout<<designate_point_order_break.back()<<" ";
			designate_point_order_break.pop_back();
		}
		for(int i=0;i<designate_point_order_break_pos.size();i++)
		{
			for(int j=0;j<num_module;j++)
			{
				if(designate_point_order_break_pos[i]==add_point[j])
				{
					add_path_location[n_order_add++]=i;
					
					add_order.push_back(j);
					break;
				}
			}
		}
		//cout<<n_order_add<<" ";
		for(int i=0;i<num_module-1;i++)
		{
			for(int j=add_path_location[i]+1;j<add_path_location[i+1];j++)
			{
				combine_path[i].push_back(designate_point_order_break_pos[j]);
			}
			//cout<<endl;
		}
		vector<int>	designate_point_order_temp;
		for(int i=0;i<num_module;i++)
		{
			for(int j=0;j<module_path[add_order[i]].size();j++)
				designate_point_order_temp.push_back(module_path[add_order[i]][j]);
			for(int j=0;j<combine_path[i].size();j++)
				designate_point_order_temp.push_back(combine_path[i][j]);
		}
		while(!designate_point_order_temp.empty())
		{
			designate_point_order.push_back(designate_point_order_temp.back());
			//cout<<designate_point_order_temp.back()<<" ";
			designate_point_order_temp.pop_back();
		}
	}	
	for(int i=0;i<num_path;i++)
	{
		delete[] change_temp[i];
	}
	delete[] change_temp;
	delete[] break_point;
	delete[] r;
	delete[] add_point;
	delete[] add_path_location;
	if(flag_combine==1)
		return 1;
	else/**/
	 	return 0;
	
}
int GA::dijkstra_search_combine(vector<int> designated_point_para,vector<int> limit_vec,int largest_temp)
{
	//cout<<largest_temp<<" ";
	int i=0,j=0,k=0,n=0;
	//新建动态二维数组
	int** prev = new int*[designated_point_para.size()];
	int** dist = new int*[designated_point_para.size()];
	int** flag = new int*[designated_point_para.size()];
	for(i = 0; i < designated_point_para.size(); i++)
	{
		prev[i]=new int[MAX];
		dist[i]=new int[MAX];
		flag[i]=new int[MAX];
	}
	//新建一维动态数组
	int* min = new int[designated_point_para.size()];
	int* tmp = new int[designated_point_para.size()];
	int* point_found = new int[designated_point_para.size()];
	//初始化访问标志，必经点作为其他所有点的前驱顶点，最短长度为必经点到其他所有点的长度
	//所有必经点并行进行
	for(i = 0; i < designated_point_para.size(); i++)
		for(j=0;j<largest_temp+1;j++)
		{
			flag[i][j]=0;
			prev[i][j]=designated_point_para[i];		//顶点i的前驱顶点为0；
			dist[i][j]=matrix[designated_point_para[i]][j].weight;
		}
	//对“顶点vs”自身进行初始化，起点和终点不能访问
	for(i = 0; i < designated_point_para.size(); i++)
	{
		flag[i][designated_point_para[i]]=1;
		dist[i][designated_point_para[i]]=0;
		flag[i][source]=1;
		flag[i][dest]=1;
		for(j=0;j<limit_vec.size();j++)
			flag[i][j]=1;
	}
	int path_all[MAX];
	//路径中已经访问过的点
	vector<int> pre_point;
	//初始化整体路径记录前驱数组
	for(i=0;i<MAX;i++)
		path_all[i]=-1;
	//整体路径中已经找到伙伴的必经点集合
	vector<int> visited_origin;
	//循环中的标记
	int flag1,flag2,flag4=0,flag_j=0,flag3=0;
	//存储上一个必经点和下一个必经点位置，为了寻找整体路径的头和尾
	int pre_path_point[MAX],next_path_point[MAX];
	for(i=0;i< designated_point_para.size(); i++)
	{
		pre_path_point[i]=-1;
		next_path_point[i]=-1;
	}
	//最多遍历largest次；每次找到出一个顶点的最短路径
	for(i=0;i<largest_temp;i++)
	{
		//每个必经点寻找自己的小伙伴
		for(int m = 0; m < designated_point_para.size(); m++)
		{
		//奇数和偶数位置的必经点步数不同，增加多样性
			
			int n_bushu=0;
			if(m%2==0)
				n_bushu=1;
			else
				n_bushu=2;
		//cout<<m<<" ";
		//进入寻找，走一步或者两步
			while(n_bushu--)
			{
				//如果必经点已经找到自己的伙伴，将其位置存入visited_origin容器中，标记不再访问
				flag2=0;
				for(int x=0; x<visited_origin.size(); x++)
					if(m==visited_origin[x])
						flag2=1;
				//如果该必经点还没有找到自己的伙伴
				if(flag2!=1)
				{
					flag1=0;
					min[m]=INF;//初始化dj最短路径
					flag_j=0;//标志着能找到dj的下一个点
					//遍历所有点，找到此次最小距离点
					for(j=0;j<largest_temp+1;j++)
					{	
						//if(designated_point_para[m]==1923)
								//cout<<"1923的下一个值"<<dist[m][1244]<<endl;
						if(flag[m][j]==0&&dist[m][j]<min[m])
						{
							//cout<<j<<" ";
							min[m]=dist[m][j];
							point_found[m]=j;
							flag_j=1;
						}
					}
					//cout<<point_found[m]<<endl;
					//如果能找到下一点
					if(flag_j==1)
					{
						//cout<<point_found[m]<<endl;
						flag3=0;
						flag[m][point_found[m]]=1;//标记顶点point_found[m]为已经获取到的最短路径
						//判断这次的最近点是否是必经点中的点
						for(int p = 0; p < designated_point_para.size(); p++)
						{
							
							//如果是必经点
							if(point_found[m]==designated_point_para[p])
							{
								//cout<<"进来了";
								pre_path_point[p]=m;//该必经点位置p的前驱必经点位置是m
								next_path_point[m]=p;//必经点位置m的后继必经点位置是p
								//查看此段是否存在之前路径中的重复点
								//cout<<"发散点"<<designated_point_para[m]<<endl;
								//cout<<"伙伴点"<<designated_point_para[p]<<endl;
								int next=designated_point_para[p];
								flag4=0;//标志有无重复点
								while(next!=designated_point_para[m])
								{
									for(int q=0; q< largest_temp+1; q++)
										if(prev[m][next]==path_all[q])
										{  
											//cout<<"有重复点"<< prev[m][next];
											flag4=1;
											break;
										}
									if(flag4==1)
										break;
									next=prev[m][next];
								}
								//如果有重复点则进行dj修复
								int dj_flag=0;//标记这如果有重复点，是否可以dj修复
								if(flag4==1)
								{	
								//-------------------------------------------------------------------
									for(int i=0;i<designated_point_para.size();i++)//bug更改点
										if(i!=p)
											pre_point.push_back(designated_point_para[i]);
									if(dijkstra_repair(pre_point,designated_point_para[m], designated_point_para[p], pre_refresh, dist_refresh)==1)	
									{	
										dj_flag=1;
										//存储dj修复后的路径段
										int first=designated_point_para[p];
										//-------------------------------------------------------
										while(first!=designated_point_para[m])
										{
											//cout<<first<<" ";
											path_all[first]=pre_refresh[first];
											pre_point.push_back(first);
											first=pre_refresh[first];
										}
										//cout<<first<<" "<<endl;
									}
									else
										flag3=1;
								}
								//如果没有重复点，则按照正常处理
								else
								{
									dj_flag=1;
									//储存路径段
									next=designated_point_para[p];
									while(next!=designated_point_para[m])
									{
										//cout<<next<<" ";
										path_all[next]=prev[m][next];
										pre_point.push_back(next);
										next=prev[m][next];
									}
									//cout<<next<<" "<<endl;
								 }
								//如果路径段确定存储后
								if(dj_flag==1)
								{
									visited_origin.push_back(m);//将该寻找伙伴的必经点关闭
									//同理，将关闭其他必经点寻找这个伙伴（必经点）
									for(int q = 0; q < designated_point_para.size(); q++)
										flag[q][point_found[m]]=1;
									//开始寻找路径的头和尾
									//寻找路径的尾
									int last_point=p;
									while(pre_path_point[last_point]!=-1)
									{
										last_point=pre_path_point[last_point];
									}
									//寻找路径的头
									int first_point=m;
									while(next_path_point[first_point]!=-1)
									{
										first_point=next_path_point[first_point];
									}
									//找到之后，来使得头不能访问尾
									flag[first_point][designated_point_para[last_point]]=1;
									flag1=1;
									n++;//记录有多少段
								}
								//如果找不到这段路径则重置前驱和后继必经点
								else
								{
									pre_path_point[p]=-1;
									next_path_point[m]=-1;
								}
							}
						}
						
						if(flag1==1)
							continue;
						if(flag3==0)//如果最短距离点是必经点，有充分点，但是没有dj修复成功//或者最短距离点不是必经点
						{
							for(j=0;j<largest_temp+1;j++)
							{	
								tmp[m]=(matrix[point_found[m]][j].weight==INF?INF:(min[m]+matrix[point_found[m]][j].weight));
								if(flag[m][j]==0&&(tmp[m]<dist[m][j]))
								{
									dist[m][j]=tmp[m];
									prev[m][j]=point_found[m];
								}
							}
						}
					}
				}
			}
		}
	}
	//释放动态二维数组
	//先释放第二维数组
	for(int i=0; i < designated_point_para.size(); i++)
	{
		delete[] prev[i];
		delete[] dist[i];
		delete[] flag[i];
	}
	//再释放第一维数组
		delete[] prev;
		delete[] dist;
		delete[] flag;
		delete[] min;
		delete[] tmp;
		delete[] point_found;
	//寻找终点必经点
	vector<int> path_last_node;
	int last_node=-1;
	for(i=0;i<designated_point_para.size();i++)
	{	
		int flag_last=0;
		for(j=0;j<largest_temp+1;j++)
		{
			if(path_all[j]==designated_point_para[i])
			{	
				//cout<<"进来了"<<designated_point_para[i];
				flag_last=1;
				break;
			}
		}
		if(flag_last==0)
		{
			last_node=designated_point_para[i];
			//cout<<"最后一个"<<designated_point_para[i];
			path_last_node.push_back(last_node);
			//break;
		}
	}
	//cout<<"有几个最后一点"<<path_last_node.size()<<endl;
	//将完整必经点路径段存入designate_point_order中
	if(path_last_node.size()==1)
	{
		last_node=path_last_node[0];
		while(last_node!=-1)
		{	
			designate_point_order_break.push_back(last_node);
			//cout<<last_node<<" ";
			last_node=path_all[last_node];
		}
	}
	if(path_last_node.size()==1)
	{
		return 1;
	}
	else 
		return 0;
	//cout<<endl;
}
//dj修复函数，在一定的限制条件下进行寻找dj路径
int GA::dijkstra_repair(vector<int> d,int vs, int vd,int prev[],int dist[])
{
	int i,j,k=-1;
	int min;
	int tmp;
	int flag[MAX];//flag[i]=1;表示“顶点vs到顶点i”的最短路径已经成功获取。
	for(i=0;i<largest+1;i++)
	{
		flag[i]=0;
		prev[i]=vs;		//顶点i的前驱顶点为0；
		dist[i]=matrix[vs][i].weight;
	}
	//cout<<"1985下一个距离"<<dist[47]<<endl;
	//对“顶点vs”自身进行初始化
	for(i=0;i<d.size();i++)
	{	flag[d[i]]=1;
		dist[d[i]]=INF;
	}
	flag[vs]=1;
	dist[vs]=0;
	int flag_j=0;
	//cout<<"起点"<<vs<<endl;
	//cout<<"1985下一个点"<<flag[47]<<endl;
	
	//遍历largest次；每次找到出一个顶点的最短路径
	for(i=0;i<largest;i++)
	{
		//寻找当前最小的路径
		flag_j=0;
		min=INF;
		for(j=0;j<largest+1;j++)
		{	
			
			if(flag[j]==0&&dist[j]<min)
			{
				min=dist[j];
				k=j;
				//cout<<j<<" ";
				flag_j=1;
			}
		}
		//标记顶点k为已经获取到的最短路径
		if(flag_j==1)
		{
			flag[k]=1;
			if(k==vd)
			{
				break;
			}
			for(j=0;j<largest+1;j++)
			{	
				tmp=(matrix[k][j].weight==INF?INF:(min+matrix[k][j].weight));
				if(flag[j]==0&&(tmp<dist[j]))
				{
					dist[j]=tmp;
					prev[j]=k;
				}
			}
		}
	}
	if(flag[vd]==1&&dist[vd]!=INF)
	{	
		return 1;
	}
	else 
		return 0;
}

//初始化一个群体
void GA::Init_Groups(int k) //初始一个群体
{
	srand(time(NULL));
	//单源随机路径生成算法
	//if(largest>500&&designated_point_num!=30)
	designated_point_num= designate_len[k];
	if(largest<=20)
	{
		while(1)
		{
			//whole_end=clock();
			
			//if((double)(whole_end - whole_begin)/CLOCKS_PER_SEC >8)	//时间达到9s就跳出
				//break;
			if(groups[k].size()==1)
				break;
			int path[MAX];
			for(int i=0;i<=largest;i++)
			{
				visited[i]=0;
				path[i]=-1;
			}
			vector<int> route;
			stack<int> s;
			s.push(source);
			visited[source]=1;
			while(!s.empty())
			{
				int v=s.top();
				int flag=0;
				vector<int> temp;
				for(int i=0;i<=largest;i++)
					if(!visited[i]&&matrix[v][i].weight!=INF)
							temp.push_back(i);
				
				if(!temp.empty())
				{
					int i=rand()%temp.size();
					path[temp[i]]=v;
					if(temp[i]==dest)
					{
						int temp1=dest;
						while(temp1!=-1)
						{
							route.push_back(temp1);
							temp1=path[temp1];
						}
						break;
					}
					s.push(temp[i]);
					visited[temp[i]]=1;
				}
				if(temp.empty())
				{
					s.pop();
				}
			}
			individual c;
			int pre_designed_num=0;
			for(int i=0;i<designate_len[k];i++)
			{
				for(int j=0;j<route.size();j++)
				{
					if(designated_point[i]==route[j])				
						pre_designed_num++;
				}
			}
			if(pre_designed_num==designate_len[k])
			{	
				
				vector<int> route_pos;
				while(!route.empty())
				{
					route_pos.push_back(route.back());
					route.pop_back();
				}
				c.path=route_pos;
				c.gene_designed_num=pre_designed_num;
				CalCulate_length(c,k);//计算代价
				vector<int> num_individual_rand;
				bool flag_repass_rand=0;
				if(groups[k].size()>0)
				{
					for(int i=0;i<groups[k].size();i++)
						if(c.length==groups[k][i].length)
						{
							num_individual_rand.push_back(i);
							flag_repass_rand=1;
						}
				}
				bool flag_diff_rand=0;
				int num_diff_rand=0;
				if(flag_repass_rand==1)
				{
					for(int i=0;i<num_individual_rand.size();i++)
						for(int j=0;j<groups[k][num_individual_rand[i]].path.size();j++)
						{
							if(groups[k][num_individual_rand[i]].path[j]!=c.path[j])
							{
								flag_diff_rand=1;
								break;
							}
						}
						if(flag_diff_rand==1)
							num_diff_rand++;
				}
				if(num_diff_rand==num_individual_rand.size())
				{
					//begin_stop=clock();
					cout<<"长度"<<c.length;
					cout<<"必经点个数"<<c.gene_designed_num<<endl;
					groups[k].push_back(c);
				}
			}
			/*end_stop=clock();
			duration=(double)(end_stop-begin_stop)/CLOCKS_PER_SEC;
			if(duration > 1)	//两个解之间时间间隔达到1s就跳出
			break;*/
		}
	}

	else
	{
		//int num_x=1;//调试时用，限制循环次数
		while(1)
		{
			//time(&begin_stop);
			//whole_end=clock();
			//if((double)(whole_end-whole_begin)/CLOCKS_PER_SEC >9)
			//	break;
			if(groups[k].size()==1)	
				break;
			//随机化必经点顺序排列
			random_shuffle(designated_point.begin(),designated_point.end());
			//进入核心算法，dj发散搜索
			dijkstra_search(designated_point);
			//得到的必经点路径段正序排列
			vector<int> designate_point_order_pos;
			while(!designate_point_order.empty())
			{
				designate_point_order_pos.push_back(designate_point_order.back());
				designate_point_order.pop_back();
			}
			//记录得到的路径中含有的必经点个数
			int num_designated_path=0;
			
			for(int i=0;i<designate_point_order_pos.size();i++)
			{
				for(int j=0;j<designate_len[k];j++)
				{
					if(designate_point_order_pos[i]==designated_point[j])
					{	
						num_designated_path++;
						break;
					}
				}
			}/**/
			//当满足含有所有的必经点个数在进入连接头和尾
			if(num_designated_path==designate_len[k])
			{	
				vector<int> first_dj(designate_point_order_pos.size()-1,-1);
				//cout<<"进来了";
				copy(designate_point_order_pos.begin()+1,designate_point_order_pos.end(),first_dj.begin());//除了必经点段的首节点，其余都放在dj访问过的容器中，以备dj修复
				//判断dj修复是否可以成功
				int prev_first[MAX]={0},dist_first[MAX]={0};
				if(dijkstra_repair(first_dj,source,designate_point_order_pos[0],prev_first,dist_first)==1)
				{
					vector<int> first_repair_temp,first_repair,last_repair_temp,last_repair;
					int next=designate_point_order_pos[0];
					while(next!=source)
					{
					
						first_repair_temp.push_back(next);
						next=prev_first[next];
					}
					first_repair_temp.push_back(source);
		
					while(!first_repair_temp.empty())
					{
						first_repair.push_back(first_repair_temp.back());
						first_repair_temp.pop_back();
					}
					//源节点与必经点段连接后存入first_combine里
					int prev_last[MAX]={0},dist_last[MAX]={0};
					for(int i=0; i<	first_dj.size();i++)
						first_repair.push_back(first_dj[i]);
					vector<int> first_combine(first_repair);
					if(dijkstra_repair(first_combine,designate_point_order_pos.back(),dest,prev_last,dist_last)==1)//如果能 连到最后一个
					{
						//cout<<"进来了";
						next=dest;
						while(next!=designate_point_order_pos.back())
						{
							last_repair_temp.push_back(next);
							next=prev_last[next];
						}
						while(!last_repair_temp.empty())
						{
							last_repair.push_back(last_repair_temp.back());
							last_repair_temp.pop_back();
						}
						//将必经点段最后一个节点连接尾节点后的段和之前源节点与必经点段的首节点连接段拼接成path_full
						vector<int> path_full(first_combine.size()+last_repair.size(),-1);
						copy(first_combine.begin(),first_combine.end(),path_full.begin());
						copy(last_repair.begin(),last_repair.end(),path_full.begin()+first_combine.size());
						//初始化个体p
						individual p;
						p.path=path_full;//个体的路径赋值
						CalCulate_length(p,k);//计算个体的长度
						//判断个体的长度是否和其他的重复
						bool flag_num_groups=0;
						vector<int> num_individual_repass_1;//存储相同长度个体的个数
						if(groups[k].size()>0)
						{
							int num_groups=groups[k].size();
						
							for(int i=0;i<num_groups;i++)
							{	
								if(p.length==(groups[k][i]).length)
								{
									flag_num_groups=1;
									num_individual_repass_1.push_back(i);
								}
							}
							//判断该个体是否和其他个体长度都不同
							int flag_diff_break_1=0,diff_num_individual_1=0;	
							if(flag_num_groups==1)
							{
								for(int i=0;i<num_individual_repass_1.size();i++)
								{
									flag_diff_break_1=0;
									for(int j=0;j<groups[k][num_individual_repass_1[i]].path.size();j++)
									{
										if(groups[k][num_individual_repass_1[i]].path[j]!=p.path[j])
										{
											flag_diff_break_1=1;
											break;	
										}
									}
									if(flag_diff_break_1==1)
										diff_num_individual_1++;
								}
							}
							//如果都不同，且该个体必经点个数也是要求的必经点个数	
							if(diff_num_individual_1==num_individual_repass_1.size()&&p.gene_designed_num==designate_len[k])//不重复则存入
							{	
								//begin_stop=clock();	
								groups[k].push_back(p);
								cout<<"长度"<<p.length;
								cout<<"必经点数"<<p.gene_designed_num<<endl;
							}
						}
						else
						{
							if(p.gene_designed_num==designate_len[k])
							{
								//begin_stop=clock();
								groups[k].push_back(p);
								cout<<"长度"<<p.length;
								cout<<"必经点数"<<p.gene_designed_num<<endl;
							}
						}
						
						
					}
					else																		//如果不能就进行切割
					{
						//概率切割
						vector<int> s_num_first;//存必经点在first_combine里的位置和源节点的位置
						s_num_first.push_back(0);//首先将源节点放在其中
						for(int i=0; i< first_combine.size(); i++)
						{
							for(int j=0;j<designated_point.size();j++)
							{
								if(first_combine[i]==designated_point[j])
									s_num_first.push_back(i);//每次遇到一个必经点则把其放在其中
							}
						}
						int *length_piece=new int[s_num_first.size()-1];//存储路径段的长度，比位置容器元素个数少一
						for(int i=0; i<s_num_first.size()-1;i++)		//对路段长度数组初始化为0
						 length_piece[i]=0;
						int n=0;				//路径段序列号
						for(int p=0; p < s_num_first.size()-1; p++)
						{
							for(int i=s_num_first[p]; i < s_num_first[p+1]; i++)
								length_piece[n]+=matrix[first_combine[i]][first_combine[i+1]].weight;
							n++;
						}
						
						//新建数组，计算不同路径段长度概率，长度大的被切割的概率大
						double* change_temp=new double[s_num_first.size()-1];
						double* change=new double[s_num_first.size()-1];
						for(int i=0;i<s_num_first.size()-1;i++)
						{
							change_temp[i]=0;
							change[i]=0;
						}
						//计算路径段总长度
						double all_length_piece=0;
						for(int i=0; i<s_num_first.size()-1;i++)
						{
							//cout<<length_piece[i]<<" ";
							all_length_piece+=length_piece[i];
						}
						//路径段切割的概率
						for(int i=0; i<s_num_first.size()-1;i++)
							change_temp[i]=length_piece[i]/all_length_piece;
						//路径段被切割的累计概率
						for(int i=0; i<s_num_first.size()-1;i++)
							for(int j=0;j<=i;j++)
								change[i]+=change_temp[j];
						delete[] length_piece;
						int dj_repair_flag1=0,dj_repair_flag2=0,dj_repair_flag3=0;
						individual p;
						int num_cut=20;
						while(num_cut--)
						{
							int n1=0,n2=0;
							while(1)
							{
								double mm=0,nn=0;
								double  r1=rand()/(RAND_MAX+1.0),r2=rand()/(RAND_MAX+1.0);
								for(int i=0;i<s_num_first.size()-1;i++)
								{
									mm=change[i];
									if(r1<=mm) 
									{
										n1=i;
										break;
									}
								}
								for(int i=0;i<s_num_first.size()-1;i++)
								{
									nn=change[i];
									if(r2<=nn)
									{ 
										n2=i;
										break;
									}
								}
								if(n1!=n2)
									break;
							}
							int num_piece=s_num_first.size()-2;
							int max1=0,max2=0,num_max1=0,num_max2=0;
							num_max1=n1;//order[n1];
							num_max2=n2;//order[n2];
							int first_break,second_break;//第一第二切割点
							//对第一和第二切割点进行排序
							if(num_max1>num_max2)
							{
								first_break=num_max2;
								second_break=num_max1;
							}
							else
							{
								first_break=num_max1;
								second_break=num_max2;
							}
							//两个切割点将期分成三段，出来两个自由断点和尾节点总共三个点
							vector<int> t1(s_num_first[first_break]+1,-1),t2(s_num_first[second_break]-s_num_first[first_break+1],-1),t3(s_num_first.back()-s_num_first[second_break+1],-1);
							copy(first_combine.begin(),first_combine.begin()+s_num_first[first_break]+1,t1.begin());
							copy(first_combine.begin()+s_num_first[first_break+1]+1,first_combine.begin()+s_num_first[second_break]+1,t2.begin());
							copy(first_combine.begin()+s_num_first[second_break+1]+1,first_combine.end(),t3.begin());
							//进行第一段dj修复
							vector<int> dj1(t1.size()+t2.size()+t3.size(),-1);
							copy(t1.begin(),t1.end(),dj1.begin());
							copy(t2.begin(),t2.end(),dj1.begin()+t1.size());
							copy(t3.begin(),t3.end(),dj1.begin()+t1.size()+t2.size());
							dj1.push_back(first_combine[s_num_first[first_break]+1]);
							dj1.push_back(dest);
							//dj1是作为第一段修复限制访问元素容器
							int prev1[MAX]={0},dist1[MAX]={0};
							dj_repair_flag1=dijkstra_repair(dj1,first_combine[s_num_first[first_break]],first_combine[s_num_first[second_break+1]],prev1,dist1);
							if(dj_repair_flag1)
							{
								int pre1=first_combine[s_num_first[second_break+1]];
								vector<int> t_first,t_first_pos;
								while(pre1!=first_combine[s_num_first[first_break]])
								{
									t_first.push_back(pre1);
									pre1=prev1[pre1];
								}
								while(!t_first.empty())
								{
									t_first_pos.push_back(t_first.back());
									t_first.pop_back();
								}
								//进行第二段修复
								vector<int> dj2(t1.size()+t2.size()+t3.size(),-1);
								copy(t1.begin(),t1.end(),dj2.begin());
								copy(t2.begin(),t2.end(),dj2.begin()+t1.size());
								copy(t3.begin(),t3.end(),dj2.begin()+t1.size()+t2.size());
								dj2.push_back(dest);
								for(int i=0;i<t_first_pos.size();i++)
									dj2.push_back(t_first_pos[i]);
								//dj2是第二段修复限制访问元素容器
								int prev2[MAX]={0},dist2[MAX]={0};
								dj_repair_flag2=dijkstra_repair(dj2,first_combine.back(),first_combine[s_num_first[first_break+1]],prev2,dist2);
								if(dj_repair_flag2)
								{
									int pre2=first_combine[s_num_first[first_break+1]];
									vector<int> t_second,t_second_pos;
									while(pre2!=first_combine.back())
									{
										t_second.push_back(pre2);
										pre2=prev2[pre2];
									}
									while(!t_second.empty())
									{
										t_second_pos.push_back(t_second.back());
										t_second.pop_back();
									}
									//进行第三段修复
									vector<int> dj3(t1.size()+t2.size()+t3.size(),-1);
									copy(t1.begin(),t1.end(),dj3.begin());
									copy(t2.begin(),t2.end(),dj3.begin()+t1.size());
									copy(t3.begin(),t3.end(),dj3.begin()+t1.size()+t2.size());
									for(int i=0;i<t_first_pos.size();i++)
									{
										dj3.push_back(t_first_pos[i]);
									}
									for(int i=0;i<t_second_pos.size();i++)
									{
										dj3.push_back(t_second_pos[i]);
									}
									//dj3是第三段修复的限制访问元素容器
									int prev3[MAX]={0},dist3[MAX]={0};
									dj_repair_flag3=dijkstra_repair(dj3,first_combine[s_num_first[second_break]],dest,prev3,dist3);
									if(dj_repair_flag3)
									{
										int pre3=dest;
										vector<int> t_third,t_third_pos;
										while(pre3!=first_combine[s_num_first[second_break]])
										{
											t_third.push_back(pre3);
											pre3=prev3[pre3];
										}
										while(!t_third.empty())
										{
											t_third_pos.push_back(t_third.back());
											t_third.pop_back();
										}
										//三段修复成功即可拼接6段路径
										vector<int> path_full1(t1.size()+t_first_pos.size()+t2.size()+t_second_pos.size()+t3.size()+t_third_pos.size(),-1);
										copy(t1.begin(),t1.end(),path_full1.begin());
										copy(t_first_pos.begin(),t_first_pos.end(),path_full1.begin()+t1.size());
										copy(t3.begin(),t3.end(),path_full1.begin()+t1.size()+t_first_pos.size());
										copy(t_second_pos.begin(),t_second_pos.end(),path_full1.begin()+t1.size()+t_first_pos.size()+t3.size());
										copy(t2.begin(),t2.end(),path_full1.begin()+t1.size()+t_first_pos.size()+t3.size()+t_second_pos.size());
										copy(t_third_pos.begin(),t_third_pos.end(),path_full1.begin()+t1.size()+t_first_pos.size()+t3.size()+t_second_pos.size()+t2.size());
							
										//拼接成功即初始化个体p
										p.path=path_full1;
										CalCulate_length(p,k);
										bool flag_groups_num_1=0;	//如果长度和之前的相同标志
										vector<int> num_individual_repass;//存重复个体的位置
										if(groups[k].size()>0)
										{
											//寻找之前的相同的个体
											for(int i=0;i<groups.size();i++)
											{
												if(p.length==(groups[k][i]).length)
												{
													flag_groups_num_1=1;
													num_individual_repass.push_back(i);
												}
											}
											//如果存在的话，来确定路径是否和之前的都不同
											int flag_diff_break=0,diff_num_individual=0;	
											if(flag_groups_num_1==1)
											{
												for(int i=0;i<num_individual_repass.size();i++)
												{
													flag_diff_break=0;
													for(int j=0;j<groups[k][num_individual_repass[i]].path.size();j++)
													{
														if(groups[k][num_individual_repass[i]].path[j]!=p.path[j])
														{
															flag_diff_break=1;
															break;
														}
													}
													if(flag_diff_break==1)
														diff_num_individual++;
												}
											}
											//如果都不同，且个体的必经点个数和要求的必经点个数相同
											if(diff_num_individual==num_individual_repass.size()&&p.gene_designed_num==designate_len[k])
											{
												//begin_stop=clock();
												groups[k].push_back(p);
												cout<<"切割后长度"<<p.length;
												cout<<"切割后必经点数"<<p.gene_designed_num<<endl;
												//cout<<groups.back();
												break;
											}
										
										}
										else
										{
										
											if(p.gene_designed_num==designate_len[k])
											{
												//begin_stop=clock();
												cout<<"切割后长度"<<p.length;
												cout<<"切割后必经点数"<<p.gene_designed_num<<endl;
												groups[k].push_back(p);
												break;
											}
										}
									
									}
								}
							}
						}/**/
						delete[] change_temp;
						delete[] change;
					}
				}
			}
			/*end_stop=clock();
			duration=(double)(end_stop-begin_stop)/CLOCKS_PER_SEC;
			if(duration > 2)	//两个解之间时间间隔达到2s就跳出
				break;*/
		}
	}
}
int GA::repass_edge(vector<int> path_1,vector<int> path_2)
{
	vector<int> path_edge_1,path_edge_2;
	for(int i=0;i<path_1.size()-1;i++)
		path_edge_1.push_back(matrix[path_1[i]][path_1[i+1]].num);
	for(int i=0;i<path_2.size()-1;i++)
		path_edge_2.push_back(matrix[path_2[i]][path_2[i+1]].num);
	int num_repass_edge=0;
	for(int i=0;i<path_edge_1.size();i++)
	{
		for(int j=0;j<path_edge_2.size();j++)
		{
			if(path_edge_1[i]==path_edge_2[j])
			{
				num_repass_edge++;
				cout<<"重复的边的编号："<<path_edge_1[i]<<" ";
				break;
			}
		}
	}
	return num_repass_edge;
}
/**
	进化过程
	选择 交叉 变异
**/
//抖动变异
void GA::Shake_var(int k)
{
	//int n1=groups.size()*5;
	//进行种群数目个分割
	begin_stop=clock();
	while(1)
	{	
		//cout<<"jin";
		var_end=clock();
		if((double)(var_end - var_begin)/CLOCKS_PER_SEC > 1)	//时间达到9s就跳出
			break;
		int num_individual=0;
		if(groups[k].size()>1)
			num_individual=rand()%groups[k].size();
		//cout<<"参与变异的个体"<<num_individual<<endl;
		int n2=designate_len[k];	//变异段的次数
		vector<int> path_v(groups[k][num_individual].path);
		//cout<<groups[k][num_individual]<<endl;
		map<int,vector<int> > path_visited_point;//存储每一个必经点变异时限制访问的点集合
		//先初始化这些必经点限制访问集合
		for(int i=0;i<designate_len[k];i++)
		{
			for(int j=0;j<path_v.size();j++)
			{
				if(path_v[j]!=designated_point[i])
				{
					path_visited_point[i].push_back(path_v[j]);
				}
			}
		}
		//进行抖动变异
		while(n2--)
		{
			int designated_point_var_order=rand()%designate_len[k];//随机确定抖动变异的必经点//有优化的空间（将source和dest加入）
			//cout<<designated_point[designated_point_var_order]<<" ";
			int break_first=-1,break_second=-1,break_second_designated=-1;
			for(int i=0;i<path_v.size();i++)
			{
				if(path_v[i]==designated_point[designated_point_var_order])
				{
					break_first=i;//确定变异的必经点在路径中的位置即变异边起端
					break;
				}
			}
			//cout<<break_first<<" ";
			//寻找下一个必经点
			bool find_break_second=0;
			for(int i=break_first+1;i<path_v.size();i++)
			{	for(int j=0;j<designate_len[k];j++)
				{	if(path_v[i]==designated_point[j])
					{
						break_second=i;
						break_second_designated=j;//确定变异边的末端
						find_break_second=1;
						break;
					
					}
				}
				if(find_break_second==1)
					break;
			}
			//如果找到下一个必经点
			if(find_break_second==1)
			{
				vector<int> path_var,path_second,path_full_dj;;
				for(int i=0;i<break_first+1;i++)
				{
					path_full_dj.push_back(path_v[i]);
				}
				for(int i=break_second+1;i<path_v.size();i++)
				{
					path_second.push_back(path_v[i]);
				}
				int prev_dj[MAX]={0},dist_dj[MAX]={0};
				vector<int> var_jj;
				for(auto &cc:path_visited_point[break_second_designated])
				{
					var_jj.push_back(cc);
				}
				
				int temp=matrix[path_v[break_first]][path_v[break_second]].weight;
				if(break_second-break_first==1)//如果两个必经点是紧邻的话，需要初始其距离为无穷大，即不联通
				{
					matrix[path_v[break_first]][path_v[break_second]].weight=INF;
				}
				//判断是否能修复成功
				int flag=dijkstra_repair(var_jj,path_v[break_first],path_v[break_second],prev_dj,dist_dj);
				if(flag==1)
				{
					//dj路径记录
					int next_dj=path_v[break_second];
					while(next_dj!=path_v[break_first])
					{
						path_var.push_back(next_dj);
						next_dj=prev_dj[next_dj];
					}
					//dj路径正向
					vector<int> path_var_pos;
					while(!path_var.empty())
					{
						path_var_pos.push_back(path_var.back());
						path_var.pop_back();
					}
					//三段路径拼接
					for(int i=0;i<path_var_pos.size();i++)
					{
						path_full_dj.push_back(path_var_pos[i]);
					}
					for(int i=0;i<path_second.size();i++)
					{
						path_full_dj.push_back(path_second[i]);
					}
					//初始化遗传个体
					individual p_var;
					p_var.path=path_full_dj;
					//还原原始权值
					if(break_second-break_first==1)
					{
						matrix[path_v[break_first]][path_v[break_second]].weight=temp;
					}
					CalCulate_length(p_var,k);//计算个体的路径和必经点个数
					//-------------------------------------------------------------------------
					//cout<<"断开的边前端"<<path_full_dj[break_first]<<" ";
					//cout<<"断开的边后端"<<designated_point[break_second_designated]<<endl;
					//for(int i=0;i<path_full_dj.size();i++)
						//cout<<path_full_dj[i]<<" ";
						//cout<<"路经长度"<<p_var.length<<endl;
					//-------------------------------------------------------------------------------
					for(int i=break_first+1;i<p_var.path.size();i++)
					{	
						if(p_var.path[i]==designated_point[break_second_designated])
						{
							break;
						}
						path_visited_point[break_second_designated].push_back(p_var.path[i]);
					}
					//存重复个体的位置
					bool flag_repass_length=0;
					vector<int> num_induvidual_repass;
					for(int i=0;i<groups[k].size();i++)
						if(groups[k][i].length==p_var.length)
						{
							num_induvidual_repass.push_back(i);
							flag_repass_length=1;
						}
					//判断是否有重复长度个体
					int flag_dif=0,diff_num_individual_var=0;
					if(flag_repass_length==1)
					{
						for(int i=0;i<num_induvidual_repass.size();i++)
						{
							flag_dif=0;
							for(int j=0;j<groups[k][num_induvidual_repass[i]].path.size();j++)
							{
								if(groups[k][num_induvidual_repass[i]].path[j]!=p_var.path[j])
								{
									flag_dif=1;
									break;
								}
							}
							if(flag_dif==1)
								diff_num_individual_var++;
						}	
					}
					//如果都不同，且个体必经点个数等于必经点个数
					if(diff_num_individual_var==num_induvidual_repass.size()&&p_var.gene_designed_num==designate_len[k])
					{	
						begin_stop=clock();
						groups[k].push_back(p_var);
						cout<<"变异后长度"<<p_var.length;
						cout<<"变异后必经点数"<<p_var.gene_designed_num<<endl;
					}
					
				}
				else
					matrix[path_v[break_first]][path_v[break_second]].weight=temp;
			}
		}
		end_stop=clock();
		duration=(double)(end_stop-begin_stop)/CLOCKS_PER_SEC;
		if(duration > 1)	//两个解之间时间间隔达到1s就跳出
			break;
	}
}	
void GA::CalCulate_length(individual &p,int num_path) //计算某个个体的代价 计算路径长
{
	p.length=0;
	int path_designed_num_temp=0;
	for(int i=0;i<p.path.size()-1;i++)
	{
		p.length+=matrix[p.path[i]][p.path[i+1]].weight;
	}
	for(int j=0;j<p.path.size();j++)
		for(int k=0;k<designate_len[num_path];k++)
			if(p.path[j]==designated_point_tran[num_path][k])
				path_designed_num_temp++;
	
	p.gene_designed_num=path_designed_num_temp;
	//cout<<"必经点个数"<<p.gene_designed_num;
	//cout<<"长度"<<p.length;
}
	

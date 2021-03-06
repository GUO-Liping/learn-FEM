#!/usr/bin/python3
# -*- coding:UTF-8 -*-
'''
Author: GUO Liping
Created by June, 2020
Function: Finite Element Method of 2D Euler beam analysis
由于FEM-2Dbar-v0.0.1.py中的单元刚度矩阵难以集成，
故在FEM-2Dbar-v0.0.2.py版本中将局部坐标系下的K矩阵以节点为单位进行分块，
尝试进一步高效集成为稀疏矩阵，但无法获得较高的集成效率
故在FEM-2Dbar-v0.0.2.py版本中直接从最小的节点自由度出发对单元刚度矩阵分块，
再进一步尝试高效集成稀疏矩阵...
程序中应用了scipy.sparse.coo_matrix稀疏矩阵生成的重要性质，即某一特定行列位置的值重复输入时，默认对输入值进行叠加
FEM-2Dbar-v0.0.4.py在FEM-2Dbar-v0.0.3.py的基础上删除了冗余代码，添加了注释信息
但仍然存在两个问题
（1）不能分析三维梁单元问题——难度系数：3星
（3）不能分析铁摩辛柯梁单元问题——难度系数：4星
（2）不能分析梁的动力学问题——难度系数：4星
（2）不能根据节点处为铰接时的自由度释放（静力凝聚）问题——难度系数：2星
（2）等效节点力的自动计算问题——难度系数：1星
（3）截面特性的自动计算问题——难度系数：1星
'''
import numpy as np
from scipy import sparse, matrix
from scipy.sparse.linalg import spsolve

def func_K_matrix_ii(E,I,A,l):
	K_matrix_row = np.array([0,1,2,2])
	K_matrix_column = np.array([0,1,1,2])
	K_matrix_value = np.array([E*A/l, 12*E*I/l**3, 6*E*I/l**2, 4*E*I/l])
	K_matrix_down = sparse.coo_matrix((K_matrix_value, (K_matrix_row,K_matrix_column)),shape=(3,3))
	K_matrix_ii = sparse.tril(K_matrix_down,k=-1) + K_matrix_down.T
	return K_matrix_ii

def func_K_matrix_ji(E,I,A,l):
	K_matrix_row = np.array([0,1,1,2,2])
	K_matrix_column = np.array([0,1,2,1,2])
	K_matrix_value = np.array([-E*A/l, -12*E*I/l**3, -6*E*I/l**2, 6*E*I/l**2, 2*E*I/l])
	K_matrix_ji = sparse.coo_matrix((K_matrix_value, (K_matrix_row,K_matrix_column)),shape=(3,3))
	return K_matrix_ji

def func_K_matrix_jj(E,I,A,l):
	K_matrix_row = np.array([0,1,2,2])
	K_matrix_column = np.array([0,1,1,2])
	K_matrix_value = np.array([E*A/l, 12*E*I/l**3, -6*E*I/l**2, 4*E*I/l])
	K_matrix_down = sparse.coo_matrix((K_matrix_value, (K_matrix_row,K_matrix_column)),shape=(3,3))
	K_matrix_jj = sparse.tril(K_matrix_down,k=-1) + K_matrix_down.T
	return K_matrix_jj

def func_freeM1_K_matrix_ii(E,I,A,l):
	K_matrix_row = np.array([0,1])
	K_matrix_column = np.array([0,1])
	K_matrix_value = np.array([E*A/l, 3*E*I/l**3])
	K_matrix_down = sparse.coo_matrix((K_matrix_value, (K_matrix_row,K_matrix_column)),shape=(3,3))
	K_matrix_ii = sparse.tril(K_matrix_down,k=-1) + K_matrix_down.T
	return K_matrix_ii

def func_freeM1_K_matrix_ji(E,I,A,l):
	K_matrix_row = np.array([0,1,2])
	K_matrix_column = np.array([0,1,1])
	K_matrix_value = np.array([-E*A/l, -3*E*I/l**3, 3*E*I/l**2])
	K_matrix_ji = sparse.coo_matrix((K_matrix_value, (K_matrix_row,K_matrix_column)),shape=(3,3))
	return K_matrix_ji

def func_freeM1_K_matrix_jj(E,I,A,l):
	K_matrix_row = np.array([0,1,1,2,2])
	K_matrix_column = np.array([0,1,2,1,2])
	K_matrix_value = np.array([E*A/l, 3*E*I/l**3,-3*E*I/l**2,-3*E*I/l**2,3*E*I/l])
	K_matrix_jj = sparse.coo_matrix((K_matrix_value, (K_matrix_row,K_matrix_column)),shape=(3,3))
	return K_matrix_jj
	
def func_freeM2_K_matrix_ii(E,I,A,l):
	K_matrix_row = np.array([0,1,2,2])
	K_matrix_column = np.array([0,1,1,2])
	K_matrix_value = np.array([E*A/l, 3*E*I/l**3, 3*E*I/l**2, 3*E*I/l])
	K_matrix_down = sparse.coo_matrix((K_matrix_value, (K_matrix_row,K_matrix_column)),shape=(3,3))
	K_matrix_ii = sparse.tril(K_matrix_down,k=-1) + K_matrix_down.T
	return K_matrix_ii

def func_freeM2_K_matrix_ji(E,I,A,l):
	K_matrix_row = np.array([0,1,1])
	K_matrix_column = np.array([0,1,2])
	K_matrix_value = np.array([-E*A/l, -3*E*I/l**3, -3*E*I/l**2])
	K_matrix_ji = sparse.coo_matrix((K_matrix_value, (K_matrix_row,K_matrix_column)),shape=(3,3))
	return K_matrix_ji

def func_freeM2_K_matrix_jj(E,I,A,l):
	K_matrix_row = np.array([0,1])
	K_matrix_column = np.array([0,1])
	K_matrix_value = np.array([E*A/l, 3*E*I/l**3])
	K_matrix_jj = sparse.coo_matrix((K_matrix_value, (K_matrix_row,K_matrix_column)),shape=(3,3))
	return K_matrix_jj

def func_T_matrix_block(theta):
	T_matrix_block_Theta = sparse.coo_matrix([[np.cos(theta), -np.sin(theta)],[np.sin(theta),np.cos(theta)]])
	T_matrix_block_Eye = sparse.eye(1)
	T_matrix_block = sparse.block_diag((T_matrix_block_Theta, T_matrix_block_Eye))
	return T_matrix_block

def func_rank_element_nodes(para_array):
	for i in range(len(para_array)):
		if para_array[i,0]<para_array[i,1]:
			pass
		else:
			max_value = para_array[i,0]
			para_array[i,0]=para_array[i,1]
			para_array[i,1]=max_value
	return para_array

if __name__ == "__main__":

	# 截面惯性矩、截面面积、截面弹性模量
	sec_I = 8.33e-6
	sec_A = 1e-2
	mat_E = 210e9
	EA = mat_E*sec_A
	EI = mat_E*sec_I

	# 问题描述：二维、三维、自由度、单元类型
	diamension = '2D'  # diamension='2D'表示二维问题
	node_degree = 3
	eleType = 2  # eleType=2表示两节点单元

	# 节点坐标
	global_coo_node = np.array([[0,0],[0,4],[4,0],[4,5],[4,10],[9,10],[9,5],[9,0]])

	# 单元编号——从1开始
	element_nodes = np.array([[1,2],[2,4],[4,5],[3,4],[5,6],[4,7],[6,7],[7,8]])

	# 已知的节点荷载条件（若为均布荷载，需转化为节点对应不同自由度的直接荷载）
	# 计算等效节点力,若为跨中某位置的荷载，按两端固结的梁计算支座反力（为了统一考虑铰接情况）
	force_known_node_ID = np.array([1,5])
	force_known_global = np.array([[0,0,0],[1e4,0, 0]])
	
	# 已知的位移边界条件（固定边界）
	delta_known_node_ID = np.array([1,3,8])
	delta_known_global = np.array([[0,0,0],[0,0,0],[0,0,0]])
	
	# 静力凝聚：单元的节点自由度释放（调整）
	x_pin, y_pin, xy_pin = 0, 1, 2

	# element_node_BC的参数分别为([单元号，节点号，铰接/刚接])，2表示释放转动约束：铰接
	element_node_BC = np.array([[2,4,xy_pin]])

	# 节点数量、单元数量，与整体刚度矩阵相关
	nodeNum = len(global_coo_node)  # 单个节点自由度，（未知量个数）,len()返回第0轴的数量-行数
	eleNum = len(element_nodes)

	# 每个单元的节点号从小到大排序，方便形成刚度矩阵
	element_nodes_rank = func_rank_element_nodes(element_nodes)

	# 新建numpy array用于存储单元长度，仅考虑局部坐标与整体坐标间的旋转角度
	# 注意：这里并没有考虑整体坐标系与局部坐标系之间的平移关系
	ele_length = np.zeros(eleNum)
	ele_angle = np.zeros(eleNum)

	# 根据单元两端节点坐标求解单元在局部坐标系中的倾斜角度与单元长度
	for i_angle in range(eleNum):
		sub_two_nodes = global_coo_node[element_nodes_rank[i_angle,1]-1,:] - global_coo_node[element_nodes_rank[i_angle,0]-1,:]
		ele_length[i_angle] = np.sqrt(sub_two_nodes.dot(sub_two_nodes))
		# 对于求解arctan()时，分母为零的情况讨论
		if sub_two_nodes[0] == 0:
			if sub_two_nodes[1] > 0:
				ele_angle[i_angle] = np.pi/2
			else:
				ele_angle[i_angle] = -np.pi/2
		else:
			# np.arctan求解角度的结果，实数部分取值范围为[-pi/2, pi/2]
			ele_angle[i_angle] = np.arctan(sub_two_nodes[1]/sub_two_nodes[0])

	# 形成整体刚度稀疏矩阵
	# 列表用于存储整体刚度矩阵的行信息（测试结果表明，列表list形式的append函数效率高于numpy append函数100倍以上）
	row_list = []
	# 列表用于存储整体刚度矩阵的列信息（测试结果表明，列表list形式的append函数效率高于numpy append函数100倍以上）
	col_list = []
	# 列表用于存储整体刚度矩阵对应行列的值信息（测试结果表明，列表list形式的append函数效率高于numpy append函数100倍以上）
	value_list = []

	# 稀疏刚度矩阵（实对称阵）只存储了整体刚度矩阵的分块矩阵下三角部分，减少了形成刚度矩阵1/4的计算量
	# 将全局坐标系单元刚度矩阵按节点分成子块3×3，（这里要非常注意子块索引从单元的开始节点编号指向结尾单元编号，与角度对应）
	# 以下循环中仅仅对于ii，jj，ji中默认i>j，使得指定为形成整体刚度矩阵的分块矩阵下三角部分
	for i_K in range(eleNum):


		# 第i_K个单元的第一个节点i节点对应3个（node_freedom）自由度的单元刚度矩阵子块
		K_matrix_local_ele_b_ii = func_K_matrix_ii(mat_E, sec_I, sec_A, ele_length[i_K])
		# 第i_K个单元的第二个节点j到第一个节点i对应3个（node_freedom）自由度的单元刚度矩阵子块
		K_matrix_local_ele_b_ji = func_K_matrix_ji(mat_E, sec_I, sec_A, ele_length[i_K])
		# 第i_K个单元的第二个节点j节点对应3个（node_freedom）自由度的单元刚度矩阵子块
		K_matrix_local_ele_b_jj = func_K_matrix_jj(mat_E, sec_I, sec_A, ele_length[i_K])

		#考虑节点铰接的情况——静力凝聚法
		for j_nodeBC in range(len(element_node_BC)):
			index_row_element_node_BC = element_node_BC[j_nodeBC,0] - 1
			if index_row_element_node_BC == i_K:
				# K_cc_wang为静力凝聚法（王勖成）需释放自由度对应位置的刚度矩阵子块
				K_cc_wang = (4*mat_E * sec_I)/ele_length[index_row_element_node_BC]
				if element_node_BC[j_nodeBC,1] == element_nodes_rank[i_K,0]:
					K_matrix_local_ele_b_ii=func_freeM1_K_matrix_ii(mat_E, sec_I, sec_A, ele_length[i_K])
					K_matrix_local_ele_b_ji=func_freeM1_K_matrix_ji(mat_E, sec_I, sec_A, ele_length[i_K])
					K_matrix_local_ele_b_jj=func_freeM1_K_matrix_jj(mat_E, sec_I, sec_A, ele_length[i_K])
					
					K0_star_11 = K_matrix_local_ele_b_ii.tocsr()[np.arange(0,2),:].tocsc()[:,np.arange(0,2)]
					K0_star_21 = K_matrix_local_ele_b_ji.tocsr()[np.arange(0,3),:].tocsc()[:,np.arange(0,2)]
					K0_star_12 = K0_star_21.T
					K0_star_22 = K_matrix_local_ele_b_jj
					
					K0_star = sparse.bmat([[K0_star_11,K0_star_12],[K0_star_21,K0_star_22]])
					K0_star_diag = K0_star.diagonal()
					BigNumber_wang = 36854775807
					K0_star_diag[2:5] = BigNumber_wang*K0_star_diag[2:5]
					K0_star.setdiag(K0_star_diag)

					K_c0_wang = matrix([0,(6*mat_E * sec_I)/(ele_length[index_row_element_node_BC]**2)])
					K_c1_wang = matrix([0,-(6*mat_E * sec_I)/(ele_length[index_row_element_node_BC]**2),(2*mat_E * sec_I)/ele_length[index_row_element_node_BC]])

					P0_star = matrix([0,0,0,0,0])
					a0_a1_wang = spsolve(K0_star.tocsc(),P0_star.T)
					a0_wang = a0_a1_wang[0:2]
					a1_wang = a0_a1_wang[2:5]
					Pc_wang = 0

					ac_wang = 1/K_cc_wang * (Pc_wang-K_c0_wang.dot(a0_wang)-K_c1_wang.dot(a1_wang))
					print('ac_wang_nodei=',ac_wang)

				elif element_node_BC[j_nodeBC,1] == element_nodes_rank[i_K,1]:
					K_matrix_local_ele_b_ii=func_freeM2_K_matrix_ii(mat_E, sec_I, sec_A, ele_length[i_K])
					K_matrix_local_ele_b_ji=func_freeM2_K_matrix_ji(mat_E, sec_I, sec_A, ele_length[i_K])
					K_matrix_local_ele_b_jj=func_freeM2_K_matrix_jj(mat_E, sec_I, sec_A, ele_length[i_K])

					K0_star_11 = K_matrix_local_ele_b_ii.tocsr()[np.arange(0,3),:].tocsc()[:,np.arange(0,3)]
					K0_star_21 = K_matrix_local_ele_b_ji.tocsr()[np.arange(0,2),:].tocsc()[:,np.arange(0,3)]
					K0_star_12 = K0_star_21.T
					K0_star_22 = K_matrix_local_ele_b_jj.tocsr()[np.arange(0,2),:].tocsc()[:,np.arange(0,2)]
					
					K0_star = sparse.bmat([[K0_star_11,K0_star_12],[K0_star_21,K0_star_22]])
					K0_star_diag = K0_star.diagonal()
					BigNumber_wang = 36854775807
					K0_star_diag[0:3] = BigNumber_wang*K0_star_diag[0:3]
					K0_star.setdiag(K0_star_diag)
					
					K_c0_wang = np.array([0, (6*mat_E * sec_I)/(ele_length[index_row_element_node_BC]**2), (2*mat_E * sec_I)/ele_length[index_row_element_node_BC], 0, -(6*mat_E * sec_I)/(ele_length[index_row_element_node_BC]**2)])
					P0_star = matrix([0,0,0,0,0])
					a0_wang = spsolve(K0_star.tocsc(),P0_star.T)
					Pc_wang = 0

					ac_wang = 1/K_cc_wang * (Pc_wang-K_c0_wang.dot(a0_wang))
					print('ac_wang_nodej=',ac_wang)

				else:
					pass
			else:
				pass
		# print('K_matrix_local_ele_b_ii = ',K_matrix_local_ele_b_ii)
		# 坐标旋转矩阵一个节点所有自由度对应的子矩阵
		T_matrix_block_ele = func_T_matrix_block(ele_angle[i_K])
		T_matrix_block_ele_T = T_matrix_block_ele.T
		np.set_printoptions(precision=6, suppress=True)
		print('T_matrix_block_ele=',T_matrix_block_ele.toarray())

		# 单元刚度矩阵各分块经变换转换为整体坐标系下的总体刚度矩阵子块值
		K_matrix_global_ele_b_ii = (T_matrix_block_ele.dot(K_matrix_local_ele_b_ii)).dot(T_matrix_block_ele_T)
		K_matrix_global_ele_b_ji = (T_matrix_block_ele.dot(K_matrix_local_ele_b_ji)).dot(T_matrix_block_ele_T)
		K_matrix_global_ele_b_jj = (T_matrix_block_ele.dot(K_matrix_local_ele_b_jj)).dot(T_matrix_block_ele_T)
		#print('K_matrix_global_ele_b_ii=',K_matrix_global_ele_b_ii.toarray())
		#print('K_matrix_global_ele_b_ji=',K_matrix_global_ele_b_ji.toarray())
		#print('K_matrix_global_ele_b_jj=',K_matrix_global_ele_b_jj.toarray())
		# 定位整体坐标系下的总体刚度矩阵子块对应的节点号，默认i节点号小于j节点号
		node_ID_i = element_nodes_rank[i_K,0]
		node_ID_j = element_nodes_rank[i_K,1]

		# 利用各单元两端的节点号信息找到该节点自由度对应的刚度矩阵子块在整体刚度矩阵中的行列索引
		# 循环检索不同单元的i节点到i节点刚度矩阵子块信息，存储行列信息
		# 考虑到python从0开始索引，而节点编号从1开始，故“-1”
		index_row_ii = node_degree*(node_ID_i - 1)
		# 考虑到python从0开始索引，而节点编号从1开始，故“-1”
		index_col_ii = node_degree*(node_ID_i - 1)
		# 整体刚度矩阵中节点自由度ii对应单元刚度矩阵的行索引
		row_list_ii = [index_row_ii,index_row_ii,index_row_ii,index_row_ii+1,index_row_ii+1,index_row_ii+1,index_row_ii+2,index_row_ii+2,index_row_ii+2]
		# 整体刚度矩阵中节点自由度ii对应单元刚度矩阵的列索引
		col_list_ii = [index_col_ii,index_col_ii+1,index_col_ii+2,index_col_ii,index_col_ii+1,index_col_ii+2,index_col_ii,index_col_ii+1,index_col_ii+2]
		# 循环检索不同单元的ii节点自由度刚度矩阵子块信息，存储行列信息
		for j_rc_ii in range(node_degree**2):
			row_list.append(row_list_ii[j_rc_ii])
			col_list.append(col_list_ii[j_rc_ii])
		for j_va_ii in range(node_degree):
			for k_va_ii in range(node_degree):
				value_list.append(K_matrix_global_ele_b_ii.toarray()[j_va_ii, k_va_ii])
		
		# 循环检索不同单元的j节点到i节点刚度矩阵子块信息，存储行列信息
		# 考虑到python从0开始索引，而节点编号从1开始，故“-1”
		index_row_ji = node_degree*(node_ID_j - 1)
		# 考虑到python从0开始索引，而节点编号从1开始，故“-1”
		index_col_ji = node_degree*(node_ID_i - 1)
		# 整体刚度矩阵中节点自由度ji对应单元刚度矩阵的行索引
		row_list_ji = [index_row_ji,index_row_ji,index_row_ji,index_row_ji+1,index_row_ji+1,index_row_ji+1,index_row_ji+2,index_row_ji+2,index_row_ji+2]
		# 整体刚度矩阵中节点自由度ji对应单元刚度矩阵的列索引
		col_list_ji = [index_col_ji,index_col_ji+1,index_col_ji+2,index_col_ji,index_col_ji+1,index_col_ji+2,index_col_ji,index_col_ji+1,index_col_ji+2]
		# 循环检索不同单元的ji节点自由度刚度矩阵子块信息，存储行列信息
		for j_rc_ji in range(node_degree**2):
			row_list.append(row_list_ji[j_rc_ji])
			col_list.append(col_list_ji[j_rc_ji])
		for j_va_ji in range(node_degree):
			for k_va_ji in range(node_degree):
				value_list.append(K_matrix_global_ele_b_ji.toarray()[j_va_ji, k_va_ji])

		# 循环检索不同单元的j节点到j节点刚度矩阵子块信息，存储行列信息
		# 考虑到python从0开始索引，而节点编号从1开始，故“-1”
		index_row_jj = node_degree*(node_ID_j - 1)
		# 考虑到python从0开始索引，而节点编号从1开始，故“-1”
		index_col_jj = node_degree*(node_ID_j - 1)
		# 整体刚度矩阵中节点自由度jj对应单元刚度矩阵的行索引
		row_list_jj = [index_row_jj,index_row_jj,index_row_jj,index_row_jj+1,index_row_jj+1,index_row_jj+1,index_row_jj+2,index_row_jj+2,index_row_jj+2]
		# 整体刚度矩阵中节点自由度jj对应单元刚度矩阵的列索引
		col_list_jj = [index_col_jj,index_col_jj+1,index_col_jj+2,index_col_jj,index_col_jj+1,index_col_jj+2,index_col_jj,index_col_jj+1,index_col_jj+2]
		# 循环检索不同单元的jj节点自由度刚度矩阵子块信息，存储行列信息	
		for j_rc_jj in range(node_degree**2):
			row_list.append(row_list_jj[j_rc_jj])
			col_list.append(col_list_jj[j_rc_jj])
		for j_va_jj in range(node_degree):
			for k_va_jj in range(node_degree):
				value_list.append(K_matrix_global_ele_b_jj.toarray()[j_va_jj, k_va_jj])

	K_matrix_global_all = sparse.coo_matrix((value_list, (row_list,col_list)), shape = (node_degree*nodeNum, node_degree*nodeNum), dtype=np.float64)

	F_row_list = []
	F_col_list = []
	F_val_list = []

	# 根据已知的节点荷载形成稀疏矩阵
	for i_rowF in range(len(force_known_node_ID)):
		F_rowID = node_degree*(force_known_node_ID[i_rowF] - 1)
		for j_Fr in range(node_degree):
			F_row_list.append(F_rowID+j_Fr)
			# print('F_row_list=',F_row_list)
			F_col_list.append(0)
			F_val_list.append(force_known_global[i_rowF,j_Fr])
	
	# '''
	# ----------------------------------------------------------------------------------------------------------------
	# 位移边界条件处理方法——对角元素乘大数法BigNumber = 36854775807
	# import sys
	# max = sys.maxsize
	# print (max)
	# max = 9223372036854775807
	BigNumber = 36854775807
	K_matrix_global_all_dia = K_matrix_global_all.diagonal()
	# print('K_matrix_global_all_dia',K_matrix_global_all_dia)
	for i_dia in range(len(delta_known_node_ID)):
		diaID = node_degree*(delta_known_node_ID[i_dia] - 1)
		K_matrix_global_all_dia[diaID:diaID+node_degree] = BigNumber*K_matrix_global_all_dia[diaID:diaID+node_degree]

	# 根据位移边界条件对应修改节点荷载稀疏矩阵-采用乘大数法
	for i_rowFB in range(len(delta_known_node_ID)):
		FB_rowID = node_degree*(delta_known_node_ID[i_rowFB] - 1)
		for j_FBr in range(node_degree):
			F_row_list.append(FB_rowID+j_FBr)
			# print('F_row_list=',F_row_list)
			F_col_list.append(0)
			# print('F_col_list',F_col_list)
			F_val_list.append(BigNumber*(K_matrix_global_all_dia[FB_rowID+j_FBr])*(delta_known_global[i_rowFB,j_FBr]))
			# print('F_val_list',F_val_list)
	# 求解刚度矩阵
	# 对于已释放节点自由度对应的刚度矩阵与荷载向量，采用“加”大数法进行处理，这样可避免刚度矩阵奇异的问题
	# K_matrix_global_all_dia[]
	K_matrix_global_all.setdiag(K_matrix_global_all_dia)  # 整体刚度矩阵的下三角分块矩阵部分
	# ----------------------------------------------------------------------------------------------------------------
	#'''
	trial_K_matrix = sparse.tril(K_matrix_global_all)  # 整体刚度矩阵的下三角矩阵部分
	K_matrix_global_all_final = sparse.tril(K_matrix_global_all) + (sparse.tril(K_matrix_global_all,k=-1)).T  # 还原整体刚度对称矩阵

	K_matrix_global_all_11 = K_matrix_global_all_final.tocsr()[np.arange(0,8),:].tocsc()[:,np.arange(0,8)]
	K_matrix_global_all_12 = K_matrix_global_all_final.tocsr()[np.arange(0,8),:].tocsc()[:,np.arange(9,24)]
	K_matrix_global_all_21 = K_matrix_global_all_final.tocsr()[np.arange(9,24),:].tocsc()[:,np.arange(0,8)]
	K_matrix_global_all_22 = K_matrix_global_all_final.tocsr()[np.arange(9,24),:].tocsc()[:,np.arange(9,24)]
	
	K_matrix_global_all_final_reduced = sparse.bmat([[K_matrix_global_all_11,K_matrix_global_all_12],[K_matrix_global_all_21, K_matrix_global_all_22]])

	np.set_printoptions(precision=6, suppress=True)
	print('K_matrix_global_all_final=',K_matrix_global_all_final_reduced.toarray())
	# 通过整体刚度矩阵子块形成完整的整体刚度矩阵(已经考虑边界约束)
	F_matrix_global_all = sparse.coo_matrix((F_val_list,(F_row_list,F_col_list)),dtype=np.float64)
	F_matrix_global_all_reduced = np.delete(F_matrix_global_all.toarray(),5)
	print('F_matrix_global_all',F_matrix_global_all_reduced)
	delta_array_node_all = spsolve(K_matrix_global_all_final_reduced.tocsc(), F_matrix_global_all_reduced)
	# '''
	print('delta_array_node_all=', delta_array_node_all)
	# print('delta_array_node_all1=', 1/K_cc_wang * Pc_wang)

'''
结点,1,0,0
结点,2,0,4
结点,3,4,0
结点,4,4,5
结点,5,4,10
结点,6,9,10
结点,7,9,5
结点,8,9,0
单元,1,2,1,1,1,1,1,1
单元,2,4,1,1,1,1,1,0
单元,4,5,1,1,1,1,1,1
单元,3,4,1,1,1,1,1,1
单元,5,6,1,1,1,1,1,1
单元,6,7,1,1,1,1,1,1
单元,7,8,1,1,1,1,1,1
单元,4,7,1,1,1,1,1,1
结点支承,1,6,0,0,0,0
结点支承,3,6,0,0,0,0
结点支承,8,6,0,0,0,0
结点荷载,5,1,10000,0
单元材料性质,1,8,210e7,1749.3e3,0,0,-1
'''
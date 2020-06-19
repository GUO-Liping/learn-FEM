'''
由于FEM-2Dbar-v0.0.1.py中的单元刚度矩阵难以集成，
故在FEM-2Dbar-v0.0.2.py版本中将局部坐标系下的K矩阵以节点为单位进行分块，
尝试进一步高效集成为稀疏矩阵，但无法获得较高的集成效率
故在FEM-2Dbar-v0.0.2.py版本中直接从最小的节点自由度出发对单元刚度矩阵分块，
再进一步尝试高效集成稀疏矩阵...
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
from scipy import sparse
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
	sec_I = 1e-2
	sec_A = 1e-2
	mat_E = 1e6
	EA = mat_E*sec_A
	EI = mat_E*sec_I

	# 问题描述：二维、三维、自由度、单元类型
	diamension = '2D'  # diamension='2D'表示二维问题
	node_degree = 3
	eleType = 2  # eleType=2表示两节点单元

	# 节点坐标
	global_coo_node = np.array([[0,0],[1,0],[2,0]])

	# 单元编号
	element_nodes = np.array([[1,2],[2,3]])

	# 已知的节点荷载条件（若为均布荷载，需转化为节点对应不同自由度的直接荷载）
	# 计算等效节点力
	force_known_node_ID = np.array([1,2])
	force_known_global = np.array([[0,0,0],[0,1e4, 0]])
	
	# 已知的位移边界条件（固定边界）
	delta_known_node_ID = np.array([1,3])
	delta_known_global = np.array([[0,0,0],[0,0,0]])
	
	# 静力凝聚：单元的节点自由度释放（调整）
	x_pin, y_pin, xy_pin = 1, 2, 3
	# element_node_BC的参数分别为([单元号，节点i自由度，节点j自由度])，默认i>j，0表示铰接，1表示刚接
	# element_node_BC = np.array([[1,fixed,pinned],[2,pinned,fixed]])
	# element_node_BC的参数分别为([单元号，节点号，铰接/刚接])，0表示铰接
	element_node_BC = np.array([[2,3,xy_pin],[3,3,xy_pin]])

	# 节点数量、单元数量，与整体刚度矩阵相关
	nodeNum = len(global_coo_node)  # 单个节点自由度，（未知量个数）,len()返回第0轴的数量-行数
	eleNum = len(element_nodes)

	# 每个单元的节点号从小到大排序，方便形成刚度矩阵
	element_nodes_rank = func_rank_element_nodes(element_nodes)

	# 新建numpy array用于存储单元长度，仅考虑局部坐标与整体坐标间的旋转角度
	# 注意：这里并没有考虑整体坐标系与局部坐标系之间的平移关系
	ele_length = np.empty(eleNum)
	ele_angle = np.empty(eleNum)

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
			#考虑节点铰接的情况——静力凝聚法
		for j_nodeBC in range(len(element_node_BC)):
			if i_K == element_node_BC[j_nodeBC,0]:
				if element_nodes_rank[i_K,0] == element_node_BC[j_nodeBC,1]:
					K_matrix_local_ele_b_ii=func_freeM1_K_matrix_ii(mat_E, sec_I, sec_A, ele_length[i_K])
					K_matrix_local_ele_b_ii=func_freeM1_K_matrix_ji(mat_E, sec_I, sec_A, ele_length[i_K])
					K_matrix_local_ele_b_ji=func_freeM1_K_matrix_jj(mat_E, sec_I, sec_A, ele_length[i_K])
					
				elif elment_nodes_rank[i_K,1] == element_node_BC[j_nodeBC,1]:
					K_matrix_local_ele_b_ii=func_freeM2_K_matrix_ii(mat_E, sec_I, sec_A, ele_length[i_K])
					K_matrix_local_ele_b_ji=func_freeM2_K_matrix_ji(mat_E, sec_I, sec_A, ele_length[i_K])
					K_matrix_local_ele_b_jj=func_freeM2_K_matrix_jj(mat_E, sec_I, sec_A, ele_length[i_K])
				else:
					pass
			else:
				pass

		# 第i_K个单元的第一个节点i节点对应3个（node_freedom）自由度的单元刚度矩阵子块
		K_matrix_local_ele_b_ii = func_K_matrix_ii(mat_E, sec_I, sec_A, ele_length[i_K])
		# 第i_K个单元的第二个节点j到第一个节点i对应3个（node_freedom）自由度的单元刚度矩阵子块
		K_matrix_local_ele_b_ji = func_K_matrix_ji(mat_E, sec_I, sec_A, ele_length[i_K])
		# 第i_K个单元的第二个节点j节点对应3个（node_freedom）自由度的单元刚度矩阵子块
		K_matrix_local_ele_b_jj = func_K_matrix_jj(mat_E, sec_I, sec_A, ele_length[i_K])
		
		# 坐标旋转矩阵一个节点所有自由度对应的子矩阵
		T_matrix_block_ele = func_T_matrix_block(ele_angle[i_K])
		T_matrix_block_ele_T = T_matrix_block_ele.T

		# 单元刚度矩阵各分块经变换转换为整体坐标系下的总体刚度矩阵子块值
		K_matrix_global_ele_b_ii = (T_matrix_block_ele.dot(K_matrix_local_ele_b_ii)).dot(T_matrix_block_ele_T)
		K_matrix_global_ele_b_ji = (T_matrix_block_ele.dot(K_matrix_local_ele_b_ji)).dot(T_matrix_block_ele_T)
		K_matrix_global_ele_b_jj = (T_matrix_block_ele.dot(K_matrix_local_ele_b_jj)).dot(T_matrix_block_ele_T)

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
	# ----------------------------------------------------------------------------------------------------------------

	# 求解刚度矩阵
	K_matrix_global_all.setdiag(K_matrix_global_all_dia)  # 整体刚度矩阵的下三角分块矩阵部分
	trial_K_matrix = sparse.tril(K_matrix_global_all)  # 整体刚度矩阵的下三角矩阵部分
	K_matrix_global_all_final = sparse.tril(K_matrix_global_all) + (sparse.tril(K_matrix_global_all,k=-1)).T  # 还原整体刚度对称矩阵
	#print('K_matrix_global_all_final=',K_matrix_global_all_final)
	# 通过整体刚度矩阵子块形成完整的整体刚度矩阵(已经考虑边界约束)
	F_matrix_global_all = sparse.coo_matrix((F_val_list,(F_row_list,F_col_list)),dtype=np.float64)
	#print('F_matrix_global_all',F_matrix_global_all)
	delta_array_node_all = spsolve(K_matrix_global_all_final.tocsc(), F_matrix_global_all)
	# '''
	print('delta_array_node_all=', delta_array_node_all)
	
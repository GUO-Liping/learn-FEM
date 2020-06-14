'''
由于FEM-2Dbar-v0.0.1.py中的单元刚度矩阵难以集成，
故在FEM-2Dbar-v0.0.2.py版本中将局部坐标系下的K矩阵以节点为单位进行分块，
尝试进一步高效集成为稀疏矩阵，但无法获得较高的集成效率
故在FEM-2Dbar-v0.0.2.py版本中直接从最小的节点自由度出发对单元刚度矩阵分块，
再进一步尝试高效集成稀疏矩阵...
'''
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve

def calcu_I(para, width=0, height=0, d=1):
	if para == 'rectangle':
		return width*height**3/12
	if para == 'triangle':
		return width*height**3/36
	if para == 'circle':
		return np.pi*d**4/64
	if para =='ring':
		return np.pi*(D**4-d**4)/64

def func_get_node(para_dia, para_num):
	if para_dia == '2D':
		nodeNum = para_num
		global_coo_node = np.empty((nodeNum, 2))
		for i_co2 in range(nodeNum):
			print('请输入节点',i_co2+1,"的x坐标")
			x_coo = float(input())
			print('请输入节点',i_co2+1,"的y坐标")
			y_coo = float(input())
			global_coo_node[i_co2,:] = np.array([x_coo, y_coo])
	elif para_dia == '3D':
		nodeNum = para_num
		global_coo_node = np.empty((nodeNum, 3))
		for i_co3 in range(nodeNum):
			print('请输入节点',i_co3+1,"的x坐标")
			x_coo = float(input())
			print('请输入节点',i_co3+1,"的y坐标")
			y_coo = float(input())
			print('请输入节点',i_co3+1,"的z坐标")
			z_coo = float(input())
			global_coo_node[i_co3,:] = np.array([x_coo, y_coo, z_coo])
	return global_coo_node

def func_get_element(para_type, para_eleNum):
	eleNum = para_eleNum
	eleType = para_type
	global_ele_num = np.empty((eleNum, eleType), dtype=int)
	for i_in in range(eleNum):
		print('请输入单元',i_in+1,"的首节点号")
		i_num = float(input())
		print('请输入单元',i_in+1,"的末节点号")
		j_num = float(input())
		global_ele_num[i_in,:] = np.array([i_num, j_num])
	return global_ele_num

def func_K_matrix(E,I,A,l):
	K_matrix_row = np.array([0,1,2,2,3,3,4,4,4,5,5,5,5])
	K_matrix_column = np.array([0,1,1,2,0,3,1,2,4,1,2,4,5])
	K_matrix_value = np.array([E*A/l, 12*E*I/l**3, 6*E*I/l**2, 4*E*I/l, -E*A/l, E*A/l, -12*E*I/l**3, -6*E*I/l**2, 12*E*I/l**3, 6*E*I/l**2, 2*E*I/l, -6*E*I/l**2, 4*E*I/l])
	K_matrix_down = sparse.coo_matrix((K_matrix_value, (K_matrix_row,K_matrix_column)))
	dia_data = K_matrix_down.diagonal()
	K_matrix_dia = sparse.dia_matrix((dia_data, np.array([0])),shape=(6,6))
	K_matrix = K_matrix_down + K_matrix_down.T - K_matrix_dia
	return K_matrix

def func_T_matrix(theta):
	block_matrix = np.array([[np.cos(theta), -np.sin(theta), 0],[np.sin(theta),np.cos(theta),0],[0,0,1]])
	index_row = np.array([0,1])
	count_row = np.array([0,1,2])
	data = np.array([block_matrix,block_matrix])
	T_matrix = sparse.bsr_matrix((data, index_row, count_row), shape=(6,6))
	return T_matrix

def func_K_matrix_ii(E,I,A,l):
	K_matrix_row = np.array([0,1,2,2])
	K_matrix_column = np.array([0,1,1,2])
	K_matrix_value = np.array([E*A/l, 12*E*I/l**3, 6*E*I/l**2, 4*E*I/l])
	K_matrix_down = sparse.coo_matrix((K_matrix_value, (K_matrix_row,K_matrix_column)))
	dia_data = K_matrix_down.diagonal()  # 提取对角线元素
	K_matrix_dia = sparse.dia_matrix((dia_data, np.array([0])),shape=(3,3))
	K_matrix_ii = K_matrix_down + K_matrix_down.T - K_matrix_dia
	return K_matrix_ii

def func_K_matrix_ji(E,I,A,l):
	K_matrix_row = np.array([0,1,1,2,2])
	K_matrix_column = np.array([0,1,2,1,2])
	K_matrix_value = np.array([-E*A/l, -12*E*I/l**3, -6*E*I/l**2, 6*E*I/l**2, 2*E*I/l])
	K_matrix_ji = sparse.coo_matrix((K_matrix_value, (K_matrix_row,K_matrix_column)))
	return K_matrix_ji

def func_K_matrix_jj(E,I,A,l):
	K_matrix_row = np.array([0,1,2,2])
	K_matrix_column = np.array([0,1,1,2])
	K_matrix_value = np.array([E*A/l, 12*E*I/l**3, -6*E*I/l**2, 4*E*I/l])
	K_matrix_down = sparse.coo_matrix((K_matrix_value, (K_matrix_row,K_matrix_column)))
	dia_data = K_matrix_down.diagonal()  # 提取对角线元素
	K_matrix_dia = sparse.dia_matrix((dia_data, np.array([0])),shape=(3,3))  # 从对角线元素建立对角阵
	K_matrix_jj = K_matrix_down + K_matrix_down.T - K_matrix_dia
	return K_matrix_jj

def func_T_matrix_block(theta):
	T_matrix_block_Theta = sparse.coo_matrix([[np.cos(theta), -np.sin(theta)],[np.sin(theta),np.cos(theta)]])
	T_matrix_block_Eye = sparse.eye(1)
	T_matrix_block = sparse.block_diag((T_matrix_block_Theta, T_matrix_block_Eye))
	return T_matrix_block

if __name__ == "__main__":

	# 截面惯性矩、截面面积、截面弹性模量
	sec_I = 15760e-8
	sec_A = 76.3e-4
	mat_E = 2e11

	# 节点编号
	diamension, nodeNum = '2D', 4	# diamension='2D'表示二维问题
	# global_coo_node = func_get_node(diamension, nodeNum)
	global_coo_node = np.array([[0,5],[6.4,5],[0,0],[9.6,0]])
	nodeNum, node_degree = 4, 3  # 单个节点自由度，（未知量个数）

	# 单元编号
	eleType, eleNum = 2, 3  # eleType=2表示两节点单元
	# element_nodes = func_get_element(eleType, eleNum)
	element_nodes = np.array([[1,2],[3,1],[2,4]])
	ele_length = np.empty(eleNum)
	ele_angle = np.empty(eleNum)

	# 根据单元两端节点坐标求解单元在局部坐标系中的倾斜角度与单元长度
	for i_angle in range(eleNum):
		sub_two_nodes = global_coo_node[element_nodes[i_angle,1]-1,:] - global_coo_node[element_nodes[i_angle,0]-1,:]
		ele_length[i_angle] = np.sqrt(sub_two_nodes.dot(sub_two_nodes))
		if sub_two_nodes[0] == 0:
			if sub_two_nodes[1] > 0:
				ele_angle[i_angle] = np.pi/2
			else:
				ele_angle[i_angle] = -np.pi/2
		else:
			# np.arctan其实数部分的取值范围为[-pi/2, pi/2]
			ele_angle[i_angle] = np.arctan(sub_two_nodes[1]/sub_two_nodes[0])

	print('ele_length=\n', ele_length)

	# 单元长度
	ele_l1 = 6.4
	ele_l2 = 5.0
	ele_l3 = 5.936

	print('ele_angle=\n', ele_angle)

	# 单元局部坐标系相当于整体坐标系中的旋转角度
	theta1 = 0
	theta2 = np.pi/2
	theta3 = -np.arctan(5/3.2)

	# 局部坐标系单元刚度矩阵6×6
	K_matrix_local_ele = func_K_matrix(mat_E, sec_I, sec_A, ele_length[0])
	T_matrix_theta_ele = func_T_matrix(ele_angle[0])
	K_matrix_global_ele = (T_matrix_theta_ele.dot(K_matrix_local_ele)).dot(T_matrix_theta_ele.T)

	#'''
	# 形成整体刚度稀疏矩阵
	#'''
	row_list = []
	col_list = []
	value_list = []
	for i_K in range(eleNum):
		K_matrix_local_ele_b_ii = func_K_matrix_ii(mat_E, sec_I, sec_A, ele_length[i_K])
		K_matrix_local_ele_b_ji = func_K_matrix_ji(mat_E, sec_I, sec_A, ele_length[i_K])
		K_matrix_local_ele_b_jj = func_K_matrix_jj(mat_E, sec_I, sec_A, ele_length[i_K])
		
		T_matrix_block_ele = func_T_matrix_block(ele_angle[i_K])

		T_matrix_block_ele_T = T_matrix_block_ele.T
		K_matrix_global_ele_b_ii = (T_matrix_block_ele.dot(K_matrix_local_ele_b_ii)).dot(T_matrix_block_ele_T)
		K_matrix_global_ele_b_ji = (T_matrix_block_ele.dot(K_matrix_local_ele_b_ji)).dot(T_matrix_block_ele_T)
		K_matrix_global_ele_b_jj = (T_matrix_block_ele.dot(K_matrix_local_ele_b_jj)).dot(T_matrix_block_ele_T)

		node_ID_i = min(element_nodes[i_K,0],element_nodes[i_K,1])
		node_ID_j = max(element_nodes[i_K,0],element_nodes[i_K,1])

		index_row_ii = node_degree*(node_ID_i - 1)  # 考虑到python从0开始索引，而节点编号从1开始，故“-1”
		index_col_ii = node_degree*(node_ID_i - 1)  # 考虑到python从0开始索引，而节点编号从1开始，故“-1”
		row_list_ii = [index_row_ii,index_row_ii,index_row_ii,index_row_ii+1,index_row_ii+1,index_row_ii+1,index_row_ii+2,index_row_ii+2,index_row_ii+2]
		col_list_ii = [index_col_ii,index_col_ii+1,index_col_ii+2,index_col_ii,index_col_ii+1,index_col_ii+2,index_col_ii,index_col_ii+1,index_col_ii+2]
		
		for j_rc_ii in range(node_degree**2):
			row_list.append(row_list_ii[j_rc_ii])
			col_list.append(col_list_ii[j_rc_ii])
		for j_va_ii in range(node_degree):
			for k_va_ii in range(node_degree):
				value_list.append(K_matrix_global_ele_b_ii.toarray()[j_va_ii, k_va_ii])

		index_row_ji = node_degree*(node_ID_j - 1)  # 考虑到python从0开始索引，而节点编号从1开始，故“-1”
		index_col_ji = node_degree*(node_ID_i - 1)  # 考虑到python从0开始索引，而节点编号从1开始，故“-1”
		row_list_ji = [index_row_ji,index_row_ji,index_row_ji,index_row_ji+1,index_row_ji+1,index_row_ji+1,index_row_ji+2,index_row_ji+2,index_row_ji+2]
		col_list_ji = [index_col_ji,index_col_ji+1,index_col_ji+2,index_col_ji,index_col_ji+1,index_col_ji+2,index_col_ji,index_col_ji+1,index_col_ji+2]
		
		for j_rc_ji in range(node_degree**2):
			row_list.append(row_list_ji[j_rc_ji])
			col_list.append(col_list_ji[j_rc_ji])
		for j_va_ji in range(node_degree):
			for k_va_ji in range(node_degree):
				value_list.append(K_matrix_global_ele_b_ji.toarray()[j_va_ji, k_va_ji])
		
		index_row_jj = node_degree*(node_ID_j - 1)  # 考虑到python从0开始索引，而节点编号从1开始，故“-1”
		index_col_jj = node_degree*(node_ID_j - 1)  # 考虑到python从0开始索引，而节点编号从1开始，故“-1”
		row_list_jj = [index_row_jj,index_row_jj,index_row_jj,index_row_jj+1,index_row_jj+1,index_row_jj+1,index_row_jj+2,index_row_jj+2,index_row_jj+2]
		col_list_jj = [index_col_jj,index_col_jj+1,index_col_jj+2,index_col_jj,index_col_jj+1,index_col_jj+2,index_col_jj,index_col_jj+1,index_col_jj+2]
		
		for j_rc_jj in range(node_degree**2):
			row_list.append(row_list_jj[j_rc_jj])
			col_list.append(col_list_jj[j_rc_jj])
		for j_va_jj in range(node_degree):
			for k_va_jj in range(node_degree):
				value_list.append(K_matrix_global_ele_b_jj.toarray()[j_va_jj, k_va_jj])

	print('row_list=',row_list)
	print('col_list=',col_list)
	K_matrix_global_all = sparse.coo_matrix((value_list, (row_list,col_list)), shape = (node_degree*nodeNum, node_degree*nodeNum), dtype=np.float64)
	#'''
	
	K_matrix_local_ele1 = func_K_matrix(mat_E, sec_I, sec_A, ele_l1)
	K_matrix_local_ele2 = func_K_matrix(mat_E, sec_I, sec_A, ele_l2)
	K_matrix_local_ele3 = func_K_matrix(mat_E, sec_I, sec_A, ele_l3)

	# print('calibration of local Stiffness matrix \n', K_matrix_local_ele1==K_matrix_local_ele[0])

	# 局部坐标系-整体坐标系转换矩阵6×6
	T_matrix_theta1 = func_T_matrix(theta1)
	T_matrix_theta2 = func_T_matrix(theta2)
	T_matrix_theta3 = func_T_matrix(theta3)

	# 全局坐标系单元刚度矩阵6×6
	K_matrix_global_ele1 = (T_matrix_theta1.dot(K_matrix_local_ele1)).dot(T_matrix_theta1.T)
	K_matrix_global_ele2 = (T_matrix_theta2.dot(K_matrix_local_ele2)).dot(T_matrix_theta2.T)
	K_matrix_global_ele3 = (T_matrix_theta3.dot(K_matrix_local_ele3)).dot(T_matrix_theta3.T)

	# print('calibration of global Stiffness matrix \n', K_matrix_global_ele1==K_matrix_global_ele[0])

	# 将全局坐标系单元刚度矩阵按节点分成子块3×3，（这里要非常注意子块索引从单元的开始节点编号指向结尾单元编号，与角度对应）
	K_matrix_global_ele1_11 = K_matrix_global_ele1.tocsr()[np.arange(3),:].tocsc()[:,np.arange(0,3)]
	K_matrix_global_ele1_12 = K_matrix_global_ele1.tocsr()[np.arange(3),:].tocsc()[:,np.arange(3,6)]
	K_matrix_global_ele1_21 = K_matrix_global_ele1.tocsr()[np.arange(3,6),:].tocsc()[:,np.arange(0,3)]
	K_matrix_global_ele1_22 = K_matrix_global_ele1.tocsr()[np.arange(3,6),:].tocsc()[:,np.arange(3,6)]

	K_matrix_global_ele2_33 = K_matrix_global_ele2.tocsr()[np.arange(3),:].tocsc()[:,np.arange(0,3)]
	K_matrix_global_ele2_31 = K_matrix_global_ele2.tocsr()[np.arange(3),:].tocsc()[:,np.arange(3,6)]
	K_matrix_global_ele2_13 = K_matrix_global_ele2.tocsr()[np.arange(3,6),:].tocsc()[:,np.arange(0,3)]
	K_matrix_global_ele2_11 = K_matrix_global_ele2.tocsr()[np.arange(3,6),:].tocsc()[:,np.arange(3,6)]

	K_matrix_global_ele3_22 = K_matrix_global_ele3.tocsr()[np.arange(3),:].tocsc()[:,np.arange(0,3)]
	K_matrix_global_ele3_24 = K_matrix_global_ele3.tocsr()[np.arange(3),:].tocsc()[:,np.arange(3,6)]
	K_matrix_global_ele3_42 = K_matrix_global_ele3.tocsr()[np.arange(3,6),:].tocsc()[:,np.arange(0,3)]
	K_matrix_global_ele3_44 = K_matrix_global_ele3.tocsr()[np.arange(3,6),:].tocsc()[:,np.arange(3,6)]

	# 通过整体刚度矩阵子块形成完整的整体刚度矩阵(已经考虑边界约束)
	K_11 = K_matrix_global_ele1_11 + K_matrix_global_ele2_11
	K_12 = K_matrix_global_ele1_12
	K_21 = K_matrix_global_ele1_21
	K_22 = K_matrix_global_ele1_22 + K_matrix_global_ele3_22
	K_all_element = sparse.bmat([[K_11,K_12],[K_21, K_22]])
	
	# set_printoptions(threshold='nan')
	print('K_all_element=\n',K_all_element)
	print('K_matrix_global_all=\n',K_matrix_global_all.todense())
	
	# 计算等效节点力
	F_node1 = np.array([0,-192000,-204800])
	F_node2 = np.array([0,-192000, 204800])
	F_all_node = np.concatenate((F_node1, F_node2), axis=0)

	# '''
	# 已知的边界条件（固定边界）
	delta_known_node_ID = np.array([3,4])
	delta_known_global = np.array([[0,0,0],[0,0,0]])

	# 位移边界条件处理方法——对角元素乘大数法BigNumber = 36854775807
	# import sys
	# max = sys.maxsize
	# print (max)
	# max = 9223372036854775807
	BigNumber = 36854775807
	K_matrix_global_all_dia = K_matrix_global_all.diagonal()
	for i_dia in range(len(delta_known_node_ID)):
		diaID = node_degree*(delta_known_node_ID[i_dia] - 1)
		K_matrix_global_all_dia[diaID:diaID+node_degree] = BigNumber*K_matrix_global_all_dia[diaID:diaID+node_degree]

	K_matrix_global_all_dia_BC = K_matrix_global_all.setdiag(K_matrix_global_all_dia)
	
	# 已知的节点荷载（若为均布荷载，需转化为节点对应不同自由度的直接荷载）
	force_known_node_ID = np.array([1,2])
	force_known_global = np.array([[0,-192000,-204800],[0,-192000, 204800]])

	F_row_list = []
	F_col_list = []
	F_val_list = []
	for i_rowF in range(len(force_known_node_ID)):
		F_rowID = node_degree*(force_known_node_ID[i_rowF] - 1)
		for j_Fr in range(node_degree):
			F_row_list.append(F_rowID+j_Fr)
			F_col_list.append(0)
			F_val_list.append(force_known_global[i_rowF,j_Fr])

	F_matrix_global_all = sparse.coo_matrix((F_val_list,(F_row_list,F_col_list)),dtype=np.float64)
	# '''

	delta_all_node = spsolve(K_all_element.tocsc(), F_all_node)

	print('delta_all_node = \n', delta_all_node)
	print(K_all_element.dot(delta_all_node)-F_all_node)

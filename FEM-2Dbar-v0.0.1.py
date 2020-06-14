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
		for i in range(nodeNum):
			print('请输入节点',i+1,"的x坐标")
			x_coo = float(input())
			print('请输入节点',i+1,"的y坐标")
			y_coo = float(input())
			global_coo_node[i,:] = np.array([x_coo, y_coo])
	elif para_dia == '3D':
		nodeNum = para_num
		global_coo_node = np.empty((nodeNum, 3))
		for i in range(nodeNum):
			print('请输入节点',i+1,"的x坐标")
			x_coo = float(input())
			print('请输入节点',i+1,"的y坐标")
			y_coo = float(input())
			print('请输入节点',i+1,"的z坐标")
			z_coo = float(input())
			global_coo_node[i,:] = np.array([x_coo, y_coo, z_coo])
	return global_coo_node

def func_get_element(para_type, para_eleNum):
	eleNum = para_eleNum
	eleType = para_type
	global_ele_num = np.empty((eleNum, eleType), dtype=int)
	for i in range(eleNum):
		print('请输入单元',i+1,"的首节点号")
		i_num = float(input())
		print('请输入单元',i+1,"的末节点号")
		j_num = float(input())
		global_ele_num[i,:] = np.array([i_num, j_num])
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

if __name__ == "__main__":

	# 截面惯性矩、截面面积、截面弹性模量
	sec_I = 15760e-8
	sec_A = 76.3e-4
	mat_E = 2e11

	# 节点编号
	diamension, nodeNum = '2D', 4	# diamension='2D'表示二维问题
	# global_coo_node = func_get_node(diamension, nodeNum)
	global_coo_node = np.array([[0,5],[6.4,5],[0,0],[9.6,0]])

	# 单元编号
	eleType, eleNum = 2, 3  # eleType=2表示两节点单元
	# element_nodes = func_get_element(eleType, eleNum)
	element_nodes = np.array([[1,2],[3,1],[2,4]])
	ele_length = np.empty(eleNum)
	ele_angle = np.empty(eleNum)

	for i in range(eleNum):
		sub_two_nodes = global_coo_node[element_nodes[i,1]-1,:] - global_coo_node[element_nodes[i,0]-1,:]
		ele_length[i] = np.sqrt(sub_two_nodes.dot(sub_two_nodes))
		if sub_two_nodes[0] == 0:
			if sub_two_nodes[1] > 0:
				ele_angle[i] = np.pi/2
			else:
				ele_angle[i] = -np.pi/2
		else:
			# np.arctan其实数部分的取值范围为[-pi/2, pi/2]
			ele_angle[i] = np.arctan(sub_two_nodes[1]/sub_two_nodes[0])

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

	for j in range(eleNum-1):
		K_matrix_local_ele_j1 = func_K_matrix(mat_E, sec_I, sec_A, ele_length[j+1])
		T_matrix_theta_ele_j1 = func_T_matrix(ele_angle[j+1])
		K_matrix_global_ele_j1 = (T_matrix_theta_ele_j1.dot(K_matrix_local_ele_j1)).dot(T_matrix_theta_ele_j1.T)

		K_matrix_local_ele = np.hstack((K_matrix_local_ele, K_matrix_local_ele_j1))
		T_matrix_theta_ele = np.hstack((T_matrix_theta_ele, T_matrix_theta_ele_j1))
		K_matrix_global_ele = np.hstack((K_matrix_global_ele, K_matrix_global_ele_j1))


	K_matrix_local_ele1 = func_K_matrix(mat_E, sec_I, sec_A, ele_l1)
	K_matrix_local_ele2 = func_K_matrix(mat_E, sec_I, sec_A, ele_l2)
	K_matrix_local_ele3 = func_K_matrix(mat_E, sec_I, sec_A, ele_l3)

	print('calibration of local Stiffness matrix \n', K_matrix_local_ele1==K_matrix_local_ele[0])

	# 局部坐标系-整体坐标系转换矩阵6×6
	T_matrix_theta1 = func_T_matrix(theta1)
	T_matrix_theta2 = func_T_matrix(theta2)
	T_matrix_theta3 = func_T_matrix(theta3)

	# 全局坐标系单元刚度矩阵6×6
	K_matrix_global_ele1 = (T_matrix_theta1.dot(K_matrix_local_ele1)).dot(T_matrix_theta1.T)
	K_matrix_global_ele2 = (T_matrix_theta2.dot(K_matrix_local_ele2)).dot(T_matrix_theta2.T)
	K_matrix_global_ele3 = (T_matrix_theta3.dot(K_matrix_local_ele3)).dot(T_matrix_theta3.T)

	print('calibration of global Stiffness matrix \n', K_matrix_global_ele1==K_matrix_global_ele[0])

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

	# 计算等效节点力
	F_node1 = np.array([0,-192000,-204800])
	F_node2 = np.array([0,-192000, 204800])
	F_all_node = np.concatenate((F_node1, F_node2), axis=0)

	delta_all_node = spsolve(K_all_element.tocsc(), F_all_node)

	print('delta_all_node = \n', delta_all_node)
	print(K_all_element.dot(delta_all_node)-F_all_node)

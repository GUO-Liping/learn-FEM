import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve
# new
def calcu_I(para, width=0, height=0, d=1):
	if para == 'rectangle':
		return width*height**3/12
	if para == 'triangle':
		return width*height**3/36
	if para == 'circle':
		return np.pi*d**4/64
	if para =='ring':
		return np.pi*(D**4-d**4)/64

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

	# 单元长度
	ele_l1 = 6.4
	ele_l2 = 5.0
	ele_l3 = 5.936

	# 单元局部坐标系相当于整体坐标系中的旋转角度
	theta1 = 0
	theta2 = np.pi/2
	theta3 = -np.arctan(5/3.2)

	# 局部坐标系单元刚度矩阵6×6
	K_matrix_local_ele1 = func_K_matrix(mat_E, sec_I, sec_A, ele_l1)
	K_matrix_local_ele2 = func_K_matrix(mat_E, sec_I, sec_A, ele_l2)
	K_matrix_local_ele3 = func_K_matrix(mat_E, sec_I, sec_A, ele_l3)

	# 局部坐标系-整体坐标系转换矩阵6×6
	T_matrix_theta1 = func_T_matrix(theta1)
	T_matrix_theta2 = func_T_matrix(theta2)
	T_matrix_theta3 = func_T_matrix(theta3)

	# 全局坐标系单元刚度矩阵6×6
	K_matrix_global_ele1 = (T_matrix_theta1.dot(K_matrix_local_ele1)).dot(T_matrix_theta1.T)
	K_matrix_global_ele2 = (T_matrix_theta2.dot(K_matrix_local_ele2)).dot(T_matrix_theta2.T)
	K_matrix_global_ele3 = (T_matrix_theta3.dot(K_matrix_local_ele3)).dot(T_matrix_theta3.T)

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

	delta_all_node = spsolve(K_all_element, F_all_node)

	print('delta_all_node = \n', delta_all_node)
	print(K_all_element.dot(delta_all_node)-F_all_node)

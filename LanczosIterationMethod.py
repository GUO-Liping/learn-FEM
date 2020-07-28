#!/usr/bin/python
# -*- coding:utf-8 -*-

'''
本程序用于实现有限单元法中的Lanczos迭代法（求解特征值问题）
This program is created to solve the eigenvalue problem in the finite element method by Lanczos iteration method.
'''

from scipy import linalg
import numpy as np
K_matrix = np.array([[2,-1,0,0,0],[-1,2,-1,0,0],[0,-1,2,-1,0],[0,0,-1,2,-1],[0,0,0,-1,1]])
M_matrix = np.array([[1,0,0,0,0],[0,1,0,0,0],[0,0,1,0,0],[0,0,0,1,0],[0,0,0,0,1/2]])

def func_Lanczos_method(K_matrix,M_matrix):
	# x_i = np.random.rand(4)
	x_start = np.array([1,1,1,1,1])
	gamma_start = np.sqrt(np.dot(x_start.dot(M_matrix), x_start.T))
	alpha_satrt = 0
	beta_start = 0
	x_1 = x_start / gamma_start
	
	x_1_ba = (np.dot(linalg.inv(K_matrix), M_matrix)).dot(x_1.T)
	alpha_1 = (np.dot(x_1_ba, M_matrix)).dot(x_1.T)
	x_1_bo = x_1_ba - alpha_1*x_1 - beta_start*x_start
	beta_1 = np.sqrt(np.dot(x_1_bo.dot(M_matrix),x_1_bo.T))
	
	loop_num = 2
	alpha_i = np.zeros(loop_num + 2)
	beta_i = np.zeros(loop_num + 2)
	x_i = np.zeros([loop_num + 3,5])
	T_matrix = np.zeros([loop_num+1,loop_num+1])

	alpha_i[0] = alpha_satrt
	alpha_i[1] = alpha_1

	beta_i[0] = beta_start
	beta_i[1] = beta_1

	x_i[0,:] = x_start
	x_i[1,:] = x_1
	x_i[2,:] = x_1_bo/beta_i[1]
	
	for i in range(loop_num ):

		x_i_ba = (np.dot(linalg.inv(K_matrix), M_matrix)).dot((x_i[i+2,:]).T)
		alpha_i[i+2] = (np.dot(x_i_ba, M_matrix)).dot((x_i[i+2,:]).T)
		x_i_bo = x_i_ba - alpha_i[i+2]*(x_i[i+2,:]) - beta_i[i+1]*(x_i[i+1,:])
		beta_i[i+2] = np.sqrt(np.dot(x_i_bo.dot(M_matrix),x_i_bo.T))
		x_i[i+3,:] = x_i_bo/beta_i[i+2]

		T_matrix[i,i] = alpha_i[i+1]
		T_matrix[i,i+1] = beta_i[i+1]
		T_matrix[i+1,i] = beta_i[i+1]
		T_matrix[i+1,i+1] = alpha_i[i+2]
	return T_matrix

if __name__ == "__main__":
	T_matrix = func_Lanczos_method(K_matrix,M_matrix)
	print(T_matrix)

#!/usr/bin/python
# -*- coding:utf-8 -*-

'''
本程序用于实现有限单元法中的逆迭代法（求解特征值问题）
This program is created to solve the eigenvalue problem in the finite element method by inverse iteration method.
'''

from scipy import linalg
import numpy as np
K_matrix = np.array([[2,-1,0,0],[-1,2,-1,0],[0,-1,2,-1],[0,0,-1,1]])
M_matrix = np.array([[0,0,0,0],[0,2,0,0],[0,0,0,0],[0,0,0,1]])

def func_inverse_humam(K_matrix,M_matrix):
	# x_i = np.random.rand(4)
	x_i = np.array([1,1,1,1])
	x_i_norm = x_i/((np.dot(x_i.T,M_matrix).dot(x_i))**0.5)
	
	for i in range(3):
		y_i = M_matrix.dot(x_i_norm)
		x_i = linalg.inv(K_matrix).dot(y_i)
		x_i_norm = x_i/((np.dot(x_i.T,M_matrix).dot(x_i))**0.5)
	
	lambda_1 = (((x_i.T).dot(K_matrix)).dot(x_i)) / (np.dot(x_i.T,M_matrix.dot(x_i)))
	phi_1 = x_i_norm
	return lambda_1,phi_1

def func_inverse_computer(K_matrix,M_matrix):
	# x_i = np.random.rand(4)
	x_i = np.array([1,1,1,1])
	y_i = M_matrix.dot(x_i)

	for i in range(3):
		y_i_norm = y_i / (x_i.dot(y_i))**0.5
		x_i = linalg.inv(K_matrix).dot(y_i_norm)
		y_i = M_matrix.dot(x_i)

	lambda_1 = x_i.dot(y_i_norm) / (x_i.dot(y_i))
	phi_1 = x_i/(x_i.dot(y_i))**0.5
	return lambda_1,phi_1

if __name__ == "__main__":
	lambda_1,phi_1 = func_inverse_computer(K_matrix,M_matrix)
	print(' lambda_1 = ',lambda_1, '\n','phi_1 = ',phi_1)

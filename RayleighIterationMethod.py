#!/usr/bin/python
# -*- coding:utf-8 -*-

'''
本程序用于实现有限单元法中的Rayleigh商迭代法（求解特征值问题）
This program is created to solve the eigenvalue problem in the finite element method by Rayleigh iteration method.
'''

from scipy import linalg
import numpy as np
Lambda_matrix = np.array([[2,0],[0,6]])
I_matrix = np.eye(2)

def func_Rayleigh_human(Lambda_matrix, I_matrix):
	# x_i = np.random.rand(4)
	x_i = np.array([1,0.1])
	rho_i = 0
	
	for i in range(2):
		x_i_before = x_i
		rho_i_before = rho_i
		Lambda_matrix_subs = Lambda_matrix - rho_i_before*I_matrix
		x_i = linalg.inv(Lambda_matrix_subs).dot(x_i_before)
		rho_i = (x_i.dot(x_i_before.T)) / (x_i.dot(x_i.T)) + rho_i_before
	
	return rho_i

if __name__ == "__main__":
	rho = func_Rayleigh_human(Lambda_matrix,I_matrix)
	print(' rho = ',rho)

#!/usr/bin/python
# -*- coding:utf-8 -*-

'''
本程序用于实现有限单元法中的子空间迭代法（求解特征值问题）
This program is created to solve the eigenvalue problem in the finite element method by Subspace iteration method.
'''

from scipy import linalg
import numpy as np

def func_Subspace_method(K_matrix,M_matrix):
	x_i = np.array([[0,0],[2,0],[0,0],[0,1]])
	for i in range(1):
		x_i_ba = np.dot(linalg.inv(K_matrix),x_i)
		K_i = (np.dot(x_i_ba.T, K_matrix)).dot(x_i_ba)
		M_i = (np.dot(x_i_ba.T, M_matrix)).dot(x_i_ba)
		Lambda_i,phi = linalg.eig(K_i,M_i)
		# Lambda_matrix = np.diag(Lambda_i.real[::-1])  # 若有需要则进行倒序排列
		# Phi_matrix = phi[:,::-1]  # 若有需要则进行倒序排列
		Lambda_matrix = np.diag(Lambda_i.real)
		Q_matrix = phi
		x_i = x_i_ba.dot(Q_matrix)

	return Lambda_matrix, Q_matrix

if __name__ == "__main__":

	K_matrix = np.array([[2,-1,0,0],[-1,2,-1,0],[0,-1,2,-1],[0,0,-1,1]])
	M_matrix = np.array([[0,0,0,0],[0,2,0,0],[0,0,0,0],[0,0,0,1]])

	Lambda_matrix, Q_matrix = func_Subspace_method(K_matrix,M_matrix)
	print('Lambda_matrix=',Lambda_matrix,'\n', 'Q_matrix=',Q_matrix)

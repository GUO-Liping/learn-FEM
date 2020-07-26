#!/usr/bin/python
# -*- coding:utf-8 -*-

'''
本程序用于实现有限单元法中的Jacobi迭代法（求解特征值问题）
This program is created to solve the eigenvalue problem in the finite element method by Jacobi iteration method.
'''

from scipy import linalg
import numpy as np

def func_Jacobi_Iteration(K_matrix):
	Phi_matrix = np.eye(4)
	for k in range(10):
		for i in range(4-1):
			for j in range(i+1,4):
				P_matrix = np.eye(4)
				theta = np.arctan(2*K_matrix[i,j]/(K_matrix[i,i]-K_matrix[j,j]))/2

				P_matrix[i,i] = np.cos(theta)
				P_matrix[i,j] =-np.sin(theta)
				P_matrix[j,i] = np.sin(theta)
				P_matrix[j,j] = np.cos(theta)
	
				K_matrix = np.dot((P_matrix.T).dot(K_matrix), P_matrix)
				Phi_matrix = Phi_matrix.dot(P_matrix)
		K_matrix = K_matrix
	return K_matrix, Phi_matrix

if __name__ == "__main__":
	K_matrix = np.array([[5,-4,1,0],[-4,6,-4,1],[1,-4,6,-4],[0,1,-4,5]])
	M_matrix = np.eye(4)
	K_matrix, Phi_matrix = func_Jacobi_Iteration(K_matrix)
	print(' Phi_matrix, = ',Phi_matrix)
	print(' K_matrix, = ',K_matrix)

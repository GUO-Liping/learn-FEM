import sympy as sp

def triangleArea(x_i,y_i,x_j,y_j,x_m,y_m):
	vector_A = sp.Matrix([[1,x_i,y_i],[1,x_j,y_j],[1,x_m,y_m]])
	return 0.5*(vector_A.det())

if __name__ == '__main__':
	x_i,y_i,x_j,y_j,x_m,y_m = sp.symbols('x_i,y_i,x_j,y_j,x_m,y_m')
	E, miu, t = sp.symbols('E, miu, t')

	'''
	type_number = int(input('请选择三角形单元类型对应的数字：1平面应力；2平面应变'))
	
	if type_number == 2:
		E = E/(1-miu**2)
		miu = miu/(1-miu)
	else:
		pass
	'''

	G = E/(2*(1+miu))
	vector_D = sp.Matrix([[E/(1-miu**2),E*miu/(1-miu**2),0],[E*miu/(1-miu**2),E/(1-miu**2),0],[0,0,G]])

	a_i = x_j*y_m - x_m*y_j
	b_i = -(y_m-y_j)
	c_i = x_m-x_j

	a_j = x_m*y_i - x_i*y_m
	b_j = -(y_i-y_m)
	c_j = x_i-x_m

	a_m = x_i*y_j - x_j*y_i
	b_m = -(y_j-y_i)
	c_m = x_j-x_i

	Area = triangleArea(x_i,y_i,x_j,y_j,x_m,y_m)
	vector_B = 1/(2*Area)*sp.Matrix([[b_i,0,b_j,0,b_m,0],[0,c_i,0,c_j,0,c_m],[c_i,b_i,c_j,b_j,c_m,b_m]])
	vector_S = vector_D * (vector_B)
	vector_Ke = (vector_B.T) * (vector_D * (vector_B))*t*Area
	# print('vector_Ke=', vector_Ke)

	my_dict = {x_i:0, y_i:0, x_j:3, y_j:0, x_m:0, y_m:4, E:2e5, miu:0.2, t:1}

	for i in range(len(vector_Ke.row(0))):
	    for j in range(len(vector_Ke.col(0))):
	    	vector_Ke[i,j] = vector_Ke[i,j].evalf(subs=my_dict)
	    	print(vector_Ke[i,j])

	'''
	x_i,y_i = map(float,input('请输入三角形单元第一个顶点xi,yi坐标，数字间逗号隔开：').split(','))
	x_j,y_j = map(float,input('请输入三角形单元第二个顶点xj,yj坐标，数字间逗号隔开：').split(','))
	x_m,y_m = map(float,input('请输入三角形单元第一个顶点xm,ym坐标，数字间逗号隔开：').split(','))
	E = np.float(input('请输入弹性模量值(MPa)：'))
	miu = np.float(input('请输入泊松比：'))
	t = np.float(input('请输入单元厚度：'))
	'''
	print('vector_Ke=', vector_Ke/8581)
    


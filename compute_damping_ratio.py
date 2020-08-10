# 该程序用于求解阻尼反映
# 数据源为落石冲击试验球体运动轨迹捕捉数据中的振幅最大值与平均周期
import numpy as np

v1 = np.array([-0.95,-0.66,-0.49,-0.24,-0.19])
T1 = 2*(0.29+0.34+0.34+0.32+0.34)/5
w_D1 = 2*np.pi/T1
v3 = np.array([-0.45241, -0.32577, -0.21952, -0.13367, -0.09905])
T3 = 2*((0.344-0.186)+(1.924-1.726)+(2.992-2.764)+(3.784-3.544)+(4.46-4.176))/5
w_D3 = 2*np.pi/T3

vi = 0
vj = 0
m1 = 144*7850*np.pi*0.3*19*np.pi*0.003**2/4
m2=2370
mass_all = m1 + m2

ksi = 0
for i in range(4):
	vi = v3[i]
	for j in range(4-i):
		vj = v3[i+j+1]
		n=i
		m=j+1
		s = vi/vj
		ksi_sim = np.log(s)/(2*np.pi*m)
		ksi_com = np.sqrt(1/(1+(np.log(s)/(2*np.pi*m))**-2))

		# print('n=',n+1,'m+n=',m+n+1,'ksi_sim ',ksi_sim)
		print('n=',n+1,'m+n=',m+n+1,'ksi_com=',ksi_com)



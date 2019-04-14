# -*- coding: utf-8 -*-
"""
Created on Wed Apr 10 16:49:00 2019

@author: GX
"""
import math
import matplotlib.pyplot as plt
pi=3.1416

dict1={9.1:0.0005, 9.2:0.0009, 9.3:0.0050, 9.4:0.0357, 9.5:0.0335, 9.6:0.0197, 9.7:0.0273,
       9.8:0.0153,9.9:0.0109,10.0:0.0011,10.1:0.0020,10.2:0.0117,10.3:0.0127,10.4:0.0009,10.5:0.0005,
       10.6:0.0005,10.7:0.0003,10.8:0.0005,10.9:0.0014,11.0:0.0014,11.1:0.0016,11.2:0.0022,
       11.3:0.0029,11.4:0.0035,11.5:0.0041,11.6:0.0046,11.7:0.0046,11.8:0.0035,11.9:0.0022,
       12.0:0.0007,12.1:0.0005,12.2:0.0005,12.3:0.0098,12.4:0.0304,12.5:0.0660,12.6:0.1027,
       12.7:0.0212,12.8:0.0524,12.9:0.0781,13.0:0.0462,13.1:0.0524,13.2:0.1109,13.3:0.2465,
       13.4:0.3352,13.5:0.5305} #二氧化碳海平面水平路程上光谱吸收系数

dict2={9.1:0.1035, 9.2:0.0996, 9.3:0.0965, 9.4:0.0941, 9.5:0.0918, 9.6:0.0895, 9.7:0.0880,
       9.8:0.0864,9.9:0.0849,10.0:0.0841,10.1:0.0841,10.2:0.0841,10.3:0.0833,10.4:0.0826,10.5:0.0818,
       10.6:0.0811,10.7:0.0811,10.8:0.0818,10.9:0.0833,11.0:0.0841,11.1:0.0849,11.2:0.0965,
       11.3:0.1027,11.4:0.0980,11.5:0.0957,11.6:0.0903,11.7:0.1342,11.8:0.0996,11.9:0.0949,
       12.0:0.0880,12.1:0.0864,12.2:0.0864,12.3:0.0880,12.4:0.0910,12.5:0.0934,12.6:0.0957,
       12.7:0.0996,12.8:0.1035,12.9:0.1075,13.0:0.1131,13.1:0.1155,13.2:0.1203,13.3:0.1251,
       13.4:0.1300,13.5:0.1350,
       7.5:1.8208,7.6:2.7326,7.7:0.7791,7.8:0.8798,7.9:0.5660,8.0:0.3419,8.1:0.1909,8.2:0.2450,
       8.3:0.1628,8.4:0.1732,8.5:0.1945,8.6:0.1205,8.7:0.1187,8.8:0.1171,8.9:0.1155,9.0:0.1115}
       #气温5摄氏度，100%湿度时海平面水平路程大气水汽光谱吸收系数
      
dict3={-20:0.89,-21:0.81,-22:0.74,-23:0.67,-24:0.61,-25:0.56,
       -10:2.15,-11:1.98,-12:1.81,-13:1.66,-14:1.52,-15:1.40,-16:1.28,-17:1.18,-18:1.08,-19:0.98,
       0:4.84,-1:4.47,-2:4.13,-3:3.81,-4:3.52,-5:3.24,-6:2.99,-7:2.75,-8:2.54,-9:2.34,
       1:5.18,2:5.54,3:5.92,4:6.33,5:6.76,6:7.22,7:7.70,8:8.22,9:8.76,10:9.33,11:9.94,
       12:10.57,13:11.25,14:11.96,15:12.0,16:13.50,17:14.34,18:15.22,19:16.14,20:17.22,
       21:18.14,22:19.22,23:20.36,24:21.55,25:22.80,26:24.11,27:25.49,28:27.00,29:28.45,
       30:30.04,31:31.70,32:33.45,33:35.28,34:37.19,35:39.19
       }#空气湿度100%和不同温度时的可凝结水毫米数 mm/km


def tau4(Jr,Js,R):      #气象(雨、雪)衰减
    J1=math.pow(Jr,0.66)
    J2=math.pow(Js,0.7)
    x=math.exp(-0.66*J1*R)
    y=math.exp(-6.5*J2*R)
    return x*y

def tau3(lamda,Dv,R):   #大气散射衰减
    if Dv<6:
        q=0.585*math.pow(Dv,0.33)
    elif Dv>50:
        q=1.6
    else:
        q=1.3
    return math.exp((-3.91/Dv)*math.pow((0.55/lamda),q)*R)

def tau2_h(lamda,R):    #水平路径二氧化碳吸收衰减
    if lamda<=9.0:
        miu=0
    else:
        miu=dict1[lamda]
    return math.exp(-miu*R)

def tau1_h(lamda,R,t,f):  #水平路径水汽吸收衰减
    w=dict3[t]
    miu0=dict2[lamda]
    miu=(w/6.76)*f*miu0
    return math.exp(-miu*R)

def average_tau_h(dlamda,delta,i):  #计算倾斜路径平均大气透射率
    lamda_1=7.5
    lamda_2=13.5
    tau_1=tau1_h(lamda,R,t,f)*tau2_h(lamda,R)*tau3(lamda_1,Dv,R)*tau4(Jr,Js,R)
    tau_2=tau1_h(lamda,R,t,f)*tau2_h(lamda,R)*tau3(lamda_2,Dv,R)*tau4(Jr,Js,R)
    n=dlamda/delta
    return (n*(0.5*(tau_1+tau_2)+i))

def tau2_v(lamda,H1,H2,alpha):    #倾斜路径二氧化碳吸收衰减
    alpha=alpha/180*pi
    if lamda<=9.0:
        miu=0
    else:
        miu=dict1[lamda]
    R=(math.exp(-0.313*H1)-math.exp(-0.313*H2))/(0.313*math.cos(alpha))
    return math.exp(-miu*R)

def tau1_v(lamda,H1,H2,t,f,alpha):  #倾斜路径水汽吸收衰减
    w=dict3[t]
    miu0=dict2[lamda]
    miu=(w/6.76)*f*miu0
    alpha=alpha/180*pi
    R=(math.exp(-0.5154*H1)-math.exp(-0.5154*H2))/(0.5154*math.cos(alpha))
    return math.exp(-miu*R)

  
def average_tau_v(dlamda,delta,i):  #计算倾斜路径平均大气透射率
    lamda_1=7.5
    lamda_2=13.5
    tau_1=tau1_v(lamda_1,H1,H2,t,f,alpha)*tau2_v(lamda_1,H1,H2,alpha)*tau3(lamda_1,Dv,R)*tau4(Jr,Js,R)
    tau_2=tau1_v(lamda_2,H1,H2,t,f,alpha)*tau2_v(lamda_2,H1,H2,alpha)*tau3(lamda_2,Dv,R)*tau4(Jr,Js,R)
    n=dlamda/delta
    return (n*(0.5*(tau_1+tau_2)+i))
    
if __name__ == "__main__":
    mode=int(input('请输入工作模式 1 辐射传输路径水平；2 辐射传输路径倾斜；0 退出: '))
    #R=float(input('请输入辐射传输距离R(km): '))
    #lamda=float(input('请输入波长λ(μm): '))   
    if mode==2:       
        Dv=float(input('请输入能见度Dv(km): '))
        Jr=float(input('请输入降雨量Jr(mm/h): '))
        Js=float(input('请输入降雪量Js(mm/h): '))
        t=float(input('请输入气温t(℃): '))
        f=float(input('请输入相对湿度f: '))
        H1=float(input('请输入起始点海拔高度H1(km): '))
        H2=float(input('请输入终点海拔高度H2(km): '))
        alpha=float(input('请输入海平面法线与辐射传输方向之间的夹角α(°): '))
        R=(H2-H1)/math.cos(alpha/180*pi)
        lamdalist=[i/10.0 for i in range(75,136)]
        taulist=[]
        i=0
        for lamda in lamdalist:
            #a=tau4(Jr,Js,R)
            #b=tau3(lamda,Dv,R)
            #c=tau2_v(lamda,H1,H2,alpha)
            #d=tau1_v(lamda,H1,H2,t,f,alpha)
            tau=tau1_v(lamda,H1,H2,t,f,alpha)*tau2_v(lamda,H1,H2,alpha)*tau3(lamda,Dv,R)*tau4(Jr,Js,R)
            i+=tau
            taulist.append(tau)
            print('经计算，大气透射率为{}'.format(tau))
        dlamda=0.1
        delta=6
        average_tau_v=average_tau_v(dlamda,delta,i)
        print('大气平均透射率为{}'.format(average_tau_v))
        plt.xlabel("λ /μm")
        plt.ylabel("Atmospheric transmittance τ")
        plt.title("7.5-13.5μm Atmospheric transmittance")
        plt.grid(True)
        plt.plot(lamdalist,taulist)
        
        
    if mode==1:
        Dv=float(input('请输入能见度Dv(km): '))
        Jr=float(input('请输入降雨量Jr(mm/h): '))
        Js=float(input('请输入降雪量Js(mm/h): '))
        t=float(input('请输入气温t(℃): '))
        f=float(input('请输入相对湿度f: '))
        R=float(input('请输入辐射传输距离R(km): '))
        lamdalist=[i/10.0 for i in range(75,136)]
        taulist=[]
        i=0
        for lamda in lamdalist:
            #a=tau4(Jr,Js,R)
            #b=tau3(lamda,Dv,R)
            #c=tau2_h(lamda,R)
            #d=tau1_h(lamda,R,t,f)
            tau=tau1_h(lamda,R,t,f)*tau2_h(lamda,R)*tau3(lamda,Dv,R)*tau4(Jr,Js,R)
            i+=tau
            taulist.append(tau)
            print('经计算，大气透射率为{}'.format(tau))
        dlamda=0.1
        delta=6
        average_tau_h=average_tau_h(dlamda,delta,i)
        print('大气平均透射率为{}'.format(average_tau_h))
        plt.xlabel("λ /μm")
        plt.ylabel("Atmospheric transmittance τ")
        plt.title("7.5-13.5μm Atmospheric transmittance")
        plt.grid(True)
        plt.plot(lamdalist,taulist)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
        
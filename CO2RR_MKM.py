from numpy import exp
from sympy import solve,symbols,nsolve,Eq
import numpy as np
from math import sqrt,log10
import xlrd,xlwt
def K_gas_dis(G):
    kb = 1.38e-23 * 6.24e18  # 玻尔兹曼常数 将J转化为eV
    A=kb*T/h #指前因子 6.2e12
    q=7.148e-29
    B=2416.7*0.0833
    K=B*exp(-G/(kb*T))
    return K
def K_CO2_dis(G):
    kb = 1.38e-23 * 6.24e18  # 玻尔兹曼常数 将J转化为eV
    B=0.478/18
    K=B*exp(-G/(kb*T))
    return K
def K_gas(G):
    kb = 1.38e-23 * 6.24e18  # 玻尔兹曼常数 将J转化为eV
    K=exp(-G/(kb*T))
    return K
def K_e(Ea):
    kb=1.38e-23*6.24e18 #玻尔兹曼常数 将J转化为eV
    h=4.136e-15 #普朗克常数
    A=kb*T/h #指前因子 6.2e12
    K=A*exp((-Ea)/(kb*T))
    return K
def judge(a):
    b=[]
    for i in a:
        m = True
        for j in range(len(i)):
            if i[j]>1 or i[j]<0:
                m=False
                break
        if m:
            b.append(i)
    return b

def over(a):
    m=False
    for i in range(len(a)):
        if a[i]>1:
            m=True
            break
    return m
def sigmoid(x):
    a=1/(1+exp(-x))
    return a

if __name__=='__main__':

    T=210+273.15 #210度
    kb = 1.38e-23 * 6.24e18  # 玻尔兹曼常数 将J转化为eV
    h = 4.136e-15  # 普朗克常数
    A = kb * T / h  # 指前因子 6.2e12
    workbook = xlrd.open_workbook(r'C:\Users\KD\Desktop\read9.xlsx')
    content = workbook.sheet_by_index(10)  # choose different index
    rownum = content.nrows
    colnum = content.ncols
    name = content.col_values(0)
    # print(name)


    #单个点的测试
    # E_H=[0.98,1.18,1.30,1.34,0.72,1.25]
    # Ea=[0.47,0.08,0.07,0.23,0.93,0.84] #能垒
    P_H2 = 24e5
    S = 205.1 * 1.03e-5  # 将1J/mol变成0.0103eV 这里取得是CO的熵值
    E_ads = 0.00076
    P_CO2 = 32e5 - P_H2
    conversion = 0.02
    P_CH3OH = conversion * P_CO2
    P_H2O = P_CH3OH
    # K_H2=[K_e(x) for x in E_H]
    # K = [K_e(x) for x in Ea]
    # a=0
    # for i in K_H2:
    #     a=a+sqrt(i*P_H2)
    # o_0=1/(a+1)
    # o_H=[]
    # for i in K_H2:
    #     o_H.append(sqrt(i*P_H2)*o_0)
    # o_COOH,o_COH2,o_CHOH2,o_CHOH,o_CH2OH,o_1 = \
    #     symbols('o_COOH o_COH2 o_CHOH2 o_CHOH o_CH2OH o_1')
    # Equation1 = [
    #     Eq(K[0]*P_CO2*o_H[0]*o_1-o_COOH*K[1]*o_H[1], 0),
    #     Eq(K[1] * o_COOH*o_H[1] - K[2]*o_COH2*o_H[2], 0),
    #     Eq(K[2] * o_COH2 * o_H[2] - K[3] * o_CHOH2 * o_H[3], 0),
    #     Eq(K[3] * o_CHOH2 * o_H[3] - K[4] * o_CHOH * o_H[4], 0),
    #     Eq(K[4] * o_CHOH * o_H[4] - K[5] * o_CH2OH * o_H[5], 0),
    #     Eq(o_1+o_CHOH+o_CHOH2+o_CH2OH+o_COH2+o_COOH, 1)
    # ]
    # solution_eq1 = solve(Equation1, [o_COOH,o_COH2,o_CHOH2,o_CHOH,o_CH2OH,o_1])
    # rate=solution_eq1[o_CH2OH]*K[5]*o_H[5]


# 循环读写
    outpath = r"C:\Users\KD\Desktop\n-5.dat"
    wp = open(outpath, 'w')
    for n in range(1, rownum):

        E_H = []
        Ea = []
        Ea_re=[]
        E_H_b = []
        Ea_b = []
        Ea_re_b = []
        E_H2 = content.cell(n, colnum - 2).value #倒数第二列是H2吸附
        E_H2_dis=content.cell(n, colnum - 1).value #最后一列是H2脱附

        #获得path：a-a-a-a-a-a
        for j in range(1, 13, 2):
            E_H.append(content.cell(n, j).value)
            Ea.append(content.cell(n, j + 1).value)
        for i in range(13,19):
            Ea_re.append(content.cell(n,i).value)
        #获得path：b-b-b-b-b-b
        for j in range(19, 31, 2):
            E_H_b.append(content.cell(n, j).value)
            Ea_b.append(content.cell(n, j + 1).value)
        for i in range(31,37):
            Ea_re_b.append(content.cell(n,i).value)
        K_H2 = [K_gas(x) for x in E_H]
        K = [K_e(x) for x in Ea]
        K_re=[K_e(x) for x in Ea_re]
        K_dis_H2=K_gas_dis(E_H2)
        K_H2_disorption=K_e(E_H2_dis)

        K_H2_b=[K_gas(x) for x in E_H_b]
        K_b=[K_e(x) for x in Ea_b]
        K_re_b = [K_e(x) for x in Ea_re_b]
        print(name[n])
        #对于a-a-a-a-a-a path
        o_COOH, o_COH2, o_CHOH2, o_CHOH, o_CH2OH, o_1 = \
            symbols('o_COOH o_COH2 o_CHOH2 o_CHOH o_CH2OH o_1')
        o_alpha_0=1
        o_dis_H2=K_H2[0]/(1+sqrt(K_H2_disorption/(K_dis_H2*P_H2))) #alpha H
        # print(K_H2[0])
        o_0=[1,1 - o_dis_H2]
        K_CO2 = K_CO2_dis(Ea[0])
        K_CO2_re = K_e(Ea_re[0])
        COOH = ((K_CO2 * P_CO2) / K_CO2_re) * (o_dis_H2 / (1 - o_dis_H2))
        o_H = [o_alpha_0]
        for m in K_H2:
            o_H.append(m * o_dis_H2) #boltzman H coverage
        o_H.append(o_alpha_0)
        o_H_b = [o_alpha_0]
        for m in K_H2_b:
            o_H_b.append(m * o_dis_H2)
        o_H_b.append(o_alpha_0)
        Equation1 = [
            Eq(K[0] * P_CO2 * o_H[1] * o_1 - o_COOH * K[1] * o_H[2] -K_b[1] * o_COOH * o_H_b[2]- o_COOH * K_re[0]* o_0[0] + o_COH2 * K_re[1]* o_0[0]+o_COH2 * K_re_b[1]* o_0[1], 0), #*COOH
            Eq(K[1] * o_COOH * o_H[2] + K_b[1] * o_COOH * o_H_b[2] - K[2] * o_COH2 * o_H[3] - K_b[2] * o_COH2 *o_H_b[3] - o_COH2 * K_re[1] * o_0[0] - o_COH2 * K_re_b[1] * o_0[1] + o_CHOH2 * K_re[2] * o_0[0] + o_CHOH2 * K_re_b[2] * o_0[1], 0),  # *COHOH的 即*COH2
            Eq(K[2] * o_COH2 * o_H[3] + K_b[2] * o_COH2 * o_H_b[3] - K[3] * o_CHOH2 * o_H[4] - K_b[3] * o_CHOH2 * o_H_b[4] - o_CHOH2 * K_re[2] * o_0[0] - o_CHOH2 * K_re_b[2] * o_0[1], 0),# *CHOHOH的 即*CHOH2
            Eq(K[3] * o_CHOH2 * o_H[4] + K_b[3] * o_CHOH2 * o_H_b[4] - K[4] * o_CHOH * o_H[5] - K_b[4] * o_CHOH * o_H_b[5] + o_CH2OH * K_re[4] * o_0[0] + o_CH2OH * K_re_b[4] * o_0[1], 0),  # *CHOH的
            Eq(K[4] * o_CHOH * o_H[5] + K_b[4] * o_CHOH * o_H_b[5] - K[5] * o_CH2OH * o_H[6] - K_b[5] * o_CH2OH * o_H_b[6] - o_CH2OH * K_re[4] * o_0[0] - o_CH2OH * K_re_b[4] * o_0[1], 0),# *CH2OH的
            Eq(o_1+o_CHOH + o_CHOH2 + o_CH2OH + o_COH2 + o_COOH, 1)
        ]
        lambda_a=lambda a,b:a/b
        K_equil=[lambda_a(K[i],K_re[i]) for i in range(len(K))]
        K_mod=[sigmoid(log10(i)) for i in K_equil]
        print(K_mod)
        # print(K_equil)

        # solution_eq1 = nsolve(Equation1, [o_COOH, o_COH2, o_CHOH2, o_CHOH, o_CH2OH,o_1],[0,0,0,0,0,0])

        solution_eq1 = solve(Equation1, [o_COOH, o_COH2, o_CHOH2, o_CHOH, o_CH2OH,o_1])
        rate = solution_eq1[o_CH2OH] * K[5] * o_H[6]+solution_eq1[o_CH2OH] * K_b[5] * o_H_b[6]
        # print(solution_eq1)
        # if rate<0:
        #     rate=-1*rate
        print("alpha rate: %.2e" %rate)

        # print(solution_eq1)

        #对于b-b-b-b-b-b path
        o_COOH_b, o_COH2_b, o_CHOH2_b, o_CHOH_b, o_CH2OH_b, o_1_b = \
            symbols('o_COOH_b o_COH2_b o_CHOH2_b o_CHOH_b o_CH2OH_b o_1_b')
        o_dis_H2_belta = K_H2_b[0] / (1 + sqrt(K_H2_disorption / (K_dis_H2 * P_H2)))  # belta H
        K_CO2_b = K_CO2_dis(Ea_b[0])
        K_CO2_re_b = K_e(Ea_re_b[0])
        COOH_b = ((K_CO2_b * P_CO2) / K_CO2_re_b) * o_dis_H2_belta * (1 - o_dis_H2)

        Equation2 = [
            Eq( K_b[0] * P_CO2 * o_H_b[1] * o_1_b - o_COOH_b * K_b[1] * o_H_b[2]-o_COOH_b * K[1] * o_H[2]- o_COOH_b * K_re_b[0] * o_0[1] + o_COH2_b * K_re_b[1] *o_0[1]+o_COH2_b * K_re[1] *o_0[0], 0), #*COOH
            Eq(K[1] * o_COOH_b * o_H[2] + K_b[1] * o_COOH_b * o_H_b[2] - K[2] * o_COH2_b * o_H[3] - K_b[2] * o_COH2_b * o_H_b[3] - o_COH2_b * K_re[1] * o_0[0] - o_COH2_b * K_re_b[1] * o_0[1] + o_CHOH2_b * K_re[2] * o_0[0] + o_CHOH2_b * K_re_b[2] * o_0[1], 0),  # *COHOH的 即*COH2
            Eq(K[2] * o_COH2_b * o_H[3] + K_b[2] * o_COH2_b * o_H_b[3] - K[3] * o_CHOH2_b * o_H[4] - K_b[3] * o_CHOH2_b * o_H_b[4] - o_CHOH2_b * K_re[2] * o_0[0] - o_CHOH2_b * K_re_b[2] * o_0[1], 0),# *CHOHOH的 即*CHOH2
            Eq(K[3] * o_CHOH2_b * o_H[4] + K_b[3] * o_CHOH2_b * o_H_b[4] - K[4] * o_CHOH_b * o_H[5] - K_b[4] * o_CHOH_b * o_H_b[5] + o_CH2OH_b * K_re[4] * o_0[0] + o_CH2OH_b * K_re_b[4] * o_0[1], 0),  # *CHOH的
            Eq(K[4] * o_CHOH_b * o_H[5] + K_b[4] * o_CHOH_b * o_H_b[5] - K[5] * o_CH2OH_b * o_H[6] - K_b[5] * o_CH2OH_b * o_H_b[6] - o_CH2OH_b * K_re[4] * o_0[0] - o_CH2OH_b * K_re_b[4] * o_0[1], 0),
            # *CH2OH的
            Eq(o_1_b+o_CHOH_b + o_CHOH2_b + o_CH2OH_b + o_COH2_b + o_COOH_b,1-o_dis_H2)
        ]
        solution_eq2 = solve(Equation2, [o_1_b,o_COOH_b, o_COH2_b, o_CHOH2_b, o_CHOH_b, o_CH2OH_b])
        rate_b=solution_eq2[o_CH2OH_b] * K[5] * o_H[6]+solution_eq2[o_CH2OH_b] * K_b[5] * o_H_b[6]
        print("belta rate: %.2e" %rate_b)

        # 对于 a-a-a-a-a-a path and b-b-b-b-b-b path
        o_COOH_s, o_COH2_s, o_CHOH2_s, o_CHOH_s, o_CH2OH_s, o_1_s = \
            symbols('o_COOH_s o_COH2_s o_CHOH2_s o_CHOH_s o_CH2OH_s o_1_s')
        o_dis_H2_belta = K_H2_b[0] / (1 + sqrt(K_H2_disorption / (K_dis_H2 * P_H2)))  # belta H
        K_CO2_b = K_CO2_dis(Ea_b[0])
        K_CO2_re_b = K_e(Ea_re_b[0])
        COOH_b = ((K_CO2_b * P_CO2) / K_CO2_re_b) * o_dis_H2_belta * (1 - o_dis_H2)
        o_H_b = [o_alpha_0]
        for m in K_H2_b:
            o_H_b.append(m * o_dis_H2)
        o_H_b.append(o_alpha_0)
        Equation3 = [
            Eq(K[0] * P_CO2 * o_H[1] * o_1_s+K_b[0] * P_CO2 * o_H_b[1] * o_1_s - o_COOH_s * K[1] * o_H[2]-o_COOH_s * K_b[1] * o_H_b[2] - o_COOH_s * K_re[0] * o_0[0]-o_COOH_s * K_re_b[0] * o_0[1] + o_COH2_s * K_re[1] * o_0[0]+o_COH2_s * K_re_b[1] * o_0[1], 0),  # *COOH
            Eq(K[1] * o_COOH_s * o_H[2] + K_b[1] * o_COOH_s * o_H_b[2] - K[2] * o_COH2_s * o_H[3] - K_b[2] * o_COH2_s * o_H_b[3] - o_COH2_s * K_re[1] * o_0[0] - o_COH2_s * K_re_b[1] * o_0[1] + o_CHOH2_s * K_re[2] * o_0[0] + o_CHOH2_s *K_re_b[2] * o_0[1], 0),  # *COHOH的 即*COH2
            Eq(K[2] * o_COH2_s * o_H[3] + K_b[2] * o_COH2_s * o_H_b[3] - K[3] * o_CHOH2_s * o_H[4] - K_b[3] * o_CHOH2_s * o_H_b[4] - o_CHOH2_s * K_re[2] * o_0[0] - o_CHOH2_s * K_re_b[2] * o_0[1], 0),  # *CHOHOH的 即*CHOH2
            Eq(K[3] * o_CHOH2_s * o_H[4] + K_b[3] * o_CHOH2_s * o_H_b[4] - K[4] * o_CHOH_s * o_H[5] - K_b[4] * o_CHOH_s * o_H_b[5] + o_CH2OH_s * K_re[4] * o_0[0] + o_CH2OH_s * K_re_b[4] * o_0[1], 0),  # *CHOH的
            Eq(K[4] * o_CHOH_s * o_H[5] + K_b[4] * o_CHOH_s * o_H_b[5] - K[5] * o_CH2OH_s * o_H[6] - K_b[5] * o_CH2OH_s * o_H_b[6] - o_CH2OH_s * K_re[4] * o_0[0] - o_CH2OH_s * K_re_b[4] * o_0[1], 0),  # *CH2OH的
            Eq(o_1_s+o_CHOH_s + o_CHOH2_s + o_CH2OH_s + o_COH2_s + o_COOH_s, 1-o_dis_H2)
        ]

        solution_eq3 = solve(Equation3, [o_1_s, o_COOH_s, o_COH2_s, o_CHOH2_s, o_CHOH_s, o_CH2OH_s])
        rate_s = solution_eq3[o_CH2OH_s] * K[5] * o_H[6]+solution_eq3[o_CH2OH_s] * K_b[5] * o_H_b[6]
        print("sum rate: %.2e" % rate_s)


        #写入变量的准备
        o_H_write=[format(x,'.2e') for x in o_H]
        K_dis_H2_write=format(K_dis_H2*P_H2,'.2e')
        K_H2_disorption_write=format(K_H2_disorption,'.2e')
        o_dis_H2_write=format(o_dis_H2,'.2e')
        K_write=[format(x,'.2e') for x in K]
        K_H2_write=[format(x,'.2e') for x in K_H2]
        K_re_write=[format(x,'.2e') for x in K_re]
        K_CO2_write=format(K_CO2*P_CO2,'.2e')
        K_CO2_re_write=format(K_CO2_re,'.2e')
        K_CO2_b_write = format(K_CO2_b*P_CO2, '.2e')
        K_CO2_re_b_write = format(K_CO2_re_b, '.2e')

        o_H_b_write = [format(x, '.2e') for x in o_H_b]
        o_dis_H2_b_write = format(o_dis_H2_belta, '.2e')
        K_b_write = [format(x, '.2e') for x in K_b]
        K_H2_b_write = [format(x, '.2e') for x in K_H2_b]
        K_re_b_write = [format(x, '.2e') for x in K_re_b]

        o_H=np.array(o_H)
        C_H=np.sum(o_H)-o_H[-1]-o_H[0]
        Pt_rate= 5.71e-10

        # 写入一个.dat文件
        wp.write("Element: %s\n" %name[n])
        wp.write("K_dis_H2:\n ")
        wp.write(str(K_dis_H2_write) + '\n')
        wp.write("K_H2_desorption:\n ")
        wp.write(str(K_H2_disorption_write) + '\n')
        wp.write("Path: a-a-a-a-a-a\n")
        wp.write("K_CO2:\n ")
        wp.write(str(K_CO2_write) + '\n')
        wp.write("K_CO2_re:\n ")
        wp.write(str(K_CO2_re_write) + '\n')
        wp.write("K_H:\n ")
        wp.write(str(K_H2_write)+'\n')
        wp.write("o_dis_H2:\n ")
        wp.write(str(o_dis_H2_write) + '\n')
        wp.write("o_H:\n ")
        wp.write(str(o_H_write)+'\n')
        wp.write("K:\n ")
        wp.write(str(K_write)+'\n')
        wp.write("K_re:\n ")
        wp.write(str(K_re_write) + '\n\n')
        wp.write("Path: b-b-b-b-b-b\n")
        wp.write("K_CO2:\n ")
        wp.write(str(K_CO2_b_write) + '\n')
        wp.write("K_CO2_re:\n ")
        wp.write(str(K_CO2_re_b_write) + '\n')
        wp.write("K_H:\n ")
        wp.write(str(K_H2_b_write) + '\n')
        wp.write("o_dis_H2:\n ")
        wp.write(str(o_dis_H2_b_write) + '\n')
        wp.write("o_H:\n ")
        wp.write(str(o_H_b_write) + '\n')
        wp.write("K:\n ")
        wp.write(str(K_b_write) + '\n')
        wp.write("K_re:\n ")
        wp.write(str(K_re_b_write) + '\n\n')
        wp.write("alpha path\n")
        wp.write(
            "*: %.2e\t*COOH: %.2e\t*COH2: %.2e\t*CHOH2: %.2e\t*CHOH: %.2e\t*CH2OH: %.2e\n"
            % (solution_eq1[o_1],solution_eq1[o_COOH],solution_eq1[o_COH2],solution_eq1[o_CHOH2], solution_eq1[o_CHOH],solution_eq1[o_CH2OH]))
        wp.write("Rate: %.2e\n" % rate)
        wp.write("Relative Rate: %.2e\n\n" % (rate/Pt_rate))
        wp.write("belta path\n")
        wp.write(
            "*: %.2e\t*COOH: %.2e\t*COH2: %.2e\t*CHOH2: %.2e\t*CHOH: %.2e\t*CH2OH: %.2e\n"
            % (solution_eq2[o_1_b],solution_eq2[o_COOH_b],solution_eq2[o_COH2_b],solution_eq2[o_CHOH2_b], solution_eq2[o_CHOH_b],solution_eq2[o_CH2OH_b]))
        wp.write("Rate: %.2e\n\n" % rate_b)
        wp.write("sum alpha and belta path\n")
        wp.write(
            "*: %.2e\t*COOH: %.2e\t*COH2: %.2e\t*CHOH2: %.2e\t*CHOH: %.2e\t*CH2OH: %.2e\n"
            % (solution_eq3[o_1_s], solution_eq3[o_COOH_s], solution_eq3[o_COH2_s], solution_eq3[o_CHOH2_s],
               solution_eq3[o_CHOH_s], solution_eq3[o_CH2OH_s]))
        wp.write("Rate: %.2e\n\n" % rate_s)









#%%
#***********************************************************#
#####작성자 :  임현선 (한국) 
#####동기   :연구에 관련된 계산 프로그램 작성을 목적으로 시작
#####목적   :기둥부재의 파괴모드 결정
#####        GENERATE EVELOPE CURVE (SHORT, INTERMEDIATE, TYPICAL)
#####모멘트 산정 프로그램
#####데이터 관리 프로그램
#####힘-변위 관계도 작성 프로그램
#***********************************************************#

print("Envelope curve of Column")

import pandas as pd
import math
import numpy as np 
import matplotlib.pyplot as plt 

#Read Data***************************************************
df=pd.read_excel('test.xlsx',sheet_name='list',header=1)

df_ne=df.loc[:,['Test ID','Reference','Title','Specimen Name','First Author','Type of Deformation','h [mm]','b [mm]','Cc [mm]','lc [mm]','Long. bar (1st face)','Long. bar (2nd face)','ρL','Trans. reinf. legs perp. to load','s [mm]','ρt','ρv','fc [MPa]','Long. bar dia. (1st face) [mm]','fy (long.) [MPa]','Trans.bar dia. [mm]','fy (trans.) [MPa]','P [kN]','Axial load ratio']]
print(df_ne)
#num=input('Enter the test ID number : ')
#print(df_ne.loc[df_ne['Test ID']==int(num)])

#def Test_inf(Title, Writer ):
 #  print("The " + Title + " written by " + Writer )
#Test_inf(df_ne["Title"], df_ne["First Author"])

r_asp=df_ne["lc [mm]"]/df_ne["h [mm]"]
#print(r_asp)
count_t, count_i, count_s = 0, 0, 0
index  = 0 
for r in r_asp:
   if r > 4.0:
      print(str(index) + " " + str(math.floor(r)) + " Typical Column")
      count_t = count_t + 1
   elif r>2 and r<=4:
      print(str(index) + " " + str(math.floor(r)) + " Intermediate Column")
      count_i = count_i + 1
   elif r<=2:
      print(str(index) + " " + str(math.floor(r)) + " Short Column")
      count_s = count_s + 1
   index = index + 1
print("Num of t : " + str(count_t) + " Num of i : " + str(count_i) + " Num of s : " + str(count_s) )

#num=input('Enter the test ID number : ')
#print(df_ne.loc[df_ne['Test ID']==int(num)])

print("***********Section Analysis***********")

#Section Analysis*************************************************
Area=df_ne["h [mm]"]*df_ne["b [mm]"]
#print(str(Area) + " mm^2 ")

d_1=df_ne["h [mm]"]-df_ne["Cc [mm]"]-df_ne["Trans.bar dia. [mm]"]-0.5*df_ne["Long. bar dia. (1st face) [mm]"]
#print(round(d_1))
kd=(0.25+0.85*df_ne["P [kN]"]*1000/Area/df_ne["fc [MPa]"])*df_ne["h [mm]"]
#print(round(kd))
M_y=df_ne["Long. bar (1st face)"]*math.pi*((0.5*df_ne["Long. bar dia. (1st face) [mm]"])**2)*df_ne["fy (long.) [MPa]"]*(d_1-kd/3)+df_ne["P [kN]"]*1000*(df_ne["h [mm]"]*0.5-kd/3)
#print(M_y)
V_y=[]
for i in range(len(r_asp)):
   if df_ne['Type of Deformation'][i] == 1:
    V_y.append(M_y[i]/df_ne['lc [mm]'][i])
   if df_ne['Type of Deformation'][i] == 2:  
    V_y.append(2*M_y[i]/df_ne['lc [mm]'][i])



#######FAILURE MODE**********************************************************

#######Calculate Mn**********************************************************

##R SECTION******************************************************************
#1.각 철근의 위치 정하기. Loc_list
#2.각 열에 위치한 철근의 면적 산정하기 As_list
#3.각 열에 위치한 철근의 탄성계수 E_list
#4.각 열의 변형률 epsilon 산정 eps_list
#5.각 열의 인장력 산정 (T=E*A*eps)
#6.P, Cc, T 값의 합이 0이 되도록 설정
#7.0이되는 c값을 찾은 후 Mn 계산값 출력
##***************************************************************************
print("V_Mn:")

Y=df_ne["Long. bar (1st face)"]
#print(Y)
N=df_ne["Long. bar (2nd face)"] ## Number of bar
#print(N)

#******************Determine Location of Rebar*******************************

##1~N개
a=np.array(df_ne["Cc [mm]"])+np.array(df_ne["Trans.bar dia. [mm]"])+np.array(df_ne["Long. bar dia. (1st face) [mm]"]/2)
b=(df_ne["h [mm]"]-2*a)/(N-1)

####1번의 위치 :  a
####2번의 위치 :  a + (H-2*a/(n-1))
####3번의 위치 :  a + (H-2*a/(n-1)) + (H-2*a/(n-1))
####4번의 위치 :  a + (H-2*a/(n-1)) + (H-2*a/(n-1)) + (H-2*a/(n-1)) 
####~~~~~N번의 위치
####도달하면 멈춤

Loc_list = []
for val1, val2, count in zip(a, b, N):
   small = [val1]
   for rep in range(1, count):
      small.append(val1 + val2 * rep)
   Loc_list.append(small)

#print(Loc_list)


####Determine A for steel****************************************************

Num = [n for n in df_ne['Long. bar (1st face)']]
Dia = [d for d in df_ne['Long. bar dia. (1st face) [mm]']]
N_list = [k for k in N]

Total_As=[]
def calc_A(N, r):
   return N * math.pi * r ** 2/4 
for i in range(len(N_list)):
   As_list = []
   for val in range(1, N_list[i]+1):
      if val == 1:
         As_list.append(calc_A(Num[i], Dia[i]))
      elif val == N_list[i]:
         As_list.append(calc_A(Num[i], Dia[i]))
      else:
         As_list.append(calc_A(2, Dia[i]))
   #print(As_list)
   Total_As.append(As_list)




####Determine E for steel****************************************************
E_s=210000
c_list=[1 for i in range(len(r_asp))]
E_slist=[E_s for i in range(len(Loc_list))]
fck = [a for a in df_ne['fc [MPa]']]
h = [a for a in df_ne['h [mm]']]
b=[b for b in df_ne['b [mm]']]
T=[]
Cc=[]
c_value=[]
half=0.5


beta=[]
for i in range(len(fck)):
   if fck[i] <= 27:
      beta.append(0.85)
   elif fck[i] > 27 and fck[i] <= 34:
      beta.append(0.8)
   elif fck[i] > 34 and fck[i] < 41:
      beta.append(0.75)

P = list(df_ne["P [kN]"]*1000)

Cc_all=[]
T_all=[]
realc_list=[]

for i in range(len(fck)):
   eps_s = []
   current_sum = []
   newEs_list=[E_s]*N[i]
   eps_y_upper=df_ne["fy (long.) [MPa]"][i]/210000
   eps_y_lower=-1*df_ne["fy (long.) [MPa]"][i]/210000
   while c_list[i] <= h[i]:
      newc_list = [] 
      for j in range(N[i]):
         newc_list.append(c_list[i])

      eps_s = list(0.003*(np.array(Loc_list[i]) - np.array(newc_list))/np.array(newc_list))
      for k in range(len(eps_s)):
         if eps_s[k] >= eps_y_upper:
            eps_s[k] = eps_y_upper
         elif eps_s[k] <= eps_y_lower:
            eps_s[k] = eps_y_lower

      T=np.array(eps_s)*np.array(newEs_list)*np.array(Total_As[i])      
      T_sum = np.sum(np.array(T))

      Cc= fck[i]*c_list[i]*0.85*beta[i]*b[i]
      if P[i]-Cc+T_sum <= 0.000001:
         realc_list.append(c_list[i])
         Cc_all.append(Cc)
         T_all.append(T_sum)
         break

      c_list[i]=c_list[i]+0.001  

Mn_list=[]
Final_T=[]
New_T=[]
epsilon=[]
for i in range(len(fck)):
   fakec_list = [realc_list[i]]*N[i]
   d1_list=[d_1[i]]*N[i]
   h_list=[h[i]]*N[i]
   nEs_list=[E_s]*N[i]
   neweps_s=[]
   eps_y_upper=df_ne["fy (long.) [MPa]"][i]/210000
   eps_y_lower=-1*df_ne["fy (long.) [MPa]"][i]/210000
   neps_s = list(0.003*(np.array(Loc_list[i]) - np.array(fakec_list))/np.array(fakec_list))
   for k in range(len(neps_s)):
      if neps_s[k] >= eps_y_upper:
         neps_s[k] = eps_y_upper
      elif neps_s[k] <= eps_y_lower:
         neps_s[k] = eps_y_lower
      neweps_s.append(neps_s[k])   
    
   epsilon.append(neweps_s)
   NewT=np.array(neweps_s)*np.array(nEs_list)*np.array(Total_As[i])
   New_T.append(NewT)
   
   difference=[]
   zip_object = zip(Loc_list[i], h_list)
   for list1_i, list2_i in zip_object:
     difference.append(list1_i-0.5*list2_i)

   NewT_sum=NewT*difference
   FinalT_sum=np.sum(np.array(NewT_sum))
   Final_T.append(FinalT_sum)
   
   NewCc= fck[i]*realc_list[i]*0.85*beta[i]*b[i]
   a=realc_list[i]*beta[i]/2
   difference2=h[i]/2-a

   FinalC=NewCc*difference2

   Mn=FinalC+FinalT_sum
   #print(round(Mn/1000000))
   Mn_list.append(Mn)

#**********V_mn(Shear Strength)**********

l_c=[l for l in df_ne['lc [mm]']]

hVmn_list = []
for val1, val2 in zip(Mn_list, l_c):
   hVmn_list.append(val1 / val2)

V_mn=[]
for i in range(len(hVmn_list)):
  if df_ne['Type of Deformation'][i] == 1:
    V_mn.append(hVmn_list[i]/1000)
  if df_ne['Type of Deformation'][i] == 2:  
    V_mn.append(2*hVmn_list[i]/1000)
print(np.round(V_mn,0))





#****************************************************************
#********SHEARSTRENGTH*******************************************
#****************************************************************

print("Cracking point : Vcr")
##Cracking point##

V_cr_w=0.27*((df_ne["fc [MPa]"])**(1/2))*df_ne["b [mm]"]*d_1+df_ne["P [kN]"]*1000*d_1/4/df_ne["h [mm]"]
V_cr_f=df_ne["b [mm]"]*d_1*(0.05*((df_ne["fc [MPa]"])**(1/2))+df_ne["h [mm]"]*(0.1*((df_ne["fc [MPa]"])**(1/2))+0.2*df_ne["P [kN]"]*1000/df_ne["b [mm]"]/df_ne["h [mm]"])/((M_y/V_y)-df_ne["h [mm]"]/2))

V_cr=[]
for i in range(len(r_asp)):
 if M_y[i]/V_y[i]-h[i]/2 <=0:
     V_cr.append(round(V_cr_w[i]/1000))
 else:
     V_cr.append(round(min(V_cr_f[i],V_cr_w[i])/1000))
#V_cr = list(V_cr_f)
#for i in range(len(V_cr_w)):
   #V_cr[i] = round(min(V_cr_f[i],V_cr_w[i])/1000)

print(V_cr)


##Strength point##

A_u=12*df_ne["ρt"]*df_ne["fy (trans.) [MPa]"]/df_ne["fc [MPa]"]
B_o=30*df_ne["ρt"]*df_ne["fy (trans.) [MPa]"]/df_ne["fc [MPa]"]

for i in range(len(A_u)):
   if A_u[i] > 1:
      A_u[i] = 1

for i in range(len(B_o)):
   if B_o[i] > 1:
      B_o[i] = 1

theta=list(r_asp)
K_o=list(r_asp)
for i in range(len(r_asp)):
   if r_asp[i]>2:
      theta[i]=math.radians(65) 
      K_o[i]=(min(math.tan(theta[i])**A_u[i]+(1/math.tan(theta[i]))**A_u[i]-1+0.14*B_o[i],1.64))
   elif r_asp[i]<=4 and r_asp[i]>2:
      theta[i]=math.radians(65) 
      K_o[i]=(min(math.tan(theta[i])**A_u[i]+(1/math.tan(theta[i]))**A_u[i]-1+0.14*B_o[i],1.64))
   else:
      theta[i]=np.arctan(df_ne["lc [mm]"][i]/(df_ne["h [mm]"][i]-2*kd[i]/3))
      K_o[i]=(min(math.tan(theta[i])**A_u[i]+(1/math.tan(theta[i]))**A_u[i]-1+0.14*B_o[i],1.64))


Soften=3.35/df_ne["fc [MPa]"]**(1/2)

for i in range(len(Soften)):
   if Soften[i] > 0.52:
      Soften[i] = 0.52

print("Strength point : Vn")

A_str=df_ne["b [mm]"]*kd


M_vd=[]

for i in range(len(Soften)):
   a=M_y[i]/V_y[i]/d_1[i]
   if a < 2:
      M_vd.append(2)
   elif a > 4:
      M_vd.append(4)
   elif a >=2 and a <=4:
      M_vd.append(a)

V_nc=K_o*Soften*df_ne["fc [MPa]"]*A_str*np.cos(np.array(theta))
V_nt=df_ne["ρt"]*df_ne["b [mm]"]*df_ne["s [mm]"]*df_ne["fy (trans.) [MPa]"]*d_1/df_ne["s [mm]"]+(0.5*((df_ne["fc [MPa]"])**(1/2)))*((1+1000*df_ne["P [kN]"]/(0.5*((df_ne["fc [MPa]"])**(1/2)))//df_ne["b [mm]"]/df_ne["h [mm]"])**(1/2))*0.8*df_ne["b [mm]"]*df_ne["h [mm]"]/M_vd



V_n = []
d_theta=[math.degrees(i) for i in theta]
for i in range(len(V_nc)):
    if r_asp[i] > 2:
        if V_nc[i] > V_nt[i]:
            V_n.append(np.round(V_nt[i]/1000))
        if V_nc[i] <= V_nt[i]:
            n_theta=[]
            n_theta.append(d_theta[i])
            nV_nc=K_o[i]*Soften[i]*df_ne["fc [MPa]"][i]*A_str[i]*np.cos(math.radians(n_theta[0]))
            if nV_nc/V_nt[i] >=0.00000001 and d_theta[i] > 45:               
               V_n.append(np.round(V_nt[i]/1000))
            if nV_nc/V_nt[i] <0.00000001 and n_theta[0] <= 45 :
               V_n.append(np.round(nV_nc/1000))                     
            d_theta[i]=d_theta[i]-0.001                       
    if r_asp[i] <= 2:
        V_n.append(round(V_nc[i]/1000))

print(V_n)


#Failure Mode
print("***********  Failure Mode  ***********")

for i in range(len(V_nc)):
 if V_mn[i]/V_n[i] > 1:
    if r_asp[i] > 4:
       print(df_ne['Specimen Name'][i] + ": Shear Failure " + "(Typical Column)")
    elif r_asp[i] <= 4 and r_asp[i] > 2:
       print(df_ne['Specimen Name'][i] + ": Shear Failure " + "(Intermediate Column)")
    elif r_asp[i] <= 4 and r_asp[i] <= 2:
       print(df_ne['Specimen Name'][i] + ": Shear Failure " + "(Short Column)")

 elif V_mn[i]/V_n[i] > 0.6 and V_mn[i]/V_n[i] <= 1:
    if r_asp[i] > 4:
         print(df_ne['Specimen Name'][i] + ": Flexural-Shear Failure " + "(Typical Column)")
    elif r_asp[i] <= 4 and r_asp[i] > 2:
         print(df_ne['Specimen Name'][i] + ": Flexural-Shear Failure " + "(Intermediate Column)")
    elif r_asp[i] <= 2:
         print(df_ne['Specimen Name'][i] + ": Flexural-Shear Failure " + "(Short Column)")
 elif V_mn[i]/V_n[i] < 0.6:
    if r_asp[i] > 4:
         print(df_ne['Specimen Name'][i] + ": Flexural Failure " + "(Typical Column)")
    elif r_asp[i] <= 4 and r_asp[i] > 2:
         print(df_ne['Specimen Name'][i] + ": Flexural Failure " + "(Intermediate Column)")
    elif r_asp[i] <= 2:
         print(df_ne['Specimen Name'][i] + ": Flexural Failure " + "(Short Column)") 

#****************************************************************
#********DISPLACEMENT********************************************
#****************************************************************
E_c=4700*(df_ne["fc [MPa]"])**(1/2)

I_g=df_ne["b [mm]"]*(df_ne["h [mm]"]**3)/12

L_r=df_ne["P [kN]"]*1000/df_ne["fc [MPa]"]/Area


#==============Graph for each failure mode==============

Ast=df_ne['Trans. reinf. legs perp. to load']*math.pi*(0.5*df_ne["Trans.bar dia. [mm]"])**2
r=[]
Ar=[i for i in df_ne['Axial load ratio']]
Roh=[i for i in df_ne['ρt']]
EI_eff = []
L_D = []
L_B = []
r = []

EI_eff_short=[]
for i in range(len(V_nc)):
   if V_mn[i]/V_n[i] > 1 and r_asp[i] > 2:
      I_eff_lower=0.3*I_g[i]*E_c[i]
      I_eff_upper=0.7*I_g[i]*E_c[i]
      if L_r[i] <= 0.1:
         EI_eff.append(I_eff_lower)
         EI_eff_short.append(0)
      elif L_r[i] >= 0.5:
         EI_eff.append(I_eff_upper)
         EI_eff_short.append(0)
      elif L_r[i] > 0.1 and L_r[i] < 0.5:
         EI_eff.append((I_eff_upper-I_eff_lower)*(L_r[i]-0.1)/(0.5-0.1)+I_eff_lower)
         EI_eff_short.append(0)

   if V_mn[i]/V_n[i] > 1 and r_asp[i] <= 2:
      I_eff_lower=0.7*I_g[i]*E_c[i]
      I_eff_upper=1.0*I_g[i]*E_c[i]
      I_eff_lower_short=0.5*I_g[i]*E_c[i]
      I_eff_upper_short=0.7*I_g[i]*E_c[i]
      if L_r[i] <= 0.1:
         EI_eff.append(I_eff_lower)
         EI_eff_short.append(I_eff_lower_short)
      elif L_r[i] >= 0.5:
         EI_eff.append(I_eff_upper)
         EI_eff_short.append(I_eff_upper_short)
      elif L_r[i] > 0.1 and L_r[i] < 0.5:
         EI_eff.append((I_eff_upper-I_eff_lower)*(L_r[i]-0.1)/(0.5-0.1)+I_eff_lower)
         EI_eff_short.append((I_eff_upper_short-I_eff_lower_short)*(L_r[i]-0.1)/(0.5-0.1)+I_eff_lower_short)
   
   if V_mn[i]/V_n[i] <= 1:
      I_eff_lower=0.3*I_g[i]*E_c[i]
      I_eff_upper=0.7*I_g[i]*E_c[i]
      if L_r[i] <= 0.1:
         EI_eff.append(I_eff_lower)
         EI_eff_short.append(0)
      elif L_r[i] >= 0.5:
         EI_eff.append(I_eff_upper)
         EI_eff_short.append(0)
      elif L_r[i] > 0.1 and L_r[i] < 0.5:
         EI_eff.append((I_eff_upper-I_eff_lower)*(L_r[i]-0.1)/(0.5-0.1)+I_eff_lower)
         EI_eff_short.append(0)
      
for i in range(len(V_nc)):
   if V_mn[i]/V_n[i] > 1:
      if Ar[i] <= 0.1 and Roh[i] >= 0.006:
       r.append(0.06)
      elif Ar[i] >= 0.6 and Roh[i] >= 0.006:
       r.append(0.008)
      elif Ar[i] <= 0.1 and Roh[i] <= 0.0005:
       r.append(0.006)
      elif Ar[i] >= 0.6 and Roh[i] <= 0.0005:
       r.append(0)
      elif Ar[i] > 0.1 and Ar[i] < 0.6 and Roh[i] >= 0.006:
       r.append((0.06-0.008)*(Ar[i]-0.6)/(0.1-0.6)+0.008)
      elif Ar[i] > 0.1 and Ar[i] < 0.6 and Roh[i] <= 0.0005:
       r.append(0.006*(0.6-Ar[i])/(0.6-0.1))
      elif Ar[i] <=0.1 and  Roh[i] < 0.006 and Roh[i] > 0.0005:
       r.append((0.06-0.006)*(Roh[i]-0.0005)/(0.006-0.0005)+0.006)
      elif Ar[i] >= 0.6 and Roh[i] < 0.006 and Roh[i] >0.0005:
       r.append((0.008-0)*(Roh[i]-0.0005)/(0.006-0.0005)+0)
      elif Ar[i] > 0.1 and Ar[i] < 0.6 and Roh[i] < 0.006 and Roh[i] > 0.0005:
       r.append((0.006*(0.6-Ar[i])*(0.006-Roh[i])+0*(Ar[i]-0.1)*(0.006-Roh[i])+0.06*(0.6-Ar[i])*(Roh[i]-0.0005)+0.008*(Ar[i]-0.1)*(Roh[i]-0.0005))/(0.6-0.1)/(0.006-0.0005))
   
   if V_mn[i]/V_n[i] <= 1 and V_mn[i]/V_n[i] > 0.6:
       if V_y[i]/b[i]/d_1[i]/math.sqrt(fck[i]) <=0.25:
           if Ar[i] <= 0.1 and Roh[i] >= 0.006:
                r.append(0.06)
           elif Ar[i] >= 0.6 and Roh[i] >= 0.006:
                r.append(0.01)
           elif Ar[i] <= 0.1 and Roh[i] <= 0.0005:
                r.append(0.012)
           elif Ar[i] >= 0.6 and Roh[i] <= 0.0005:
                r.append(0.004)
           elif Ar[i] > 0.1 and Ar[i] < 0.6 and Roh[i] >= 0.006:
                r.append((0.06-0.01)*(Ar[i]-0.6)/(0.1-0.6)+0.01)
           elif Ar[i] > 0.1 and Ar[i] < 0.6 and Roh[i] <= 0.0005:
                r.append((0.012-0.004)*(0.6-Ar[i])/(0.6-0.1)+0.004)
           elif Ar[i] <=0.1 and  Roh[i] < 0.006 and Roh[i] > 0.0005:
                r.append((0.06-0.012)*(Roh[i]-0.0005)/(0.006-0.0005)+0.012)
           elif Ar[i] >= 0.6 and Roh[i] < 0.006 and Roh[i] >0.0005:
                r.append((0.01-0.004)*(Roh[i]-0.0005)/(0.006-0.0005)+0.004)
           elif Ar[i] > 0.1 and Ar[i] < 0.6 and Roh[i] < 0.006 and Roh[i] > 0.0005:
                r.append((0.012*(0.6-Ar[i])*(0.006-Roh[i])+0.004*(Ar[i]-0.1)*(0.006-Roh[i])+0.06*(0.6-Ar[i])*(Roh[i]-0.0005)+0.01*(Ar[i]-0.1)*(Roh[i]-0.0005))/(0.6-0.1)/(0.006-0.0005))
       if V_y[i]/b[i]/d_1[i]/math.sqrt(fck[i]) >=0.5:
           if Ar[i] <= 0.1 and Roh[i] >= 0.006:
                r.append(0.06)
           elif Ar[i] >= 0.6 and Roh[i] >= 0.006:
                r.append(0.008)
           elif Ar[i] <= 0.1 and Roh[i] <= 0.0005:
                r.append(0.006)
           elif Ar[i] >= 0.6 and Roh[i] <= 0.0005:
                r.append(0)
           elif Ar[i] > 0.1 and Ar[i] < 0.6 and Roh[i] >= 0.006:
                r.append((0.06-0.008)*(Ar[i]-0.6)/(0.1-0.6)+0.008)
           elif Ar[i] > 0.1 and Ar[i] < 0.6 and Roh[i] <= 0.0005:
                r.append(0.006*(0.6-Ar[i])/(0.6-0.1))
           elif Ar[i] <=0.1 and  Roh[i] < 0.006 and Roh[i] > 0.0005:
                r.append((0.06-0.006)*(Roh[i]-0.0005)/(0.006-0.0005)+0.006)
           elif Ar[i] >= 0.6 and Roh[i] < 0.006 and Roh[i] >0.0005:
                r.append((0.008-0)*(Roh[i]-0.0005)/(0.006-0.0005)+0)
           elif Ar[i] > 0.1 and Ar[i] < 0.6 and Roh[i] < 0.006 and Roh[i] > 0.0005:
                r.append((0.006*(0.6-Ar[i])*(0.006-Roh[i])+0*(Ar[i]-0.1)*(0.006-Roh[i])+0.06*(0.6-Ar[i])*(Roh[i]-0.0005)+0.008*(Ar[i]-0.1)*(Roh[i]-0.0005))/(0.6-0.1)/(0.006-0.0005))

   if V_mn[i]/V_n[i] <= 0.6:
       r.append(0)



for i in range(len(V_nc)):
   if V_mn[i]/V_n[i] > 0.6 and r_asp[i] > 2:
      L_D.append(d_1[i]*math.tan(theta[i]))
      L_B.append(df_ne['lc [mm]'][i]-2*L_D[i])
   else:
      L_D.append(0)
      L_B.append(0)

d_yeild=[]
d_Vmn=[]
d_ns=[]
d_Vs=[]
d_Vc=[]
d_c=[]

for i in range(len(V_nc)):
   if V_mn[i]/V_n[i] > 1 and r_asp[i] > 4:     
      d_ns.append((2*0.006*math.sin(2*theta[i])*L_D[i])+(V_n[i]*L_B[i]/(0.4*E_c[i]*b[i]*d_1[i]))+(V_n[i]*df_ne['lc [mm]'][i]**3*1000/12/EI_eff[i]))
   elif V_mn[i]/V_n[i] > 1 and r_asp[i] <=4 and r_asp[i] >2:                 
      d_ns.append((0.006*math.sin(2*theta[i])*df_ne['lc [mm]'][i])+(V_n[i]*df_ne['lc [mm]'][i]**3*1000/12/EI_eff[i]))
   elif V_mn[i]/V_n[i] > 1 and r_asp[i] <=4 and r_asp[i] <=2:
      d_ns.append((0.006*math.sin(2*theta[i])*df_ne['lc [mm]'][i])+(V_n[i]*df_ne['lc [mm]'][i]**3*1000/12/EI_eff_short[i])) 
   else:
      d_ns.append(0)


      
for i in range(len(V_nc)):
   if V_mn[i]/V_n[i] > 1:
      d_Vc.append(d_ns[i]+r[i]*l_c[i])
   if V_mn[i]/V_n[i] <= 1:
      d_Vc.append(0)

d_cr=[]     

for i in range(len(V_nc)):     
   if V_mn[i]/V_n[i] > 1:
      d_cr.append(V_cr[i]*df_ne["lc [mm]"][i]**3*1000/12/EI_eff[i]+V_cr[i]*df_ne["lc [mm]"][i]*1000/0.4/E_c[i]/df_ne["b [mm]"][i]/d_1[i])
   if V_mn[i]/V_n[i] <= 1:
      d_cr.append(0)


for i in range(len(V_nc)):
   if V_mn[i]/V_n[i] > 1:
      x = [0,d_cr[i],d_ns[i],d_Vc[i]]
      y = [0,V_cr[i],V_n[i],0]

      plt.plot(x, y) 
  
      # naming the x axis 
      plt.xlabel('Displacement(mm)') 
      # naming the y axis 
      plt.ylabel('Shear Strength(kN)') 
  
      # giving a title to my graph 
      plt.title(df_ne['Specimen Name'][i]) 
  
      # function to show the plot 
      plt.show()

d_flc=[]
d_fssc=[]

for i in range(len(V_nc)): 
   if V_mn[i]/V_n[i] > 0.6 and V_mn[i]/V_n[i] <= 1:
     d_y=V_y[i]*l_c[i]**3/0.7/12/EI_eff[i]+V_y[i]*l_c[i]/0.4/E_c[i]/b[i]/d_1[i]
     d_mn=V_mn[i]*1000*l_c[i]**3/0.35/12/EI_eff[i]+V_mn[i]*1000*l_c[i]/0.4/E_c[i]/b[i]/d_1[i]
     d_s=(l_c[i]*(3/100+4*df_ne["ρt"][i]-(1/133)*(V_mn[i]/b[i]/d_1[i])/math.sqrt(fck[i])-df_ne["P [kN]"][i]*1000/40/h[i]/b[i]/fck[i]))
     #d_s=2*0.006*math.sin(2*theta[i])*L_D[i]*V_mn[i]/V_n[i]+V_mn[i]*1000*L_B[i]/0.4/E_c[i]/b[i]/d_1[i]+d_mn
     d_fsc=(l_c[i]*0.04*(1+math.tan(theta[i])**2)/((math.tan(theta[i])**2)+df_ne["P [kN]"][i]*1000*df_ne["s [mm]"][i]/Ast[i]/df_ne['fy (trans.) [MPa]'][i]/h[i]/math.tan(theta[i])))
     #d_fsc=d_s+r[i]*h[i]
     d_yeild.append(d_y)
     d_Vmn.append(d_mn)
     d_Vs.append(d_s)
     d_fssc.append(d_fsc)
     #if d_c[i] < d_s:
        #d_c[i]=d_s
     x = [0,d_y,d_mn,d_s,d_fsc]
     y = [0,V_y[i]/1000,V_mn[i],V_mn[i],0] 
     plt.plot(x, y) 
  
     # naming the x axis 
     plt.xlabel('Displacement(mm)') 
     # naming the y axis 
     plt.ylabel('Shear Strength(kN)') 
  
     # giving a title to my graph 
     plt.title(df_ne['Specimen Name'][i]) 
  
     # function to show the plot 
     plt.show() 

   if V_mn[i]/V_n[i] <= 0.6: 
     d_y=V_y[i]*l_c[i]**3/0.7/12/EI_eff[i]+V_y[i]*l_c[i]/0.4/E_c[i]/b[i]/d_1[i]
     d_mn=V_mn[i]*1000*l_c[i]**3/0.35/12/EI_eff[i]+V_mn[i]*1000*l_c[i]/0.4/E_c[i]/b[i]/d_1[i]
     d_fc=(l_c[i]*3.25*(1+40*df_ne["ρv"][i]*df_ne["fy (trans.) [MPa]"][i]*df_ne["Trans.bar dia. [mm]"][i]/fck[i]/h[i])*(1-Ar[i])*(1+l_c[i]/10/h[i]))/100
     d_yeild.append(d_y)
     d_Vmn.append(d_mn)
     d_flc.append(d_fc)
     x = [0,d_y,d_mn,d_fc]
     y = [0,V_y[i]/1000,V_mn[i],V_mn[i]] 
     plt.plot(x, y) 
  
     # naming the x axis 
     plt.xlabel('Displacement(mm)') 
     # naming the y axis 
     plt.ylabel('Shear Strength(kN)') 
  
     # giving a title to my graph 
     plt.title(df_ne['Specimen Name'][i]) 
  
     # function to show the plot 
     plt.show() 


#df = pd.DataFrame.from_records(V_cr, V_n, d_s,)




# %%

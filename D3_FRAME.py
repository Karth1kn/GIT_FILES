#INSTALL MAYAVI AND MATPLOTLIB USING PIP
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from mpl_toolkits import mplot3d
import scipy.linalg as slin
import math 
import csv
import time
import os
import sys
import scipy
from scipy.linalg import solveh_banded
from mayavi import mlab

script_dir = os.path.dirname(os.path.abspath(__file__))
coord_path = os.path.join(script_dir, 'coord.csv')
con_path = os.path.join(script_dir, 'connect.csv')
prop_path = os.path.join(script_dir, 'properties.csv')
bc_path = os.path.join(script_dir, 'fr_bound_con.csv')
force_path = os.path.join(script_dir, 'fr_force_list.csv')
option_path = os.path.join(script_dir, 'stress_options.csv')
disp_path = os.path.join(script_dir,'displacements.csv')
rot_path = os.path.join(script_dir,'rotates.csv')



def timer(func):
    def wrapper(*args, **kwargs):
        t1 = time.time()
        mat = func(*args, **kwargs)
        t2 = time.time()
        print(f'{func.__name__} finished in: {t2-t1}')
        return mat
    return wrapper

with open(coord_path, 'r') as file:
    try:
        val_list = []
        for i in csv.reader(file):
            val_list.append([float(j) for j in i])
        val_array= np.array(val_list)
        with open(con_path,'r') as con:
            con_list= []

            for c in csv.reader(con):
                con_list.append([int(b) for b in c])
            con_array= np.array(con_list)

            for k in con_array:
                x_ele = [val_array[k[0],0],val_array[k[1],0]]
                y_ele = [val_array[k[0],1],val_array[k[1],1]]
                z_ele = [val_array[k[0],2],val_array[k[1],2]]
                #ax.plot3D(x_ele, y_ele, z_ele, 'blue')

    except Exception:
        print('error')

cord=val_array
m= con_array

#### BOUNDARY CONDITIIONS
with open(bc_path, 'r') as bc:
        bc_lst = []
        c = 0
        for c1 in csv.reader(bc):
            pos = int(c1[0])*6
            if c1[1] == 'a':
                for r1 in range(6):
                    bc_lst.append(pos+r1)
                    c= c1[0]
                continue
            elif c1[1]=='x':
                if c1[0]==c:
                    continue
                bc_lst.append(pos)
            elif c1[1] == 'y':
                if c1[0]==c:
                    continue
                bc_lst.append(pos+1)
            elif c1[1] == 'z':
                if c1[0]==c:
                    continue
                bc_lst.append(pos+2)
            elif c1[1]=='rx':
                if c1[0]==c:
                    continue
                bc_lst.append(pos+3)
            elif c1[1] == 'ry':
                if c1[0]==c:
                    continue
                bc_lst.append(pos+4)
            elif c1[1] == 'rz':
                if c1[0]==c:
                    continue
                bc_lst.append(pos+5)


######## FORCES
with open(force_path, 'r') as force:
        force_list =[]
        lst = csv.reader(force)
        for row in lst:
            if row == []: continue
            for i in row[6:]:
                force_list.append([int(i),[float(j) for j in row[:6]]])

with open(prop_path,'r')as fl:
    props = []
    for i in csv.reader(fl):
        for j in i:
            props.append(float(j))
    plop = props[9]
    density = props[8]
    do = props[7]
    area =props[4]   #math.pi*(do**2-di**2)/4  #
    youngs_modulus = props[0] #5e11
    ks = props[6]
    Iy = props[1]  #im
    Iz = props[2]    #im
    shear_modulus = props[3]    #2e11
    polar_inertia = Iy+Iz
    yield_stress = props[5]     #1550


ln = lambda n1,n2: (( ( cord[(n2),0] - cord[(n1),0] )**2 + ( cord[(n2),1] - cord[(n1),1] )**2 + ( cord[(n2),2] - cord[(n1),2] )**2 )**0.5)

op = lambda node1,node2,i: ((cord[node2,i]-cord[node1,i])/ln(node1,node2))
cosx = lambda node1,node2: op(node1,node2,0)
cosy = lambda node1,node2: op(node1,node2,1)
cosz = lambda node1,node2: op(node1,node2,2)

np.set_printoptions(linewidth=np.inf)



with open(rot_path, 'r') as rot:
    jt = []
    for i in csv.reader(rot):
        jt.append(i)
if len(jt)!=0:
    elem = [int(k) for k in jt[0][1:-1]]
else : elem=[]

trans_mats = []
def lsm(element):
    al = 0
    if element[0] in elem and element[1] in elem:
        #print('yes')
        al = float(jt[0][0])*math.pi/180

    rotate = np.array([[1,0,0],
                    [0, np.cos(al), np.sin(al) ],
                    [0,-np.sin(al), np.cos(al)]])

    c1,c2,c3,c4,l = area*youngs_modulus, shear_modulus*polar_inertia, youngs_modulus*Iy, youngs_modulus*Iz,ln(*element)
    oy = (12*youngs_modulus*Iy)/(ks*area*shear_modulus*l)
    oz = (12*youngs_modulus*Iz)/(ks*area*shear_modulus*l)
    stiff = np.zeros((12,12),dtype = float)
    main_diag= 2*[c1/l , (12*c4/(l**3))/(1+oz) , (12*c3/(l**3))/(1+oy) , c2/l , (c3/l)*(4+oy)/(1+oy) , (c4/l)*(4+oz)/(1+oz)]
    main_off_diag = [-c1/l, (-12*c4/l**3)/(1+oz) , (-12*c3/l**3)/(1+oy) , -c2/l , (c3/l)*(2-oy)/(1+oy) , (c4/l)*(2-oz)/(1+oz)]
    minor_diag = [(6*c4/l**2)/(1+oz) , (-6*c3/l**2)/(1+oy) , -c2/l , (6*c3/l**2)/(1+oy) , (-6*c4/l**2)/(1+oz) , c1/l , (-6*c4/l**2)/(1+oz) , (6*c3/l**2)/(1+oy) , -c2/l , (-6*c3/l**2)/(1+oy) , (6*c4/l**2)/(1+oz)]
    minor_up = [-c1/l , (6*c4/l**2)/(1+oz) , (-6*c3/l**2)/(1+oy) , c2/l , (-6*c3/l**2)/(1+oy) , (6*c4/l**2)/(1+oz) , -c1/l]
    minor_down = [(-6*c4/l**2)/(1+oz) , (6*c3/l**2)/(1+oy) , c2/l , (6*c3/l**2)/(1+oy) , (-6*c4/l**2)/(1+oz)]


    def filler(lst,start_row, end_row,a,b):
        for i, val in enumerate(lst):
            stiff[start_row+a*i, end_row+b*i] = val #CSR STORAGE

    filler(main_off_diag,0,6,1,1)
    filler(main_off_diag,6,0,1,1)
    filler(minor_diag,11,1,-1,1)
    filler(minor_up,6,0,-1,1)
    filler(minor_down,11,7,-1,1)
    np.fill_diagonal(stiff,main_diag)


    li,mi,ni = cosx(*element), cosy(*element), cosz(*element)


    if mi==0 and ni==0 :
        if li<0:
            delta = np.array([[-1,0,0],
                            [0,-1,0],
                            [0,0,1]])
            
        elif li>0:
            delta = np.array([[1,0,0],
                            [0,1,0],
                            [0,0,1]])  
    else:
        di = (mi**2 + ni**2)**0.5
        delta = np.array([[li,        mi,        ni],      
                          [0 ,    -ni/di,     mi/di],
                          [di, -mi*li/di, -ni*li/di]])

        
        dlta = np.array([[li,mi,di],
                          [-mi/di,li/di,0],
                          [-li*ni/di,-mi*ni/di,di]])

    trans_mats.append(delta)
    trans = np.zeros((12,12))
    for i in range(4):
        trans[3*i:3*i+3 , 3*i:3*i+3] = rotate@delta
    
    return np.transpose(trans)@stiff@trans
    #print(local_stiff)

def global_stiff_mat(cord_list,dof):
    ran = int((dof/2))
    matrix= np.zeros([len(cord)*ran,len(cord)*ran])
    global e_len
    e_len= []
    for i in cord_list:
        k = lsm(i)
        t1 = np.array([[i[a],i[b]]for a in range(2) for b in range(2)])
        t2 = np.array([[r,s]for r in range(0,ran+1,ran) for s in range(0,ran+1,ran)])
        for j in range(4):
           matrix[6*t1[j,0]:6*t1[j,0]+6 , 6*t1[j,1]:6*t1[j,1]+6] += k[t2[j,0]:t2[j,0]+6 , t2[j,1]:t2[j,1]+6]
        e_len.append(ln(*i))
    np.set_printoptions(precision=2)
    return matrix

def force_matrix(node,force):
    f_matrix= np.zeros([6*len(cord),1])
    k = 6*node
    for i in range(6):
        f_matrix[k+i] += force[i]

    return f_matrix

force_mat=force_matrix(0,[0.0,0.0,0.0,0.0,0.0,0.0])

for ac in force_list:
    force_mat += force_matrix(*ac)

######SELF_WEIGHT
for i,wt in enumerate(m):
    weight = density*area*ln(*wt)/2
    force_mat[6*wt[0]+2] += -weight
    force_mat[6*wt[1]+2] += -weight

def boundary_conditions(matrix,bc_list):
    #matrix= np.delete(matrix,bc_list,0)
    matrix[bc_list , :] =0

    if len(matrix[0])!=1:
        #matrix= np.delete(matrix,bc_list,1)
        matrix[: , bc_list] =0
        for i in bc_list:
            matrix[i,i]=1
    np.set_printoptions(precision=2)
    return matrix

global_stiff_matrix= global_stiff_mat(m,12)

#APPLYING BOUNDARY CONDITIONS
bc_stiff_mat= boundary_conditions(global_stiff_matrix,bc_lst)
#print(len(bc_stiff_mat))
bc_force_matrix= boundary_conditions(force_mat,bc_lst)

def rnge():
    t5 = time.time()
    max_arr = [abs(i[0]-i[1]) for i in con_array ]
    t6 = time.time()
    return max(max_arr)*6+6

def sparser():
    oy = []
    for q in range(len(cord)):
        for w in range(len(cord)):
            if not np.array_equal(global_stiff_matrix[6*q:6*q+6 , 6*w:6*w+6] , np.zeros([6,6])):
                oy.append([w,q])

    print(f'Sparsity = {1- (len(oy)/len(cord)**2)}')
    plt.scatter(np.array(oy)[:,0] , np.array(oy)[:,1],s=1)
    jk = 0
    plt.xlim(0-jk,len(cord)+jk)
    plt.ylim(0-jk,len(cord)+jk)
    plt.gca().invert_yaxis()
    plt.show()

if __name__ == "__main__":
    if len(sys.argv) > 1:
        function_name = sys.argv[1]
        if function_name == "sparser":
            sparser()
t5 = time.time()
def main_sequence():
    global t3
    t3 = time.time()
    bandwidth = rnge()         
    diagonals = np.array([np.append(np.diag(bc_stiff_mat, -i), np.zeros([1,i])[0]) for i in range(bandwidth)  ])
    nodal_displacement = solveh_banded(diagonals, bc_force_matrix, lower= True)

    global t4
    t4 = time.time()
    # Uncomment for solution using LU decomposition
    #nodal_displacement = slin.solve(bc_stiff_mat,bc_force_matrix )
    global t5
    t5 = time.time()
    nodal_displacement = [item for sublist in nodal_displacement for item in sublist]
    with open(disp_path,'w') as disp:
        writer = csv.writer(disp, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        writer.writerow(nodal_displacement)

if __name__ == "__main__":
    # Check if the script is run as the main program

    if len(sys.argv) > 1:
        function_name = sys.argv[1]
        if function_name == "main_sequence":
            t1 = time.time()
            main_sequence()
            t2 = time.time()
            print(f'Solution time only Thomas: {t4-t3}')
            print(f'Solution time only LU: {t5-t4}')
            #print('Solution time plus write: ', t2-t1, '\n')
            #print(f'Total time: {(l2-l1) + (t2-t1)}', '\n')


nodal_disp = []
with open(disp_path, 'r') as disp1:
    reader = csv.reader(disp1)
    for row in reader:
        nodal_disp.extend(row)

nodal_disp = np.array([float(i) for i in nodal_disp])
nodal_disp1 = nodal_disp.copy()


shape_nodal_disp= nodal_disp1.reshape(len(cord),6)
for i,val in enumerate(cord):
    shape_nodal_disp[i,0:3]+=val
cord = shape_nodal_disp[:,0:3]
def_e_len= [ln(*i) for i in m]

def transfmat(element):
    li,mi,ni = cosx(*element), cosy(*element), cosz(*element)


    if mi==0 and ni==0 :
        if li<0:
            delta = np.array([[-1,0,0],
                            [0,-1,0],
                            [0,0,1]])
            
        elif li>0:
            delta = np.array([[1,0,0],
                            [0,1,0],
                            [0,0,1]])  
    else:
        if li>0:
            li= -li
            mi = -mi
            ni = -ni
        di = (mi**2 + ni**2)**0.5
        delta = np.array([[li,        mi,        ni],      
                            [0 ,    -ni/di,     mi/di],
                            [di, -mi*li/di, -ni*li/di]])

        
        dlta = np.array([[li,mi,di],
                            [-mi/di,li/di,0],
                            [-li*ni/di,-mi*ni/di,di]])

    trans_mats.append(delta)

loc_disp = []
for ind,val in enumerate(m):
    #transfmat(val)
    trans = np.zeros((12,12))
    o = np.array(trans_mats)[ind]
    for l in range(4):
        trans[3*l:3*l+3 , 3*l:3*l+3] = o
 
    t = np.insert(nodal_disp[6*val[0]:6*val[0]+6],6,nodal_disp[6*val[1]:6*val[1]+6])
    
    loc_disp.append(trans@t)

    np.set_printoptions(precision=2)
loc_disp = np.array(loc_disp)

axial_stress = []
zb_strs  =[]
yb_strs = []
x_tors = []
fac = 10**6
for a,ij in enumerate(loc_disp):

    lt = abs(ln(*m[a]))
    #AXIAL TENSION
    axial_stress.append(((ij[6] - ij[0])/lt)*youngs_modulus/fac)
    #BENDING STRESS ALONG z
    d1 = lambda a: ij[7]*((12*a-6)/lt**2) + ij[1]*((-12*a+6)/lt**2) + ij[5]*((6*a-4)/lt) + ij[11]*((6*a-2)/lt)
    avg1 = (d1(0.25)+d1(0.5)+d1(0.75))/3
    #av1 = (ij[11]-ij[5])/lt
    zb_strs.append(avg1*(do/2)*youngs_modulus/fac)
    #BENDING STRESS ALONG y
    d2 = lambda a: ij[8]*((12*a-6)/lt**2) + ij[2]*((-12*a+6)/lt**2) - ij[10]*((6*a-4)/lt) - ij[4]*((6*a-2)/lt)
    avg2 = (d2(0.25)+d2(0.5)+d2(0.75))/3
    #av2 = (ij[8] - ij[2])/lt
    yb_strs.append(avg2*(do/2)*youngs_modulus/fac)

    x_tors.append(((ij[9] - ij[3])*(do/2)/lt)*shear_modulus/fac)
def color_plotter(coordinate):
    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(111, projection= '3d')
    global t3
    t3 = time.time()
    with open(option_path,'r') as file:
        for i in csv.reader(file):
            if i[0] == 'tensile':
                value = axial_stress
            elif i[0] == 'bendingY':
                value = yb_strs
            elif i[0] == 'bendingZ':
                value = zb_strs
            elif i[0] == 'torsion':
                value = x_tors
            for j,k in enumerate(m):
                ks = np.array([coordinate[k[0]],coordinate[k[1]]])
                c = np.array([ks[0,0],ks[1,0]])
                d = np.array([ks[0,1],ks[1,1]])
                e = np.array([ks[0,2],ks[1,2]])

                fos= value[j]/yield_stress
                #print(fos)
                if fos**2>=1:
                    r,g,b,a = 0 ,1 ,0 ,1
                    #plt.annotate('material will fail at this load',(3000,-5000),color= 'white')
                    
                    if plop ==0:
                        ax.plot3D(c,d,e,color=(r,g,b,a))
                    else: 
                        line1 = mlab.plot3d(c,d,e,color=(r,g,b))
                elif fos<=0:
                    r,g,b,a = 1+fos, 1+fos, 1,  1
                    if plop ==0:
                        ax.plot3D(c,d,e,color=(r,g,b,a))
                    else: 
                        line1 = mlab.plot3d(c,d,e,color=(r,g,b))
                elif fos>0:
                    r,g,b,a= 1, 1-fos, 1-fos,  1         
                    if plop ==0:
                        ax.plot3D(c,d,e,color=(r,g,b,a))
                    else: 
                        line1 = mlab.plot3d(c,d,e,color=(r,g,b))


    ax.set_facecolor('black')
    ax.grid(visible=False)
    sc = 0
    x = 5
    ax.set_xlim(-x-sc,x+sc)
    ax.set_ylim(-x-sc,x+sc)
    ax.set_zlim(-x-sc,x+sc)
    global t4
    t4 = time.time()
    if plop ==0:
        plt.show()
    else: 
        mlab.show()
def ploter():
    color_plotter(shape_nodal_disp[:,:3])


if __name__ == "__main__":
    # Check if the script is run as the main program
    if len(sys.argv) > 1:
        function_name = sys.argv[1]
        if function_name == "ploter":
            ploter()
            print(f'Plotting time: {t4-t3}')


            

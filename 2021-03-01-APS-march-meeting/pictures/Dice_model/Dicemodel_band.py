#!/usr/bin/env python

# one dimensional chain

# Copyright under GNU General Public License 2010, 2012, 2016
# by Sinisa Coh and David Vanderbilt (see gpl-pythtb.txt)
from pythtb import * # import TB model class
import matplotlib.pyplot as plt
import matplotlib as mpl
# mpl.use('agg')
mpl.rcParams['figure.dpi']= 300
# font stuff
font = {'family' : 'DejaVu Sans',
        'weight' : 'light',
        'size'   : 18}

font1 = {'family' : 'DejaVu Sans',
        'weight' : 'light',
        'size'   : 14}

font2 = {'family' : 'DejaVu Sans',
        'weight' : 'light',
        'size'   : 8}

mpl.rc('font', **font)



# def cal_occ(mu,evals,evacs,temperature=0.001):
#     occ = np.zeros(3)
#     for iorb in range(3):
#         for ikpt in range(evals.shape[1]):
#             for ibnd in range(2):
#                 # occ[iorb] += 1./evals.shape[1] * np.abs(evacs[ibnd,ikpt,iorb])**2 * 1./(np.exp((evals[ibnd,ikpt]-mu)/temperature)+1.)
#                 # print np.conj(evacs[ibnd,ikpt,iorb])
#                 occ[iorb] += 1./evals.shape[1] * np.abs(evacs[ibnd,ikpt,iorb]**2)
#
#     return occ

#%%
# LaBr FM down
# t12=0.3
# onsite = [-2.677544,-1.963672,-1.963672]
# hoppings = [t01,t01,t12]#[-0.092875,-0.092875,0.045314]

# LaBr AFM down
# onsite = [-2.845791,-2.016203,-2.911754]
# hoppings = [t01,-0.16694,t12]#[-0.085837,-0.16694,,0.047878]

# LaBr AFM down test
onsite = [0,0,0]
hoppings = [-0.3, -0.3, 0]


lat=[[2.0,0.0],[-1.0,np.sqrt(3)]]
orb=[[0.0,0.0],[1./3.,2./3.],[2./3.,1./3.]]
at3_model=tb_model(2,2,lat,orb)
at3_model1=tb_model(2,2,lat,orb)
at3_model2=tb_model(2,2,lat,orb)

#--
at3_model.set_onsite(onsite)

at3_model.set_hop( hoppings[0], 0, 1, [0,0])
at3_model.set_hop( hoppings[0], 1, 0, [0,1])
at3_model.set_hop( hoppings[0], 1, 0, [1,1])

at3_model.set_hop( hoppings[1], 0, 2, [0,0])
at3_model.set_hop( hoppings[1], 2, 0, [1,0])
at3_model.set_hop( hoppings[1], 2, 0, [1,1])

at3_model.set_hop( hoppings[2], 1, 2, [0,0])
at3_model.set_hop( hoppings[2], 1, 2, [0,1])
at3_model.set_hop( hoppings[2], 1, 2, [-1,0])
#1--
at3_model1.set_onsite(onsite)

at3_model1.set_hop( hoppings[0], 0, 1, [0,0])
at3_model1.set_hop( hoppings[0], 1, 0, [0,1])
at3_model1.set_hop( hoppings[0], 1, 0, [1,1])

at3_model1.set_hop( 0.1, 0, 2, [0,0])
at3_model1.set_hop( 0.1, 2, 0, [1,0])
at3_model1.set_hop( 0.1, 2, 0, [1,1])

at3_model1.set_hop( hoppings[2], 1, 2, [0,0])
at3_model1.set_hop( hoppings[2], 1, 2, [0,1])
at3_model1.set_hop( hoppings[2], 1, 2, [-1,0])
#2--
at3_model2.set_onsite(onsite)

at3_model2.set_hop( hoppings[0], 0, 1, [0,0])
at3_model2.set_hop( hoppings[0], 1, 0, [0,1])
at3_model2.set_hop( hoppings[0], 1, 0, [1,1])

at3_model2.set_hop( hoppings[1], 0, 2, [0,0])
at3_model2.set_hop( hoppings[1], 2, 0, [1,0])
at3_model2.set_hop( hoppings[1], 2, 0, [1,1])

at3_model2.set_hop( -0.1, 1, 2, [0,0])
at3_model2.set_hop( -0.1, 1, 2, [0,1])
at3_model2.set_hop( -0.1, 1, 2, [-1,0])


plot = True

if plot == True:
    path=[[0.0,0.0],[0.5,0.0],[1./3.,1./3.],[0.0,0.0]]
    k_label=( r'$\Gamma $',r'$\mathrm{X}$',r'$\mathrm{M}$',r'$\Gamma $')
    (k_vec,k_dist,k_node)=at3_model.k_path(path,301,report=False)
else:
    mesh_size = [60,60]
    k_vec = at3_model.k_uniform_mesh(mesh_size)

# solve model
(evals,evacs)=at3_model.solve_all(k_vec,eig_vectors=True)
(evals1,evacs)=at3_model1.solve_all(k_vec,eig_vectors=True)
(evals2,evacs)=at3_model2.solve_all(k_vec,eig_vectors=True)

# evals=at3_model.solve_all(k_vec)

# plot band structure
if plot == True:
    fig, ax = plt.subplots(figsize=(3.6,4.4))
    #--
    ax.plot(k_dist,evals[0],color='#1f77b4',label='Dice')
    ax.plot(k_dist,evals[1],color='#1f77b4')
    ax.plot(k_dist,evals[2],color='#1f77b4')
    #--
    # ax.plot(k_dist,evals1[0],color='#ff7f0e',label=r'$\mathrm{\alpha-t_3}$')
    # ax.plot(k_dist,evals1[1],color='#ff7f0e')
    # ax.plot(k_dist,evals1[2],color='#ff7f0e')
    # #--
    # ax.plot(k_dist,evals2[0],color='#2ca02c',label='Triangle')
    # ax.plot(k_dist,evals2[1],color='#2ca02c')
    # ax.plot(k_dist,evals2[2],color='#2ca02c')

    # ax.set_title(str(0.16))
    # ax.set_xlabel("Path in k-space",**font)
    ax.set_ylabel("Energy [eV]",**font)
    ax.set_xticks(k_node)
    ax.set_xticklabels(k_label,**font)
    ax.set_xlim(k_node[0],k_node[-1])
    # ax.set_ylim(-4,-1.2)
    for n in range(len(k_node)):
      ax.axvline(x=k_node[n], linewidth=0.5, color='k')

    # ax.legend(fontsize='small')

    fig.tight_layout()
    # fig.show()
#%
# data = cal_occ(0,evals,evacs,temperature=0.001)
# # print '%01.2f %s' % (t12,data[1]*data[2])
# print '%s' % (data)

"""
aut_lib.py

This module contains various utility functions for performing transformations and automated FEMM simulations.

Functions:
    T_dq_abc(phys_in, ar):
        Transforms dq-axis currents to abc-axis currents.
    T_abc_dq(phys_in, ar):
        Transforms abc-axis currents to dq-axis currents.
    calc_FEMM_faulted(id, iq, i0, ang, pp):
        Performs FEMM simulation for given current inputs and rotor angle, and returns various electromagnetic properties for the faulted machine.
    calc_FEMM_healthy(id, iq, ang, pp):
        Performs FEMM simulation for given current inputs and rotor angle, and returns various electromagnetic properties for the healthy machine.

Dependencies:
    numpy
    femm
    scipy.io
    os
"""

import numpy                as np           
import femm
from scipy.io import savemat
import os

pi = np.pi

def T_dq_abc(phys_in,ar):
    """
    Transforms dq-axis currents to abc-axis currents.
    Parameters:
        phys_in (array-like): Input currents in dq-axis.
        ar (float): Rotor angle in radians.
    Returns:
        array: Transformed currents in abc-axis.
    """
    transformation = np.matmul(np.array([[np.cos(ar),-np.sin(ar),1],[np.cos(ar-2*pi/3),-np.sin(ar-2*pi/3),1] ,[np.cos(ar+2*pi/3),-np.sin(ar+2*pi/3),1]]),phys_in)
    return transformation

def T_abc_dq(phys_in,ar):
    """
    Transforms  abc-axis  currents to dq-axis currents.
    Parameters:
        phys_in (array-like): Input currents in abc-axis.
        ar (float): Rotor angle in radians.
    Returns:
        array: Transformed currents in dq-axis.
    """
    transformation = (2/3)*np.matmul(np.array([[np.cos(ar),np.cos(ar-2*pi/3),np.cos(ar+2*pi/3)],[-np.sin(ar),-np.sin(ar-2*pi/3),-np.sin(ar+2*pi/3)] ,[1/2,1/2,1/2]]),np.transpose(phys_in))
    return transformation

def calc_FEMM_faulted(id,iq,i0,ang,pp):
    """
    Performs FEMM simulation for given current inputs and rotor angle, and returns various electromagnetic properties.
    Parameters:
        id (float): Direct-axis current.
        iq (float): Quadrature-axis current.
        i0 (float): Zero-sequence current.
        ang (float): Rotor angle in degrees.
        pp (int): Number of pole pairs.
    Returns:
        tuple: Contains the following electromagnetic properties:
            - psi_d (float): Direct-axis flux linkage.
            - psi_q (float): Quadrature-axis flux linkage.
            - psi_f (float): Flux linkage of the faulted turn.
            - i_af (float): Adjusted field current.
            - b_rad (float): Radial component of the air gap flux density.
            - b_tan (float): Tangential component of the air gap flux density.
            - m_tot (float): Torque map).
            - psi_d2 (float): Direct-axis flux linkage with faulted turn.
            - psi_q2 (float): Quadrature-axis flux linkage with faulted turn.
    """
    femm.openfemm(True,femmpath='C:/femm42')
    theta_act_mech_deg = ang
    theta_act_el_rad = (theta_act_mech_deg*pp) * pi/180

    femm.opendocument('prius_ITSC.FEM')
    femm.mi_saveas(f'temp_{i0}.FEM')

    femm.mi_selectgroup(1)
    femm.mi_modifyboundprop('APairgap', 10, theta_act_mech_deg)


    i_dq    =   [id,iq,0]
    i_abc   =   T_dq_abc(i_dq,theta_act_el_rad)
    i_af    =   i_abc[2] - i0 

    femm.mi_modifycircprop('U',1,i_abc[0])
    femm.mi_modifycircprop('V',1,i_abc[1])
    femm.mi_modifycircprop('W',1,i_abc[2])
    femm.mi_modifycircprop('F',1,i_af)

    # starts femm analysis
    femm.mi_analyze(1)
    femm.mi_loadsolution()

    # get results of circuits and flux linkage
    prop_a              =       femm.mo_getcircuitproperties('U')
    prop_b              =       femm.mo_getcircuitproperties('V')
    prop_c              =       femm.mo_getcircuitproperties('W')
    prop_f              =       femm.mo_getcircuitproperties('F')
    gap_B               =       femm.mo_getgapb('APairgap',ang)
    gap_M               =       femm.mo_gapintegral("APairgap",0)


    psi_abc = [prop_a[2], prop_b[2], prop_c[2]]
    psi_abc_2 = [prop_a[2], prop_b[2], prop_c[2]+prop_f[2]]
    psi_dq = T_abc_dq(psi_abc,theta_act_el_rad)
    psi_dq2 = T_abc_dq(psi_abc_2,theta_act_el_rad)
    psi_d = psi_dq[0]
    psi_q = psi_dq[1]
    psi_d2 = psi_dq2[0]
    psi_q2 = psi_dq2[1]
    psi_f = prop_f[2]
    b_rad = gap_B[0]
    b_tan = gap_B[1]
    m_tot = gap_M

    femm.mi_selectgroup(1)

    femm.mi_close
    femm.closefemm
    os.remove(f'temp_{i0}.FEM')
    os.remove(f'temp_{i0}.ans')
    return psi_d,psi_q,psi_f,i_af,b_rad,b_tan,m_tot,psi_d2,psi_q2

def calc_FEMM_healthy(id,iq,ang,pp):
    """
    Performs FEMM simulation for given current inputs and rotor angle, and returns various electromagnetic properties.
    Parameters:
        id (float): Direct-axis current.
        iq (float): Quadrature-axis current.
        ang (float): Rotor angle in degrees.
        pp (int): Number of pole pairs.
    Returns:
        tuple: Contains the following electromagnetic properties:
            - psi_d (float): Direct-axis flux linkage.
            - psi_q (float): Quadrature-axis flux linkage.
            - b_rad (float): Radial component of the air gap flux density.
            - b_tan (float): Tangential component of the air gap flux density.
            - m_tot (float): Torque map.
    """
    femm.openfemm(True,femmpath='C:/femm42')
    theta_act_mech_deg = (-90/pp + ang)*0
    theta_act_el_rad = (theta_act_mech_deg*pp) * pi/180

    femm.opendocument('prius.FEM')
    femm.mi_saveas(f'temp_{iq}.FEM')

    femm.mi_selectgroup(1)
    femm.mi_moverotate(0,0,theta_act_mech_deg)
    femm.mi_modifyboundprop('APairgap', 10, theta_act_mech_deg)



    i_dq=[id,iq,0]
    i_abc=T_dq_abc(i_dq,theta_act_el_rad)
    femm.mi_modifycircprop('U',1,i_abc[0])
    femm.mi_modifycircprop('V',1,i_abc[1])
    femm.mi_modifycircprop('W',1,i_abc[2])

    # starts femm analysis
    femm.mi_analyze(1)
    femm.mi_loadsolution()

    # get results of circuits and flux linkage
    prop_a              =       femm.mo_getcircuitproperties('U')
    prop_b              =       femm.mo_getcircuitproperties('V')
    prop_c              =       femm.mo_getcircuitproperties('W')
    gap_B               =       femm.mo_getgapb('APairgap',ang)
    gap_M               =       femm.mo_gapintegral("APairgap",0)

    psi_abc = [prop_a[2], prop_b[2], prop_c[2]]
    psi_dq = T_abc_dq(psi_abc,theta_act_el_rad)
    psi_d = psi_dq[0]
    psi_q = psi_dq[1]
    b_rad = gap_B[0]
    b_tan = gap_B[1]
    m_tot = gap_M

    femm.mi_selectgroup(1)

    os.remove(f'temp_{iq}.FEM')
    os.remove(f'temp_{iq}.ans')
    return psi_d,psi_q,b_rad,b_tan,m_tot
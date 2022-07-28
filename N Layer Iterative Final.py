import numpy as np

import matplotlib.pyplot as plt

''' 
    This program contains functions to compute the reflectivity of N layers for a multilayer system. It contains the functions
    calculate_rp_comp and calculate_rs_comp. The details of each are provided below. 

    Author: Kyle Sprague

    Date Written: 13 June 2022
'''


def calculate_rp_comp(n, d, theta_one, lamb):

    '''
        Function to sum the reflection vectors stemming from N layers of a multilayer material for the P spin case /
        It begins with the "base case", which is to say the angle of incidence and refraction for the last layer of the material. These /
        are used to determine the reflectivity of the last layer, indicated by r_base. We then drop into a for loop that calculates
        the rest of the reflectivities by starting at the second layer going up, calculating the angle of refraction and angle of incidence, then
        finding the reflectivity for this particular layer r_part; subsequently, k for the given layer is calculated and then r_part is used to find
        r_base, with r_base representing the reflectivity of the given layer and the previous layer. The for loop moves backwards through all layers until/
        it reaches the top, as indicated by the negative indexation .

        list -- n : refractive indexes list
        list -- d : list of separations between each layer of material
        theta_one -- float : the first angle of incidence
        lamb -- float : the wavelength of the incident light

        returns: r_base, the value of the reflectivity associated with P spin for all layers of the material

    '''

    theta_one=theta_one+0.*1j #converts theta one to a complex number such that all subsequent values of theta are complex and the graphed RP/RS value is then also complex
    cos_th_inc = np.sqrt(1-(n[0]/n[-2])**2*np.sin(theta_one)**2)#base case
    cos_th_ref = np.sqrt((1-(n[0]/n[-1])**2*np.sin(theta_one)**2))#base case
    r_base = ((n[-2]*cos_th_ref)-(n[-1]*cos_th_inc))/((n[-2]*cos_th_ref)+(n[-1]*cos_th_inc))#base case

    for layer_num in range(2,len(n)): #complex case
        cos_th_ref = np.sqrt(1-((n[0]/n[-layer_num])**2)*np.sin(theta_one)**2)
        cos_th_inc = np.sqrt(1-(n[0]/n[-layer_num-1])**2*np.sin(theta_one)**2) 
        r_part = ((n[-layer_num-1]*cos_th_ref)-(n[-layer_num]*cos_th_inc))/((n[-layer_num-1]*cos_th_ref)+(n[-layer_num]*cos_th_inc))
        k = ((2*(np.pi)/lamb)*cos_th_ref)*n[-layer_num] 
        r_base = (r_part + r_base*np.exp(1j*2*k*d[-layer_num])) / (1+r_base*r_part*np.exp(1j*2*k*d[-layer_num])) 

    return r_base

def calculate_rs_comp(n,d,theta_one,lamb):

    '''
        Function to sum the reflection vectors stemming from N layers of a multilayer material for the S spin case /
        It begins with the "base case", which is to say the angle of incidence and refraction for the last layer of the material. These /
        are used to determine the reflectivity of the last layer, indicated by r_base. We then drop into a for loop that calculates
        the rest of the reflectivities by starting at the second layer going up, calculating the angle of refraction and angle of incidence, then
        finding the reflectivity for this particular layer r_part; subsequently, k for the given layer is calculated and then r_part is used to find
        r_base, with r_base representing the reflectivity of the given layer and the previous layer. The for loop moves backwards through all layers until/
        it reaches the top, as indicated by the negative indexation .

        list -- n : refractive indexes list
        list -- d : list of separations between each layer of material
        theta_one -- float : the first angle of incidence
        lamb -- float : the wavelength of the incident light

        returns: r_base, the value of the reflectivity associated with S spin for all layers of the material

    '''

    theta_one=theta_one+0.*1j
    cos_th_inc = np.sqrt(1-(n[0]/n[-2])**2*np.sin(theta_one)**2) #base case
    cos_th_ref = np.sqrt((1-(n[0]/n[-1])**2*np.sin(theta_one)**2))#base case
    r_base = ((n[-2]*cos_th_inc)-(n[-1]*cos_th_ref))/((n[-2]*cos_th_inc)+(n[-1]*cos_th_ref)) #base case

    for layer_num in range(2,len(n)): #complex case
        cos_th_ref = np.sqrt(1-((n[0]/n[-layer_num])**2)*np.sin(theta_one)**2)
        cos_th_inc = np.sqrt(1-(n[0]/n[-layer_num-1])**2*np.sin(theta_one)**2) 
        r_part = ((n[-layer_num-1]*cos_th_inc)-(n[-layer_num]*cos_th_ref))/((n[-layer_num-1]*cos_th_inc)+(n[-layer_num]*cos_th_ref))
        k = ((2*(np.pi)/lamb)*cos_th_ref)*n[-layer_num] 
        r_base = (r_part + r_base*np.exp(1j*2*k*d[-layer_num])) / (1+r_base*r_part*np.exp(1j*2*k*d[-layer_num]))

    return r_base

 

def main():

    th_list = np.zeros(1000, dtype = complex)
    RPN_list = np.zeros(1000, dtype = complex)
    RSN_list = np.zeros(1000, dtype = complex) 
    lamb = 650.
    n = [1.5,0.08 + 4.24j,1.3, 1.0]
    d = [0,50., 5., 0]
    dth=(np.pi)/2000

    for m in range(1000):
        theta_one=m*dth
        RPN_list[m]=calculate_rp_comp(n, d, theta_one, lamb)
        th_list[m]=theta_one
    plt.plot(np.real(th_list*180/np.pi), np.abs(RPN_list)**2)
    plt.show()

 
    for m in range(1000):
        theta_one=m*dth
        RSN_list[m]=calculate_rs_comp(n, d, theta_one, lamb)
        th_list[m]=theta_one
    plt.plot(th_list*180/np.pi, np.abs(RSN_list)**2)
    plt.show()

main()

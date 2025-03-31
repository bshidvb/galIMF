import numpy as np
import scipy.integrate
import csv  # csv and izip/zip are used to create output files

def Random_sampling(Z_H, M_ecl):
    delta_alpha = 63
    Z_sun  = 0.0142
    rho = 10 ** (0.61 * np.log10(M_ecl) + 2.85)
    y = -0.14*Z_H + 0.99*np.log10(rho)
    
    m_min = 0.08
    m_1 = 0.5
    m_2 = 1
    m_max = 150
    
    alpha_1 = 1.3 + delta_alpha*(10**Z_H - 1)*Z_sun
    alpha_2 = 2.3 + delta_alpha*(10**Z_H - 1)*Z_sun
    
    if y < -0.87:
        alpha_3 = 2.3
    elif y >= -0.87:
        alpha_3 = -0.41*y + 1.94
    
    ## determining the normalization constant for the power law
    ## the normalization constant is chosen such that the integral of the power law from m_0 to m_max is unity
    Term_1=1/(1-alpha_1)*(m_1**(1-alpha_1)-m_min**(1-alpha_1))
    Term_2=m_1**(-alpha_1+alpha_2)*1/(1-alpha_2)*(m_2**(1-alpha_2)-m_1**(1-alpha_2))
    
    if alpha_3 == 1:
        Term_3=m_1**(-alpha_1+alpha_2)*m_2**(-alpha_2+alpha_3)*(np.log(m_max) - np.log(m_2))
    else:
        Term_3=m_1**(-alpha_1+alpha_2)*m_2**(-alpha_2+alpha_3)*1/(1-alpha_3)*(m_max**(1-alpha_3)-m_2**(1-alpha_3))

    ## calculate normalization constant
    k=1/(Term_1+Term_2+Term_3)
    ## calculate normalization constants for each part of the power law
    k_1=k
    k_2=k_1*m_1**(-alpha_1 + alpha_2)
    k_3=k_2*m_2**(-alpha_2 + alpha_3)  
    
    
    ## IMF function
    def xi(m):
        
        if m < m_min:
            return 0
        if m >= m_min and m < m_1:
            return k_1*m**(-alpha_1)
        if m >= m_1 and m < m_2:
            return k_2*m**(-alpha_2)
        if m >= m_2 and m < m_max:
            return k_3*m**(-alpha_3)
        if m >= m_max:
            return 0
    
    # X_1, representing the random variable at which the first part of the power law ends
    X_1= scipy.integrate.quad(xi, m_min, m_1)[0]

    # X_2, here being taken as the integral of the power law from m_1 to m_2
    # The second part of the power law ends at the integral from m_0 to m_2, so at X_1+X_2
    X_2= scipy.integrate.quad(xi, m_1, m_2)[0]

    # X_3, here being taken as the integral of the power law from m_2 to m_max
    # The third part of the power law ends at the integral from m_0 to m_max, so at X_1+X_2+X_3
    X_3= scipy.integrate.quad(xi, m_2, m_max)[0]

    if np.round((X_1 + X_2 + X_3), 2) != 1.0:
        print('X_2+ X_2 + X_3 are not 1 but' ,X_1 + X_2 + X_3, '-> wrong normalization')
    m_array=([])
    ### drawing while the sum of masses is smaller than the embedded cluster mass
    while np.sum(m_array)<M_ecl:


        ###draw new random number between zero and one
        X=np.random.uniform(0,1)

        ## when X<X_1 the mass is drawn from the first power law
        if X < X_1:
            m=((1-alpha_1)*X/k_1+m_min**(1-alpha_1))**(1/(1-alpha_1))
            m_array=np.append(m_array,m)
        
        ## when X>=X_1 and X<X_1+X_2 the mass is drawn from the second power law
        elif X >= X_1 and X < X_1+X_2:
            m=((1-alpha_2)*(X-X_1)/k_2+m_1**(1-alpha_2))**(1/(1-alpha_2))
            m_array=np.append(m_array,m)
        
        ## when X>=X_1+X_2 and X<X_1+X_2+X_3 the mass is drawn from the third power law
        elif X >= X_1+X_2 and X < X_1+X_2+X_3:
            if alpha_3 == 1:
                m = np.exp((X-X_1-X_2)/k_3+np.log(m_2))
            else:
                m=((1-alpha_3)*(X-X_1-X_2)/k_3+m_2**(1-alpha_3))**(1/(1-alpha_3))
            m_array=np.append(m_array,m)

        ## case that total sum eventually gets larger than M_ecl
        if np.sum(m_array)>M_ecl:
            m_array=np.delete(m_array,-1)
            #m_array = np.append(m_array, (M_ecl-np.sum(m_array)))
            break

    return m_array

#creating random sampled clusters
StarClusterMasses = np.logspace(2.5, 5.0, num=20, base=10)
M_over_Hs = [0, -0.27576, -0.58336]#[0, -0.1, -0.3, -0.4, -0.5, -0.6, -0.7] #np.arange(-0.5, 0, 0.1)
ages = np.arange(1e6,11.e6, 1.e6)

for mass in StarClusterMasses:
    for metallicity in M_over_Hs:
        for age in ages:
            m_array = Random_sampling(metallicity, mass)
            m_array = np.sort(m_array)
            numbers_array = np.full((len(m_array)), 1,  dtype=int)
        
            with open(f'metallicity_{metallicity}/Random_sampled_star_cluster_{metallicity}_{np.round(age, 1)}_{np.round(mass, 1)}_ini.txt', 'w') as file:
                    writer = csv.writer(file, delimiter=' ')
                    writer.writerow(np.array([np.round(metallicity, 3), age,  np.round(mass, 2)]))
                    writer.writerows(zip(m_array, numbers_array))   

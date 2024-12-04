# ZEPHYR Cell
import numpy as np
import matplotlib.pyplot as plt

class Zephyr:
    def __init__(self) -> None:
        self.bulk_youngs_modulus = 20e6 # Pa
        self.bulk_pratio = 0.3
        self.l = 0.008          # m
        self.h = 0.022           # m
        self.z = 0.01           # m
        self.t1 = 0.001         # m
        self.t2 = 0.002         # m
        self.theta_1 = 1.047    # 60 degrees
        self.theta_2 = 1.047    # 60 degrees
        self.a1_lattice_modulus = 1
        self.a2_lattice_modulus = 1
        self.a12_lattice_poissons_ratio = 0
        self.a21_lattice_poissons_ratio = 0
        self.m_of_inertia_1 = 1
        self.m_of_inertia_2 = 1
    
    def set_m_of_inertia(self):
        self.m_of_inertia_1 = (self.z*self.t1**3)/(12)
        self.m_of_inertia_2 = (self.z*self.t2**3)/(12)
        return 0

    def set_lattice_modulus(self):
        t_over_l_ratio = (1/2)*((self.t1/self.l)**3 + (self.t2/self.l)**3)
        # print(t_over_l_ratio)
        E1 = self.bulk_youngs_modulus*t_over_l_ratio*(self.l/self.h)*((np.cos(self.theta_1)**2)/(np.sin(self.theta_1)))
        self.a1_lattice_modulus = E1
        return 0

    def get_lattice_pratio(self):
        return (1 - self.theta_1/self.theta_2)

    def get_lattice_modulus(self):
        return self.a1_lattice_modulus

    def set_lattice_poissons_ratio(self):
        return 0
    
    def get_a1_force_displacement(self, F1):
        """
        For an F1 force, get resulting a1 displacement
        """
        a1_sigma = F1/(self.z*self.h)
        a1_delta = (a1_sigma*self.h*self.z*(np.sin(self.theta_1)**2)*(self.l**3))/(12*self.bulk_youngs_modulus*self.m_of_inertia_1)
        return a1_delta

    def plot_stress_and_strain(self, F_list, t1_str):
        """
        For a given list of Forces, find stress and strain state for a Zephyr cell
        plot stress vs strain
        """
        # Testing results
        F_test = [0,1,2,3,4,5] # N approximate from spring scale
        epsilon_test = [0,0.072, 0.102, 0.131, 0.155, 0.168]
        sigma_test = []
        for i in F_test:
            # force applied over limited area
            sigma_value = i/(0.0008*self.z)
            sigma_test.append(sigma_value)
        
        # For each force, use Young's mod and buckling equation to estimate stress and strain
        bulk_epsilon = []
        bulk_sigma = []
        test_epsilon = np.linspace(0,0.25,100)
        for i in test_epsilon:
            bulk_epsilon.append(i)
            bulk_sigma.append(self.bulk_youngs_modulus*i)
        sigma = []
        epsilon = []

        for i in F_list:
            sigma_i = i/(self.h*self.z)
            epsilon_i = self.get_a1_force_displacement(i)
            epsilon.append(sigma_i/self.a1_lattice_modulus)
            sigma.append(100*sigma_i)

        # Get linear lattice Young's Modulus
        delta_sigma = sigma[2] - sigma[1]
        delta_epsilon = epsilon[2] - epsilon[1]
        youngs_modulus = (delta_sigma)/(delta_epsilon)
        print(youngs_modulus)
        plt.plot(epsilon, sigma, label="Lattice")
        plt.plot(bulk_epsilon, bulk_sigma, label="Solid")

        order = [1,2,3,4,5,6,7,8]
        growth_factor = 1
        for count, o in enumerate(order, 1):
            sigma_o = [i*growth_factor for i in sigma]
            plt.plot(epsilon, sigma_o, linestyle="dashed", label="HZ {}".format(count))
            growth_factor = growth_factor+0.7**o

        plt.scatter(epsilon_test, sigma_test, label="Experiment")
        plt.xlabel("Strain Along a1")
        plt.ylabel("Stress [Pa]")
        plt.xlim([0, 0.25])
        plt.ylim([0, 5e6])
        plt.title("Stress and Strain in a1 Loading of H-ZEPHYR (t1 = {} m, t2 = 0.002 m)".format(t1_str))
        plt.legend()
        plt.savefig("./stress_strain1_Hzephyr_{}.png".format(t1_str))
        plt.show()
        return 0
    
    def plot_strains_and_poissons_ratio(self):
        """
        Plot a1 and a2 strains and resulting P ratio
        """
        L_a1 = 0.050 # m
        L_a2 = 0.176 # m
        test_epsilon = np.linspace(0,0.25,100)
        bulk_epsilon2 = []
        epsilon1 = test_epsilon
        for i in test_epsilon:
            bulk_epsilon2.append(self.bulk_pratio*i)
            epsilon2 = epsilon1*self.get_lattice_pratio()

        epsilon1_test = [0,0.072, 0.102, 0.131, 0.155, 0.168]
        epsilon2_test = [0,0.001, 0.002, 0.01, 0.014, 0.035]
        plt.plot(epsilon1, epsilon2, label="Lattice")
        plt.plot(test_epsilon, bulk_epsilon2, label="Solid")

        theta_diff = [10,20,30,40]
        for count, o in enumerate(theta_diff, 1):
            bulk_eps_o = [i**(1.1+o/90) for i in bulk_epsilon2]
            plt.plot(test_epsilon, bulk_eps_o, linestyle="dashed", label="Theta Diff. {} Deg".format(str(o)))

        plt.scatter(epsilon1_test, epsilon2_test, label="Experiment")
        plt.xlabel("Epsilon 1")
        plt.ylabel("Epsilon 2")
        plt.title("Poisson's Ratio of ZEPHYR Lattice (t1 = 1 mm, t2 = 2 mm)")
        plt.legend()
        plt.savefig("./a1_a2_PR_Hzephyr.png")
        plt.show()
        
        return 0


if __name__ == "__main__":
    test_device = Zephyr()
    test_device.set_lattice_modulus()
    test_device.set_lattice_poissons_ratio()
    test_device.set_m_of_inertia()

    test_device.t1 = 0.001
    test_device.set_lattice_modulus()
    test_device.set_lattice_poissons_ratio()
    test_device.set_m_of_inertia()
    test_device.plot_stress_and_strain(F_list=[0,0.1,0.2,1,2,4,5,10,20], t1_str=str(test_device.t1))
    test_device.plot_strains_and_poissons_ratio()

    test_theta = list(np.linspace(0.5,1.4, 100))
    modulus = []
    for i in test_theta:
        test_device.theta_1 = i
        test_device.set_lattice_modulus()
        modulus.append(80*test_device.get_lattice_modulus())
    plt.plot(np.rad2deg(test_theta), modulus, label="ZEPHR Lattice")
    order = [1,2,3,4,5,6,7,8]
    growth_factor = 1
    for count, o in enumerate(order, 1):
        modulus_o = [i*(growth_factor) for i in modulus]
        plt.plot(np.rad2deg(test_theta), modulus_o, linestyle='dashed', label="HZ {}".format(count))
        growth_factor = growth_factor+0.7**o

    plt.xlabel("Theta (Degrees)")
    plt.ylabel("Young's Modulus")
    plt.title("Comparing Theta to Young's Modulus for Standard ZEPHYR")
    plt.legend()
    plt.savefig("./HZ_theta_and_modulus.png")
    plt.show()

    test_device.theta_1 = 1

    test_h = list(np.linspace(0.001,0.030,100))
    modulus = []
    for i in test_theta:
        test_device.h = i
        test_device.set_lattice_modulus()
        res = modulus.append(4000*test_device.get_lattice_modulus())
    
    plt.plot(test_h, modulus, label="ZEPHR Lattice")
    order = [1,2,3,4,5,6,7,8]
    growth_factor = 1
    for count, o in enumerate(order, 1):
        modulus_o = [i*(growth_factor) for i in modulus]
        plt.plot(test_h, modulus_o, linestyle='dashed', label="HZ {}".format(count))
        growth_factor = growth_factor+0.7**o
    
    plt.xlabel("h (mm)")
    plt.ylabel("Young's Modulus")
    plt.title("Comparing h to Young's Modulus for Standard ZEPHYR")
    plt.legend()
    plt.savefig("./HZ_h_and_modulus.png")
    plt.show()
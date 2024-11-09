# Shear Lag Model
import numpy as np
import matplotlib.pyplot as plt

class ShearLagComposite:
    def __init__(self) -> None:
        self.matrix_modulus = 1e6 # Pa
        self.matrix_poissonsratio = 0.4 
        self.fiber_modulus = 1e9
        self.fiber_poissonsratio = 0.3
        self.applied_strain = 0.01
        self.volume_fraction = 0.3
        self.n_factor = 1
        self.s_factor = 2

    def set_n_factor(self):
        E_m = self.matrix_modulus
        E_f = self.fiber_modulus
        nu_m = self.matrix_poissonsratio
        nu_f = self.fiber_poissonsratio
        f = self.volume_fraction
        n = np.sqrt((2*E_m)/(E_f*(1+nu_m)*np.log(1/f)))
        self.n_factor = n
        return 0

    def get_ddx_axial_stress(self):
        return 0

    def get_interfacial_shear_stress(self):
        return 0

    def get_ddr_shear_strain(self):
        return 0

    def get_axial_stress(self):
        """
        From composite object parameters, get list of x positions and axial stress
        """
        E_f = self.fiber_modulus
        eps_1 = self.applied_strain
        n = self.n_factor
        s = self.s_factor

        x = np.linspace(-1,1,1000)
        axial_stress = []
        for i in x:
            axial_stress_at_x = E_f*eps_1*(1 - (np.cosh(n*s*i))/(np.cosh(n*s)))
            axial_stress.append(axial_stress_at_x)
        return x, axial_stress
    
    def get_interfacial_shear_stress(self):
        """
        From composite object parameters, get list of x positions and IS stress
        """
        E_f = self.fiber_modulus
        eps_1 = self.applied_strain
        n = self.n_factor
        s = self.s_factor

        x = np.linspace(-2,2,1000)
        interfacial_shear_stress = []
        for i in x:
            interfacial_shear_stress_at_x = ((E_f*eps_1*n)/(2))*((np.sinh(n*s*i))/(np.cosh(n*s)))
            interfacial_shear_stress.append(interfacial_shear_stress_at_x)
        return x, interfacial_shear_stress


def plot_axial_stress(x, y, x_label, y_label):
    """
    """
    plt.plot(x, y)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    #plt.show()
    return 0


p1a_material = ShearLagComposite()
p1a_material.set_n_factor()

p1a_material.s_factor = 2
p1a_x, p1a_ax_stress_2 = p1a_material.get_axial_stress()
p1a_material.s_factor = 5
p1a_x, p1a_ax_stress_5 = p1a_material.get_axial_stress()
p1a_material.s_factor = 10
p1a_x, p1a_ax_stress_10 = p1a_material.get_axial_stress()
p1a_material.s_factor = 20
p1a_x, p1a_ax_stress_20 = p1a_material.get_axial_stress()
p1a_material.s_factor = 30
p1a_x, p1a_ax_stress_30 = p1a_material.get_axial_stress()
p1a_material.s_factor = 40
p1a_x, p1a_ax_stress_40 = p1a_material.get_axial_stress()
p1a_material.s_factor = 50
p1a_x, p1a_ax_stress_50 = p1a_material.get_axial_stress()

plt.plot(p1a_x, p1a_ax_stress_2, label="s=2")
plt.plot(p1a_x, p1a_ax_stress_5, label="s=5")
plt.plot(p1a_x, p1a_ax_stress_10, label="s=10")
plt.plot(p1a_x, p1a_ax_stress_20, label="s=20")
plt.plot(p1a_x, p1a_ax_stress_30, label="s=30")
plt.plot(p1a_x, p1a_ax_stress_40, label="s=40")
plt.plot(p1a_x, p1a_ax_stress_50, label="s=50")
plt.xlabel("x")
plt.ylabel("Axial Stress")
plt.legend()
plt.savefig("./axial_stress_p1a.png")
plt.show()


p1a_material.s_factor = 2
p1a_x, p1a_is_stress_2 = p1a_material.get_interfacial_shear_stress()
p1a_material.s_factor = 5
p1a_x, p1a_is_stress_5 = p1a_material.get_interfacial_shear_stress()
p1a_material.s_factor = 10
p1a_x, p1a_is_stress_10 = p1a_material.get_interfacial_shear_stress()
p1a_material.s_factor = 20
p1a_x, p1a_is_stress_20 = p1a_material.get_interfacial_shear_stress()
p1a_material.s_factor = 30
p1a_x, p1a_is_stress_30 = p1a_material.get_interfacial_shear_stress()
p1a_material.s_factor = 40
p1a_x, p1a_is_stress_40 = p1a_material.get_interfacial_shear_stress()
p1a_material.s_factor = 50
p1a_x, p1a_is_stress_50 = p1a_material.get_interfacial_shear_stress()

plt.plot(p1a_x, p1a_is_stress_2, label="s=2")
plt.plot(p1a_x, p1a_is_stress_5, label="s=5")
plt.plot(p1a_x, p1a_is_stress_10, label="s=10")
plt.plot(p1a_x, p1a_is_stress_20, label="s=20")
plt.plot(p1a_x, p1a_is_stress_30, label="s=30")
plt.plot(p1a_x, p1a_is_stress_40, label="s=40")
plt.plot(p1a_x, p1a_is_stress_50, label="s=50")
plt.xlabel("x")
plt.ylabel("Interfacial Shear Stress")
plt.legend()
plt.savefig("./is_stress_p1a.png")
plt.show()



plt.plot(p1a_x, p1a_is_stress_30, label="s=30")
plt.plot(p1a_x, p1a_is_stress_40, label="s=40")
plt.plot(p1a_x, p1a_is_stress_50, label="s=50")
plt.xlabel("x")
plt.ylabel("Interfacial Shear Stress")
plt.legend()
plt.savefig("./is_stress_p1a.png")
plt.show()
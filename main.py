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
            interfacial_shear_stress_at_x = ((E_f*eps_1*n)/(1))*((np.sinh(n*s*i))/(np.cosh(n*s)))
            interfacial_shear_stress.append(interfacial_shear_stress_at_x)
        return x, interfacial_shear_stress

    def get_modulus_of_composite(self):
        """
        """
        E1 = self.volume_fraction*self.fiber_modulus*(1-(np.tanh(self.n_factor*self.s_factor)/(self.n_factor*self.s_factor)))+(1-self.volume_fraction)*self.matrix_modulus
        return E1

"""
def plot_axial_stress(x, y, x_label, y_label):
    plt.plot(x, y)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    #plt.show()
    return 0
"""

p1a_material = ShearLagComposite()
p1a_material.set_n_factor()


# P1a
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
plt.xlabel("x/L")
plt.ylabel("Axial Stress")
plt.title("P1A")
plt.legend()
plt.savefig("./axial_stress_p1a.png")
plt.show()


# P1b
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
plt.xlabel("x/L")
plt.ylabel("Interfacial Shear Stress")
plt.title("P1B")
plt.legend()
plt.savefig("./is_stress_p1b.png")
plt.show()


# P1c
# examine range of fiber modulus from 1 to 1000
fiber_modulus = np.linspace(1e6, 1e9, 100)
E_f_over_E_m = np.linspace(1,1000,100)
p1a_material.s_factor = 5
p1_mod_0p2 = []
p1_mod_0p5 = []
p1_mod_0p8 = []

p1a_material.volume_fraction = 0.2
for i in fiber_modulus:
    p1a_material.fiber_modulus = i
    p1_mod_0p2.append( (p1a_material.get_modulus_of_composite())/(1000000) )

p1a_material.volume_fraction = 0.5
for i in fiber_modulus:
    p1a_material.fiber_modulus = i
    p1_mod_0p5.append( (p1a_material.get_modulus_of_composite())/(1000000) )

p1a_material.volume_fraction = 0.8
for i in fiber_modulus:
    p1a_material.fiber_modulus = i
    p1_mod_0p8.append( (p1a_material.get_modulus_of_composite())/(1000000) )

plt.scatter(E_f_over_E_m, p1_mod_0p2, label="f=0.2")
plt.scatter(E_f_over_E_m, p1_mod_0p5, label="f=0.5")
plt.scatter(E_f_over_E_m, p1_mod_0p8, label="f=0.8")
plt.xlabel("Fiber/Matrix Modulus Ratio")
plt.ylabel("Composite/Matrix Modulus Ratio")
plt.title("P1C")
plt.legend()
plt.savefig("./mod_ratio_p1c.png")
plt.show()

#P1d
# examine range of fiber modulus from 1 to 1000
fiber_modulus = np.linspace(1e6, 1e9, 100)
E_f_over_E_m = np.linspace(1,1000,100)
p1a_material.s_factor = 5
p1_mod_0p2 = []
p1_mod_0p5 = []
p1_mod_0p8 = []

p1a_material.volume_fraction = 0.2
for i in fiber_modulus:
    p1a_material.fiber_modulus = i
    p1_mod_0p2.append( (p1a_material.get_modulus_of_composite())/(i) )

p1a_material.volume_fraction = 0.5
for i in fiber_modulus:
    p1a_material.fiber_modulus = i
    p1_mod_0p5.append( (p1a_material.get_modulus_of_composite())/(i) )

p1a_material.volume_fraction = 0.8
for i in fiber_modulus:
    p1a_material.fiber_modulus = i
    p1_mod_0p8.append( (p1a_material.get_modulus_of_composite())/(i) )

ax = plt.gca()
ax.scatter(E_f_over_E_m, p1_mod_0p2, label="f=0.2")
ax.scatter(E_f_over_E_m, p1_mod_0p5, label="f=0.5")
ax.scatter(E_f_over_E_m, p1_mod_0p8, label="f=0.8")
ax.set_yscale('log')
ax.set_xscale('log')
plt.xlabel("Fiber/Matrix Modulus Ratio")
plt.ylabel("Composite/Fiber Modulus Ratio")
plt.title("P1D")
plt.legend()
plt.savefig("./mod_ratio_p1d.png")
plt.show()

# P1E
# examine range of fiber modulus from 1 to 1000
fiber_modulus = np.linspace(1e6, 1e9, 100)
E_f_over_E_m = np.linspace(1,1000,100)
p1a_material.volume_fraction = 0.5
p1_s_5 = []
p1_s_25 = []
p1_s_50 = []

p1a_material.s_factor = 5
for i in fiber_modulus:
    p1a_material.fiber_modulus = i
    p1_s_5.append( (p1a_material.get_modulus_of_composite())/(i) )

p1a_material.s_factor = 25
for i in fiber_modulus:
    p1a_material.fiber_modulus = i
    p1_s_25.append( (p1a_material.get_modulus_of_composite())/(i) )

p1a_material.s_factor = 50
for i in fiber_modulus:
    p1a_material.fiber_modulus = i
    p1_s_50.append( (p1a_material.get_modulus_of_composite())/(i) )

ax = plt.gca()
ax.scatter(E_f_over_E_m, p1_s_5, label="s=5")
ax.scatter(E_f_over_E_m, p1_s_25, label="s=25")
ax.scatter(E_f_over_E_m, p1_s_50, label="s=50")
ax.set_yscale('log')
ax.set_xscale('log')
plt.xlabel("Fiber/Matrix Modulus Ratio")
plt.ylabel("Composite/Fiber Modulus Ratio")
plt.title("P1E")
plt.legend()
plt.savefig("./mod_ratio_p1e.png")
plt.show()





# P2
p2_material = ShearLagComposite()
p2_material.matrix_modulus = 1e6
p2_material.matrix_poissonsratio = 0.4
p2_material.fiber_poissonsratio = 0.3

# P2a
p2a_s_3 = []
p2a_s_30 = []
p2_material.volume_fraction = 0.3

p2_material.s_factor = 3
for i in fiber_modulus:
    p2_material.fiber_modulus = i
    p2a_s_3.append( (p2_material.get_modulus_of_composite())/(i) )

p2_material.s_factor = 30
for i in fiber_modulus:
    p2_material.fiber_modulus = i
    p2a_s_30.append( (p2_material.get_modulus_of_composite())/(i) )

# rule of mixture result
rule_of_mix_modulus = []
for i in fiber_modulus:
    rule_of_mix_modulus.append((0.3*i + (1-0.3)*p2_material.matrix_modulus)/(i))

# modulus ratio
fiber_modulus = np.linspace(1e6, 1e9, 100)
E_f_over_E_m = np.linspace(1,1000,100)

ax = plt.gca()
ax.scatter(E_f_over_E_m, p1_s_5, label="s=3")
ax.scatter(E_f_over_E_m, p1_s_25, label="s=30")
ax.scatter(E_f_over_E_m, rule_of_mix_modulus, label="Rule of Mixture")
ax.set_yscale('log')
ax.set_xscale('log')
plt.xlabel("Fiber/Matrix Modulus Ratio")
plt.ylabel("Composite/Fiber Modulus Ratio")
plt.title("P2A")
plt.legend()
plt.savefig("./mod_ratio_p2a.png")
plt.show()

# P2b
p2a_s_3 = []
p2a_s_30 = []
p2_material.volume_fraction = 0.9

p2_material.s_factor = 3
for i in fiber_modulus:
    p2_material.fiber_modulus = i
    p2a_s_3.append( (p2_material.get_modulus_of_composite())/(i) )

p2_material.s_factor = 30
for i in fiber_modulus:
    p2_material.fiber_modulus = i
    p2a_s_30.append( (p2_material.get_modulus_of_composite())/(i) )

# rule of mixture result
rule_of_mix_modulus = []
for i in fiber_modulus:
    rule_of_mix_modulus.append((0.9*i + (1-0.9)*p2_material.matrix_modulus)/(i))

# modulus ratio
fiber_modulus = np.linspace(1e6, 1e9, 100)
E_f_over_E_m = np.linspace(1,1000,100)

ax = plt.gca()
ax.scatter(E_f_over_E_m, p1_s_5, label="s=3")
ax.scatter(E_f_over_E_m, p1_s_25, label="s=30")
ax.scatter(E_f_over_E_m, rule_of_mix_modulus, label="Rule of Mixture")
ax.set_yscale('log')
ax.set_xscale('log')
plt.xlabel("Fiber/Matrix Modulus Ratio")
plt.ylabel("Composite/Fiber Modulus Ratio")
plt.title("P2B")
plt.legend()
plt.savefig("./mod_ratio_p2b.png")
plt.show()

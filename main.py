from ase.calculators.espresso import Espresso
from ase.atoms import Atoms
from ase.visualize import view
from ase.build import make_supercell
from ase.io.trajectory import Trajectory
from ase.io import read as aread
from ase.units import kJ
from ase.eos import EquationOfState
import numpy as np
import matplotlib.pylab as plt

a0_expt = 5.45
cell_expt = a0_expt * np.array([
    (0, 0.5, 0.5),
    (0.5, 0, 0.5),
    (0.5, 0.5, 0)
])
atoms = Atoms(
    ["Ga", "P"],
    scaled_positions=[(0, 0, 0), (0.75, 0.25, 0.75)],
    cell=cell_expt,
    pbc=True
)
view(atoms)
P = np.array([
    (-1, 1, 1),
    (1, -1, 1),
    (1, 1, -1)
])
np.dot(P, atoms.cell)
np.array([[5.45, 0.0, 0.0],
          [0.0, 5.45, 0.0],
          [0.0, 0.0, 5.45]])
sc = make_supercell(atoms, P)
view(sc)
view(sc * (4, 4, 4))

nk = 6
ecut = 40

prefix = 'GaP'
input_data = {
    'control': {
        'prefix': prefix,
        'calculation': 'scf',
        'outdir': './tmp',
        'tprnfor': True, 'tstress': True
    },
    'system': {'ecutwfc': ecut, 'occupations': 'tetrahedra'},
    'electrons': {'conv_thr': 1e-9}
}

calc = Espresso(
    input_data=input_data,
    pseudopotentials={"Ga": "Ga.pbe-hgh.UPF", "P": "P.pbe-hgh.UPF"}, pseudo_dir=".",
    kpts=(nk, nk, nk)
)

atoms.calc = calc
atoms.get_potential_energy()

ecut_list = list(range(40, 110, 10))
etot_at_list = []

for ecut in ecut_list:
    input_data['system']['ecutwfc'] = ecut
    calc.set(input_data=input_data)
    calc.calculate(atoms)
    etot = atoms.get_potential_energy()
    etot_at = etot / atoms.get_global_number_of_atoms()
    etot_at_list.append(etot_at)
    print(f"ecut = {ecut:3} Ry; etot = {etot:9.6f} eV; etot_at = {etot_at:9.6f} eV")

etot_at_list = np.array(etot_at_list)
de_list = etot_at_list[0:-1] - etot_at_list[1:]

plt.figure()
plt.plot(ecut_list[0:-1], de_list, "b.-")
plt.plot(ecut_list[0:-1], 1e-3 * np.ones(np.shape(ecut_list[0:-1])), "k:")
plt.yscale("log")
plt.xlabel(r"$\epsilon_{cut}$, Ry")
plt.ylabel(r"$\Delta \epsilon_{tot}$, eV per atom")
plt.show()

ecut = 70
input_data['system']['ecutwfc'] = ecut
nk_list = list(range(4, 15, 1))
etot_at_list = []

for nk in nk_list:
    calc.set(kpts=(nk, nk, nk))
    calc.calculate(atoms)
    etot = atoms.get_potential_energy()
    etot_at = etot / atoms.get_global_number_of_atoms()
    etot_at_list.append(etot_at)
    print(f"nk = {nk:3}; etot = {etot:9.6f} eV; etot_at = {etot_at:9.6f} eV")

etot_at_list = np.array(etot_at_list)
de_list = np.abs(etot_at_list[0:-1] - etot_at_list[1:])

plt.figure()
plt.plot(nk_list[0:-1], de_list, "b.-")
plt.plot(nk_list[0:-1], 1e-3 * np.ones(np.shape(nk_list[0:-1])), "k:")
plt.yscale("log")
plt.xlabel(r"$n_k$, Ry")
plt.ylabel(r"$\Delta \epsilon_{tot}$, eV per atom")
plt.show()

cell = atoms.get_cell()
traj_file = f'{prefix}_eos.traj'
traj = Trajectory(traj_file, 'w')
n_cfg = 11
for x in np.linspace(0.95, 1.06, n_cfg):
    atoms.set_cell(cell_expt * x, scale_atoms=True)
    atoms.get_potential_energy()
    traj.write(atoms)

configs = aread(traj_file + "@0:%d" % n_cfg)  # read n_cfg configurations
# Extract volumes and energies:
volumes = [atoms.get_volume() for atoms in configs]
energies = [atoms.get_potential_energy() for atoms in configs]
eos = EquationOfState(volumes, energies)
v0, e0, B = eos.fit()
print(B / kJ * 1.0e24, 'GPa')
eos.plot(f'{prefix}-eos.png')

a0 = (4 * v0) ** (1 / 3)
print("Lattice parameters:")
print(f"expt: a0 = {a0_expt:6.3f}A, theor: a0 = {a0:6.3f}A, rerr: {(a0 - a0_expt) / a0_expt * 100:4.2f}%")

atoms.set_cell(cell_expt * a0 / a0_expt, scale_atoms=True)
atoms.cell[0, 1] * 2

"""
Зонная структура
"""

input_data['control'].update({
    'calculation': 'scf',
    'verbosity': 'high'
})
calc.set(kpts=(nk, nk, nk), input_data=input_data)
atoms.get_potential_energy()

Ef = calc.get_fermi_level()
print(Ef)

lat = atoms.cell.get_bravais_lattice()
print(lat.description())
lat.plot_bz(show=True)

path = atoms.cell.bandpath(npoints=200)
print(path)

input_data['control'].update({
    'calculation': 'bands',
    'verbosity': 'high'
})
input_data['system'].update({'nbnd': 12})
calc.set(kpts=path, input_data=input_data)
calc.calculate(atoms)

bs = calc.band_structure()

bs.plot(emin=-7, emax=Ef + 4, filename='bands_ase.png')
(kpt, k_lab, lab) = bs.get_labels()
bands = bs.energies
(n_spin, n_kpt, n_bands) = bands.shape
print(n_spin, n_kpt, n_bands)

Emin = -3
Emax = 3
kmin = np.min(kpt)
kmax = np.max(kpt)
plt.figure()
for i_s in range(n_spin):
    for i_b in range(n_bands):
        plt.plot(kpt, bands[i_s, :, i_b] - Ef, "k-")
#
lab_tex = [r"$\Gamma$" if l == "G" else r"$%s$" % l for l in lab]
plt.xticks(k_lab, lab_tex)
for k in k_lab:
    plt.plot([k, k], [Emin, Emax], "k:")
#
plt.ylim(Emin, Emax)
plt.xlim(kmin, kmax)
#
plt.ylabel(r"$\epsilon-\epsilon_F$, eV")
plt.savefig("bands.png", dpi=300)

"""

"""

input_data['control'].update({
    'calculation': 'scf',
    'verbosity': 'high'
})
calc.set(kpts=(nk, nk, nk), input_data=input_data)
atoms.get_potential_energy()
Ef = calc.get_fermi_level()
print(Ef)
nk_dos = 32

input_data['control'].update({
    'calculation': 'nscf',
    'verbosity': 'high'
})
input_data['system'].update({'nbnd': 12})
calc.set(kpts=(nk_dos, nk_dos, nk_dos), input_data=input_data)
calc.calculate(atoms)

open("dos.in", "w").write(
    f"""&dos
  prefix  = '{prefix}'
  outdir = './tmp'
  fildos = '{prefix}.dos'
  Emin = {Ef + Emin}, Emax = {Ef + Emax}, DeltaE = 0.005
/ 
"""
)

dos = np.loadtxt(f"{prefix}.dos")
plt.figure()
plt.plot(dos[:, 0], dos[:, 1])
plt.xlabel('energy [eV]')
plt.ylabel('DOS')
plt.show()
plt.show()

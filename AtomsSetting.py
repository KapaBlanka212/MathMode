"""

"""
from ase.calculators.espresso import Espresso
from ase.atoms import Atoms
from ase.visualize import view
from ase.build import make_supercell
import numpy as np

from Constants import ConstantPATHS, ConstantNAMES


class AtomControl:
    def __init__(self, latticeConstant_exp: float = None):
        """

        :param latticeConstant_exp: the experimental value of material lattice constant
        """
        if latticeConstant_exp is None:
            self.latticeConstant_exp = 5.45
        else:
            self.latticeConstant_exp = latticeConstant_exp
        self.scaledPositionGaP = [(0, 0, 0), (0.75, 0.25, 0.75)]
        self.cellMatrices = np.array([
            (0, 0.5, 0.5),
            (0.5, 0, 0.5),
            (0.5, 0.5, 0)
        ])
        self.transformationMatrix = np.array([
            (-1, 1, 1),
            (1, -1, 1),
            (1, 1, -1)
        ])
        self.experimentalCell = self.latticeConstant_exp * self.cellMatrices
        self.atom: Atoms = np.eye(3)

    def createAtom(self, isView: bool = False):
        """

        :param isView:
        :return:
        """
        self.atom = Atoms(
            [ConstantNAMES().Ga, ConstantNAMES().P],
            scaled_positions=self.scaledPositionGaP,
            cell=self.experimentalCell,
            pbc=True
        )
        if isView:
            view(self.atom)
        return self.atom

    def makeSuperCell(self, isView: bool = False):
        transformedCell = self.experimentalCell @ self.transformationMatrix
        superCell = make_supercell(self.atom, transformedCell)
        if isView:
            view(superCell)
        return superCell


if __name__ == "__main__":
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
    atom = AtomControl().createAtom()
    atom.calc = calc
    print(atom.get_potential_energy())

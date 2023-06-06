"""

"""
import os


class ConstantNAMES:
    def __init__(self):
        """

        """
        # Material names constant
        self.Ga = "Ga"
        self.P = "P"
        self.GaP = "GaP"
        # Solver constants name for input data
        self.control = "control"
        self.prefix = "prefix"
        self.calculation = "calculation"
        self.scf = "scf"
        self.outDir = "outdir"
        self.tprnFor = "tprnfor"
        self.tstress = "tstress"
        self.system = "system"
        self.EnergyCut = "ecutwfc"
        self.occupations = "occupations"
        self.tetrahedra = "tetrahedra"
        self.electrons = "electrons"
        self.conv_thr = "conv_thr"


class ConstantPATHS:
    def __init__(self):
        self.resultsPATH = "Results"
        self.dirPseudoPOT = "pseudoPot"
        self.pseudoPOT_Ga = "Ga.pbe-hgh.UPF"
        self.pseudoPOT_P = "P.pbe-hgh.UPF"
        self.dirPseudoPot_Ga = os.path.join(self.dirPseudoPOT, self.pseudoPOT_Ga)
        self.dirPseudoPot_P = os.path.join(self.dirPseudoPOT, self.pseudoPOT_P)
        self.tmp = r"tmp"
        self.dirTMP = os.path.join(self.resultsPATH, self.tmp)


if __name__ == "__main__":
    print(os.path.isfile(ConstantPATHS().dirPseudoPot_P))
    print(ConstantPATHS().dirPseudoPOT)
    print(os.path.isdir(ConstantPATHS().dirPseudoPOT))
    print(os.path.dirname(ConstantPATHS().dirPseudoPOT))
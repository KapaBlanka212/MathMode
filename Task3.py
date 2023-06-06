"""

"""
from typing import Union, NoReturn, List
import json
import csv

from ase.calculators.espresso import Espresso
from ase.atoms import Atoms

import numpy as np
import matplotlib.pylab as plt

from AtomsSetting import AtomControl
from Constants import ConstantPATHS, ConstantNAMES


class PlaneWave:
    def __init__(self, initialEnergyCut: Union[float, int],
                 initialNk: int):
        """

        :param initialNk:
        :param initialEnergyCut:
        """
        self.inputDataDict = None
        """
        UNIT: SET ATTRIBUTES:
        ~~~~~~~~~~~~~~~~~~~~~
        Creating Attributes from Objects of Constant and Atom Classes
        """
        self.atoms: Atoms = AtomControl().createAtom(False)
        self.constantNAMES: ConstantNAMES = ConstantNAMES()
        self.constantPATHS: ConstantPATHS = ConstantPATHS()

        self.__prefix = ConstantNAMES().GaP
        self.__EnergyCut = initialEnergyCut
        self.__initialNk = initialNk
        self.calculationConfig = None

    def createInputDict(self) -> dict:
        controlDict = {self.constantNAMES.prefix: self.__prefix,
                       self.constantNAMES.calculation: self.constantNAMES.scf,
                       self.constantNAMES.outDir: self.constantPATHS.dirTMP,
                       self.constantNAMES.tprnFor: True,
                       self.constantNAMES.tstress: True}
        systemDict = {self.constantNAMES.EnergyCut: self.__EnergyCut,
                      self.constantNAMES.occupations: self.constantNAMES.tetrahedra}
        electronDict = {self.constantNAMES.conv_thr: 1e-9}
        inputDataDict = {self.constantNAMES.control: controlDict,
                         self.constantNAMES.system: systemDict,
                         self.constantNAMES.electrons: electronDict}
        PlaneWave(self.__EnergyCut, self.__initialNk).exportToJSON(inputDataDict, "ConfigGaP")
        return inputDataDict

    @staticmethod
    def exportToJSON(parametersDict: dict, name) -> NoReturn:
        """
        Saved the dictionary with material parameters to JSON format
        :param parametersDict: The dict that contained the solver config
        :param name: The material names ("Header" of dict)
        :return: NoReturn
        """
        with open(f"{name}.json", 'w') as output:
            json.dump(parametersDict, output, indent=4)

    @staticmethod
    def exportToCSV(rows: List, name):
        """

        :param rows:
        :param name:
        :return:
        """
        with open(f"{name}.csv", "w") as output:
            writer = csv.writer(output, delimiter=",")
            writer.writerows(rows)

    def setConfigCalculation(self, inputData):
        """

        :param inputData:
        :return:
        """
        # TODO: add comments
        self.calculationConfig: Espresso = Espresso(
            input_data=inputData,
            pseudopotentials={"Ga": "Ga.pbe-hgh.UPF", "P": "P.pbe-hgh.UPF"}, pseudo_dir=".",
            kpts=(self.__initialNk, self.__initialNk, self.__initialNk)
        )
        self.atoms.calc = self.calculationConfig

    def findEnergyCut(self, amountOfPoint: int):
        """

        :param amountOfPoint:
        :return:
        """
        inputData = self.createInputDict()
        self.setConfigCalculation(inputData)
        energyCut_linspace = np.linspace(self.__EnergyCut, self.__EnergyCut * 3, amountOfPoint)
        energyTotal_list = []
        for energy in energyCut_linspace:
            inputData[self.constantNAMES.system][self.constantNAMES.EnergyCut] = energy
            self.setConfigCalculation(inputData)
            self.calculationConfig.calculate(self.atoms)
            energyTotal = self.atoms.get_potential_energy()
            energyTotal_eV = energyTotal / self.atoms.get_global_number_of_atoms()
            energyTotal_list.append([energy, energyTotal, energyTotal_eV])
        print(f"ecut = {energy:3} Ry; etot = {energyTotal:9.6f} eV; etot_at = {energyTotal_eV:9.6f} eV")
        print(energyTotal_list)

    def findNk(self, amountOfPoint: int):
        """

        :param amountOfPoint:
        :return:
        """
        inputData = self.createInputDict()
        self.setConfigCalculation(inputData)
        nk_list = np.arange(self.__initialNk, self.__initialNk * 3)
        energyTotal_at_list = []
        for nk in nk_list:
            self.calculationConfig.set(kpts=(nk, nk, nk))
            self.calculationConfig.calculate(self.atoms)
            energyTotal = self.atoms.get_potential_energy()
            energyTotal_at = energyTotal / self.atoms.get_global_number_of_atoms()
            energyTotal_at_list.append(energyTotal_at)
            print(f"nk = {nk:3}; etot = {energyTotal:9.6f} eV; etot_at = {energyTotal_at:9.6f} eV")

    @staticmethod
    def showPlot(self, plotDataX: Union[List, np.ndarray], plotDataY: Union[List, np.ndarray],
                 *args, **kwargs):
        """

        :param self:
        :param plotDataX:
        :param plotDataY:
        :param args:
        :param kwargs:
        :return:
        """
        plt.figure()
        plt.plot(plotDataX[0:-1], plotDataY, "b.-")
        plt.plot(plotDataX[0:-1], 1e-3 * np.ones(np.shape(plotDataX[0:-1])), "k:")
        if kwargs:
            plotSetting = {key: kwargs[key] for key in kwargs.keys()}
            plt.yscale("log")
            plt.xlabel(r"$\epsilon_{cut}$, Ry")
            plt.ylabel(r"$\Delta \epsilon_{tot}$, eV per atom")
        plt.show()


if __name__ == "__main__":
    PlaneWave(40, 6).findEnergyCut(3)
    PlaneWave(40, 6).findNk(10)

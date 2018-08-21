import PhylloxeraPy
PhylloxeraPy.loadLibraries(True)
import ROOT

from morpho.utilities import morphologging, reader
logger = morphologging.getLogger(__name__)

from morpho.processors import BaseProcessor
from mermithid.misc import Constants

from morpho.processors.plots import RootCanvas, RootHistogram
increase_range = 10.  # energy increase required for the convolution product to work


class DistortedGaussianGenerator(BaseProcessor):
    '''
    Generate a smeared tritium spectrum.
    '''

    def InternalConfigure(self, config_dict={}):
        '''
        Required class attributes:
        - volume [m3]
        - density [1/m3]
        - experiment duration [s]
        - neutrino mass [eV]
        - energy window [KEmin,KEmax]
        - background [counts/eV/s]
        - energy resolution [eV]
        '''

        return True

    def _PrepareWorkspace(self):
        # for key in config_dict:
        #     setattr(self, key, config_dict[key])

        # We have to define the range here: the setRange() methods are only useful for the calculating integrals
        # The energy window is increased to accommodate the convolution
        KE = ROOT.RooRealVar("KE", "KE", 0, 100)
        # logger.debug(self.KEmin-self.increase_range,self.KEmax+self.increase_range)
        # KE.setRange("FullRange",0,Constants.tritium_endpoint()*2)
        # KE.setRange("Window",self.KEmin-10,self.KEmax+10)
        mean = ROOT.RooRealVar("mean", "mean", 50, 40, 60)
        width = ROOT.RooRealVar("width", "width", 10, 5, 20)
        # Define the standard tritium spectrum
        gaussian = ROOT.RooGaussian(
            "gaussian", "gaussian", KE, mean, width)

        # Define PdfFactory to add background, smearing and efficiency multiplication
        pdffactory = ROOT.PdfFactory("myPdfFactory")

        coefficients_values = [100,-1]
        cname = "coeff0"
        c0 = ROOT.RooRealVar(cname, cname, coefficients_values[0])
        #c1 = ROOT.RooRealVar("coeff1", "coeff1", 0., self.efficiency_coefficients[1])
        polynomialCoefficients = ROOT.RooArgList(c0)
        logger.info("HERERERERE")
        polynomialCoefficients.Print()
        polynomialCoefficients.removeAll()
        # polynomialCoefficients.add(KE)

        # for i in range(1, polynomialEfficiencyOrder):
        c = {}
        for i in range(0, len(coefficients_values)):
            cname = "coeff{}".format(i)
            c[i] = ROOT.RooRealVar(
                cname, cname, coefficients_values[i])
            c[i].Print()

            polynomialCoefficients.add(c[i])
        polyn = ROOT.RooPolynomial(
            "abc", "abc", KE, polynomialCoefficients, int(0))
        #polyn = ROOT.RooFormulaVar("f", "f","polynomial", polynomialCoefficients)

        logger.info("HERERERERE")
        polynomialCoefficients.Print()

        NEvents = ROOT.RooRealVar(
            "NEvents", "NEvents", 100)
        NBkgd = ROOT.RooRealVar(
            "NBkgd", "NBkgd", 100)


        backedSpectrum = pdffactory.AddBackground(ROOT.RooAbsPdf)(
            "backedSpectrum", KE, gaussian, NEvents, NBkgd)

        totalSpectrum = pdffactory.MultiplyPolynomialEfficiency(ROOT.RooAbsPdf)(
            "totalSpectrum", KE, backedSpectrum, polynomialCoefficients, int(0))
            #totalSpectrum = ROOT.RooClassFactory.makePdf("totalSpectrum", "x", "trueSpectrum(x) * f(x)")
        
        # Save things in a Workspace
        self.workspace = ROOT.RooWorkspace()

        # if self.doMultiplication:
        #    getattr(self.workspace, 'import')(polyn)
        #    getattr(self.workspace, 'import')(trueSpectrum)
        #    self.workspace.factory("PROD::totalSpectrum(trueSpectrum, abc)")
        # else:
        # getattr(self.workspace, 'import')(totalSpectrum)
        getattr(self.workspace, 'import')(totalSpectrum)
        # getattr(self.workspace, 'import')(background)

        self.workspace.Print()

        self.canvas = ROOT.TCanvas('a', 'a', 600, 400)
        self.canvas.Divide(2, 1)
        xframe = ROOT.RooPlot(KE, 0, 100, 100)
        xframe2 = ROOT.RooPlot(KE, 0, 100, 100)

        self.canvas.cd(1)

        #data = ROOT.RooDataSet(totalSpectrum.generate(ROOT.RooArgSet(KE),1000))
        # totalSpectrum.fitTo(data)
        gaussian.plotOn(xframe)
        backedSpectrum.plotOn(xframe)
        # data.plotOn(xframe)
        xframe.Draw()

        self.canvas.cd(2)
        polyn.plotOn(xframe2)
        xframe2.Draw()

        self.canvas.SaveAs('a.pdf')

    def InternalRun(self):
        self._PrepareWorkspace()
        self.results = self._GenerateData()
        return True

    def _GenerateData(self):

        logger.debug("Generate data")
        KE = self.workspace.var("KE")
        KE.setRange("window", 0, 100)

        #KE1 = self.workspace.var("KE1")
        #KE1.setRange("window", self.KEmin, self.KEmax)

        totalSpectrum = self.workspace.pdf("totalSpectrum")
        # background = self.workspace.pdf("background")

        dataLarge = totalSpectrum.generate(ROOT.RooArgSet(
            KE), 100000, ROOT.RooFit.Range("window"))
        dataLarge.Print()
        dataList = []
        for i in range(100000):
            dataList.append(dataLarge.get(i).getRealValue("KE"))
        # Have to delete the workspace to precent some nasty errors at the end...
        del self.workspace
        return {"KE": dataList}

if __name__ == "__main__":

    from morpho.processors.plots import Histogram
    proc = DistortedGaussianGenerator("a_name")
    proc.Configure({})
    proc.Run()

    hist = Histogram("a_histogram")
    hist.Configure({
        "data": "KE",
        "n_bins_x": 300,
        "output_path": "."
    })
    hist.data = proc.results
    hist.Run()
    # print(proc.results)
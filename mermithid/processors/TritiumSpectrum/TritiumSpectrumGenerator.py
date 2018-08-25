import PhylloxeraPy
PhylloxeraPy.loadLibraries(False)
import ROOT

from morpho.utilities import morphologging, reader
logger = morphologging.getLogger(__name__)

from morpho.processors import BaseProcessor
from mermithid.misc import Constants

from morpho.processors.plots import RootCanvas, RootHistogram
increase_range = 10.  # energy increase required for the convolution product to work


class TritiumSpectrumGenerator(BaseProcessor):
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

        self.KEmin, self.KEmax = reader.read_param(config_dict, "energy_window", [
                                                   Constants.tritium_endpoint()-1e3, Constants.tritium_endpoint()+1e3])
        self.volume = reader.read_param(config_dict, "volume", 1e-6)
        self.density = reader.read_param(config_dict, "density", 1e18)
        self.duration = reader.read_param(config_dict, "duration", 1e18)
        self.neutrinomass = reader.read_param(
            config_dict, "neutrino_mass", 1e18)
        self.background = reader.read_param(config_dict, "background", 1e-6)
        self.poisson_fluctuations = reader.read_param(
            config_dict, "poisson_fluctuations", False)
        self.energy_resolution = reader.read_param(
            config_dict, "energy_resolution", 0)
        self.numberDecays = reader.read_param(config_dict, "number_decays", -1)
        self.efficiency_coefficients = reader.read_param(config_dict, "efficiency_coefficients", [0])
        return True

    def _GetNEvents_Window(self, KE, spectrum):
        '''
        Calculate the number of decays events generated in the energy window
        '''
        KE = ROOT.RooRealVar("KE_tmp", "KE_tmp", 0,
                             Constants.tritium_endpoint()*2)
        KE.setRange("FullRange", 0, Constants.tritium_endpoint()*2)
        KE.setRange("Window", self.KEmin-self.increase_range,
                    self.KEmax+self.increase_range)

        m_nu = ROOT.RooRealVar("m_nu_tmp", "m_nu_tmp",
                               self.neutrinomass, -200, 200)
        endpoint = ROOT.RooRealVar("endpoint_tmp", "endpoint_tmp", Constants.tritium_endpoint(),
                                   Constants.tritium_endpoint()-10.,
                                   Constants.tritium_endpoint()+10.)
        # Define the standard tritium spectrum
        spectrum = ROOT.RealTritiumSpectrum(
            "spectrum_tmp", "spectrum_tmp", KE, endpoint, m_nu)

        # Define PdfFactory to add background and smearing
        pdffactory = ROOT.PdfFactory("myPdfFactory_tmp")
        if self.doSmearing:
            meanSmearing = ROOT.RooRealVar(
                "meanSmearing_tmp", "meanSmearing_tmp", 0)
            widthSmearing = ROOT.RooRealVar(
                "widthSmearing_tmp", "widthSmearing_tmp", self.energy_resolution)
            smearedspectrum = pdffactory.GetSmearedPdf(ROOT.RealTritiumSpectrum)(
                "smearedspectrum_tmp", 2, KE, spectrum, meanSmearing, widthSmearing, 1000000)
        fullSpectrumIntegral = spectrum.createIntegral(
            ROOT.RooArgSet(KE), ROOT.RooFit.Range("FullRange"))
        windowSpectrumIntegral = spectrum.createIntegral(
            ROOT.RooArgSet(KE), ROOT.RooFit.Range("Window"))
        ratio = windowSpectrumIntegral.getVal()/fullSpectrumIntegral.getVal()
        logger.debug("Fraction in window [{},{}]: {}".format(
            self.KEmin-self.increase_range, self.KEmax+self.increase_range, ratio))
        return ratio

    def _GetNBackground_Window(self):
        '''
        Calculate the number of background events generated in the energy window
        '''
        return self.background * (self.KEmax - self.KEmin + 2*self.increase_range) * self.duration

    def _PrepareWorkspace(self):
        # for key in config_dict:
        #     setattr(self, key, config_dict[key])
        if hasattr(self, "energy_resolution") and self.energy_resolution > 0.:
            logger.debug("Will use a smeared spectrum with {} eV energy res.".format(
                self.energy_resolution))
            self.increase_range = 10*self.energy_resolution
            self.doSmearing = True
        else:
            logger.debug("Will use a normal spectrum")
            self.increase_range = 0
            self.doSmearing = False

        if hasattr(self, "efficiency_coefficients") and len(self.efficiency_coefficients)>1:
            logger.debug("Will multiply pdf with efficiency")
            self.doMultiplication = True
            logger.debug(self.efficiency_coefficients)
        else:
            self.doMultiplication = False

        # We have to define the range here: the setRange() methods are only useful for the calculating integrals
        # The energy window is increased to accommodate the convolution
        KE = ROOT.RooRealVar("KE", "KE", self.KEmin -
                             self.increase_range, self.KEmax+self.increase_range)
        X = ROOT.RooRealVar("x", "x", self.KEmin -
                             self.increase_range, self.KEmax+self.increase_range)
        # logger.debug(self.KEmin-self.increase_range,self.KEmax+self.increase_range)
        # KE.setRange("FullRange",0,Constants.tritium_endpoint()*2)
        # KE.setRange("Window",self.KEmin-10,self.KEmax+10)
        m_nu = ROOT.RooRealVar("m_nu", "m_nu", self.neutrinomass)
        endpoint = ROOT.RooRealVar("endpoint", "endpoint", Constants.tritium_endpoint())
        b = ROOT.RooRealVar("background", "background", -
                            self.background)

        # Define the standard tritium spectrum
        spectrum = ROOT.RealTritiumSpectrum(
            "spectrum", "spectrum", KE, endpoint, m_nu)

        # Define gaussian
        mean = ROOT.RooRealVar("mean", "mean", 18600)
        width = ROOT.RooRealVar("width", "width", 500)
        # Define the standard tritium spectrum
        gaussian = ROOT.RooGaussian(
            "gaussian", "gaussian", KE, mean, width)

        # Define PdfFactory to add background, smearing and efficiency multiplication
        pdffactory = ROOT.PdfFactory("myPdfFactory")
        # Smearing of the spectrum
        if self.doSmearing:
            meanSmearing = ROOT.RooRealVar("meanSmearing", "meanSmearing", 0.)
            widthSmearing = ROOT.RooRealVar(
                "widthSmearing", "widthSmearing", self.energy_resolution, 0., 10*self.energy_resolution)
            # KE.setBins(100000, "cache")
            smearedSpectrum = pdffactory.GetSmearedPdf(ROOT.RooAbsPdf)(
                "smearedSpectrum", 2, KE, spectrum, meanSmearing, widthSmearing, 100000)



        # Background
        background = ROOT.RooUniform(
            "background", "background", ROOT.RooArgSet(KE))

        # Calculate number of events and background
        if self.numberDecays <= 0:
            number_atoms = self.volume*self.density
            total_number_decays = number_atoms*self.duration/Constants.tritium_lifetime()
            if self.doSmearing:
                number_decays_window = total_number_decays * \
                    self._GetNEvents_Window(KE, smearedSpectrum)
            else:
                number_decays_window = total_number_decays * \
                    self._GetNEvents_Window(KE, spectrum)
        else:
            number_decays_window = self.numberDecays
        logger.debug("Number decays in window: {}".format(
            number_decays_window))
        self.number_bkgd_window = self._GetNBackground_Window()
        logger.debug("Number bkgd in window: {}".format(
            self.number_bkgd_window))

        # Calculate the number of events to generate:
        # - If the poisson_fluctuations is True, it uses a Poisson process to get the number of events to generate
        #   (the mean being the number of events in the window)
        # - Else use the value calculated
        if self.poisson_fluctuations:
            ran = ROOT.TRandom3()
            self.number_decays_window_to_generate = int(
                ran.Poisson(number_decays_window))
            self.number_bkgd_window_to_generate = int(
                ran.Poisson(self.number_bkgd_window))
        else:
            self.number_decays_window_to_generate = int(number_decays_window)
            self.number_bkgd_window_to_generate = int(self.number_bkgd_window)

        NEvents = ROOT.RooRealVar(
            "NEvents", "NEvents", self.number_decays_window_to_generate)
        NBkgd = ROOT.RooRealVar(
            "NBkgd", "NBkgd", self.number_bkgd_window_to_generate)

        if self.doSmearing:
            backgroundSpectrum = pdffactory.AddBackground(ROOT.RooAbsPdf)(
                    "backgroundSpectrum", KE, smearedSpectrum, NEvents, NBkgd)
        else:
            backgroundSpectrum = pdffactory.AddBackground(ROOT.RooAbsPdf)(
                    "backgroundSpectrum", KE, spectrum, NEvents, NBkgd)

        # Efficiency
        if self.doMultiplication:
            cname = "coeff0"
            c0 = ROOT.RooRealVar(cname, cname, self.efficiency_coefficients[0])
            #c1 = ROOT.RooRealVar("coeff1", "coeff1", 0., self.efficiency_coefficients[1])
            polynomialCoefficients = ROOT.RooArgList(c0)
            polynomialCoefficients.removeAll()
            polynomialCoefficients.add(KE)
            polynomialCoefficients.add(c0)
            formula = '{}'.format(cname)


            #for i in range(1, polynomialEfficiencyOrder):
            c = {}
            for i in range(1, len(self.efficiency_coefficients)):
                cname = "coeff{}".format(i)
                c[i] = ROOT.RooRealVar(cname, cname, self.efficiency_coefficients[i])
                c[i].Print()
                formula+='+{}*pow(KE,{})'.format(cname, i)

                polynomialCoefficients.add(c[i])

            #polyn = ROOT.RooPolynomial("abc", "abc", KE, polynomialCoefficients, int(0))
            polyn = ROOT.RooFormulaVar("f", "f","@1+@2*@0",polynomialCoefficients)
            polyn.Print()

            polynomialCoefficients.Print()

            totalSpectrum = ROOT.RooEffProd("totalSpectrum", "totalSpectrum", backgroundSpectrum, polyn)

            #totalSpectrum = pdffactory.AddBackground(ROOT.RooAbsPdf)(
            #        "totalSpectrum", KE, spectrum, NEvents, NBkgd)


        # Save things in a Workspace
        self.workspace = ROOT.RooWorkspace()
        getattr(self.workspace, 'import')(totalSpectrum)
        """if self.doMultiplication:
            getattr(self.workspace, 'import')(polyn)
            getattr(self.workspace, 'import')(backgroundSpectrum)
            self.workspace.factory("PROD::totalSpectrum(backgroundSpectrum, abc)")"""

        self.workspace.Print()

        self.canvas = ROOT.TCanvas('a','a',800,400)
        self.canvas.Divide(2,1);
        xframe = ROOT.RooPlot(KE, self.KEmin, self.KEmax, 100)
        yframe = ROOT.RooPlot(KE, self.KEmin, self.KEmax, 100)

        self.canvas.cd(1)

        #data = ROOT.RooDataSet(totalSpectrum.generate(ROOT.RooArgSet(KE),1000))
        #totalSpectrum.fitTo(data)
        #self.workspace.pdf("totalSpectrum").plotOn(xframe)
        totalSpectrum.plotOn(xframe)
        #backgroundSpectrum.plotOn(xframe)
        #gaussian.Print()
        #data.plotOn(xframe)
        xframe.Draw()

        self.canvas.cd(2)
        polyn.plotOn(yframe)
        yframe.Draw()

        self.canvas.SaveAs('a.pdf')

    def InternalRun(self):
        self._PrepareWorkspace()
        self.results = self._GenerateData()
        return True

    def _GenerateData(self):

        logger.debug("Generate data")
        KE = self.workspace.var("KE")
        KE.setRange("window", self.KEmin, self.KEmax)

        #KE1 = self.workspace.var("KE1")
        #KE1.setRange("window", self.KEmin, self.KEmax)

        #background = self.workspace.pdf("background")
        totalSpectrum = self.workspace.pdf("totalSpectrum")
        print('events in spectrum: {}'.format(self.number_decays_window_to_generate))
        print('background counts in spectrum: {}'.format(self.number_bkgd_window_to_generate))
        totalEvents = self.number_decays_window_to_generate + \
            self.number_bkgd_window_to_generate
        dataLarge = totalSpectrum.generate(ROOT.RooArgSet(
            KE), totalEvents, ROOT.RooFit.Range("window"))
        dataLarge.Print()
        data = dataLarge.reduce(ROOT.RooFit.CutRange("window"))
        dataList = []
        for i in range(data.numEntries()):
            dataList.append(data.get(i).getRealValue("KE"))
        # Have to delete the workspace to precent some nasty errors at the end...
        del self.workspace
        return {"KE": dataList}

if __name__ == "__main__":

    from morpho.processors.plots import Histogram
    from mermithid.misc.Constants import seconds_per_year, tritium_endpoint

    specGen_config = {
            "volume": 7e-6*1e-2, # [m3]
            "density": 3e17, # [1/m3]
            "duration": 1.*seconds_per_year()/12., # [s]
            "neutrino_mass" :0, # [eV]
            "energy_window": [tritium_endpoint()-1e3,tritium_endpoint()+1e3], # [KEmin,KEmax]
            # "energy_window": [0.,tritium_endpoint()+1e3], # [KEmin,KEmax]
            "background": 10e-6, # [counts/eV/s]
            #"energy_resolution": 1,# [eV]
            "efficiency_coefficients": [1, -1./20000]
        }
    proc = TritiumSpectrumGenerator("specGen")
    proc.Configure(specGen_config)
    proc.Run()

    hist = Histogram("a_histogram")
    hist.Configure({
        "data": "KE",
        "n_bins_x": 300,
        "output_path": ".",
        "title": "spectrum"
    })
    hist.data = proc.results
    hist.Run()
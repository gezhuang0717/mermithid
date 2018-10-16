import PhylloxeraPy
PhylloxeraPy.loadLibraries(False)
import ROOT

from morpho.utilities import morphologging, reader
logger = morphologging.getLogger(__name__)

from morpho.processors.sampling import RooFitLikelihoodSampler
from mermithid.misc import Constants

class TritiumSpectrumLikelihoodSampler(RooFitLikelihoodSampler):

    def InternalConfigure(self,config_dict):
        super().InternalConfigure(config_dict)
        self.null_m_nu = reader.read_param(config_dict,"null_neutrino_mass",True)
        self.neutrinomass = reader.read_param(
            config_dict, "neutrino_mass", 1e18)
        self.energy_resolution = reader.read_param(
            config_dict, "energy_resolution", 0)
        self.efficiency_coefficients = reader.read_param(config_dict, "efficiency_coefficients", [0])
        self.background_coefficients = reader.read_param(config_dict, "background_shape_coefficients", [0])
        self.NTotal = reader.read_param(config_dict, "NTotal", 0)
        self.TritiumRatio= reader.read_param(config_dict, "TritiumRatio", 0.5)

        self.NEvents = int(self.TritiumRatio*self.NTotal)
        self.NBkgd = int((1-self.TritiumRatio)*self.NTotal)

        print(self.NTotal, self.NEvents, self.NBkgd)
        print(self.binned)
        return True


    def definePdf(self,wspace):
        '''
        Defines the Pdf that RooFit will sample and add it to the workspace.
        Users should edit this function.
        '''
        self.l = []
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

        if hasattr(self, "background_coefficients") and len(self.background_coefficients) > 1:
            logger.debug("Will use non uniform background.")
            self.background_shape = 1
            logger.debug(self.background_coefficients)
        else:
            logger.debug("Will use a unifrom background")
            self.background_shape = 0

        if hasattr(self, "efficiency_coefficients") and len(self.efficiency_coefficients) > 1:
            logger.debug("Will multiply pdf with efficiency")
            self.doMultiplication = True
            logger.debug(self.efficiency_coefficients)
        else:
            self.doMultiplication = False


        var = wspace.var(self.varName)

        m_nu = ROOT.RooRealVar("m_nu", "m_nu", self.neutrinomass)
        endpoint = ROOT.RooRealVar("endpoint", "endpoint", Constants.tritium_endpoint())
                                  # Constants.tritium_endpoint()-10.,
                                  # Constants.tritium_endpoint()+10.)


        #eventRatio = ROOT.RooRealVar("TritiumRatio", "TritiumRatio", self.TritiumRatio, 0, 1)
        NBkgd = ROOT.RooRealVar("NBkgd", "NBkgd", 1)
        NEvents = ROOT.RooRealVar("NEvents", "NEvents", self.NEvents, 0, self.NTotal)


        # Define the standard tritium spectrum
        spectrum = ROOT.RealTritiumSpectrum(
            "spectrum", "spectrum", var, endpoint, m_nu)

        # Define gaussian
        mean = ROOT.RooRealVar("mean", "mean", 18600)
        width = ROOT.RooRealVar("width", "width", 500)
        gaussian = ROOT.RooGaussian(
            "gaussian", "gaussian", var, mean, width)
        self.l.append(spectrum)
        self.l.append(gaussian)

        # Define PdfFactory to add background, smearing and efficiency multiplication
        pdffactory = ROOT.PdfFactory("myPdfFactory")
        # Smearing of the spectrum
        if self.doSmearing:
            meanSmearing = ROOT.RooRealVar("meanSmearing", "meanSmearing", 0.)
            widthSmearing = ROOT.RooRealVar(
                "widthSmearing", "widthSmearing", self.energy_resolution, 0., 10*self.energy_resolution)
            # KE.setBins(100000, "cache")
            smearedSpectrum = pdffactory.GetSmearedPdf(ROOT.RooAbsPdf)(
                "smearedSpectrum", 2, var, spectrum, meanSmearing, widthSmearing, 100000)
            smearedSpectrum.Print()

        # Efficiency
        if self.doMultiplication:
            while len(self.efficiency_coefficients) < 5:
                self.efficiency_coefficients.extend([0])
            print(self.efficiency_coefficients)

            cname = "coeff0"
            c0 = ROOT.RooRealVar(cname, cname, self.efficiency_coefficients[0])
            polynomialCoefficients = ROOT.RooArgList(c0)
            polynomialCoefficients.removeAll() # this is ugly but I don't know what else to do
            polynomialCoefficients.add(var) # for RooFormulaVar
            polynomialCoefficients.add(c0)

            c = {}
            for i in range(1, len(self.efficiency_coefficients)):
                cname = "coeff{}".format(i)
                c[i] = ROOT.RooRealVar(cname, cname, self.efficiency_coefficients[i])
                #self.l.append(c[i])
                polynomialCoefficients.add(c[i])

            self.l.append(polynomialCoefficients)
            polyn = ROOT.RooFormulaVar("f", "f","@1+@2*TMath::Power(@0,1)+@3*TMath::Power(@0,2)+@4*TMath::Power(@0,3)+@5*TMath::Power(@0,4)",polynomialCoefficients)
            self.l.append(polyn)

            if self.doSmearing:
                multipliedSpectrum = pdffactory.MultiplyPolynom(ROOT.RooAbsPdf)(
                        "ps", "multipliedSpectrum", polynomialCoefficients, smearedSpectrum)
                #multipliedSpectrum = ROOT.RooEffProd("multipliedSpectrum", "multipliedSpectrum", smearedSpectrum, polyn)
            else:
                multipliedSpectrum = pdffactory.MultiplyPolynom(ROOT.RooAbsPdf)(
                        "ps", "multipliedSpectrum", polynomialCoefficients, spectrum)
                #multipliedSpectrum = ROOT.RooEffProd("multipliedSpectrum", "multipliedSpectrum", gaussian, polyn)

            self.l.append(multipliedSpectrum)

        # Background
        if self.background_shape == 1:
            print("extending background coefficient list")
            while len(self.background_coefficients) < 5:
                self.background_coefficients.extend([0])
            cname = "coeffb0"
            b0 = ROOT.RooRealVar(cname, cname, self.background_coefficients[0])
            backgroundCoefficients = ROOT.RooArgList(b0)
            backgroundCoefficients.removeAll() # this is ugly but I don't know what else to do
            backgroundCoefficients.add(var) # for RooFormulaVar
            backgroundCoefficients.add(b0)


            b = {}
            for i in range(1, len(self.background_coefficients)):
                cname = "coeffb{}".format(i)
                b[i] = ROOT.RooRealVar(cname, cname, self.background_coefficients[i])
                backgroundCoefficients.add(b[i])
                #self.l.append(b[i])
            backgroundCoefficients.Print()
            self.l.append(backgroundCoefficients)
            polyn1 = ROOT.RooFormulaVar("f", "f","@1+@2*TMath::Power(@0,1)+@3*TMath::Power(@0,2)+@4*TMath::Power(@0,3)+@5*TMath::Power(@0,4)", backgroundCoefficients)
            self.l.append(polyn1)


        else:
            backgroundCoefficients = ROOT.RooArgList()



        if self.doMultiplication:
            pdf = pdffactory.AddBackground(ROOT.RooAbsPdf)(
                        "pdf", self.background_shape, var, multipliedSpectrum, NEvents, NBkgd, backgroundCoefficients)
        elif self.doSmearing:
            pdf = pdffactory.AddBackground(ROOT.RooAbsPdf)(
                    "pdf", self.background_shape, var, smearedSpectrum, NEvents, NBkgd, backgroundCoefficients)
        else:
            pdf = pdffactory.AddBackground(ROOT.RooAbsPdf)(
                    "pdf", self.background_shape, var, spectrum, NEvents, NBkgd, backgroundCoefficients)

        self.l.append(pdf)

        # Save pdf: this will save all required variables and functions
        getattr(wspace,'import')(pdf)
        return wspace

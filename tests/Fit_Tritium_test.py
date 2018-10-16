'''
This scripts aims at testing Tritium specific processors.
Author: M. Guigue
Date: Apr 1 2018
'''

import unittest
import yaml


from morpho.utilities import morphologging
logger = morphologging.getLogger(__name__)

class FitTritiumTests(unittest.TestCase):

    def test_FitSpectrum(self):
        from mermithid.processors.TritiumSpectrum import TritiumSpectrumLikelihoodSampler
        from morpho.processors.plots import TimeSeries, APosterioriDistribution
        from mermithid.processors.plots import KuriePlotGeneratorProcessor
        from morpho.processors.plots import Histogram
        from mermithid.misc.Constants import seconds_per_year, tritium_endpoint

        with open('tritium_data.yaml', 'r') as infile:
            result = yaml.load(infile)

        spectrumSampler_config = {
            "binned": False,
            "n_jobs": 2,
            "warmup": 0,
            "varName": "KE",
            "nuisanceParams":[],
            "interestParams": ["NEvents"],
            "neutrinomass": 0, # [eV]
            "energy_window": [tritium_endpoint()-1.3e3,tritium_endpoint()+0.4e3], # [KEmin,KEmax]
            # "energy_window": [0.,tritium_endpoint()+1e3], # [KEmin,KEmax]
            "background_shape_coefficients": [0, 1/20000],
            #"energy_resolution": 5,# [eV]
            "efficiency_coefficients": [-1.66135862e+05, 3.66101578e+01, -3.02446439e-03, 1.11017923e-07, -1.52774133e-12],
            "NTotal": len(result['KE']),
            "TritiumRatio": 0.5,
            "NEvents": 50,
            "NBkgd": 50
        }
        histo_plot = {
            "data": "KE",
            "n_bins_x": 100,
            "title": "spectrum"
        }
        kurie_plot = {
            "data": "KE",
            "n_bins_x": 1000,
            "title": "kurie_plot"
        }
        aposteriori_config = {
            "n_bins_x": 100,
            "n_bins_y": 100,
            "data": ["NEvents"],
            "title": "aposteriori_distribution",
            "output_path": "."
        }
        timeSeries_config = {
            "data": ["NEvents"],
            "height": 1200,
            "title": "timeseries",
            "output_path": "."
        }

        with open('tritium_data.yaml', 'r') as infile:
            result = yaml.load(infile)

        spectrumSampler = TritiumSpectrumLikelihoodSampler("spectrumSampler")
        aposterioriPlotter = APosterioriDistribution("posterioriDistrib")
        timeSeriesPlotter = TimeSeries("timeSeries")
        histo = Histogram("histo")
        kurieHisto = KuriePlotGeneratorProcessor("kurieHisto")

        spectrumSampler.Configure(spectrumSampler_config)
        aposterioriPlotter.Configure(aposteriori_config)
        timeSeriesPlotter.Configure(timeSeries_config)
        histo.Configure(histo_plot)
        kurieHisto.Configure(kurie_plot)

        spectrumSampler.data = result
        histo.data = result
        kurieHisto.data = result
        histo.Run()
        kurieHisto.Run()

        spectrumSampler.Run()
        print(spectrumSampler.results)
        aposterioriPlotter.data = spectrumSampler.results
        timeSeriesPlotter.data = spectrumSampler.results
        aposterioriPlotter.Run()
        timeSeriesPlotter.Run()


if __name__ == '__main__':
    unittest.main()
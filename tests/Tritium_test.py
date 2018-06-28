'''
This scripts aims at testing Tritium specific processors.
Author: M. Guigue
Date: Apr 1 2018
'''

import unittest

from morpho.utilities import morphologging
logger = morphologging.getLogger(__name__)

class TritiumTests(unittest.TestCase):
    
    def test_KuriePlot(self):
        from mermithid.processors.TritiumSpectrum import TritiumSpectrumGenerator
        from mermithid.processors.plots import KuriePlotGeneratorProcessor
        from morpho.processors.plots import Histogram
        from mermithid.misc.Constants import seconds_per_year, tritium_endpoint

        specGen_config = {
            "volume": 7e-6*1e-2, # [m3]
            "density": 3e17, # [1/m3]
            "duration": 1.*seconds_per_year()/12., # [s]
            "neutrino_mass" :0, # [eV]
            "energy_window": [tritium_endpoint()-1e3,tritium_endpoint()+1e3], # [KEmin,KEmax]
            # "energy_window": [0.,tritium_endpoint()+1e3], # [KEmin,KEmax]
            "background": 1e-6, # [counts/eV/s]
            "energy_resolution": 5# [eV]
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

        specGen = TritiumSpectrumGenerator("specGen")
        histo = Histogram("histo")
        kurieHisto = KuriePlotGeneratorProcessor("kurieHisto")

        specGen.Configure(specGen_config)
        histo.Configure(histo_plot)
        kurieHisto.Configure(kurie_plot)

        specGen.Run()
        result = specGen.results
        histo.data = result
        kurieHisto.data = result
        histo.Run()
        kurieHisto.Run()

if __name__ == '__main__':
    unittest.main()
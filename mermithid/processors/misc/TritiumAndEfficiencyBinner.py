'''
Bin tritium start frequencies and calculate efficiency for each bin
function.
Author: A. Ziegler, E. Novitski, C. Claessens
Date:1/20/2020
'''

from __future__ import absolute_import

import sys
from ROOT import TF1, TMath, TH1F
from morpho.utilities import morphologging, reader
from morpho.processors import BaseProcessor
from morpho.processors.plots import RootCanvas, RootHistogram
from mermithid.misc import Constants

logger = morphologging.getLogger(__name__)



__all__ = []
__all__.append(__name__)

def Energy(f, B=None, Theta=None):
    #print(type(F))
    if B==None:
        B = 0.95777194923080811
    emass_kg = Constants.m_electron()*Constants.e()/(Constants.c()**2)
    if isinstance(f, list):
        gamma = [(Constants.e()*B)/(2.0*TMath.Pi()*emass_kg) * 1/(F) for F in f]
        return [(g -1)*Constants.m_electron() for g in gamma]
    else:
        gamma = (Constants.e()*B)/(2.0*TMath.Pi()*emass_kg) * 1/(f)
        return (gamma -1)*Constants.m_electron()

class TritiumAndEfficiencyBinner(BaseProcessor):

    def InternalConfigure(self, params):
        '''
        Configure
        '''
        self.rootcanvas = RootCanvas(params, optStat=0)


        # Read other parameters
        self.namedata = reader.read_param(params, 'variables', "required")
        self.N = reader.read_param(params, 'N', 'N')
        self.eff_eqn = reader.read_param(params, 'efficiency', '1')
        self.mode = reader.read_param(params, 'mode', 'unbinned')
        self.histogram_or_dictionary = reader.read_param(params, 'histogram_or_dictionary',
         'histogram')
        self.n_bins_x = reader.read_param(params, 'n_bins_x', 100)
        self.range = reader.read_param(params, 'range', [0., -1.])
        self.eff_func = TF1('eff_func', self.eff_eqn, self.range[0], self.range[1])
        self.histo = TH1F('histo', 'histo', self.n_bins_x, self.range[0], self.range[1])
        self.asInteger = reader.read_param(params, 'asInteger', False)
        self.energy_or_frequency = reader.read_param(params, 'energy_or_frequency', 'energy')
        self.total_counts = 0
        self.total_weighted_counts = 0
        self.eff_norm = 0

        if self.energy_or_frequency == 'energy':
            self.corrected_data = TH1F('corrrected_data', 'corrected_data', self.n_bins_x, Energy(self.range[1]), Energy(self.range[0]))
            print(sys.getrefcount(self.corrected_data))
            self.output_bin_variable='KE'
        elif self.energy_or_frequency == 'frequency':
            self.corrected_data = TH1F('corrrected_data', 'corrected_data', self.n_bins_x, self.range[0], self.range[1])
            self.output_bin_variable='F'

        return True

    def InternalRun(self):
        print('hi')
        if self.mode == 'unbinned':

            for i in self.data.get(self.namedata):
                self.histo.Fill(i)

        elif self.mode == 'binned':

            for i in range(self.n_bins_x):
                self.histo.SetBinContent(i + 1, self.data.get('N')[i])

        self.histo.Sumw2()

        for i in range(self.n_bins_x):

            self.total_counts += self.histo.GetBinContent(i + 1)
            self.total_weighted_counts += self.histo.GetBinContent(i + 1)/self.eff_func.Eval(self.histo.GetBinCenter(i + 1))

        self.eff_norm = self.total_counts/self.total_weighted_counts

        if self.asInteger:

            if self.energy_or_frequency == 'frequency':
                for i in range(self.n_bins_x):
                    self.corrected_data.SetBinContent(i + 1,
                    int(self.eff_norm * self.histo.GetBinContent(i + 1)/self.eff_func.Eval(self.histo.GetBinCenter(i + 1))))
            elif self.energy_or_frequency == 'energy':
                for i in range(self.n_bins_x):
                    self.corrected_data.SetBinContent(self.n_bins_x - i,
                    int(self.eff_norm * self.histo.GetBinContent(i + 1)/self.eff_func.Eval(self.histo.GetBinCenter(i + 1))))

        else:

            if self.energy_or_frequency == 'frequency':
                for i in range(self.n_bins_x):
                    self.corrected_data.SetBinContent(i + 1,
                    self.eff_norm * self.histo.GetBinContent(i + 1)/self.eff_func.Eval(self.histo.GetBinCenter(i + 1)))
            elif self.energy_or_frequency == 'energy':
                for i in range(self.n_bins_x):
                    self.corrected_data.SetBinContent(self.n_bins_x - i,
                    self.eff_norm * self.histo.GetBinContent(i + 1)/(self.eff_func.Eval(self.histo.GetBinCenter(i + 1))))

        self.corrected_data.Sumw2()

        if self.histogram_or_dictionary == 'dictionary':
            temp_dictionary = {self.output_bin_variable: [], 'N': []}
            for i in range(self.n_bins_x):
                temp_dictionary[self.output_bin_variable].append(self.corrected_data.GetBinCenter(i + 1))
                temp_dictionary['N'].append(int(self.corrected_data.GetBinContent(i + 1)))
            self.corrected_data = temp_dictionary
            #print(self.corrected_data.keys())

        return True

        def EfficiencyAssignment(self, f):
            efficiency = 0.9
            efficiency_error = 0.05
            return efficiency, efficiency_error

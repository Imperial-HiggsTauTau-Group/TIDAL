# Analysis.py

from collections import defaultdict, OrderedDict
import ROOT
import glob
import os
from Draw.python import MultiDraw
from uncertainties import ufloat
import ctypes
import yaml
import json

ROOT.TH1.AddDirectory(0)


def WriteToTFile(obj, file, path):
    file.cd()
    as_vec = path.split('/')
    if len(as_vec) >= 1:
        for i in range(0, len(as_vec) - 1):
            if not ROOT.gDirectory.GetDirectory(as_vec[i]):
                ROOT.gDirectory.mkdir(as_vec[i])
            ROOT.gDirectory.cd(as_vec[i])
    if not ROOT.gDirectory.FindKey(as_vec[-1]):
        obj.SetName(as_vec[-1])
        ROOT.gDirectory.WriteTObject(obj, as_vec[-1])
    ROOT.gDirectory.cd('/')


class Shape(object):
    def __init__(self, hist, rate=None):
        self.hist = hist
        if rate is not None:
            self.rate = rate

    def _Int(self):
        if self._hist.GetDimension() == 1:
            return self._hist.Integral(
                0, self._hist.GetNbinsX() + 1
            )
        elif self._hist.GetDimension() == 2:
            return self._hist.Integral(
                0, self._hist.GetNbinsX() + 1,
                0, self._hist.GetNbinsY() + 1
            )
        elif self._hist.GetDimension() == 3:
            return self._hist.Integral(
                0, self._hist.GetNbinsX() + 1,
                0, self._hist.GetNbinsY() + 1,
                0, self._hist.GetNbinsZ() + 1
            )

    def _IntErr(self):
        err = ctypes.c_double()
        if self._hist.GetDimension() == 1:
            val = self._hist.IntegralAndError(
                0, self._hist.GetNbinsX() + 1, err
            )
        elif self._hist.GetDimension() == 2:
            val = self._hist.IntegralAndError(
                0, self._hist.GetNbinsX() + 1,
                0, self._hist.GetNbinsY() + 1, err
            )
        elif self._hist.GetDimension() == 3:
            val = self._hist.IntegralAndError(
                0, self._hist.GetNbinsX() + 1,
                0, self._hist.GetNbinsY() + 1,
                0, self._hist.GetNbinsZ() + 1, err
            )
        return ufloat(val, err.value)

    @property
    def hist(self):
        # print 'Getting hist for %s' % self.__repr__()
        return self._hist

    @hist.setter
    def hist(self, hist):
        # print 'Setting hist for %s ' % self.__repr__()
        self._hist = hist.Clone()
        self._rate = self._IntErr()

    @property
    def rate(self):
        # print 'Getting rate for %s ' % self.__repr__()
        return self._rate

    @rate.setter
    def rate(self, rate):
        # print 'Setting rate to {:.3f}'.format(rate)
        self._rate = rate
        integral = self._Int()
        if integral == 0.:
            print('Error, histogram integral is zero')
            return
        self._hist.Scale(rate.n / integral)

    def copy(self):
        return Shape(self.hist.Clone(), self.rate)

    def __iadd__(self, other):
        if isinstance(other, Shape):
            self.hist.Add(other.hist)
            self.rate += other.rate
        else:
            self.rate += other
        return self

    def __isub__(self, other):
        if isinstance(other, Shape):
            self.hist.Add(other.hist, -1)
            self.rate -= other.rate
        else:
            self.rate -= other
        return self

    def __imul__(self, other):
        if isinstance(other, Shape):
            self.hist.Multiply(other.hist)
            self.rate *= other.rate
        else:
            self.rate *= other
        return self

    def __itruediv__(self, other):
        if isinstance(other, Shape):
            self.hist.Divide(other.hist)
            self.rate /= other.rate
        else:
            self.rate /= other
        return self

    def __add__(self, other):
        cpy = self.copy()
        return cpy.__iadd__(other)

    __radd__ = __add__

    def __sub__(self, other):
        cpy = self.copy()
        return cpy.__isub__(other)

    def __mul__(self, other):
        cpy = self.copy()
        return cpy.__imul__(other)

    def __truediv__(self, other):
        cpy = self.copy()
        return cpy.__itruediv__(other)

    __rmul__ = __mul__

    def Print(self):
        print('rate={:.3f}, entries={:g}, sum={:g}'.format(
            self.rate,
            self.hist.GetEntries(),
            self.hist.GetSumOfWeights()
        ))


class TTreeEvaluator:
    def __init__(self, tree_name, filename):
        self.tree_name = tree_name
        self.filename = filename
        self.tree = None  # initially none -> we'll open it later
        self.file = None

    def PrepareTree(self):
        self.file = ROOT.TFile(self.filename)
        self.tree = self.file.Get(self.tree_name)

    def Draw(self, draw_list, compiled=False):
        self.PrepareTree()
        otree = MultiDraw.MultiDraw(self.tree, draw_list, Compiled=compiled)
        self.file.Close()
        return otree


class BaseNode:
    def __init__(self, name, WriteSubnodes=True):
        self.name = name
        self.WriteSubnodes = WriteSubnodes

    def GetNameStr(self):
        return '%s::%s' % (self.__class__.__name__, self.name)

    def GetInfoStr(self):
        return ''

    def SubNodes(self):
        return []

    def PrintTree(self, indent=0):
        print('%s%s%s' % (' ' * indent, self.GetNameStr(), self.GetInfoStr()))
        for node in self.SubNodes():
            node.PrintTree(indent=(indent + 2))

    def AddRequests(self, manifest):
        pass

    def Objects(self):
        return dict()

    def OutputPrefix(self, node=None):
        return self.name

    def Output(self, file, prefix=''):
        objects = self.Objects()
        for key, val in objects.items():
            WriteToTFile(val, file, '%s/%s' % (prefix, key))
        if self.WriteSubnodes:
            for node in self.SubNodes():
                node.Output(file, '%s/%s' % (prefix, self.OutputPrefix(node)))

    def Run(self):
        for node in self.SubNodes():
            node.Run()
        self.RunSelf()

    def RunSelf(self):
        pass


class BasicNode(BaseNode):
    def __init__(self, name, sample, variable, selection, factors=list(), WriteSubnodes=True):
        BaseNode.__init__(self, name, WriteSubnodes)
        self.sample = sample
        self.variable = variable
        self.selection = selection
        self.factors = factors
        self.shape = None

    def GetInfoStr(self):
        if self.shape is not None:
            self.shape.Print()
        factors_str = ''
        if len(self.factors) > 0:
            factor_vals = [('%.4g' % x) for x in self.factors]
            factors_str = ', factors=' + '*'.join(factor_vals)
        return '[%s, variable=%s, selection=%s%s]' % (
            self.sample,
            self.variable.split(';')[0],
            self.selection,
            factors_str
        )

    # This means processing to the point where we can hand our shape over
    def RunSelf(self):
        for factor in self.factors:
            self.shape *= factor

    def AddRequests(self, manifest):
        manifest.append((self.sample, self.variable, self.selection, self, 'shape'))

    def Objects(self):
        return {self.name: self.shape.hist}


class ListNode(BaseNode):
    def __init__(self, name):
        BaseNode.__init__(self, name)
        self.nodes = OrderedDict()
        self.add_output_prefix = True

    def __getitem__(self, key):
        return self.nodes[key]

    def AddNode(self, node):
        self.nodes[node.name] = node

    def SubNodes(self):
        return self.nodes.values()

    def AddRequests(self, manifest):
        for node in self.SubNodes():
            node.AddRequests(manifest)

    def OutputPrefix(self, node=None):
        if self.add_output_prefix:
            return self.name
        else:
            return ''


class SummedNode(ListNode):
    def __init__(self, name):
        ListNode.__init__(self, name)
        self.shape = None

    def RunSelf(self):
        self.shape = sum([node.shape for node in self.nodes.values()])

    def Objects(self):
        return {self.name: self.shape.hist}

    def OutputPrefix(self, node=None):
        if self.add_output_prefix:
            return self.name + '.subnodes'
        else:
            return ''


class SubtractNode(BaseNode):
    def __init__(self, name, initial, subtract):
        BaseNode.__init__(self, name)
        self.shape = None
        self.initial_node = initial
        self.subtract_node = subtract

    def RunSelf(self):
        self.shape = self.initial_node.shape - self.subtract_node.shape

    def Objects(self):
        return {self.name: self.shape.hist}

    def OutputPrefix(self, node=None):
        return self.name + '.subnodes'

    def SubNodes(self):
        return [self.initial_node, self.subtract_node]

    def AddRequests(self, manifest):
        for node in self.SubNodes():
            node.AddRequests(manifest)
            
class FF_Node(BaseNode):
    def __init__(self, name, QCD_node, W_node, Top_node, QCD_frac, W_frac, Top_frac):
        BaseNode.__init__(self, name)
        self.shape = None
        self.QCD_node = QCD_node
        self.W_node = W_node
        self.Top_node = Top_node
        self.QCD_frac = QCD_frac
        self.W_frac = W_frac
        self.Top_frac = Top_frac

    def RunSelf(self):
        for i in range(1, self.QCD_frac.shape.hist.GetNbinsX() + 1):
            # set errors for fractions to 0 so don't influence bbb uncertainties of final templates
            self.QCD_frac.shape.hist.SetBinError(i, 0.0)
            self.W_frac.shape.hist.SetBinError(i, 0.0)
            self.Top_frac.shape.hist.SetBinError(i, 0.0)
        total_jet = self.QCD_frac.shape + self.W_frac.shape + self.Top_frac.shape
        self.shape = self.QCD_node.shape * self.QCD_frac.shape / total_jet + self.W_node.shape * self.W_frac.shape/total_jet + self.Top_node.shape * self.Top_frac.shape / total_jet

    def Objects(self):
        return {self.name: self.shape.hist}

    def OutputPrefix(self, node=None):
        return self.name + '.subnodes'

    def SubNodes(self):
        return [self.QCD_node, self.W_node, self.Top_node, self.QCD_frac, self.W_frac, self.Top_frac]

    def AddRequests(self, manifest):
        for node in self.SubNodes():
            node.AddRequests(manifest)


class Analysis(object):
    def __init__(self):
        self.trees = {}
        self.nodes = ListNode('')
        self.info = {}
        self.remaps = {}
        self.compiled = False
        self.WriteSubnodes = True

    def Run(self):
        manifest = []
        self.nodes.AddRequests(manifest)
        if len(manifest) != len(set(manifest)):
            raise ValueError('Duplicate requests in manifest!')
        drawdict = defaultdict(list)
        outdict = defaultdict(list)
        for entry in manifest:
            drawdict[entry[0]].append(entry[1:3])
            outdict[entry[0]].append(entry[3:5])

        for sample in drawdict:
            print(sample)
            res = self.trees[sample].Draw(drawdict[sample], compiled=self.compiled)
            res = [x for x in res if isinstance(x, ROOT.TH1)]
            for i, hist in enumerate(res):
                setattr(outdict[sample][i][0], outdict[sample][i][1], Shape(hist))
        self.nodes.Run()

    def AddSamples(self, dir, tree, fallback=None,sample_name=None):
        files = glob.glob(dir)
        if fallback is not None:
            files += glob.glob(fallback)
        seen_names = set()
        for f in files:
            testf = ROOT.TFile(f)
            if testf.Get(tree) is not None:
                if sample_name is not None:
                    name = sample_name
                else:
                    name = os.path.splitext(os.path.basename(f))[0]
                if name in seen_names:
                    print('>> Skipping %s because we already loaded it' % f)
                    continue
                seen_names.add(name)
                newname = name
                if name in self.remaps:
                    newname = self.remaps[name]
                self.trees[newname] = TTreeEvaluator(tree, f)
            testf.Close()

    def writeSubnodes(self,WriteSubnodes):
        self.WriteSubnodes = WriteSubnodes

    def AddInfo(self, file, scaleTo=None,add_name=None):

        file_ext = os.path.splitext(file)[1].lower()

        if file_ext == '.json':
            with open(file) as f:
                info = json.load(f)
        elif file_ext == '.yaml' or file_ext == '.yml':
            with open(file) as f:
                info = yaml.load(f, Loader=yaml.FullLoader)
        else:
            raise ValueError("Unsupported file format: only .json and .yaml are supported")

        lumi = info["lumi"]
        for sample in info:
            name = sample
            if add_name is not None:
                name += add_name
            if sample in self.remaps:
                name = self.remaps[sample]
            if name in self.trees:
                self.info[sample] = info[sample]

        if scaleTo is not None:
            for name, params in self.info.items():
                if 'xs' in params and 'eff' in params:
                    params["sf"] = lumi * params["xs"] / params["eff"]
                elif 'xs' in params and 'evt' in params:
                    params["sf"] = lumi * params["xs"] / params["evt"]
                else:
                    params["sf"] = 1.0
                if 'filter_efficiency' in params:
                    params["sf"] *= params["filter_efficiency"]

    def BasicFactory(self, name, sample=None, var='', sel='', factors=[], scaleToLumi=True,add_name=None):
        if sample is None:
            sample = name
        if scaleToLumi:
            myfactors = factors[:]
            if sample in self.info:
                myfactors.append(self.info[sample]['sf'])
            else:
                myfactors.append(1.0)
        if add_name is not None:
            name += add_name
            sample += add_name
        return BasicNode(name, sample, var, sel, factors=myfactors,WriteSubnodes=self.WriteSubnodes)

    def SummedFactory(self, name, samples, var='', sel='', factors=[], scaleToLumi=True,add_name=None):
        res = SummedNode(name)
        for sa in samples:
            res.AddNode(self.BasicFactory(sa, sa, var, sel, factors, scaleToLumi,add_name))
        return res


class HttWOSSSNode(BaseNode):
    def __init__(self, name, data_os, subtract_os, data_ss, subtract_ss, w_control, w_signal, w_os, w_ss, w_shape, qcd_factor=1, get_os=True, btag_extrap_num_node=None, btag_extrap_den_node=None):
        BaseNode.__init__(self, name)
        self.shape = None
        self.data_os_node = data_os
        self.subtract_os_node = subtract_os
        self.data_ss_node = data_ss
        self.subtract_ss_node = subtract_ss
        self.w_control_node = w_control
        self.w_signal_node = w_signal
        self.w_os_node = w_os
        self.w_ss_node = w_ss
        self.w_shape = w_shape
        self.qcd_factor = qcd_factor
        self.get_os = get_os
        self.btag_extrap_num_node = btag_extrap_num_node
        self.btag_extrap_den_node = btag_extrap_den_node

    def RunSelf(self):
        w_factor = self.w_os_node.shape.rate / self.w_ss_node.shape.rate
        if self.w_shape is None:
            self.shape = ((self.data_os_node.shape.rate - self.subtract_os_node.shape.rate) - (self.data_ss_node.shape.rate - self.subtract_ss_node.shape.rate) * self.qcd_factor) / (w_factor - self.qcd_factor) / self.w_control_node.shape.rate * self.w_signal_node.shape
        else:
            self.shape = ((self.data_os_node.shape.rate - self.subtract_os_node.shape.rate) - (self.data_ss_node.shape.rate - self.subtract_ss_node.shape.rate) * self.qcd_factor) / (w_factor - self.qcd_factor) / self.w_control_node.shape.rate * self.w_signal_node.shape.rate / self.w_shape.shape.rate * self.w_shape.shape
        if self.get_os:
            self.shape *= w_factor
        if self.btag_extrap_num_node is not None and self.btag_extrap_den_node is not None:
            self.shape *= self.btag_extrap_num_node.shape.rate / self.btag_extrap_den_node.shape.rate

    def Objects(self):
        return {self.name: self.shape.hist}

    def OutputPrefix(self, node=None):
        return self.name + '.subnodes'

    def SubNodes(self):
        subnodes = [self.data_os_node, self.subtract_os_node, self.data_ss_node, self.subtract_ss_node, self.w_control_node, self.w_signal_node, self.w_os_node, self.w_ss_node]
        if self.btag_extrap_num_node is not None and self.btag_extrap_den_node is not None:
            subnodes.append(self.btag_extrap_num_node)
            subnodes.append(self.btag_extrap_den_node)
        if self.w_shape is not None:
            subnodes.append(self.w_shape)
        return subnodes

    def AddRequests(self, manifest):
        for node in self.SubNodes():
            node.AddRequests(manifest)


class HttQCDNode(BaseNode):
    def __init__(self, name, data, subtract, factor, qcd_shape=None, ratio_num_node=None, ratio_den_node=None,WriteSubnodes=True, add_weight=None):
        BaseNode.__init__(self, name,WriteSubnodes)
        self.shape = None
        self.data_node = data
        self.subtract_node = subtract
        self.factor = factor
        self.ratio_num_node = ratio_num_node
        self.ratio_den_node = ratio_den_node
        self.qcd_shape = qcd_shape
        self.add_weight = add_weight

    def RunSelf(self):
        if self.subtract_node is not None:
            self.shape = self.factor * (self.data_node.shape - self.subtract_node.shape)
        else:
            self.shape = self.factor * self.data_node.shape
        if self.ratio_num_node is not None and self.ratio_den_node is not None:
            self.shape *= self.ratio_num_node.shape.rate / self.ratio_den_node.shape.rate
        if self.qcd_shape is not None:
            self.shape = self.shape.rate / self.qcd_shape.shape.rate * self.qcd_shape.shape
        if self.add_weight is not None:
            self.shape *= self.add_weight

    def Objects(self):
        return {self.name: self.shape.hist}

    def OutputPrefix(self, node=None):
        return self.name + '.subnodes'

    def SubNodes(self):
        subnodes = [self.data_node]
        if self.subtract_node is not None:
            subnodes.append(self.subtract_node)
        if self.ratio_num_node is not None and self.ratio_den_node is not None:
            subnodes.append(self.ratio_num_node)
            subnodes.append(self.ratio_den_node)
        if self.qcd_shape is not None:
            subnodes.append(self.qcd_shape)
        return subnodes

    def AddRequests(self, manifest):
        for node in self.SubNodes():
            node.AddRequests(manifest)


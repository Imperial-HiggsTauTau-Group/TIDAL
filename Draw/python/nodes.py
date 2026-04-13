from Draw.python import Analysis
import json
import re
from functools import lru_cache
import ROOT

def build_tf1_from_fit(name, fit_info):
    """Build a ROOT TF1 from fit information dictionary."""
    formula = fit_info["formula"]
    xmin = fit_info["xmin"]
    xmax = fit_info["xmax"]

    tf1 = ROOT.TF1(name, formula, xmin, xmax)
    return tf1

def _replace_x_in_formula(formula, x_expr):
    """Replace bare TF1 variable x with a ROOT-valid expression."""
    return re.sub(r'(?<![A-Za-z0-9_])x(?![A-Za-z0-9_])', f'({x_expr})', formula)


def _clip_expr(var, xmin, xmax):
    """ROOT expression for clipping a variable into [xmin, xmax]."""
    return f'TMath::Min({float(xmax)}, TMath::Max({float(xmin)}, {var}))'


def _get_classical_ff_tag(systematic, process=None):
    """
    Supported examples:
      nominal
      FF_syst:aiso
      FF_syst:nominal
      FF_syst:QCD:aiso
      FF_syst:Wjets:aiso
      FF_syst:ttbarMC:aiso
    """
    if not isinstance(systematic, str) or 'FF_syst:' not in systematic:
        return 'nominal'

    payload = systematic.split('FF_syst:', 1)[1].strip()

    if payload in ['nominal', 'aiso']:
        return payload

    if process is not None:
        if payload == f'{process}:nominal':
            return 'nominal'
        if payload == f'{process}:aiso':
            return 'aiso'

    return 'nominal'


def _get_ff_category_cuts(process, pt_col='pt_1', jpt_col='jpt_1', npreb_col='n_prebjets'):
    """
    ROOT-string equivalent of the pandas category definitions.
    Returns list of (category_name, cut_string).
    """
    if process == 'ttbarMC':
        return [
            ('jet_pt_low_inclusive',  f'(({jpt_col}/{pt_col}) >= 0 && ({jpt_col}/{pt_col}) < 1.25)'),
            ('jet_pt_med_inclusive',  f'({jpt_col}/{pt_col} >= 1.25 && {jpt_col}/{pt_col} < 1.5)'),
            ('jet_pt_high_inclusive', f'({jpt_col}/{pt_col} >= 1.5)'),
        ]
    else:
        return [
            ('jet_pt_low_0jet',  f'({npreb_col} == 0 && ({jpt_col}/{pt_col}) >= 0 && ({jpt_col}/{pt_col}) < 1.25)'),
            ('jet_pt_med_0jet',  f'({npreb_col} == 0 && ({jpt_col}/{pt_col}) >= 1.25 && ({jpt_col}/{pt_col}) < 1.5)'),
            ('jet_pt_high_0jet', f'({npreb_col} == 0 && ({jpt_col}/{pt_col}) >= 1.5)'),
            ('jet_pt_low_1jet',  f'({npreb_col} > 0 && ({jpt_col}/{pt_col}) >= 0 && ({jpt_col}/{pt_col}) < 1.25)'),
            ('jet_pt_med_1jet',  f'({npreb_col} > 0 && ({jpt_col}/{pt_col}) >= 1.25 && ({jpt_col}/{pt_col}) < 1.5)'),
            ('jet_pt_high_1jet', f'({npreb_col} > 0 && ({jpt_col}/{pt_col}) >= 1.5)'),
        ]


@lru_cache(maxsize=None)
def _load_classical_ff_for_root(path, channel, tt_pt_suffix=None):
    """
    Load classical FF JSON and convert fits into parameter-substituted ROOT formulas.

    Returns:
      classical_data[cat][tag] = {
          "formula": "...",
          "xmin": ...,
          "xmax": ...,
      }

    For tt, tt_pt_suffix should usually be "1" or "2" so that entries like
    *_pt_1 and *_pt_2 do not overwrite each other.
    """
    classical_data = {}

    with open(path, 'r') as f:
        data = json.load(f)

    for full_cat, payload in data.items():
        if channel == 'tt':
            if tt_pt_suffix is not None:
                m = re.match(rf'^(.*?)(?:_aiso2)?_pt_{tt_pt_suffix}$', full_cat)
            else:
                m = re.match(r'^(.*?)(?:_aiso2)?_pt_[12]$', full_cat)
        else:
            m = re.match(r'^(.*?)(?:_aiso2_ss)?$', full_cat)

        if m is None:
            continue

        cat = m.group(1)
        tag = 'aiso' if '_aiso' in full_cat else 'nominal'

        fit_info = payload['fit']
        xmin = fit_info['xmin']
        xmax = fit_info['xmax']

        tf1_name = f'ff_{channel}_{cat}_{tag}_{abs(hash((path, full_cat))) % 10**8}'
        tf1 = build_tf1_from_fit(tf1_name, fit_info)

        exp_formula = tf1.GetExpFormula('P')
        formula_str = exp_formula.Data() if hasattr(exp_formula, 'Data') else str(exp_formula)

        if cat not in classical_data:
            classical_data[cat] = {
                'nominal': {},
                'aiso': {},
            }

        classical_data[cat][tag] = {
            'formula': formula_str,
            'xmin': xmin,
            'xmax': xmax,
        }

    return classical_data


def _build_classical_ff_expr(path, channel, process, pt_var, tag='nominal',
                             jpt_col='jpt_pt', npreb_col='n_prebjets',
                             tt_pt_suffix=None):
    """
    Build a ROOT formula string:
      sum_over_categories( category_mask * fitted_function(clipped_pt) )
    """
    classical_data = _load_classical_ff_for_root(path, channel, tt_pt_suffix)

    pieces = []
    for cat_name, cat_cut in _get_ff_category_cuts(process, pt_col=pt_var, jpt_col=jpt_col, npreb_col=npreb_col):
        if cat_name not in classical_data:
            raise KeyError(f"Category '{cat_name}' not found in classical FF json: {path}")
        if tag not in classical_data[cat_name] or 'formula' not in classical_data[cat_name][tag]:
            raise KeyError(f"Tag '{tag}' missing for category '{cat_name}' in classical FF json: {path}")

        fit_info = classical_data[cat_name][tag]
        clipped_pt = _clip_expr(pt_var, fit_info['xmin'], fit_info['xmax'])
        ff_formula = _replace_x_in_formula(fit_info['formula'], clipped_pt)
        pieces.append(f'(({cat_cut}) * ({ff_formula}))')

    return '(' + ' + '.join(pieces) + ')'


def BuildCutString(wt='', sel='', cat='', sign='os',bkg_sel=''):
    full_selection = '(1)'
    if wt != '':
        full_selection = '(' + wt + ')'
    if sel != '':
        full_selection += '* (' + sel + ')'
    if sign != '':
        full_selection += '* (' + sign + ')'
    if bkg_sel != '':
        full_selection += '* (' + bkg_sel + ')'
    if cat != '':
        full_selection += '* (' + cat + ')'
    return full_selection


# ------------------------------------------
# ZTT, ZLL, ZL, ZJ nodes
def GetZTTNode(ana, add_name='', samples=[], plot='', wt='', sel='', cat='', z_sels={}, get_os=True):
    if get_os:
        OSSS = 'os'
    else:
        OSSS = '!os'

    full_selection = BuildCutString(wt, sel, cat, OSSS, z_sels['ztt_sel'])
    return ana.SummedFactory('ZTT'+add_name, samples, plot, full_selection)


def GetZLLNode(ana, add_name='', samples=[], plot='', wt='', sel='', cat='', z_sels={}, get_os=True):
    if get_os:
        OSSS = 'os'
    else:
        OSSS = '!os'
    full_selection = BuildCutString(wt, sel, cat, OSSS, z_sels['zll_sel'])
    return ana.SummedFactory('ZLL'+add_name, samples, plot, full_selection)


def GetZLNode(ana, add_name='', samples=[], plot='', wt='', sel='', cat='', z_sels={}, get_os=True):
    if get_os:
        OSSS = 'os'
    else:
        OSSS = '!os'
    full_selection = BuildCutString(wt, sel, cat, OSSS, z_sels['zl_sel'])
    return ana.SummedFactory('ZL'+add_name, samples, plot, full_selection)


def GetZJNode(ana, add_name='', samples=[], plot='', wt='', sel='', cat='', z_sels={}, get_os=True):
    if get_os:
        OSSS = 'os'
    else:
        OSSS = '!os'
    full_selection = BuildCutString(wt, sel, cat, OSSS, z_sels['zj_sel'])
    return ana.SummedFactory('ZJ'+add_name, samples, plot, full_selection)


def GetTTJNode(ana, add_name='', samples=[], plot='', wt='', sel='', cat='', top_sels={}, get_os=True):
  if get_os:
      OSSS = 'os'
  else:
      OSSS = '!os'
  full_selection = BuildCutString(wt, sel, cat, OSSS, top_sels['ttj_sel'])
  return ana.SummedFactory('TTJ'+add_name, samples, plot, full_selection)


def GetTTTNode(ana, add_name='', samples=[], plot='', wt='', sel='', cat='', top_sels={}, get_os=True):
  if get_os:
      OSSS = 'os'
  else:
      OSSS = '!os'
  full_selection = BuildCutString(wt, sel, cat, OSSS, top_sels['ttt_sel'])
  return ana.SummedFactory('TTT'+add_name, samples, plot, full_selection)


def GetWNode(ana, add_name='', samples=[], plot='', wt='', sel='', cat='', w_sels={}, get_os=True):
  if get_os:
      OSSS = 'os'
  else:
      OSSS = '!os'
  full_selection = BuildCutString(wt, sel, cat, OSSS, '')
  return ana.SummedFactory('W'+add_name, samples, plot, full_selection)


def GetVVTNode(ana, add_name ='', samples=[], plot='', wt='', sel='', cat='', vv_sels={}, get_os=True):
  if get_os:
      OSSS = 'os'
  else:
      OSSS = '!os'
  full_selection = BuildCutString(wt, sel, cat, OSSS, vv_sels['vvt_sel'])
  return ana.SummedFactory('VVT'+add_name, samples, plot, full_selection)


def GetVVJNode(ana, add_name ='', samples=[], plot='', wt='', sel='', cat='', vv_sels={}, get_os=True):
  if get_os:
      OSSS = 'os'
  else:
      OSSS = '!os'
  full_selection = BuildCutString(wt, sel, cat, OSSS, vv_sels['vvj_sel'])
  return ana.SummedFactory('VVJ'+add_name, samples, plot, full_selection)


def GetSubtractNode(ana, add_name, plot, plot_unmodified, wt, sel, cat_name, categories, categories_unmodified, method, qcd_factor, get_os, samples_dict, gen_sels_dict, includeW=False, w_shift=None):
    cat = categories[cat_name]
    cat_data = categories_unmodified[cat_name]
    subtract_node = Analysis.SummedNode('total_bkg'+add_name)
    if includeW:
        w_wt = wt
        w_node = GetWNode(ana, "", samples_dict['wjets_samples'], plot, wt, sel, cat, "", get_os)
        subtract_node.AddNode(w_node)
    ttt_node = GetTTTNode(ana, "", samples_dict['top_samples'], plot, wt, sel, cat, gen_sels_dict['top_sels'], get_os)
    ttj_node = GetTTJNode(ana, "", samples_dict['top_samples'], plot, wt, sel, cat, gen_sels_dict['top_sels'], get_os)
    vvt_node = GetVVTNode(ana, "", samples_dict['vv_samples'], plot, wt, sel, cat, gen_sels_dict['vv_sels'], get_os)
    vvj_node = GetVVJNode(ana, "", samples_dict['vv_samples'], plot, wt, sel, cat, gen_sels_dict['vv_sels'], get_os)
    subtract_node.AddNode(ttt_node)
    subtract_node.AddNode(ttj_node)
    subtract_node.AddNode(vvt_node)
    subtract_node.AddNode(vvj_node)

    ztt_node = GetZTTNode(ana, "", samples_dict['ztt_samples'], plot, wt, sel, cat, gen_sels_dict['z_sels'], get_os)
    subtract_node.AddNode(ztt_node)

    zl_node = GetZLNode(ana, "", samples_dict['ztt_samples']+samples_dict["zll_samples"], plot, wt, sel, cat, gen_sels_dict['z_sels'], get_os)
    zj_node = GetZJNode(ana, "", samples_dict['ztt_samples']+samples_dict["zll_samples"], plot, wt, sel, cat, gen_sels_dict['z_sels'], get_os)
    subtract_node.AddNode(zl_node)
    subtract_node.AddNode(zj_node)

    return subtract_node

def GetFakeFractionNode(ana, add_name, plot, plot_unmodified, wt, sel, cat_name, categories, categories_unmodified, method, qcd_factor, get_os, samples_dict, gen_sels_dict):
    if get_os:
        OSSS = 'os'
    else:
        OSSS = '!os'

    cat = categories[cat_name]
    cat_data = categories_unmodified[cat_name]
    
    # Get QCD yield in AR
    AR_sel =  BuildCutString(("weight"), sel, categories[cat_name], OSSS) # do this to make sure correct selection passed to data
    total_mc_bkg = GetSubtractNode(ana, '', plot, plot_unmodified, wt, sel, 'jetfake_estimate', categories, categories_unmodified, method, qcd_factor, get_os, samples_dict, gen_sels_dict, includeW=True)
    data_node = ana.SummedFactory('data_AR', samples_dict['data_samples'], plot_unmodified, AR_sel)  # data
    
    # Fraction for QCD
    QCD_frac_node = Analysis.SubtractNode('ff_qcd_frac'+add_name,
            data_node,
            total_mc_bkg)
    
    # Get W yield in AR
    w_node = GetWNode(ana, "_AR", samples_dict['wjets_samples'], plot, wt, sel, cat, "", get_os)
    vvj_node = GetVVJNode(ana, "_AR", samples_dict['vv_samples'], plot, wt, sel, cat, gen_sels_dict['vv_sels'], get_os)
    zj_node = GetZJNode(ana, "_AR", samples_dict['ztt_samples']+samples_dict["zll_samples"], plot, wt, sel, cat, gen_sels_dict['z_sels'], get_os)

    # Fraction for W fakes
    W_frac_node = Analysis.SummedNode('ff_W_frac'+add_name)
    W_frac_node.AddNode(w_node)
    W_frac_node.AddNode(vvj_node)
    W_frac_node.AddNode(zj_node)
    
    # Fraction for top
    top_frac_node = Analysis.SummedNode('ff_top_frac'+add_name)
    TTJ_frac = GetTTJNode(ana, "_AR", samples_dict['top_samples'], plot, wt, sel, cat, gen_sels_dict['top_sels'], get_os)
    top_frac_node.AddNode(TTJ_frac)

    return QCD_frac_node, W_frac_node, top_frac_node


def GetWNodeHighmT(ana, name='W', samples_dict={}, gen_sels_dict={}, plot='',plot_unmodified='', wt='', sel='', cat_name='', categories={}, categories_unmodified={}, method=1, qcd_factor=1.0, get_os=True):
    cat = categories['cat']
    cat_data = categories_unmodified['cat']
    if get_os:
        OSSS = 'os'
    else:
        OSSS = '!os'
    full_selection = BuildCutString(wt, sel, cat, OSSS, '')
    if categories['w_shape'] != '':
        shape_cat = categories['w_shape']
    else:
        shape_cat = cat
    shape_selection = BuildCutString(wt, sel, shape_cat, OSSS, '')

    if method != 2:
        raise ValueError("GetWNodeHighmT only works for method 2 (high mT W control region)")

    full_selection = BuildCutString(wt, sel, cat, OSSS)
    ss_selection = BuildCutString(wt, '', cat, '!os', '')
    os_selection = BuildCutString(wt, '', cat, 'os', '')
    control_sel = categories['w_sdb']
    w_control_full_selection = BuildCutString(wt, control_sel, cat, OSSS)
    w_control_full_selection_os_data = BuildCutString("weight", control_sel, cat_data)
    w_control_full_selection_ss_data = BuildCutString("weight", control_sel, cat_data, '!os')
    btag_extrap_num_node = None
    btag_extrap_den_node = None
    subtract_node_os = GetSubtractNode(ana, '_os', plot, plot_unmodified, wt,control_sel, 'cat', categories, categories_unmodified, method, qcd_factor, True, samples_dict, gen_sels_dict, False)
    subtract_node_ss = GetSubtractNode(ana, '_ss', plot, plot_unmodified, wt,control_sel, 'cat', categories, categories_unmodified, method, qcd_factor, False, samples_dict, gen_sels_dict, False)

    if shape_selection == full_selection:
        w_shape = None
    else:
        w_shape = ana.SummedFactory('w_shape', samples_dict['wjets_samples'], plot, shape_selection)

    w_node = Analysis.HttWOSSSNode(name,
    ana.SummedFactory('data_os', samples_dict['data_samples'], plot_unmodified, w_control_full_selection_os_data),
    subtract_node_os,
    ana.SummedFactory('data_ss', samples_dict['data_samples'], plot_unmodified, w_control_full_selection_ss_data),
    subtract_node_ss,
    ana.SummedFactory('W_cr', samples_dict['wjets_samples'], plot, w_control_full_selection),
    ana.SummedFactory('W_sr', samples_dict['wjets_samples'], plot, full_selection),
    ana.SummedFactory('W_os', samples_dict['wjets_samples'], plot, os_selection),
    ana.SummedFactory('W_ss', samples_dict['wjets_samples'], plot, ss_selection),
    w_shape,
    qcd_factor,
    get_os,
    btag_extrap_num_node,
    btag_extrap_den_node)

    return w_node

# ------------------------------------------

# ------------------------------------------
# Generating Nodes
def GenerateZLL(ana, nodename, add_name='', samples=[], plot='', wt='', sel='', cat='', z_sels={}, get_os=True, doZL=True, doZJ=True):
    if doZL:
        zl_node = GetZLNode(ana, add_name, samples, plot, wt, sel, cat, z_sels, get_os)
        ana.nodes[nodename].AddNode(zl_node)
    if doZJ:
        zj_node = GetZJNode(ana, add_name, samples, plot, wt, sel, cat, z_sels, get_os)
        ana.nodes[nodename].AddNode(zj_node)


def GenerateZTT(ana, nodename, add_name='', samples=[], plot='', wt='', sel='', cat='', z_sels={}, get_os=True):
    ztt_node = GetZTTNode(ana, add_name, samples, plot, wt, sel, cat, z_sels, get_os)
    ana.nodes[nodename].AddNode(ztt_node)


def GenerateTop(ana, nodename, add_name='', samples=[], plot='', wt='', sel='', cat='', top_sels={}, get_os=True, doTTT=True, doTTJ=True):
  wt_=wt
  if doTTT:
      ttt_node = GetTTTNode(ana, add_name, samples, plot, wt_, sel, cat, top_sels, get_os)
      ana.nodes[nodename].AddNode(ttt_node)

  if doTTJ:
      ttj_node = GetTTJNode(ana, add_name, samples, plot, wt_, sel, cat, top_sels, get_os)
      ana.nodes[nodename].AddNode(ttj_node)


def GenerateVV(ana, nodename, add_name ='', samples=[], plot='', wt='', sel='', cat='', vv_sels={}, get_os=True, doVVT=True, doVVJ=True):
  if doVVT:
      vvt_node = GetVVTNode(ana, add_name, samples, plot, wt, sel, cat, vv_sels, get_os)
      ana.nodes[nodename].AddNode(vvt_node)

  if doVVJ:
      vvj_node = GetVVJNode(ana, add_name, samples, plot, wt, sel, cat, vv_sels, get_os)
      ana.nodes[nodename].AddNode(vvj_node)


def GenerateW(ana, nodename, add_name='', samples_dict={}, gen_sels_dict={}, plot='', plot_unmodified='', wt='', sel='', cat_name='', categories={}, categories_unmodified={}, method=1, qcd_factor=1.0, get_os=True):
  w_node_name = 'W'
  if method == 2: # high mt W control region
    w_node = GetWNodeHighmT(ana, w_node_name+add_name, samples_dict, gen_sels_dict, plot, plot_unmodified, wt, sel, cat_name, categories, categories_unmodified,  method, qcd_factor, get_os)
  else:
    cat = categories['cat']
    w_node = GetWNode(ana, add_name, samples_dict['wjets_samples'], plot, wt, sel, cat, "", get_os)
  ana.nodes[nodename].AddNode(w_node)

def GenerateFakes(ana, nodename, add_name='', samples_dict={}, gen_sels_dict={}, systematic='', plot='', plot_unmodified='', wt='', sel='', cat_name='', categories={}, categories_unmodified={}, method=3, qcd_factor=1.0, get_os=True, flatten_y=False, classical_ff_cfg=None):
    shape_node = None
    if get_os:
        OSSS = "os"
    else:
        OSSS = "!os"

    cat = categories['cat']
    cat_data = categories_unmodified['cat']

    # this is for implementing the uncertainty due to the real-tau subtraction
    sub_wt=''
    if 'sub_syst' in add_name:
        if 'Up' in add_name: sub_wt='*1.1'
        if 'Down' in add_name: sub_wt='*0.9'

    ## Add estimation of fake with anti-isolated (fake) leading tau
    if 'flat_fake_sub_up' in systematic:
        syst_weight = 1.2
    elif 'flat_fake_sub_down' in systematic:
        syst_weight = 0.8
    else:
        syst_weight = 1.0

    if method in [3,4]: # fake factor for tt channel are different
        
        if method == 3: # Flat fake factor method
            data_weight = '(weight)'

            categories['qcd_sdb_cat'] = categories[cat_name]+'&&'+categories['tt_qcd_norm']
            categories_unmodified['qcd_sdb_cat'] = categories_unmodified[cat_name]+'&&'+categories_unmodified['tt_qcd_norm']

            subtract_node = GetSubtractNode(ana, '', plot, plot_unmodified, wt+sub_wt, sel+'(genPartFlav_1 != 0)', 'cat', categories, categories_unmodified, method, qcd_factor, False, samples_dict, gen_sels_dict, includeW=True)
            num_selection = BuildCutString(data_weight, sel, cat_data, '!os')
            num_node = Analysis.SubtractNode('ratio_num',
                        ana.SummedFactory('data', samples_dict['data_samples'], plot_unmodified, num_selection),
                        subtract_node)

            subtract_node = GetSubtractNode(ana, '', plot, plot_unmodified, wt+sub_wt, sel+'(genPartFlav_1 != 0)', 'qcd_sdb_cat', categories, categories_unmodified, method, qcd_factor, False, samples_dict, gen_sels_dict, includeW=True)
            den_selection = BuildCutString(data_weight, sel, categories_unmodified['qcd_sdb_cat'], '!os')
            den_node = Analysis.SubtractNode('ratio_den',
                        ana.SummedFactory('data', samples_dict['data_samples'], plot_unmodified, den_selection),
                        subtract_node)

            shape_node = None
            full_selection = BuildCutString(data_weight, sel, categories_unmodified['qcd_sdb_cat'], OSSS)
            subtract_node = GetSubtractNode(ana, '', plot, plot_unmodified, wt+sub_wt, sel+'(genPartFlav_1 != 0)', 'qcd_sdb_cat', categories, categories_unmodified, method, qcd_factor, get_os, samples_dict, gen_sels_dict, includeW=True)

            ana.nodes[nodename].AddNode(Analysis.HttQCDNode('JetFakes'+add_name,
                ana.SummedFactory('data', samples_dict['data_samples'], plot_unmodified, full_selection),
                subtract_node,
                1,
                shape_node,
                num_node,
                den_node,
                add_weight=syst_weight))
            if systematic == 'nominal':
                ## Add estimation of fakes with fake subleading tau
                categories['sublead_fakes_estimate'] = categories[cat_name]+'&&'+categories['subleadfake']
                fake_sublead_selection = BuildCutString(f"(weight)*({syst_weight})", sel, categories['sublead_fakes_estimate'], OSSS)
                fake_sublead_node = ana.SummedFactory('JetFakesSublead'+add_name, samples_dict['ztt_samples']+samples_dict['zll_samples']+samples_dict['wjets_samples']+samples_dict['vv_samples']+samples_dict['top_samples'], plot_unmodified, fake_sublead_selection)
                ana.nodes[nodename].AddNode(fake_sublead_node)

        elif method == 4:  # Full Fake Factor Method for tt channel
            if systematic == 'nominal' or 'sub_syst' in add_name or 'FF_syst:' not in systematic: ff_weight = f'(weight) * (w_FakeFactor_cmb)' # apply the nominal fake factor weight
            else: ff_weight = f'(weight) * ({systematic.replace("FF_syst:", "").replace("*ff_nom", "w_FakeFactor_cmb")})'
            # application region
            categories['qcd_ff_estimate'] = categories[cat_name]+'&&'+categories['tt_ff_AR']
            categories_unmodified['qcd_ff_estimate'] = categories_unmodified[cat_name]+'&&'+categories_unmodified['tt_ff_AR']
            ff_selection = BuildCutString(ff_weight, sel, categories['qcd_ff_estimate'], OSSS)
            # Get MC background and data yields
            mc_bkg_node = GetSubtractNode(ana, '', plot, plot_unmodified, ff_weight+sub_wt, sel+'(genPartFlav_1 != 0)', 'qcd_ff_estimate', categories, categories_unmodified, method, qcd_factor, get_os, samples_dict, gen_sels_dict, includeW=True)
            data_node = ana.SummedFactory('data', samples_dict['data_samples'], plot_unmodified, ff_selection)
            # Data - MC background yield (left with only jet fakes as doing ALL fakes MINUS non-jet fakes)
            ff_estimate = Analysis.SubtractNode('JetFakes'+add_name,
                        data_node,
                        mc_bkg_node)
            # Store FF yield
            ana.nodes[nodename].AddNode(ff_estimate)

    elif method == 7:  # Full ML Fake Factor Method for tt
        # nominal branch
        lead_branch = "BDT_FF_score_QCD_lead"
        if systematic != "nominal" and "FF_syst:" in systematic:
            lead_branch = systematic.replace("FF_syst:", "")
        ff_weight_lead = f'(weight) * ({lead_branch})' # apply the nominal fake factor weight
        # application region
        categories['qcd_ff_estimate'] = categories[cat_name]+'&&'+categories['tt_ff_AR']
        categories_unmodified['qcd_ff_estimate'] = categories_unmodified[cat_name]+'&&'+categories_unmodified['tt_ff_AR']
        ff_selection_lead = BuildCutString(ff_weight_lead, sel, categories['qcd_ff_estimate'], OSSS)
        # Get MC background and data yields
        mc_sel_lead = f"({sel}) && (genPartFlav_1 != 0) && (genPartFlav_2 == 0)" if sel else "(genPartFlav_1 != 0) && (genPartFlav_2 == 0)"
        mc_bkg_node_lead = GetSubtractNode(ana, '', plot, plot_unmodified, ff_weight_lead+sub_wt, mc_sel_lead, 'qcd_ff_estimate', categories, categories_unmodified, method, qcd_factor, get_os, samples_dict, gen_sels_dict, includeW=True)
        data_node_lead = ana.SummedFactory('data_lead', samples_dict['data_samples'], plot_unmodified, ff_selection_lead)
        # Data - MC background yield (left with only jet fakes as doing ALL fakes MINUS non-jet fakes)
        ff_estimate_lead = Analysis.SubtractNode('JetFakes'+add_name,
                    data_node_lead,
                    mc_bkg_node_lead)
        # application region for subleading tau
        # ff_weight_sublead = f'(weight) * (BDT_FF_score_sublead)'
        # categories['qcd_ff_estimate_sublead'] = categories[cat_name]+'&&'+categories['tt_ff_AR_sublead']
        # categories_unmodified['qcd_ff_estimate_sublead'] = categories_unmodified[cat_name]+'&&'+categories_unmodified['tt_ff_AR_sublead']
        # ff_selection_sublead = BuildCutString(ff_weight_sublead, sel, categories['qcd_ff_estimate_sublead'], OSSS)
        # # Get MC background and data yields
        # mc_sel_sublead = f"({sel}) && (genPartFlav_2 != 0) && (genPartFlav_1 == 0)" if sel else "(genPartFlav_2 != 0) && (genPartFlav_1 == 0)"
        # mc_bkg_node_sublead = GetSubtractNode(ana, '', plot, plot_unmodified, ff_weight_sublead+sub_wt, mc_sel_sublead, 'qcd_ff_estimate_sublead', categories, categories_unmodified, method, qcd_factor, get_os, samples_dict, gen_sels_dict, includeW=True)
        # data_node_sublead = ana.SummedFactory('data_sublead', samples_dict['data_samples'], plot_unmodified, ff_selection_sublead)
        # # Data - MC background yield (left with only jet fakes as doing ALL fakes MINUS non-jet fakes)
        # ff_estimate_sublead = Analysis.SubtractNode('JetFakesSublead'+add_name,
        #             data_node_sublead,
        #             mc_bkg_node_sublead)
        # # Store FF yield
        # # Final fake estimate nodes
        # avg_fakes = Analysis.LinearCombinationNode('AvgFakes'+add_name,
        #                 nodes=[ff_estimate_lead, ff_estimate_sublead],
        #                 coeffs=[0.5, 0.5])
        # final_fakes = Analysis.SummedNode('JetFakes'+add_name)
        # final_fakes.AddNode(avg_fakes)
        # ana.nodes[nodename].AddNode(final_fakes)
        ana.nodes[nodename].AddNode(ff_estimate_lead)
        if systematic == 'nominal':
            ## Add estimation of fakes with fake subleading tau
            categories['sublead_fakes_estimate'] = categories[cat_name]+'&&'+categories['subleadfake']
            fake_sublead_selection = BuildCutString(f"(weight)*({syst_weight})", sel, categories['sublead_fakes_estimate'], OSSS)
            fake_sublead_node = ana.SummedFactory('JetFakesSublead'+add_name, samples_dict['ztt_samples']+samples_dict['zll_samples']+samples_dict['wjets_samples']+samples_dict['vv_samples']+samples_dict['top_samples'], plot_unmodified, fake_sublead_selection)
            ana.nodes[nodename].AddNode(fake_sublead_node)

    elif method == 8:  # # Full ML Fake Factor Method for mt and et
        # use nominal FF weights, which are overwritten for some systematics
        ff_qcd_wt = '(weight) * (BDT_FF_score_QCD_sublead)' # apply the nominal fake factor weight
        ff_W_wt = '(weight) * (BDT_FF_score_Wjets_sublead)' # apply the nominal fake factor weight
        ff_top_wt = '(weight) * (BDT_FF_score_ttbarMC_sublead) * (BDT_FF_score_Wjets_sublead) * (1 / BDT_FF_score_WjetsMC_sublead)' # apply the nominal fake factor weight

        # --- Apply BDT_FakeFactors systematics by swapping the relevant FF branch
        if isinstance(systematic, str) and ("FF_syst:" in systematic):
            # Extract branch name from "FF_syst:(*BRANCH)"
            ff_branch = systematic.split("FF_syst:", 1)[1].strip()

            # Swap only the affected component, keep others nominal
            if "BDT_FF_score_QCD_" in ff_branch:
                ff_qcd_wt = f"(weight) * ({ff_branch})"

            elif "BDT_FF_score_WjetsMC_" in ff_branch:
                # affects the MC-based Wjets factor in the top formula (denominator)
                ff_top_wt = (
                    ff_top_wt.replace("BDT_FF_score_WjetsMC_sublead", ff_branch)
                )

            elif "BDT_FF_score_ttbarMC_" in ff_branch:
                # affects the top factor in the top formula (numerator)
                ff_top_wt = (
                    ff_top_wt.replace("BDT_FF_score_ttbarMC_sublead", ff_branch)
                )

            elif "BDT_FF_score_Wjets_" in ff_branch:
                # affects W contribution and also the W factor inside the top formula
                ff_W_wt = f"(weight) * ({ff_branch})"
                ff_top_wt = (
                    ff_top_wt.replace("BDT_FF_score_Wjets_sublead", ff_branch)
                )


        categories['jetfake_estimate'] = categories[cat_name]+'&&'+categories['lt_ff_AR']
        categories_unmodified['jetfake_estimate'] = categories_unmodified[cat_name]+'&&'+categories_unmodified['lt_ff_AR']

        # get fractions of different contributions in the AR
        frac_QCD, frac_W, frac_top = GetFakeFractionNode(ana, '', plot, plot_unmodified, '(weight)'+sub_wt, sel, 'jetfake_estimate', categories, categories_unmodified, method, qcd_factor, get_os, samples_dict, gen_sels_dict)

        # get MC in the AR that is not a jet fake
        ff_QCD_selection = BuildCutString(ff_qcd_wt, sel, categories['jetfake_estimate'], OSSS)
        ff_W_selection = BuildCutString(ff_W_wt, sel, categories['jetfake_estimate'], OSSS)
        ff_top_selection = BuildCutString(ff_top_wt, sel, categories['jetfake_estimate'], OSSS)

        # QCD contribution
        mc_sel_qcd = f"({sel}) && (genPartFlav_2 != 0)" if sel and sel != "(1)" else "(genPartFlav_2 != 0)"
        qcd_substract_node = GetSubtractNode(ana, '', plot, plot_unmodified, ff_qcd_wt+sub_wt, mc_sel_qcd, 'jetfake_estimate', categories, categories_unmodified, method, qcd_factor, get_os, samples_dict, gen_sels_dict, includeW=True)
        qcd_data_node = ana.SummedFactory('data_AR_qcdFF', samples_dict['data_samples'], plot_unmodified, ff_QCD_selection)
        qcd_ff_estimate = Analysis.SubtractNode('qcd_jetfakes'+add_name,
                                                qcd_data_node,
                                                qcd_substract_node)

        # W contribution
        mc_sel_w = f"({sel}) && (genPartFlav_2 != 0)" if sel and sel != "(1)" else "(genPartFlav_2 != 0)"
        W_substract_node = GetSubtractNode(ana, '', plot, plot_unmodified, ff_W_wt+sub_wt, mc_sel_w, 'jetfake_estimate', categories, categories_unmodified, method, qcd_factor, get_os, samples_dict, gen_sels_dict, includeW=True)
        W_data_node = ana.SummedFactory('data_AR_wjFF', samples_dict['data_samples'], plot_unmodified, ff_W_selection)
        W_ff_estimate = Analysis.SubtractNode('wj_jetfakes'+add_name,
                                              W_data_node,
                                              W_substract_node)

        # Top contribution
        mc_sel_top = f"({sel}) && (genPartFlav_2 != 0)" if sel and sel != "(1)" else "(genPartFlav_2 != 0)"
        Top_substract_node = GetSubtractNode(ana, '', plot, plot_unmodified, ff_top_wt+sub_wt, mc_sel_top, 'jetfake_estimate', categories, categories_unmodified, method, qcd_factor, get_os, samples_dict, gen_sels_dict, includeW=True)
        Top_data_node = ana.SummedFactory('data_AR_topFF', samples_dict['data_samples'], plot_unmodified, ff_top_selection)
        Top_ff_estimate = Analysis.SubtractNode('top_jetfakes'+add_name,
                                                Top_data_node,
                                                Top_substract_node)

        weighted_jet_fake = Analysis.FF_Node("JetFakes"+add_name, qcd_ff_estimate, W_ff_estimate, Top_ff_estimate, frac_QCD, frac_W, frac_top, flatten_y=flatten_y)
        ana.nodes[nodename].AddNode(weighted_jet_fake)
            
    elif method == 6: # lt fake factor method
        
        # use nominal FF weights, which are overwritten for some systematics
        ff_qcd_wt = f'(weight) * (FF_qcd_ipcut_nom)' # apply the nominal fake factor weight
        ff_W_wt = f'(weight) * (FF_wj_ipcut_nom)' # apply the nominal fake factor weight
        ff_top_wt = f'(weight) * (FF_mc_top_ipcut_nom)' # apply the nominal fake factor weight

        if 'FF_uct_qcd_stat:' in systematic: # statistical unct on QCD
            ff_qcd_wt = f'(weight) * ({systematic.replace("FF_uct_qcd_stat:", "")})'
        if 'FF_uct_wj_stat:' in systematic: # statistical unct on W+jets
            ff_W_wt = f'(weight) * ({systematic.replace("FF_uct_wj_stat:", "")})'
        if 'FF_uct_mc_top_stat:' in systematic: # statistical unct on top
            ff_top_wt = f'(weight) * ({systematic.replace("FF_uct_mc_top_stat:", "")})'
        if 'FF_uct_qcd_syst:' in systematic: # systematic unct on QCD (placeholder)
            ff_qcd_wt = f'(weight) * ({systematic.replace("FF_uct_qcd_syst:", "")})'
        if 'FF_uct_wj_syst:' in systematic: # systematic unct on W+jets (placeholder)
            ff_W_wt = f'(weight) * ({systematic.replace("FF_uct_wj_syst:", "")})'
        if 'FF_uct_mc_top_syst:' in systematic: # systematic unct on top (placeholder)
            ff_top_wt = f'(weight) * ({systematic.replace("FF_uct_mc_top_syst:", "")})'

        categories['jetfake_estimate'] = categories[cat_name]+'&&'+categories['lt_ff_AR']
        categories_unmodified['jetfake_estimate'] = categories_unmodified[cat_name]+'&&'+categories_unmodified['lt_ff_AR']
        
        # get fractions of different contributions in the AR
        frac_QCD, frac_W, frac_top = GetFakeFractionNode(ana, '', plot, plot_unmodified, '(weight)'+sub_wt, sel, 'jetfake_estimate', categories, categories_unmodified, method, qcd_factor, get_os, samples_dict, gen_sels_dict)
        
        # get MC in the AR that is not a jet fake
        ff_QCD_selection = BuildCutString(ff_qcd_wt, sel, categories['jetfake_estimate'], OSSS)
        ff_W_selection = BuildCutString(ff_W_wt, sel, categories['jetfake_estimate'], OSSS)
        ff_top_selection = BuildCutString(ff_top_wt, sel, categories['jetfake_estimate'], OSSS)
        
        # QCD contribution
        qcd_substract_node = GetSubtractNode(ana, '', plot, plot_unmodified, ff_qcd_wt+sub_wt, sel+'*(genPartFlav_2 != 0)', 'jetfake_estimate', categories, categories_unmodified, method, qcd_factor, get_os, samples_dict, gen_sels_dict, includeW=True)
        qcd_data_node = ana.SummedFactory('data_AR_qcdFF', samples_dict['data_samples'], plot_unmodified, ff_QCD_selection)
        qcd_ff_estimate = Analysis.SubtractNode('qcd_jetfakes'+add_name,
                        qcd_data_node,
                        qcd_substract_node)
   
        # W contribution
        W_substract_node = GetSubtractNode(ana, '', plot, plot_unmodified, ff_W_wt+sub_wt, sel+'*(genPartFlav_2 != 0)', 'jetfake_estimate', categories, categories_unmodified, method, qcd_factor, get_os, samples_dict, gen_sels_dict, includeW=True)
        W_data_node = ana.SummedFactory('data_AR_wjFF', samples_dict['data_samples'], plot_unmodified, ff_W_selection)
        W_ff_estimate = Analysis.SubtractNode('wj_jetfakes'+add_name,
                        W_data_node,
                        W_substract_node)
        
        # Top contribution
        Top_substract_node = GetSubtractNode(ana, '', plot, plot_unmodified, ff_top_wt+sub_wt, sel+'*(genPartFlav_2 != 0)', 'jetfake_estimate', categories, categories_unmodified, method, qcd_factor, get_os, samples_dict, gen_sels_dict, includeW=True)
        Top_data_node = ana.SummedFactory('data_AR_topFF', samples_dict['data_samples'], plot_unmodified, ff_top_selection)
        Top_ff_estimate = Analysis.SubtractNode('top_jetfakes'+add_name,
                        Top_data_node,
                        Top_substract_node)
        
        weighted_jet_fake = Analysis.FF_Node("JetFakes"+add_name, qcd_ff_estimate, W_ff_estimate, Top_ff_estimate, frac_QCD, frac_W, frac_top, flatten_y=flatten_y)
        ana.nodes[nodename].AddNode(weighted_jet_fake)

    elif method == 9:  # Full classical Fake Factor method for tt from jsons
        if classical_ff_cfg is None or 'tt' not in classical_ff_cfg:
            raise ValueError("method == 9 requires classical_ff_cfg['tt']")

        tt_cfg = classical_ff_cfg['tt']
        ff_tag = _get_classical_ff_tag(systematic)

        ff_expr = _build_classical_ff_expr(
            path=tt_cfg['path'],
            channel='tt',
            process=tt_cfg.get('process', 'QCD'),
            pt_var=tt_cfg.get('pt_var', 'pt_1'),
            tag=ff_tag,
            jpt_col=tt_cfg.get('jpt_col', 'jpt_pt'),
            npreb_col=tt_cfg.get('npreb_col', 'n_prebjets'),
            tt_pt_suffix=str(tt_cfg.get('tt_pt_suffix', '1'))
        )

        ff_weight = f'(weight) * ({ff_expr})'

        categories['qcd_ff_estimate'] = categories[cat_name] + '&&' + categories[tt_cfg.get('ar_key', 'tt_ff_AR')]
        categories_unmodified['qcd_ff_estimate'] = categories_unmodified[cat_name] + '&&' + categories_unmodified[tt_cfg.get('ar_key', 'tt_ff_AR')]
        ff_selection = BuildCutString(ff_weight, sel, categories['qcd_ff_estimate'], OSSS)

        # match method-4 behavior
        mc_sel_lead = (sel + '*(genPartFlav_1 != 0)') if sel else '(genPartFlav_1 != 0)'

        mc_bkg_node_lead = GetSubtractNode(
            ana, '', plot, plot_unmodified,
            ff_weight + sub_wt,
            mc_sel_lead,
            'qcd_ff_estimate',
            categories, categories_unmodified,
            method, qcd_factor, get_os,
            samples_dict, gen_sels_dict,
            includeW=True
        )

        data_node_lead = ana.SummedFactory(
            'data_classicalFF_tt',
            samples_dict['data_samples'],
            plot_unmodified,
            ff_selection
        )

        ff_estimate_lead = Analysis.SubtractNode(
            'JetFakes' + add_name,
            data_node_lead,
            mc_bkg_node_lead
        )
        ana.nodes[nodename].AddNode(ff_estimate_lead)

        # keep sublead as a separate contribution (same pattern as method 3/4)
        if systematic == 'nominal':
            sublead_key = tt_cfg.get('sublead_key', 'subleadfake')
            categories['sublead_fakes_estimate'] = categories[cat_name] + '&&' + categories[sublead_key]
            fake_sublead_selection = BuildCutString(
                f"(weight)*({syst_weight})",
                sel,
                categories['sublead_fakes_estimate'],
                OSSS
            )
            fake_sublead_node = ana.SummedFactory(
                'JetFakesSublead' + add_name,
                samples_dict['ztt_samples'] + samples_dict['zll_samples'] + samples_dict['wjets_samples'] + samples_dict['vv_samples'] + samples_dict['top_samples'],
                plot_unmodified,
                fake_sublead_selection
            )
            ana.nodes[nodename].AddNode(fake_sublead_node)

    elif method == 10:  # Full classical Fake Factor method for lt from jsons
        if classical_ff_cfg is None or 'lt' not in classical_ff_cfg:
            raise ValueError("method == 10 requires classical_ff_cfg['lt']")

        lt_cfg = classical_ff_cfg['lt']
        lt_channel = lt_cfg.get('channel', 'mt')
        pt_var = lt_cfg.get('pt_var', 'pt_2')
        jpt_col = lt_cfg.get('jpt_col', 'jpt_pt')
        npreb_col = lt_cfg.get('npreb_col', 'n_prebjets')
        ar_key = lt_cfg.get('ar_key', 'lt_ff_AR')
        paths = lt_cfg['paths']

        ff_qcd_expr = _build_classical_ff_expr(
            path=paths['QCD'],
            channel=lt_channel,
            process='QCD',
            pt_var=pt_var,
            tag=_get_classical_ff_tag(systematic, process='QCD'),
            jpt_col=jpt_col,
            npreb_col=npreb_col
        )
        ff_w_expr = _build_classical_ff_expr(
            path=paths['Wjets'],
            channel=lt_channel,
            process='Wjets',
            pt_var=pt_var,
            tag=_get_classical_ff_tag(systematic, process='Wjets'),
            jpt_col=jpt_col,
            npreb_col=npreb_col
        )
        ff_top_expr = _build_classical_ff_expr(
            path=paths['ttbarMC'],
            channel=lt_channel,
            process='ttbarMC',
            pt_var=pt_var,
            tag=_get_classical_ff_tag(systematic, process='ttbarMC'),
            jpt_col=jpt_col,
            npreb_col=npreb_col
        )

        ff_qcd_wt = f'(weight) * ({ff_qcd_expr})'
        ff_W_wt = f'(weight) * ({ff_w_expr})'
        ff_top_wt = f'(weight) * ({ff_top_expr})'

        categories['jetfake_estimate'] = categories[cat_name] + '&&' + categories[ar_key]
        categories_unmodified['jetfake_estimate'] = categories_unmodified[cat_name] + '&&' + categories_unmodified[ar_key]

        frac_QCD, frac_W, frac_top = GetFakeFractionNode(
            ana, '', plot, plot_unmodified,
            '(weight)' + sub_wt,
            sel,
            'jetfake_estimate',
            categories, categories_unmodified,
            method, qcd_factor, get_os,
            samples_dict, gen_sels_dict
        )

        ff_QCD_selection = BuildCutString(ff_qcd_wt, sel, categories['jetfake_estimate'], OSSS)
        ff_W_selection = BuildCutString(ff_W_wt, sel, categories['jetfake_estimate'], OSSS)
        ff_top_selection = BuildCutString(ff_top_wt, sel, categories['jetfake_estimate'], OSSS)

        mc_sel_lt = f"({sel}) && (genPartFlav_2 != 0)" if sel and sel != "(1)" else "(genPartFlav_2 != 0)"

        qcd_substract_node = GetSubtractNode(
            ana, '', plot, plot_unmodified,
            ff_qcd_wt + sub_wt,
            mc_sel_lt,
            'jetfake_estimate',
            categories, categories_unmodified,
            method, qcd_factor, get_os,
            samples_dict, gen_sels_dict,
            includeW=True
        )
        qcd_data_node = ana.SummedFactory('data_AR_qcdClassicalFF', samples_dict['data_samples'], plot_unmodified, ff_QCD_selection)
        qcd_ff_estimate = Analysis.SubtractNode('qcd_jetfakes' + add_name, qcd_data_node, qcd_substract_node)

        W_substract_node = GetSubtractNode(
            ana, '', plot, plot_unmodified,
            ff_W_wt + sub_wt,
            mc_sel_lt,
            'jetfake_estimate',
            categories, categories_unmodified,
            method, qcd_factor, get_os,
            samples_dict, gen_sels_dict,
            includeW=True
        )
        W_data_node = ana.SummedFactory('data_AR_wjClassicalFF', samples_dict['data_samples'], plot_unmodified, ff_W_selection)
        W_ff_estimate = Analysis.SubtractNode('wj_jetfakes' + add_name, W_data_node, W_substract_node)

        Top_substract_node = GetSubtractNode(
            ana, '', plot, plot_unmodified,
            ff_top_wt + sub_wt,
            mc_sel_lt,
            'jetfake_estimate',
            categories, categories_unmodified,
            method, qcd_factor, get_os,
            samples_dict, gen_sels_dict,
            includeW=True
        )
        Top_data_node = ana.SummedFactory('data_AR_topClassicalFF', samples_dict['data_samples'], plot_unmodified, ff_top_selection)
        Top_ff_estimate = Analysis.SubtractNode('top_jetfakes' + add_name, Top_data_node, Top_substract_node)

        weighted_jet_fake = Analysis.FF_Node(
            "JetFakes" + add_name,
            qcd_ff_estimate,
            W_ff_estimate,
            Top_ff_estimate,
            frac_QCD,
            frac_W,
            frac_top,
            flatten_y=flatten_y
        )
        ana.nodes[nodename].AddNode(weighted_jet_fake)


def GenerateQCD(ana, nodename, add_name='', samples_dict={}, gen_sels_dict={}, systematic='', plot='', plot_unmodified='', wt='', sel='', cat_name='', categories={}, categories_unmodified={}, method=1, qcd_factor=1.0, get_os=True,w_shift=None):
    shape_node = None
    if get_os:
        OSSS = "os"
    else:
        OSSS = "!os"

    cat = categories['cat']
    cat_data = categories_unmodified['cat']

    if method in [1,2,5]:
        sub_shift='*1.0'
        if 'qcd_sub_up' in systematic:
            sub_shift = '*1.1'
        if 'qcd_sub_down' in systematic:
            sub_shift = '*0.9'

        # TODO: Weight for data
        data_weight = '(weight)'
        full_selection = BuildCutString(data_weight, sel, cat_data, '!os')

        if method in [5]:
            categories['qcd_loose_shape_cat'] = categories[cat_name]+'&&'+categories['qcd_loose_shape']
            categories_unmodified['qcd_loose_shape_cat'] = categories_unmodified[cat_name]+'&&'+categories_unmodified['qcd_loose_shape']
            shape_selection = BuildCutString(data_weight, sel, categories_unmodified['qcd_loose_shape_cat'], '!os')
            subtract_node = GetSubtractNode(ana , '', plot, plot_unmodified, wt+sub_shift, sel, 'qcd_loose_shape_cat', categories, categories_unmodified, method, qcd_factor, False, samples_dict, gen_sels_dict, includeW=True, w_shift=w_shift)
            shape_node = Analysis.SubtractNode('shape', ana.SummedFactory('data_ss', samples_dict['data_samples'], plot_unmodified, shape_selection), subtract_node)

        # method 5 subtract node gets overrided here as it's only used for the shape node
        subtract_node = GetSubtractNode(ana , '', plot, plot_unmodified, wt+sub_shift, sel, 'cat', categories, categories_unmodified, method, qcd_factor, False, samples_dict, gen_sels_dict, includeW=True, w_shift=w_shift)

        if get_os:
            qcd_ratio = qcd_factor
        else:
            qcd_ratio = 1.0

        ana.nodes[nodename].AddNode(Analysis.HttQCDNode('QCD'+add_name,
          ana.SummedFactory('data_ss', samples_dict['data_samples'], plot_unmodified, full_selection),
          subtract_node,
          qcd_ratio,
          shape_node))


def GenerateReweightedCPSignal(ana, nodename='', add_name='', samples={}, masses=[], plot='', wt='', sel='', cat='', get_os=True):
    #TODO: we probably want to reweight pT distributions to NNLOPS to take into account quark mass effects (or use similar tool)
    weights = {"sm": "wt_cp_sm", "ps": "wt_cp_ps", "mm": "wt_cp_mm", "flat": "1.0"}
    if get_os:
        OSSS = 'os'
    else:
        OSSS = '!os'

    for key, sample in samples.items():
        non_cp = True
        for name in weights:
            for mass in masses:
                if key.split("_")[1] == name:
                    non_cp=False
                    weight=wt+"*"+weights[name]
                    valid_spinner = '&& (' + weights[name] + '>= 0)' # avoid issues with Nan broadcast to -9999
                    # this part takes care of scaling to the LHE weight to take care of the CP in production
                    # currently reweighting the different samples is not fully implemted
                    # this is because the weights appear to change the cross sections slightly which isn't expected (need to check gridpack setup to understand why)
                    # if we want to combine SM, CPodd and MM samples then we would also need to modify the params file to ensure we don't triple count the events
                    # if 'prod_sm' in key:   weight+='*'+'LHEReweightingWeight_SM'
                    # elif 'prod_ps' in key: weight+='*'+'LHEReweightingWeight_PS'
                    # elif 'prod_mm' in key: weight+='*'+'LHEReweightingWeight_MM'
                    full_selection = BuildCutString(weight, sel + valid_spinner, cat, OSSS)
                    name = key

                    sample_names=[]
                    if isinstance(samples[key], (list,)):
                      for i in samples[key]:
                        sample_names.append(i.replace('*',mass))
                    else: sample_names = [samples[key].replace('*',mass)]
                    ana.nodes[nodename].AddNode(ana.SummedFactory(key.replace('*',mass)+add_name, sample_names, plot, full_selection))
            if non_cp:
                 full_selection = BuildCutString(wt, sel, cat, OSSS)
                 name = key

                 sample_names=[]
                 if isinstance(samples[key], (list,)):
                   for i in samples[key]:
                     sample_names.append(i.replace('*',mass))
                 else: sample_names = [samples[key].replace('*',mass)]
                 ana.nodes[nodename].AddNode(ana.SummedFactory(key.replace('*',mass)+add_name, sample_names, plot, full_selection))


def GenerateMSSMReweightedSignal(ana, nodename='', add_name='', samples={}, masses=[], plot='', wt='', sel='', cat='', get_os=True):
    """
    For MSSM signals: plain loop over masses with no reweighting (TODO: Irene).
    """
    OSSS = 'os' if get_os else '!os'

    for key, sample in samples.items():
        if '*' not in key:
            print(f"Warning: sample key '{key}' does not contain '*', skipping mass loop and using key as is.")
            full_selection = BuildCutString(wt, sel, cat, OSSS)
            print("Full selection for MSSM signal: ", full_selection)
            proc_name = key + add_name  
            sample_names = sample if isinstance(sample, list) else [sample] 
            ana.nodes[nodename].AddNode(
                ana.SummedFactory(proc_name, sample_names, plot, full_selection)
            )
        else:
            for mass in masses:
                full_selection = BuildCutString(wt, sel, cat, OSSS)
                print("Full selection for MSSM signal: ", full_selection)
                proc_name = key.replace('*', mass) + add_name

                sample_names = []
                if isinstance(sample, list):
                    for s in sample:
                        sample_names.append(s.replace('*', mass))
                else:
                    sample_names = [sample.replace('*', mass)]

                ana.nodes[nodename].AddNode(
                    ana.SummedFactory(proc_name, sample_names, plot, full_selection)
                )

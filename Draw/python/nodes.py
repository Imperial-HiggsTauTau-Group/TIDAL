from Draw.python import Analysis


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

def GenerateFakes(ana, nodename, add_name='', samples_dict={}, gen_sels_dict={}, systematic='', plot='', plot_unmodified='', wt='', sel='', cat_name='', categories={}, categories_unmodified={}, method=3, qcd_factor=1.0, get_os=True, flatten_y=False):
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

            subtract_node = GetSubtractNode(ana, '', plot, plot_unmodified, wt+sub_wt, sel+'*(genPartFlav_1 != 0)', 'cat', categories, categories_unmodified, method, qcd_factor, False, samples_dict, gen_sels_dict, includeW=True)
            num_selection = BuildCutString(data_weight, sel, cat_data, '!os')
            num_node = Analysis.SubtractNode('ratio_num',
                        ana.SummedFactory('data', samples_dict['data_samples'], plot_unmodified, num_selection),
                        subtract_node)

            subtract_node = GetSubtractNode(ana, '', plot, plot_unmodified, wt+sub_wt, sel+'*(genPartFlav_1 != 0)', 'qcd_sdb_cat', categories, categories_unmodified, method, qcd_factor, False, samples_dict, gen_sels_dict, includeW=True)
            den_selection = BuildCutString(data_weight, sel, categories_unmodified['qcd_sdb_cat'], '!os')
            den_node = Analysis.SubtractNode('ratio_den',
                        ana.SummedFactory('data', samples_dict['data_samples'], plot_unmodified, den_selection),
                        subtract_node)

            shape_node = None
            full_selection = BuildCutString(data_weight, sel, categories_unmodified['qcd_sdb_cat'], OSSS)
            subtract_node = GetSubtractNode(ana, '', plot, plot_unmodified, wt+sub_wt, sel+'*(genPartFlav_1 != 0)', 'qcd_sdb_cat', categories, categories_unmodified, method, qcd_factor, get_os, samples_dict, gen_sels_dict, includeW=True)

            ana.nodes[nodename].AddNode(Analysis.HttQCDNode('JetFakes'+add_name,
                ana.SummedFactory('data', samples_dict['data_samples'], plot_unmodified, full_selection),
                subtract_node,
                1,
                shape_node,
                num_node,
                den_node,
                add_weight=syst_weight))

        elif method == 4:  # Full Fake Factor Method for tt channel
            if systematic == 'nominal' or 'sub_syst' in add_name or 'FF_syst:' not in systematic: ff_weight = f'(weight) * (w_FakeFactor_cmb)' # apply the nominal fake factor weight
            else: ff_weight = f'(weight) * ({systematic.replace("FF_syst:", "").replace("*ff_nom", "w_FakeFactor_cmb")})'
            # application region
            categories['qcd_ff_estimate'] = categories[cat_name]+'&&'+categories['tt_ff_AR']
            categories_unmodified['qcd_ff_estimate'] = categories_unmodified[cat_name]+'&&'+categories_unmodified['tt_ff_AR']
            ff_selection = BuildCutString(ff_weight, sel, categories['qcd_ff_estimate'], OSSS)
            # Get MC background and data yields
            mc_bkg_node = GetSubtractNode(ana, '', plot, plot_unmodified, ff_weight+sub_wt, sel+'*(genPartFlav_1 != 0)', 'qcd_ff_estimate', categories, categories_unmodified, method, qcd_factor, get_os, samples_dict, gen_sels_dict, includeW=True)
            data_node = ana.SummedFactory('data', samples_dict['data_samples'], plot_unmodified, ff_selection)
            # Data - MC background yield (left with only jet fakes as doing ALL fakes MINUS non-jet fakes)
            ff_estimate = Analysis.SubtractNode('JetFakes'+add_name,
                        data_node,
                        mc_bkg_node)
            # Store FF yield
            ana.nodes[nodename].AddNode(ff_estimate)

        if systematic == 'nominal':
            ## Add estimation of fakes with fake subleading tau
            categories['sublead_fakes_estimate'] = categories[cat_name]+'&&'+categories['subleadfake']
            fake_sublead_selection = BuildCutString(f"(weight)*({syst_weight})", sel, categories['sublead_fakes_estimate'], OSSS)
            fake_sublead_node = ana.SummedFactory('JetFakesSublead'+add_name, samples_dict['ztt_samples']+samples_dict['zll_samples']+samples_dict['wjets_samples']+samples_dict['vv_samples']+samples_dict['top_samples'], plot_unmodified, fake_sublead_selection)
            ana.nodes[nodename].AddNode(fake_sublead_node)
            
    if method == 6: # lt fake factor method
        
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

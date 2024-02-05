from ppta_dr3_models import *
from scipy.signal import correlate
import pickle
import subprocess
import os
import sys

psrnames = ['J0030+0451',
            'J0125-2327',
            'J0437-4715',
            'J0613-0200',
            'J0614-3329',
            'J0711-6830', 
            'J0900-3144',
            'J1017-7156',
            'J1022+1001',
            'J1024-0719', 
            'J1045-4509',
            'J1125-6014',
            'J1446-4701', 
            'J1545-4550',
            'J1600-3053',
            'J1603-7202',
            'J1643-1224',
            'J1713+0747',
            'J1730-2304',
            #'J1741+1351', #trash
            'J1744-1134',
            #'J1824-2452A', #also trash
            'J1857+0943',
            'J1832-0836',
            'J1902-5105',
            'J1909-3744',
            'J1933-6211',
            'J1939+2134',
            'J2124-3358',
            'J2129-5721',
            'J2145-0750',
            'J2241-5236']



#commented lines below show pulsars whose group noise posteriors were constrained
#but have group noise not sufficiently favoured (log BF < 0.5) in model slection
#selected groups were those with constrained posteriors and with spectral indices
#with posterior maximum > 0 (posterior maximum gamma of 0 indicates unmodelled white
#noise power)

psr_groupecorrlist_dict_dr2 = {'J0437-4715': ['CASPSR_40CM', 'PDFB_20CM'],
                               #'J0613-0200': ['PDFB_40CM'],
                               #'J1017-7156': ['PDFB_20CM', 'PDFB_40CM']
                               #'J1600-3053': ['CASPSR_40CM', 'PDFB_20CM', 'PDFB_10CM'],
                               #'J1643-1224': ['PDFB_10CM', 'PDFB_40CM']#,
                               'J1713+0747': ['WBCORR_10CM', 'CPSR2_20CM'], #PDFB1_early_20CM,
                               'J1909-3744': ['CPSR2_50CM', 'CASPSR_40CM', 'PDFB1_1433', 'PDFB1_early_20CM']
}

psr_groupecorrlist_dict_uwl = {'J0437-4715': ['UWL_PDFB4_20CM', 'UWL_sbA', 'UWL_sbG'],
                               #'J0613-0200': ['UWL_sbG'],
                               'J1017-7156': ['UWL_sbA', 'UWL_sbD'], #UNCOMMENT ME
                               'J1022+1001': ['UWL_sbE', 'UWL_sbH'],
                               #'J1600-3053': ['UWL_sbA'],
                               #'J1643-1224': ['UWL_sbA', 'UWL_sbD', 'UWL_sbE', 'UWL_sbF', 'UWL_sbH'],
                               'J1713+0747': ['UWL_sbA', 'UWL_sbE', 'UWL_sbF']
                               #'J1744-1134': ['UWL_PDFB4_10CM',  'UWL_sbE'],
                               #'J1909-3744': ['UWL_PDFB4_10CM', 'UWL_PDFB4_20CM', 'UWL_sbG', 'UWL_sbH'],
                               #'J2241-5236': ['UWL_PDFB4_20CM', 'UWL_sbD']
}

psr_groupecorrlist_dict_all = {}
for psr, item in psr_groupecorrlist_dict_dr2.items():
    if psr in psr_groupecorrlist_dict_uwl:
        psr_groupecorrlist_dict_all[psr] = psr_groupecorrlist_dict_uwl[psr] + item
    else:
        psr_groupecorrlist_dict_all[psr] = item

for psr, item in psr_groupecorrlist_dict_uwl.items():
    if psr in psr_groupecorrlist_dict_dr2:
        continue
    else:
        psr_groupecorrlist_dict_all[psr] = item
        
psr_groupecorr_dict_dict = {'dr2': psr_groupecorrlist_dict_dr2, 'uwl': psr_groupecorrlist_dict_uwl, 'all': psr_groupecorrlist_dict_all}


psr_groupnoiselist_dict_dr2 = {'J0437-4715': ['CASPSR_40CM'],#,-->taken out after fixing number of harmonics for group noise 'PDFB_20CM'],
                               #'J0613-0200': ['PDFB_40CM'],
                               #'J1017-7156': ['PDFB_20CM', 'PDFB_40CM']
                               #'J1600-3053': ['CASPSR_40CM', 'PDFB_20CM', 'PDFB_10CM'],
                               #'J1643-1224': ['PDFB_10CM', 'PDFB_40CM']#,
                               'J1713+0747': ['WBCORR_10CM'], #PDFB1_early_20CM, -->taken out after fixing number of harmonics for group noise'CPSR2_20CM' 
                               'J1909-3744': ['CPSR2_50CM'] #-->taken out after fixing number of harmonics for group noise#'CASPSR_40CM', 'PDFB1_1433', 'PDFB1_early_20CM'
}

psr_groupnoiselist_dict_uwl = {'J0437-4715': ['UWL_PDFB4_20CM', 'UWL_sbA', 'UWL_sbG'],
                               #'J0613-0200': ['UWL_sbG'],
                               'J1017-7156': ['UWL_sbA', 'UWL_sbD'], #UNCOMMENT ME
                               'J1022+1001': ['UWL_sbE', 'UWL_sbH'],
                               #'J1600-3053': ['UWL_sbA'],
                               #'J1643-1224': ['UWL_sbA', 'UWL_sbD', 'UWL_sbE', 'UWL_sbF', 'UWL_sbH'],
                               'J1713+0747': ['UWL_sbA', 'UWL_sbE', 'UWL_sbF']
                               #'J1744-1134': ['UWL_PDFB4_10CM',  'UWL_sbE'],
                               #'J1909-3744': ['UWL_PDFB4_10CM', 'UWL_PDFB4_20CM', 'UWL_sbG', 'UWL_sbH'],
                               #'J2241-5236': ['UWL_PDFB4_20CM', 'UWL_sbD']
}

psr_groupnoiselist_dict_all = {}
for psr, item in psr_groupnoiselist_dict_dr2.items():
    if psr in psr_groupnoiselist_dict_uwl:
        psr_groupnoiselist_dict_all[psr] = psr_groupnoiselist_dict_uwl[psr] + item
    else:
        psr_groupnoiselist_dict_all[psr] = item

for psr, item in psr_groupnoiselist_dict_uwl.items():
    if psr in psr_groupnoiselist_dict_dr2:
        continue
    psr_groupnoiselist_dict_all[psr] = item

psr_groupnoise_dict_dict = {'dr2': psr_groupnoiselist_dict_dr2, 'uwl': psr_groupnoiselist_dict_uwl, 'all': psr_groupnoiselist_dict_all}

# define useful selections
by_backend = selections.Selection(selections.by_backend)
no_selection = selections.Selection(selections.no_selection)
low_freq = selections.Selection(low_frequencies)
mid_freq = selections.Selection(mid_frequencies)
high_freq = selections.Selection(high_frequencies)
high_mid_freq = selections.Selection(high_mid_frequencies)
ecorr_selection = selections.Selection(band_split)
by_group_all = selections.Selection(sel_by_group)
by_group_uwl = selections.Selection(sel_by_group_factory(['UWL_sbA',
                                                          'UWL_sbB',
                                                          'UWL_sbC',
                                                          'UWL_sbD',
                                                          'UWL_sbE',
                                                          'UWL_sbF',
                                                          'UWL_sbG',
                                                          'UWL_sbH'])._sel_by_group)
uwl_selection = selections.Selection(uwl_all)
not_uwl_selection = selections.Selection(not_uwl_all)
global_ecorr_selection = selections.Selection(global_ecorr)




def update_tidy_noisedict(noisefiles, dir, dirnoise):
    
        """
        Takes in noisefiles and constructs a noisedict,
        taking care of duplicates and differing labels
        """
        noisedict = {}
	for nf in noisefiles:
            if 'dr2' in nf:
                continue
            with open(nf, 'r') as f:
                noisedict.update(json.load(f))
            #rename 'n_earth' in single-pulsar noise dict to '<psrname>_n_earth' for common runs
            if 'n_earth' in noisedict.keys():
                    noisedict['{}_n_earth'.format(os.path.basename(nf).split('_')[0])] = noisedict['n_earth']
                    del(noisedict['n_earth'])

        #rename ecorr_all grabbed from UWL runs to ecorr_all_uwl: a UWL-specific global ecorr
        keys = np.array([str(key) for key in noisedict.keys()]).squeeze()
        for key in keys:
            if "basis_ecorr_all" in key and not 'uwl' in key and not 'dr2' in key:
               new_key = key.replace("basis_ecorr_all", "basis_ecorr_all_uwl")
               noisedict[new_key] = noisedict[key]


        ### NOW DO DR2 NOISE, TO RE-NAME ECORR_ALL
        for nf in noisefiles:
            if 'uwl' in nf:
                continue
            with open(nf, 'r') as f:
                noisedict.update(json.load(f))

            #remove n_earth parameter from noisedict
            if 'n_earth' in noisedict.keys():
                del(noisedict['n_earth'])

        #rename ecorr_all grabbed from DR2 runs to ecorr_all_uwl: a DR2-specific global ecorr
        keys = np.array([str(key) for key in noisedict.keys()]).squeeze()
        for key in keys:
            if "basis_ecorr_all" in key and not 'uwl' in key and not 'dr2' in key:
               new_key = key.replace("basis_ecorr_all", "basis_ecorr_all_dr2")
               noisedict[new_key] = noisedict[key]

        #NOW FUSS OVER GROUP TERM ECORR NAMING:
        for key in keys:
            if 'ecorr' not in key:
                continue
            for grouppsr, grouplist_ in psr_groupecorr_dict_dict[dir].items():

                if grouppsr not in key:

                    continue
                for group_ in grouplist_:
                    if group_ not in key:
                        continue
                    group_done = False
                    if f'{grouppsr}_basis_ecorr_{group_}' in key:
                        new_key = key.replace(f'basis_ecorr_{group_}', f'basis_ecorr_group_{group_}')
                        noisedict[new_key] = noisedict[key]
                        print(f'renamed {key} to {new_key}')
                        del noisedict[key]
                        group_done = True
                    #check if it's already been converted to new name
                    elif f'{grouppsr}_basis_ecorr_group_{group_}' in key:
                        group_done = True
                    #check if group isn't added in - e.g. it is not in the `all` noise files
                    if not group_done:
                        print('COULD NOT FIND {} TERM FOR {} (KEY IS {})'.format(group_, grouppsr, key))
                        raise ValueError
                        print('looking in uwl noise file for {} ecorr term'.format(group_))
                        uwl_gnf = sorted(glob.glob(dirnoise+f'/uwl/noiseFiles_maxlike/{grouppsr}*.json'))[0]
                        with open(uwl_gnf, 'r') as uwl_nf:
                            #uwl noisedict
                            uwl_nd = json.load(uwl_nf)
                            keys__ = np.array([str(key) for key in uwl_nd.keys()]).squeeze()
                            for k in keys__:
                                if f'{grouppsr}_basis_ecorr_{group_}' in k or f'{grouppsr}_basis_ecorr_group_{group_}' in k:
                                    noisedict[k.replace(f'basis_ecorr_{group_}', f'basis_ecorr_group_{group_}')] = uwl_nd[k]

        return noisedict


def get_3sig_noisedict(dir, dirnoise):
        """
        function to read in the two-sided 3-sigma bounds on the single-pulsar noise terms, for prior clipping
        """
        if dir == 'all':
            if len(glob.glob(dirnoise + '/all/noiseFiles/3sig/*p0015s.json')) > 0:
                noisedict_p0015 = {}
                noisefiles_p0015 = sorted(glob.glob(dirnoise+'/all/noiseFiles/3sig/*p0015s.json'))
                for nf in noisefiles_p0015:
                    with open(nf, 'r') as f:
                        noisedict_p0015.update(json.load(f))
                        psrname = nf.split('_')[0].split('/')[-1]
                        print('noise dict pulsar name {}'.format(psrname))
                        if "n_earth" in noisedict_p0015.keys():
                            noisedict_p0015[psrname + '_n_earth'] = noisedict_p0015['n_earth']
                        del noisedict_p0015['n_earth']
                noisedict_p9985 = {}
                noisefiles_p9985 = sorted(glob.glob(dirnoise+'/all/noiseFiles/3sig/*p9985s.json'))
                for nf in noisefiles_p9985:
                    with open(nf, 'r') as f:
                        noisedict_p9985.update(json.load(f))
                        psrname = nf.split('_')[0].split('/')[-1]
                        if 'n_earth' in noisedict_p9985.keys():
                            noisedict_p9985[psrname + '_n_earth'] = noisedict_p9985['n_earth']
                        del noisedict_p9985['n_earth']

            else:
                # take 3-sigma bounds from dr2 only
                print("Couldn't find 3-sigma noisefiles!")
                dir_3sig = os.path.abspath("./data/")
                noisefiles_p9985_dr2  = sorted(glob.glob(dir_3sig + '/dr2/noiseFiles/3sig/*p9985s.json'))
                noisefiles_p0015_dr2  = sorted(glob.glob(dir_3sig + '/dr2/noiseFiles/3sig/*p0015s.json'))
                noisefiles_p9985_uwl  = sorted(glob.glob(dir_3sig + '/uwl/noiseFiles/3sig/*p9985s.json'))
                noisefiles_p0015_uwl  = sorted(glob.glob(dir_3sig + '/uwl/noiseFiles/3sig/*p0015s.json'))    
                #print(noisefiles)

                noisedict_p9985_dr2 = {}
                for nf in noisefiles_p9985_dr2:
                    with open(nf, 'r') as f:
                        #print(nf)
                        noisedict_p9985_dr2.update(json.load(f))
                noisedict_p0015_dr2 = {}
                for nf in noisefiles_p0015_dr2:
                    with open(nf, 'r') as f:

                        #print(nf)
                        noisedict_p0015_dr2.update(json.load(f))
                noisedict_p9985_uwl = {}
                for nf in noisefiles_p9985:
                    with open(nf, 'r') as f:
                        #print(nf)
                        noisedict_p9985.update(json.load(f))
                noisedict_p0015_uwl = {}
                for nf in noisefiles_p0015:
                    with open(nf, 'r') as f:
                        #print(nf)
                        noisedict_p0015.update(json.load(f))

                #go through each noise term in DR2 and add it into the overall p0015 dict
                #if the noise term is also in UWL, add in the minimum of the DR2 or UWL values
                noisedict_p0015 = {}
                for key, item in noisedict_p0015_dr2.items():
                    if key in noisedict_p0015_uwl:
                        val_ = np.min(noisedict_p0015_uwl[key], item)
                    else:
                        val_ = item
                    noisedict_p0015[key] = val_
                #add in any leftover terms that weren't in DR2 dict:
                for key, item in noisedict_p0015_uwl.items():
                    #if they're not in the overall dict then they weren't in DR2
                    if not key in noisedict_p0015:
                        noisedict_p0015[key] = item

                #go through each noise term in DR2 and add it into the overall 3sig_max dict
                #if the noise term is also in UWL, add in the max of the DR2 or UWL values
                noisedict_p9985 = {}
                for key, item in noisedict_p9985_dr2.items():
                    if key in noisedict_p9985_uwl:
                        val_ = np.max(noisedict_p9985_uwl[key], item)
                    else:
                        val_ = item
                    noisedict_p9985[key] = val_
                #add in any leftover terms that weren't in DR2 dict:
                for key, item in noisedict_p9985_uwl.items():
                    #if they're not in the overall dict then they weren't in DR2
                    if not key in noisedict_p9985:
                        noisedict_p9985[key] = item

        else:
            dir_3sig = os.path.abspath("./data/")
            noisefiles_p9985  = sorted(glob.glob(dir_3sig + '/' + str(dir) + '/noiseFiles/3sig/*p9985s.json'))
            noisefiles_p0015  = sorted(glob.glob(dir_3sig + '/' + str(dir) + '/noiseFiles/3sig/*p0015s.json'))
            noisedict_p0015 = {}
            for nf in noisefiles_p0015:
                with open(nf, 'r') as f:
                    #print(nf)
                    noisedict_p0015.update(json.load(f))
                if 'n_earth' in noisedict_p0015.keys():
                    noisedict_p0015['{}_n_earth'.format(os.path.basename(nf).split('_')[0])] = noisedict_p0015['n_earth']
                    del(noisedict_p0015['n_earth'])
            noisedict_p9985 = {}
            for nf in noisefiles_p9985:
                with open(nf, 'r') as f:
                    #print(nf)
                    noisedict_p9985.update(json.load(f))
                if 'n_earth' in noisedict_p9985.keys():
                    noisedict_p9985['{}_n_earth'.format(os.path.basename(nf).split('_')[0])] = noisedict_p9985['n_earth']
                    del(noisedict_p9985['n_earth'])

        return noisedict_p0015, noisedict_p9985



def parse_crn_name(crn_name):
        """
        function to parse user input for common noise inference runs.
        Determines if model comparison is requested and returns a dictionary of lsits
        """
        CRN_Name = {}
        #set up model comparison betwenn crn models
        if "-" in crn_name:

            #e.g. pl_nocorr_freegam-pl_hd_freegam: no correlations vs hd correlated, free-gamma run
            model_comp = True
            print("Doing model comparison!")
            crn_name1 = crn_name.split("-")[0]
            crn_name2 = crn_name.split("-")[1]
            #allow for multiple common noise models separated by comma
            if "," in crn_name1:
                crn_names_1 = crn_name1.split(",")
            else:
                crn_names_1 = [crn_name1]
            if "," in crn_name2:
                crn_names_2 = crn_name2.split(",")
            else:
                crn_names_2 = [crn_name2]
            CRN_Name[0] = crn_names_1
            CRN_Name[1] = crn_names_2

        else:
            model_comp = False
            #allow for multiple common noise models separated by comma
            if "," in crn_name:
                crn_names = crn_name.split(",")
            else:
                crn_names = [crn_name]
            CRN_Name[0] = crn_names
            
        return(model_comp, CRN_Name)


def construct_psrs(parfiles, timfiles, islice, yearslice, chainnum, ephem, load_existing = True, min_mjd = None, max_mjd = None):
        """
        construct Pulsar objects and does time-slicing of ToAs if necessary.
        """
        psrs = []
        if os.path.exists('./psrs.pkl') and islice is None and load_existing:
            psrs = pickle.load(open('./psrs.pkl', 'rb'))
            print('found pickle with {} psrs'.format(len(psrs)))
            return psrs
        for p, t in zip(parfiles, timfiles):

            #slice the tim files if necessary
            if islice is not None:
                print(min_mjd, max_mjd)
                tnew = t.replace('.tim', '_cut_{}_{}_{}.tim'.format(islice, yearslice, chainnum))
                command = "awk '{{if (($3 < {0}) || ($3 > {1})) print \"# \", $0; else print $0}}' {2}  > {3}".format(min_mjd, max_mjd, t, tnew)
                os.system(command)
                nlines = int(subprocess.check_output(f'grep -v "#" {tnew} | wc -l', shell = True, text = True))
                if nlines <= 1:
                    print("{} Has <= 1 ToA after cutting. Continuing.".format(p.split('/')[-1].split('.')[0]))
                    os.system('rm {}'.format(tnew))
                    continue
                print("loading...", p, tnew)
                psr = Pulsar(p, tnew, ephem=ephem)
                groups = np.unique(psr.flags["group"])
                #now delete groups with less than 
                groupcounts = [np.sum(psr.flags["group"] == g) for g in groups]
                data_tspan = psr.toas.max() - psr.toas.min()
                if data_tspan < 240.0*86400.0:
                    os.system('rm {}'.format(tnew))
                    #ensure that all pulsars in the slice have enough timespan for red noise inference
                    print('{} dataspan not long enough. continuing.'.format(psr.name))
                    continue
                for group, count in zip(groups, groupcounts):
                    if count <= 10:
                        #cut out groups with less than 10 ToAs
                        print('{} has less than 10 ToAs. cutting it out.'.format(group))
                        command = "grep -v '{}' {} > {}2".format(group, tnew, tnew)
                        os.system(command)
                        command = "mv {}2 {}".format(tnew, tnew)
                        os.system(command)
                        psr = Pulsar(p, tnew, ephem=ephem)

                try:
                    print("loading...", p, tnew)
                    psr = Pulsar(p, tnew, ephem=ephem)
                    groups = np.unique(psr.flags["group"])
                    #now delete groups with less than 
                    groupcounts = [np.sum(psr.flags["group"] == g) for g in groups]
                    data_tspan = psr.toas.max() - psr.toas.min()
                    if data_tspan < 240.0*86400.0:
                        os.system('rm {}'.format(tnew))
                        #ensure that all pulsars in the slice have enough timespan for red noise inference
                        print('{} dataspan not long enough. continuing.'.format(psr.name))
                        continue
                    for group, count in zip(groups, groupcounts):
                        if count <= 10:
                            #cut out groups with less than 10 ToAs
                            print('{} has less than 10 ToAs. cutting it out.'.format(group))
                            command = "grep -v '{}' {} > {}2".format(group, tnew, tnew)
                            os.system(command)
                            command = "mv {}2 {}".format(tnew, tnew)
                            os.system(command)
                    psr = Pulsar(p, tnew, ephem=ephem)

                except ValueError as e:
                    print(e)
                    print("{} Has no data after cutting. Continuing.".format(p.split('/')[-1].split('.')[0]))
                    os.system('rm {}'.format(tnew))
                    continue
                os.system('rm {}'.format(tnew))
            else:
                #no slicing required: just load in the plain par and tim files
                print("loading...", p, t)
                psr = Pulsar(p, t, ephem=ephem)

            if psr.name in psrnames:
                psrs.append(psr)
        return(psrs)


    
def get_tspan_fundamental_freq(psrs):
    """
    gets the timespan and fundamental freq (1 / tspan)
    given a list of Pulsar objects
    """
    toamins = [p.toas.min() for p in psrs]
    toamaxs = [p.toas.max() for p in psrs]
    tspan =  np.max(toamaxs) - np.min(toamins)
    fundamental_freq = 1.0/tspan
    return tspan, fundamental_freq
    
"""
Define common red noise
"""
#components = 5

def get_crn_model_dict(tspan):
    """
    function to return a dictionary "menu" of common red noise models set up via crn_block
    with appropriate number of commponents.
    """
    fundamental_freq = 1.0/tspan
    max_freq = 1.0/240.0/86400.0 #1 / 240 days
    float_components = max_freq/fundamental_freq 
    components = int(float_components)

    crn_model_dict = {
        #powerlaw, no correlation, free spectral index
        'pl_nocorr_freegam': crn_block(psd='powerlaw', prior='log-uniform', gamma_val=None,
                                       components=components, orf=None, name='gw_pl_nocorr_freegam'),
        #powerlaw, no correlation, fixed 13/3 spectral index
        'pl_nocorr_fixgam': crn_block(psd='powerlaw', prior='log-uniform', gamma_val=4.333,
                                      components=components, orf=None, name='gw_pl_nocorr_fixgam'),
        #broken power law, no correlation
        'bpl_nocorr_freegam': crn_block(psd='broken_powerlaw', prior='log-uniform',
                                        components=60, orf=None, name='gw_bpl_nocorr_freegam'),
        #free-spectrum no correlation
        'freespec_nocorr': crn_block(psd = 'spectrum', prior = "log-uniform", components = 40,
                                    orf = None, name = 'gw_freespec_nocorr'),
        #power law red noise, ORF sampled in angular bins
        'pl_orf_bins': crn_block(psd='powerlaw', prior='log-uniform', gamma_val=4.333,
                                 components=components, orf='bin_orf', name='gw_apl_orf_bins'),
        #power law red noise, ORF sampled in angular bins
        'pl_orf_singlebin': common_red_noise_block(psd='powerlaw', prior='log-uniform', gamma_val=4.333,
                                 components=components, orf='single_bin_orf', name='gw_apl_orf_bins'),
        'pl_orf_singlebin_fixA': common_red_noise_block(psd='powerlaw', prior='log-uniform', gamma_val=4.333,
                                                   log10_A_val=-14.69, components=components, orf='single_bin_orf', name='gw_apl_orf_bins_fixA'),
        #power law red noise, spline-interpolated ORF bins
        'pl_orf_spline': crn_block(psd='powerlaw', prior='log-uniform', gamma_val=4.333,
                                   components=components, orf='spline_orf', name='gw_apl_orf_spline'),
        #power law red noise, fixed spectral index, Hellings-Downs correlations,
        'pl_hd_fixgam': crn_block(psd='powerlaw', prior='log-uniform',
                                  components=components, orf='hd', name='gw_pl_hd_fixgam', gamma_val = 4.333),
        #power law red noise, free spectral index, Hellings-Downs correlations,
        'pl_hd_freegam': crn_block(psd='powerlaw', prior='log-uniform',
                                   components=components, orf='hd', name='gw_pl_hd_freegam', gamma_val = None),    
        #power law red noise, fixed spectral index, Hellings-Downs correlations excluding auto-correlations
        'pl_hdnoauto_fixgam': crn_block(psd='powerlaw', prior='log-uniform',
                                        components=components, orf='zero_diag_hd', name='gw_pl_hdnoauto_fixgam', gamma_val = 4.333),
        #power law red noise, fixed spectral index, Hellings-Downs correlations excluding auto-correlations
        'pl_crnnoauto_fixgam': common_red_noise_block(psd='powerlaw', prior='log-uniform',
                                        components=components, orf='zero_diag_crn', name='gw_pl_hdnoauto_fixgam', gamma_val = 4.333),
        #free spectrum red noise, Hellings-Downs correlations
        'freespec_hd': crn_block(psd = 'spectrum', prior = "log-uniform", components = 40,
                                 orf = 'hd', name = 'gw_freespec_hd'),
        #power law red noise, free spectral index, Dipole correlations,
        'pl_dipole_freegam': crn_block(psd='powerlaw', prior='log-uniform',
                                       components=components, orf='dipole', name='gw_pl_dipole_freegam', gamma_val = None),
        #free spectrum red noise, Dipole correlations
        'freespec_dipole': crn_block(psd = 'spectrum', prior = "log-uniform", components = 40,
                                     orf = 'dipole', name = 'gw_freespec_dipole'),
        #power law red noise, free spectral index, Monopole correlations,
        'pl_monopole_freegam': crn_block(psd='powerlaw', prior='log-uniform',
                                         components=components, orf='monopole', name='gw_pl_monopole_freegam', gamma_val = None),
        #free spectrum red noise, Monopole correlations
        'freespec_monopole': crn_block(psd = 'spectrum', prior = "log-uniform", components = 40,
                                       orf = 'monopole', name = 'gw_freespec_monopole')    
    }

    return crn_model_dict


def get_informed_nearth_priors(psr, noisedict_3sig_min, noisedict_3sig_max, priors):

    if psr.name[1:5] in ['0030', '0125', '0437', '0613', '1022', '1024', '1545', '1600', '1643', '1713', '1730', '1744', '1824', '1832', '1909', '2145', '2241']:
        key_ = psr.name + '_n_earth'
        print('getting prior for {}'.format(key_))
        try:
            n_earth_min = np.max([0.0, noisedict_3sig_min[key_] - 1 ])
            n_earth_max = np.min([20.0, noisedict_3sig_max[key_] + 1 ])
            print(f'found n_earth_min = {n_earth_min}, n_earth_max = {n_earth_max}')
        except KeyError as e:
            print('KeyError: {}'.format(e))
            n_earth_min = 0.0
            n_earth_max = 20.0
        n_earth = parameter.Uniform(n_earth_min, n_earth_max)
        priors[key_ + "_min"] = n_earth_min
        priors[key_ + "_max"] = n_earth_max
    else:
        print('nearth constant')
        n_earth = parameter.Constant(4)

    deter_sw = solar_wind(n_earth=n_earth)
    mean_sw = deterministic_signals.Deterministic(deter_sw)#, name='n_earth')
    sw = mean_sw

    #vary gp_sw for pulsars with constrained gp_sw parameters
    if psr.name[1:5] in ['0437', '0711', '0900', '1024', '1643', '1713', '1730', '1744', '1909', '2145']:
        key_ = psr.name + '_gp_sw_log10_A'
        print('getting prior for {}'.format(key_))
        try:
            log10_A_min = np.max([-10, noisedict_3sig_min[key_] - 1 ])
            log10_A_max = np.min([-3, noisedict_3sig_max[key_] + 1 ])
            print(f'found log10_A_min = {log10_A_min}, log10_A_max = {log10_A_max}')
            key_ = psr.name + '_gp_sw_gamma'
            gamma_min = np.max([-4, noisedict_3sig_min[key_] - 0.5])
            gamma_max = np.min([4, noisedict_3sig_max[key_] + 0.5])
            print(f'found gamma_min = {gamma_min}, gamma_max = {gamma_max}')
        except KeyError as e:
            print('KeyError:', e)
            log10_A_min = -10
            log10_A_max = -3
            gamma_min = -4
            gamma_max = 4
        log10_A_sw = parameter.Uniform(log10_A_min, log10_A_max)
        gamma_sw = parameter.Uniform(gamma_min, gamma_max)
        print(f"""{psr.name}_gp_sw_noise prior:
        log10_A in [{log10_A_min}, {log10_A_max}]
        gamma in [{gamma_min}, {gamma_max}]
        """)
        key_ = psr.name + '_gp_sw_log10_A'
       	priors[key_ + "_min"] = log10_A_min
       	priors[key_ + "_max"] = log10_A_max
       	key_ = psr.name + '_gp_sw_gamma'
       	priors[key_ + "_min"] = gamma_min
       	priors[key_ + "_max"] = gamma_max        

        Tspan = psr.toas.max() - psr.toas.min()
        max_cadence = 60
        sw_components = int(Tspan / (max_cadence*86400))
        sw_prior = utils.powerlaw(log10_A=log10_A_sw, gamma=gamma_sw)
        sw_basis = createfourierdesignmatrix_solar_dm(nmodes=sw_components, Tspan=Tspan)
        sw += gp_signals.BasisGP(sw_prior, sw_basis, name='gp_sw')

    return sw, priors


def get_informed_rednoise_priors(psr, noisename, noisedict_3sig_min, noisedict_3sig_max, priors, return_priorvals = False, use_basic_priors = False, log10_A_min_basic=-20, log10_A_min_informed=-18, log10_A_min_bound=-2):
    key_ = psr.name + '_' + noisename + '_log10_A'
    print('getting prior for {}'.format(key_))
    if not use_basic_priors:
        try:
            log10_A_min = np.max([log10_A_min_informed, noisedict_3sig_min[key_] + log10_A_min_bound ])
            log10_A_max = np.min([-11, noisedict_3sig_max[key_] + 1 ])
            print(f'found log10_A_min = {log10_A_min}, log10_A_max = {log10_A_max}')
            key_ = psr.name + '_' + noisename + '_gamma'
            gamma_min = np.max([0, noisedict_3sig_min[key_] - 0.5])
            gamma_max = np.min([7, noisedict_3sig_max[key_] + 0.5])
            print(f'found gamma_min = {gamma_min}, gamma_max = {gamma_max}')
            
        except KeyError as e:
            print('KeyError:', e)
            log10_A_min = -18
            log10_A_max = -11
            gamma_min = 0
            gamma_max = 7
    else:
        log10_A_min = log10_A_min_basic
        log10_A_max = -11
        gamma_min = 0
        gamma_max = 7
    log10_A_prior = parameter.Uniform(log10_A_min, log10_A_max)
    gamma_prior = parameter.Uniform(gamma_min, gamma_max)
    print(f"""{psr.name}_{noisename}_noise prior:
    log10_A in [{log10_A_min}, {log10_A_max}]
    gamma in [{gamma_min}, {gamma_max}]
    """)
    # # powerlaw
    rednoise_model = gp_priors.powerlaw(log10_A=log10_A_prior,
                                        gamma=gamma_prior)
    key_ = psr.name + '_' + noisename + '_log10_A'
    priors[key_ + "_min"] = log10_A_min
    priors[key_ + "_max"] = log10_A_max
    key_ = psr.name + '_' + noisename + '_gamma'
    priors[key_ + "_min"] = gamma_min
    priors[key_ + "_max"] = gamma_max
    if return_priorvals:
        return rednoise_model, log10_A_prior, gamma_prior, log10_A_min, log10_A_max, gamma_min, gamma_max, priors
    
    return rednoise_model, log10_A_prior, gamma_prior, priors



def get_initial_samples_chain(chainfiles, hyper_model, priors, informed_sample = True):

    x0 = hyper_model.initial_sample()
    ndim = len(x0) 
    params = hyper_model.param_names
    if not informed_sample:
        return x0
    
    for chainfile in chainfiles:

            if chainfile.split('_')[0].split('/')[-1] not in psrnames:
                continue
            print("Loading {}".format(chainfile))

            chain = np.load(chainfile)
            parsfile = chainfile.replace('_chain.npy', '_pars.npy')
            pars = np.load(parsfile)
            print(pars)

            for i in range(0, ndim):
                psrname = params[i].split('_')[0]
                if 'J' not in psrname: continue
                if psrname not in chainfile: continue

                print("Parameter {}, prior sample {}".format(params[i], x0[i]))

                no_sample = True
                while no_sample:        
                    ind_sample = np.random.randint(0, len(chain[:, 0]))
                    param_new = params[i]
                    if 'n_earth' in params[i]:
                        param_new = 'n_earth'
                    if '_low_low' in params[i]:
                        param_new = params[i].replace('_low_low', '_low')
                    ind_par = np.argwhere([param_new == p for p in pars])
                    if '_mid_mid' in params[i] or '_high_high' in params[i]:
                        param_new = param_new.replace('_mid_mid', '_mid').replace('_high_high', '_high')
                    if len(ind_par) == 0:
                        print("Parameter {} not in chain - continuing.".format(param_new))
                        break
                    sample = chain[ind_sample, ind_par].squeeze()
                    if 'n_earth' in params[i]:
                        param_new = params[i]  # change back for prior name
                    if not 'dmexp' in params[i] and not 'dm1yr' in params[i] and not 'dmgaus' in params[i]:
                        try:
                            if (sample < priors[param_new + '_min'] or sample >  priors[param_new + '_max']):
                                print("################ SAMPLE {} OUTSIDE OF PRIOR RANGE FOR {} ###############".format(sample, params[i]))
                                continue
                        except KeyError as e:
                            print(e)
                            break
                    x0[i] = sample
                    print("Parameter {}, posterior sample {}".format(params[i], x0[i]))
                    no_sample = False

    for i in range(0, len(params)):
        if 'bin' in params[i] and not 'log10' in params[i]:
            x0[i] = 0.0

    return x0



def get_noisedict_maxlike(psrname, chain, pars):
    x = {}
    for ct, par in enumerate(pars):
        ind = np.argmax(chain[:, -4])  # find index of maximum posterior likelihood
        x[par] = chain[ind, ct] 

    #os.system('mkdir -p {}'.format(outdir))
    #with open(outdir + '/{}_{}_noise.json'.format(psrname, noise_model), 'w') as fout:
    #    json.dump(x, fout, sort_keys=True, indent=4, separators=(',', ': '))
    return x


def acor(arr):
    """
    for thinning chains
    """
    arr -= np.mean(arr)
    auto_correlation = correlate(arr, arr, mode='full')
    auto_correlation = auto_correlation[auto_correlation.size//2:]
    auto_correlation /= auto_correlation[0]
    indices = np.where(auto_correlation<0.5)[0]
    if len(indices)>0:
        if indices[0] == 0:
            indices[0] = 1
        return indices[0]
    else:
        return 1


def get_chain(noise_model, datadir, max_nchain = 500, nburnmin = 10000, plot_trace = False, load_existing = True, save_loaded = True):
        """
        Used for burning-in, thinning, and combining chains
        """
        #max_nchain = 100  # number of chains to read
        #nburnmin = 1200
        first_chain = True

        for i in range(1, max_nchain+1): 
            outdir = datadir + "/chains/{}/{}".format(noise_model, i)
            if load_existing and os.path.exists("{}/chains/{}.npy".format(datadir, noise_model)):
                pars = np.loadtxt(outdir + '/pars.txt', dtype=np.unicode_)
                chain_total = np.load("{}/chains/{}.npy".format(datadir, noise_model))
                break
            chainfile = outdir + '/chain_1.txt'
            print(chainfile)
            nlines = int(subprocess.check_output(f'cat {chainfile} | wc -l', shell = True, text = True))
            if nlines < nburnmin:
                print(nlines, "dead")
                continue
            if not os.path.exists(chainfile):
                continue
            if os.path.getsize(chainfile) == 0:
                continue
            print("loading {}".format(chainfile))
            chain_i = np.loadtxt(chainfile).squeeze()
            pars = np.loadtxt(outdir + '/pars.txt', dtype=np.unicode_)
            print('loaded pars')

            if plot_trace:
                pp = model_utils.PostProcessing(chain_i, pars)
                pp.plot_trace()
                plt.savefig(chainfile.replace('.txt', '_{}_{}_trace.png'.format(psrname, dir)))
                plt.close()
                print("Made {}".format(chainfile.replace('.txt', '_{}_{}_trace.png'.format(psrname, dir))))

            # Burn larger of first 25% or 2500
            burn = int(max([0.25*chain_i.shape[0], nburnmin]))
            inds = np.argwhere([('nmodel' not in p and 'gw' not in p and 'bin' not in p) for p in pars])
            thin = 1
            good_pars = 0
            for indi in inds:
                try:
                    thin_i = acor(chain_i[burn:, indi])
                except ValueError:
                    continue
                thin += thin_i
                good_pars += 1
                #print("Autocorrelation on {} is {}".format(pars[ind], thin_i))
            if good_pars == 0: continue
            thin /= good_pars
            thin = int(thin)
            chain = chain_i[burn::thin, :]    
            #chain = chain_i[burn:, :]
            if chain.size == 0 or np.std(chain) < 1e-2:
                continue

            if first_chain:
                chain_total = chain
                first_chain = False
            else:
                chain_total = np.concatenate((chain_total, chain)).squeeze()

        np.save("{}/chains/{}.npy".format(datadir, noise_model), chain_total)
        return pars, chain_total


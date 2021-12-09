import os
#import numpy as np
import csv
from scm.plams import Settings, Molecule, Units, AMSJob, MultiJob, CRSJob, KFFile, ADFResults, JobRunner, GridRunner
import multiprocessing

wd = os.getcwd()


"""
Before running the script there are few paths to be adjusted:

1. The paths to the folder with starting geometries, folder with 
optimised and not optimised structures and to the cosmo-rs solvent coskf file

2. The configurations - local and remote - for parralel job running

3. Version of ADF (to be loaded on remote nodes)

4. At the bottom of the script calculation settings 

"""

starting_geometries = os.path.join(wd,'start')


# Path to optimised and non optimised structures
preopt          = os.path.join(wd, 'Preoptimized')
non_preopt      = os.path.join(wd, 'Not_Preoptimized')

GeoOpt          = os.path.join(wd, "Optimized_Geo")
NotGeoOpt       = os.path.join(wd, "Not_Optimized_Geo")

# Solvent path
database_path   = os.path.join(wd,'ADFCRS', 'Dichloromethane.coskf')


# Configure parralel job running localy 
maxjobs = multiprocessing.cpu_count()
jr = JobRunner(parallel=True, maxjobs=maxjobs)

# Configure parralel job running remote
gr = GridRunner(maxjobs=10)

# Version of ADF (to be loaded on remote nodes)
module_load = "module load 2021 \n module load AMS/2021.101-intelmpi" 

# Settings for preoptimisation of the molecules that haven't been previously 
#optimised or optimised on FF level

parallel 
conv
################################################################################################################
############################################# Settings #########################################################
################################################################################################################


# Optimisation on the DFTB level
def DFTB_optimisation(model = ["DFTB3", 'DFTB.org/3ob-3-1'], nf_tolerance = -20, charge = None):
    
    settings = Settings()
    settings.input.ams.task  = 'GeometryOptimization' 
    settings.input.ams.Properties.NormalModes           = 'Yes'
    settings.input.ams.Properties.PESPointCharacter     = 'Yes'
    settings.input.ams.NormalModes.ReScanFreqRange      = '-1000 0'
    settings.input.ams.PESPointCharacter.NegativeFrequenciesTolerance = nf_tolerance
    settings.input.DFTB
    settings.input.DFTB.Model  = model[0]
    if len(model)==2:
        settings.input.DFTB.ResourcesDir = model[1]
    if charge:
        settings.input.ams.System.Charge = charge

    settings.runscript.nproc = 1

    return settings


# Single point calculation 
def SinglePoint(xc_fun=['LibXC', 'CAM-B3LYP'], basis='DZP', charge=0, tddft=False, solvent=False):
    
    settings = Settings()
    
    settings.runscript.pre  = module_load
    #settings.keep = 'all'

    settings.input.ams.task           = 'SinglePoint'
    settings.input.adf.basis.type     = basis
    settings.input.adf.basis.core     = 'None'
    settings.input.adf.xc[xc_fun[0]]  = xc_fun[1]
    settings.input.adf.xc.Dispersion  = 'GRIMME3 BJDAMP'
    settings.input.adf.Relativity.Level   = 'None'
    settings.input.adf.NumericalQuality   = 'Normal'
    settings.input.ams.UseSymmetry            = 'No'
    settings.input.adf.RIHartreeFock.UseMe    = 'Yes'
    settings.input.adf.RIHartreeFock.Quality  = 'Normal'
    if charge:
        settings.input.adf.Unrestricted     = 'Yes'
        settings.input.adf.SpinPolarization = charge
        settings.input.ams.System.Charge    = charge
    if tddft:
        settings.input.adf.Excitations.lowest     = 10
        settings.input.adf.Excitations.OnlySing   = ''
        
        if tddft == 'Davidson':
            settings.input.adf.Excitations[tddft] 
            
        else:
            settings.input.adf.Excitations[tddft] = ''   

        if tddft == 'stddft':
            settings.input.adf.ModifyExcitation.GrimmeAlpha = 0.9
            settings.input.adf.ModifyExcitation.GrimmeBeta  = 1.86
            settings.input.adf.ModifyExcitation.GrimmeAex   = 0.38
    
    if solvent:
        if solvent == 'CRS':
            solvation_block = {
                'surf': 'Delley',
                'solv': 'name=CRS cav0=0.0 cav1=0.0',
                'charged': 'method=Conj corr',
                'c-mat': 'Exact',
                'scf': 'Var All',
                'radii': {
                    'H': 1.30,
                    'C': 2.00,
                    'N': 1.83,
                    'O': 1.72,
                    'F': 1.72,
                    'Si': 2.48,
                    'P': 2.13,
                    'S': 2.16,
                    'Cl': 2.05,
                    'Br': 2.16,
                    'I': 2.32
                }
            }
        else:

            solvation_block = solvent

        settings.input.adf.solvation    = solvation_block    
        
    return settings

# COSMO-RS calculation
def COSMORSsett(solute_path, solvent_path, num_compounds =2):
    
    actcoef = Settings()
    actcoef.input.property._h = 'ACTIVITYCOEF'

    compounds = [Settings() for i in range(num_compounds)]
    compounds[0]._h = solvent_path
    compounds[1]._h = solute_path
    #set compound mole fractions
    compounds[0].frac1 = 1
    compounds[1].frac1 = 0

    actcoef.input.temperature = "298.15"
    actcoef.input.compound = compounds

    return actcoef

#################################################################################################################
####################################### PREOPTIMISATION  ########################################################
#################################################################################################################


def make_csv(file_name, header, data):
    """
    file_name : string
        name of new csv file 
    header : list 
        list of strings 
    data : dic
        key is a string with a name of the molecule
        value is list of floats
    """
    
    with open(file_name + '.csv', 'w') as f:
        writer = csv.writer(f)
        writer.writerow(header)
        for n in data:
            x = [n] + data[n]
            writer.writerow(x)

csv

def export_coords(geometries, path):
    """
    geometries: dictionary with PLAMS Molecules
        The dictionary contains geometries to be places in corresponding folder
    path : string
        Path to the folder where geometries will be written

    """
    try:
        if not os.path.exists(path):
            os.makedirs(path)
    
        for name, mol in geometries.items():
            mol.write(os.path.join(path, name + '.xyz' ))  
    
    except AttributeError:
        print("Geometries for exporting not found")



def list_AMSjobs(job_sett, prefix, molecules):
    """
    job_sett    : PLAMS Settings
    prefix      : string
        to be added to molecules name to indicate type of job
    molecules     : dictionary 
        key is name of a molecules, value is PLAMS Molecule
    """
    list_of_jobs = [ AMSJob(name        = prefix + m, 
                            molecule    = molecules[m], 
                            settings    = job_sett)
                    for m in molecules
                    ]

    return list_of_jobs




def select_opt_mol(results, prefix, path):
    """
    results  : PLAMS Result object
    prefix   : string
        to be added to molecules name to indicate type of job
    path     : string
        optimised molecules are exported to this path
    """

    success_jobs = [job for job in results if job.ok()]
    unsuccess_jobs = [job for job in results if job not in success_jobs]

    opt_mol = {res.name.replace(prefix,'') : res.get_main_molecule() for res in success_jobs}
    export_coords(opt_mol, os.path.join(path))

    return success_jobs, unsuccess_jobs


def parallel_preopt(settings, prefix, mols, maxjobs=False):
    """
    maxjobs     : integer
        numeber of jobs allowed to run in parallel
        if not specified, it will equal the numper of cpu's
    prefix      : string
        to be added to molecules name to indicate type of job
    mols     : dictionary 
        key is name of a molecules, value is PLAMS Molecule
    """
    if maxjobs:
        pass
    else:
        maxjobs = multiprocessing.cpu_count()
    
    jobs = list_AMSjobs(settings, prefix, mols)
    results = [job.run(jr) for job in jobs]

    return results



def optimization(mols, prefix, settings, path):
    """
    Runs optimization with a semiempirical models of the initial structures. 
    Once the calculation is done it checks the status of the calculation,
    extracts final geometries - sucessfuly and unsucessfuly optimised

    mols     : dictionary 
        key is name of a molecules, value is PLAMS Molecule
    prefix   : string
        to be added to molecules name to indicate type of job
    settings : Plams Settings 
    path     : string
        optimised molecules are exported to this path

    """
    
    opt_mol, nonopt_mol = {}, {}
    
    results     = parallel_preopt(settings, prefix, mols)
    
    opt_mol, fail    = select_opt_mol(results, prefix, os.path.join(path))

    if fail:
        not_converged = [job for job in fail if job.grep_output(
            'ERROR: Geometry optimization failed! (Converged, but did not find a minimum on the PES.)')]
        not_converged_mols  = {res.name : res.get_main_molecule() for res in not_converged}
        
        if not_converged:
            increase_accuracy = Settings()
            increase_accuracy.input.ams.GeometryOptimization.Convergence.Gradients   = 0.0001 
            increase_accuracy.input.ams.GeometryOptimization.Convergence.Energy      = 1e-06
            

            settings += increase_accuracy
            results     = parallel_preopt(settings, prefix, not_converged_mols)
            
            opt_mol, fail   = select_opt_mol(results, prefix, path)
    
            if fail:
                nonopt_mol = {res.name : res.get_main_molecule() for res in fail}
                export_coords(nonopt_mol, non_preopt)
        else:
            nonopt_mol = {res.name.replace(prefix,'') : res.get_main_molecule() for res in fail}
            export_coords(nonopt_mol, os.path.join(path, non_preopt))




###############################################################################################################
####################################### MAIN CALCULATIONS  ####################################################
###############################################################################################################


# TDDFT calculations
##############################################################################################################

class TDDFT_jobs(MultiJob):
    
    def postrun(self):
        """
        Reads the computational details: the number of mulecules, the time (available from AMS2021) and final 
        status of the calculation and writes it in csv file
        Reads the values of excitation energies and oscilatior strengths
        and writes it in csv file 
        """
        data = {}
        tddft_CD = {}

        for job in self.children:
            if job.ok():
                try:
                    numA        = job.results.readrkf('Molecule','nAtoms')
                    Etime       = job.results.readrkf('General', 'ElapsedTime')
                    Tstatus     = job.results.readrkf('General', 'termination status')
                    
                    cd = [Tstatus] + [numA] + [Etime]
                    a = {job.name : cd}
                    tddft_CD.update(a)
                except:
                    cd = []
                    a = {job.name : cd}
                    tddft_CD.update(a)
            else:
                a = {job.name : False}
                tddft_CD.update(a)
            

            if job.ok():
                try:
                    #NumOfExc        = job.results.readrkf('All excitations','nr excitations','adf')
                    ExcEnergies  = Units.convert(job.results.readrkf('Excitations SS A', 'excenergies','adf'), 'au', 'eV')
                    OscilatorStr = job.results.readrkf('Excitations SS A', 'oscillator strengths','adf')
                    
                    i = {job.name : ExcEnergies + OscilatorStr}
                    data.update(i)
                except: 
                    i = {job.name : False}
                    data.update(i)
            else:
                i = {job.name : False}
                data.update(i)
        

        header = ['Name'] + list(range(1, 11)) + list(range(1, 11))
        make_csv(self.name, header, data)

        header_cd = ['Name', 'TermStatus', 'NAtoms', 'Etime']
        make_csv('td_comput_details', header_cd, tddft_CD)

def getExcitations(dirname, molecules, settings):
    """
    Creates the TDDFT_jobs MultiJob of and runs them with prefered jobrunner
    """

    jobs = [ AMSJob(name        = m, 
                    molecule    = molecules[m], 
                    settings    = settings[s])
                    for m in molecules
                    for s in settings
                    ]

    ManyJobs = TDDFT_jobs(name=dirname, children=jobs)

    ManyJobs.run(jobrunner=gr, nodes=1, cores=32, walltime='00-01:00:00')
    #ManyJobs.run(jr)


# ADF with COSMO to create the starting surface for COSMO-RS
"""
# COSMO-RS 
############################################################################################################

class MultipleCRS(MultiJob):
    def prerun():
        print('starting cosmo-rs calculations')

def calculateGibbsFree(dirname, solute, solvent, settings):
    
    CRS_jobs = [CRSJob( settings = COSMORSsett(solute_path=v, solvent_path=solvent), 
                        name     = 'crs' + n) 
                for n,v in solute.items()
                ]
                   
    ManyCRSJobs = MultipleCRS(name=dirname, children=CRS_jobs)
    #ManyJobs.settings.save = ['-', '$CH/adf.rkf']

    ManyCRSJobs.run(jr)
"""
############################################################################################################        


class ADFwithCOSMO(MultiJob):


    def postrun(self):

        data = {}
        comp_data = {}

        for job in self.children:
            """
            Reads the details the time (available from AMS2021) and final status of the calculation
            and writes it in cvs file
            Reads the COSMO section and creates a file with its copy
            Makes the COSMO-RS job scritps and runs COSMO-RS calculations
            Reads the value of solution-phase Gibbs free energy and writes it in csv file 
            """
            if job.ok():
                try:
                    numA        = job.results.readrkf('Molecule','nAtoms')
                    #Etime       = job.results.readrkf('General', 'ElapsedTime')
                    Tstatus     = job.results.readrkf('General', 'termination status')
                    
                    cd = [Tstatus] + [numA] #+ [Etime]
                    a = {job.name : cd}
                    comp_data.update(a)
                except:
                    cd = []
                    a = {job.name : cd}
                    comp_data.update(a)
                try:
                    # read COSMO section of adf.rkf file for a child job
                    resfile = KFFile(job.path + '/adf.rkf')
                    cosmo   = resfile.read_section("COSMO")
        
                    # make a new file .coskf 
                    coskf_path = os.path.join(job.path, job.name + '.coskf')
                    coskf   = KFFile(coskf_path)
                    for k,v in cosmo.items():
                        coskf.write("COSMO",k,v) 
        
                    crs_sett = COSMORSsett(solute_path=coskf_path, solvent_path=database_path)
                    crs_job = CRSJob(settings = crs_sett, 
                                    name     = 'crs' + job.name)
                    crs_res = crs_job.run()
                    
                    #extract data and update data dictionary
                    extracted = crs_res.get_results()
                    

                    if extracted:
                        i = {job.name: [Units.convert( float(extracted["G solute"][1]), 'kcal/mol', 'eV')]} 
                        data.update(i)
                    else: 
                        i = {job.name : False}
                        data.update(i)

                except:
    
                    i = {job.name : False}
                    data.update(i)
            else:
                i = {job.name : False}
                data.update(i)


        header = ['Name', 'G_'+self.name]
        make_csv(self.name, header, data)

        header_cd = ['Name', 'TermStatus', 'NAtoms'] #+ [Etime]]
        make_csv('comput_details', header_cd, comp_data)


def getGibbsFreeEnergy(dirname, molecules, settings):
    """
    Creates the ADFwithCOSMO MultiJob of and runs them with prefered jobrunner
    """

    jobs = [ AMSJob(name        = m, 
                    molecule    = molecules[m], 
                    settings    = settings[s])
                    for m in molecules
                    for s in settings
                    ]

    ManyJobs = ADFwithCOSMO(name=dirname, children=jobs)


    ManyJobs.run(jobrunner=gr, nodes=1, cores=32, walltime='00-01:00:00')
    #ManyJobs.run(jr)



############################################################################################################
#############################################     Run      #################################################
############################################################################################################


# Preoptimisation

in_molecules = read_molecules(starting_geometries)
optimization(mols    = in_molecules, 
             prefix     = '', 
             settings   = DFTB_optimisation(model = ["DFTB3", 'DFTB.org/3ob-3-1']), 
             path       = preopt
             )



# Geometry optimisation 

# GNV - Geometry optimisation for Neutral molecule in Vacuum
# GOV - Geometry optimisation for Oxidised molecule in Vacuum
GO = {  'gnv' : DFTB_optimisation(model = ['GFN1-xTB']),

        'gov' : DFTB_optimisation(model = ['GFN1-xTB'], charge = 1)
    }

molecules = read_molecules(preopt)
for job_name, settings in GO.items():
    optimization(mols       = molecules, 
                 prefix     = job_name, 
                 settings   = settings, 
                 path       = os.path.join(GeoOpt,job_name)
                 )



#Calculate Gibbs free energies

# get settings 
neu = {  'ns' : SinglePoint(xc_fun=['GGA', 'BLYP'], basis='DZ', solvent='CRS') }
ox =  {  'os' : SinglePoint(xc_fun=['GGA', 'BLYP'], basis='DZ', solvent='CRS', charge=1)}

# get neutral and oxidised optimised 
molecules_neu   = read_molecules(os.path.join(GeoOpt,'gnv'))
molecules_ox    = read_molecules(os.path.join(GeoOpt,'gov'))

# intersection of 
available_structures = list(set(list(molecules_neu.keys())).intersection(list(molecules_ox.keys())))


molecules_neu   = {mol: molecules_neu[mol] for mol in available_structures}
molecules_ox    = {mol: molecules_ox[mol] for mol in available_structures}

oxi_data     = getGibbsFreeEnergy(dirname='oxi', molecules=molecules_ox, settings = ox )
neutral_data = getGibbsFreeEnergy(dirname='neu', molecules=molecules_neu, settings = neu )



# Calculate excitations

stddft = { 'exc' : SinglePoint(tddft='stddft')}
excitations = getExcitations(dirname='exc', molecules=molecules_neu, settings=stddft)







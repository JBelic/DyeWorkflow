import os
#import numpy as np
import csv
from scm.plams import Settings, Molecule, Units, AMSJob, MultiJob, CRSJob, KFFile, ADFResults, JobRunner, GridRunner
import multiprocessing

wd = os.getcwd()


"""
Before running the script adjust the following:

1. The paths to the folder with starting geometries, a folder with 
optimized and not optimized structures, and to the cosmo-rs solvent coskf file

2. The configurations - local and remote - for a parallel job running

3. Version of ADF (to be loaded on remote nodes)

4. At the bottom of the script calculation settings 

"""

starting_geometries = os.path.join(wd,'start')


#  Path to optimized and non-optimized structures
preopt          = os.path.join(wd, 'Preoptimized')
non_preopt      = os.path.join(wd, 'Not_Preoptimized')

GeoOpt          = os.path.join(wd, "Optimized_Geo")
NotGeoOpt       = os.path.join(wd, "Not_Optimized_Geo")

# Solvent path
database_path   = os.path.join(wd,'ADFCRS', 'Dichloromethane.coskf')


# Configure parallel job running locally 
maxjobs = multiprocessing.cpu_count()
jr = JobRunner(parallel=True, maxjobs=maxjobs)


# A version of ADF (to be loaded on remote nodes)
module_load = "module load 2021 \n module load AMS/2021.101-intelmpi" 

# Settings for preoptimization of the molecules that haven't been 
# previously optimized or optimized on FF level



################################################################################################################
############################################# Settings #########################################################
################################################################################################################


# Optimization on the DFTB level
def DFTB_optimization(model = ["DFTB3", 'DFTB.org/3ob-3-1'], nf_tolerance = -20, charge = None):
    
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
        a key is a string with a name of the molecule
        value is a list of floats
    """
    
    with open(file_name + '.csv', 'w') as f:
        writer = csv.writer(f)
        writer.writerow(header)
        for n in data:
            x = [n] + data[n]
            writer.writerow(x)


def export_coords(geometries, path):
    """
    geometries: dictionary with PLAMS Molecules
        The dictionary contains geometries to be placed in the corresponding folder
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
        to be added to molecules name to indicate the type of job
    molecules : dictionary 
        the key is the name of a molecule, value is PLAMS Molecule
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
        to be added to molecules name to indicate the type of job
    path     : string
        optimized molecules are exported to this path
    """

    success_jobs = [job for job in results if job.ok()]
    unsuccess_jobs = [job for job in results if job not in success_jobs]

    opt_mol = {res.name.replace(prefix,'').replace('.002','') : res.get_main_molecule() for res in success_jobs}
    export_coords(opt_mol, os.path.join(path))

    return success_jobs, unsuccess_jobs


def parallel_preopt(settings, prefix, mols, maxjobs=False):
    """
    maxjobs     : integer
        number of jobs allowed to run in parallel
        if not specified, it will equal the number of CPU's
    prefix      : string
        to be added to molecules name to indicate the type of job
    mols     : dictionary 
        the key is the name of a molecule, value is PLAMS Molecule
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
    Runs optimization with a semiempirical model of the initial structures. 
    Once the calculation is done it checks the status of the calculation,
    extracts final geometries - successfully and unsuccessfully optimised

    mols     : dictionary 
        the key is the name of a molecule, value is PLAMS Molecule
    prefix   : string
        to be added to molecules name to indicate the type of job
    settings : PLAMS Settings 
    path     : string
        optimized molecules are exported to this path

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




############################################################################################################
#############################################     Run      #################################################
############################################################################################################


# Preoptimisation

in_molecules = read_molecules(starting_geometries)
optimization(mols    = in_molecules, 
             prefix     = '', 
             settings   = DFTB_optimization(model = ["DFTB3", 'DFTB.org/3ob-3-1']), 
             path       = preopt
             )



# Geometry optimisation 

# GNV - Geometry optimization for Neutral molecule in Vacuum
# GOV - Geometry optimization for Oxidised molecule in Vacuum
GO = {  'gnv' : DFTB_optimization(model = ['GFN1-xTB']),

        'gov' : DFTB_optimization(model = ['GFN1-xTB'], charge = 1)
    }

molecules = read_molecules(preopt)
for job_name, settings in GO.items():
    optimization(mols       = molecules, 
                 prefix     = job_name, 
                 settings   = settings, 
                 path       = os.path.join(GeoOpt,job_name)
                 )










import os
import re
from optparse import OptionParser

def get_options():
    parser = OptionParser()
    parser.add_option('-N', dest='n_samples', default=1, type='int', help="Number of samples to process. If -1 process all samples")
    parser.add_option('-d', dest='dataset', default=None, help="Dataset name in DAS")
    parser.add_option('-c', dest='config', default=None, help="CMSSW config file to run")
    parser.add_option('-n', dest='name', default='submit_crab', help="Submission file name")
    parser.add_option('-o', dest='outLFNDir', default='HGCAL_L1T_emulator_Feb22', help="Output LFN directory")
    parser.add_option('--storage_site', dest='storage_site', default='T2_UK_London_IC', help="User storage site")
    return parser.parse_args()
(opt,args) = get_options()

def write_sub(options):
    with open("%s_cfg.py"%options.name, "w") as fsub:
        fsub.write("from CRABClient.UserUtilities import config\n")
        fsub.write("config = config()\n\n")
        fsub.write("config.General.requestName = \'%s\'\n"%options.name)
        fsub.write("config.General.workArea = \'crab_area\'\n")
        fsub.write("config.General.transferOutputs = True\n")
        fsub.write("config.General.transferLogs = True\n\n")
        fsub.write("config.JobType.pluginName = \'Analysis\'\n")
        fsub.write("config.JobType.psetName = \'%s\'\n"%options.config)
        fsub.write("config.JobType.maxMemoryMB = 2500\n\n")
        fsub.write("config.Data.inputDataset = \'%s\'\n"%options.dataset)
        fsub.write("config.Data.inputDBS = \'global\'\n")
        fsub.write("config.Data.splitting = \'FileBased\'\n")
        fsub.write("config.Data.unitsPerJob = 1\n")
        fsub.write("NJOBS = %g\n"%options.n_samples)
        fsub.write("config.Data.totalUnits = config.Data.unitsPerJob * NJOBS\n")
        fsub.write("config.Data.outLFNDirBase = \'/store/user/%s/%s\'\n"%(os.environ['USER'],options.outLFNDir))
        fsub.write("config.Data.publication = True\n")
        fsub.write("config.Data.outputDatasetTag = \'%s\'\n\n"%(re.sub("submit_","",options.name)))
        fsub.write("config.Site.storageSite = \'%s\'"%opt.storage_site)
        fsub.close()
    

if opt.n_samples == -1:
    from Utilities.General.cmssw_das_client import get_data as das_query
    def safe_das_query( search, cmd ):
        output = das_query( search, cmd=cmd )
        if not 'data' in output:
            raise Exception('Your das query has not worked properly - check your proxy is valid')
        return output

    files = safe_das_query('file dataset=%s'%opt.dataset, cmd='dasgoclient')['data']
    opt.n_samples = len(files)

write_sub(opt)
     



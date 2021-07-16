import sys, csv, re
import pyopenms
import pandas as pd
from pyteomics.mztab import MzTab 
from dataclasses import dataclass

@dataclass
class Feature:
    spec_id: str
    mz: str
    z: str
    rt_mean: str
    seq: str
    scan: str

    def to_list(self):
        return [self.spec_id, self.mz, self.z, self.rt_mean, self.seq, self.scan, "0.0:1.0", "1.0"]


def scan2seq(percolator: pd.DataFrame) -> dict:
    scan = percolator['scan'].str.split(",")
    scan2idx_dict = {}
    last2nd = percolator.columns[-2] # last column is scan, last second is sequence
    for i, row in percolator.iterrows():
        scan = row['scan'].split(",")
        for s in scan:
            scan2idx_dict[s] = row[last2nd] # row['sequence']
    # for i, row in scan.iteritems():
    #     for s in row:
    #         scan2idx_dict[s] = i
    return scan2idx_dict
        

def parse(msdata, perlocator: pd.DataFrame, mgf, writer, sampleID=""):
    # Iterate through all spectra, skip all MS1 spectra and then write mgf format
    scan2seq_dict = scan2seq(perlocator)
    nr_ms2_spectra = 0
    for spectrum in msdata:
        if spectrum.getMSLevel() == 1:
            continue
        nr_ms2_spectra += 1

        scan = spectrum.getNativeID().split(" ")[-1][5:]
        
        if scan in scan2seq_dict:
            seq = scan2seq_dict[scan]
        else:
            seq = ''

        # drop problematic spectrums. these spectrums will make downstream work stop.
        try:     
            mz =  spectrum.getPrecursors()[0].getMZ()
            ch = spectrum.getPrecursors()[0].getCharge()
        except IndexError:
            continue
            # mgf.write("PEPMASS=unknown\n")
        ch = spectrum.getPrecursors()[0].getCharge()
        if ch <= 0: continue

        ## write mgf output
        mgf.write("\nBEGIN IONS\n")
        mgf.write("TITLE=%s\n" % spectrum.getNativeID())
        mgf.write("PEPMASS=%s\n" % spectrum.getPrecursors()[0].getMZ())
        mgf.write("CHARGE=%s\n" % ch)

        # FIXME: makesure the scan attribute is corret
        # re-number scan id (add sample id)
        scanID = f"F{sampleID}:{scan}"
        mgf.write("SCANS=%s\n" % scanID)
        mgf.write("RTINSECONDS=%s\n" % spectrum.getRT())
        #mgf.write("SEQ=%s\n"%spectrum.getPeptideIdentifications())

        for peak in spectrum:
            mgf.write("%s %s\n" % (peak.getMZ(), peak.getIntensity() ))
        mgf.write("END IONS\n")

        feature = Feature(spec_id=scanID, 
                          mz=mz, 
                          z=ch, 
                          rt_mean=spectrum.getRT(), 
                          seq=seq, 
                          scan=scanID,
                          )

        writer.writerow(feature.to_list())

    if nr_ms2_spectra == 0:
        print("Did not find any MS2 spectra in your input, thus the output file is empty!")

if __name__ == "__main__":
    if len(sys.argv) <= 5:
        print ("Usage: mzML2pointnovo.py in_mzML in_mzTab out_mgf out_feature.csv mzML_ID")
        sys.exit()

    ## read inputs
    msdata = pyopenms.MSExperiment()
    pyopenms.FileHandler().loadExperiment(sys.argv[1], msdata)
    # percolator = pd.read_table(sys.argv[2])
    # read mzTab data, only need psm
    mztab = MzTab(sys.argv[2])
    last = mztab.spectrum_match_table.columns[-1]
    cols = ['sequence','PSM_ID', 'accession', 'unique','modifications', 
            'retention_time', 'charge', 'exp_mass_to_charge',
             'calc_mass_to_charge', 'spectra_ref', last]
    percolator = mztab.spectrum_match_table.loc[:,cols]
    
    # extract scan 
    ### NOTE: if your scan is not the same format here, modify your code please
    pat = re.compile("scan=([0-9]+)")
    percolator['scan'] = percolator['spectra_ref'].apply(lambda x: ",".join(pat.findall(x))) # join multi scans if more than 1

    # need to reindex from 0-n after drop some rows
    # percolator = percolator[percolator['percolator q-value'] < 0.01].reset_index(drop=True)
    mgf = open(sys.argv[3], "w")
    ft = open(sys.argv[4], "w")
    sampleID = sys.argv[5]

    writer = csv.writer(ft, delimiter=',')
    header = ["spec_group_id","m/z","z","rt_mean","seq","scans","profile","feature area"]
    writer.writerow(header)
    
    # write mgf header
    mgf.write(f"COM=F{sampleID}:{sys.argv[1]}\n")
    mgf.write("ITOL=1\n")
    mgf.write("ITOLU=Da\n")
    mgf.write("CLE=NO_ENZYME\n")
    mgf.write("CHARGE=1,2,3\n")

    ### parse and write data
    flag = False
    parse(msdata, percolator, mgf, writer, sampleID)
    mgf.close()
    ft.close()
    print("Job Done")
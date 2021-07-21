import sys, csv, re
import pandas as pd
from pyteomics import mgf
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
    scan2idx_dict = {}
    percolator['scan'] = percolator['scan'].astype(str)
    for i, row in percolator.iterrows():
        scan = row['scan'].split(",")
        for s in scan:
            scan2idx_dict[int(s)] = row['sequence'] 
    return scan2idx_dict
        

def parse(msdata, perlocator: pd.DataFrame, omgf, writer, sampleID=""):
    """Iterate through all spectra, skip all MS1 spectra and then write mgf format
    """
    # scan to dict[int: string]
    scan2seq_dict = scan2seq(perlocator)
    scan = 0
    for spectrum in msdata:
        ## Title format: TITLE=fileName.scanNumber.scanNumber.charge.dta. 
        # or: TITLE=fileName.scanNumber.scanNumber.charge File:"fileNameWithExtension", NativeID:"controllerType=0 controllerNumber=1 scan=scanNumber"
        ## scan = spectrum['params']['title'].split(".")[1]
        
        ## NOTE: the scan number in the percolator file is just the index of spectrums of mgf !!!
        scan +=1
        if scan in scan2seq_dict:
            seq = scan2seq_dict[scan]
        else:
            seq = ''

        chs = spectrum['params']['charge']
        if len(chs) > 1: continue
        ch = chs[0]
        pepmass = spectrum['params']['pepmass'][0]

        ## write mgf output
        omgf.write("\nBEGIN IONS\n")
        omgf.write("TITLE=%s\n" % spectrum['params']['title'])
        omgf.write("PEPMASS=%s\n" % pepmass)
        omgf.write("CHARGE=%s\n" % ch)

        # FIXME: makesure the scan attribute is corret
        # re-number scan id (add sample id)
        scanID = f"F{sampleID}:{scan}"
        omgf.write("SCANS=%s\n" % scanID)
        omgf.write("RTINSECONDS=%s\n" % str(spectrum['params']['rtinseconds']))
        omgf.write("SEQ=%s\n"%(seq if seq != '' else 'UNKNOWN'))

        for m, intensity in zip(spectrum['m/z array'], spectrum['intensity array']):
            omgf.write("%s %s\n" % (m, intensity ))
        omgf.write("END IONS\n")

        feature = Feature(spec_id=scanID, 
                          mz=pepmass, 
                          z=int(ch), 
                          rt_mean=str(spectrum['params']['rtinseconds']), 
                          seq=seq, 
                          scan=scanID,
                          )

        writer.writerow(feature.to_list())


if __name__ == "__main__":
    if len(sys.argv) <= 5:
        print ("Usage: crux2mgf.py pid.mgf percolator.target.peptide.txt out_mgf out_feature.csv SAMPLE_ID")
        sys.exit()

    ## read inputs
    msdata = mgf.read(sys.argv[1])
    percolator = pd.read_table(sys.argv[2])
    # read mzTab data, only need psm

    # extract scan 
    # ### NOTE: if your scan is not the same format here, modify your code please
    # pat = re.compile("scan=([0-9]+)")
    # percolator['scan'] = percolator['spectra_ref'].apply(lambda x: ",".join(pat.findall(x))) # join multi scans if more than 1

    # need to reindex from 0-n after drop some rows
    # percolator = percolator[percolator['percolator q-value'] < 0.01].reset_index(drop=True)
    mgf_out = open(sys.argv[3], "w")
    ft = open(sys.argv[4], "w")
    sampleID = sys.argv[5]

    writer = csv.writer(ft, delimiter=',')
    header = ["spec_group_id","m/z","z","rt_mean","seq","scans","profile","feature area"]
    writer.writerow(header)
    
    # write mgf header
    mgf_out.write(f"COM=F{sampleID}:{sys.argv[1]}\n")
    mgf_out.write("ITOL=1\n")
    mgf_out.write("ITOLU=Da\n")
    mgf_out.write("CLE=NO_ENZYME\n")
    mgf_out.write("CHARGE=1,2,3\n")

    ### parse and write data
    flag = False
    parse(msdata, percolator, mgf_out, writer, sampleID)
    mgf_out.close()
    ft.close()
    msdata.close()
    print("Job Done")